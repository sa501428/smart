/*
 * The MIT License (MIT)
 *
 * Copyright (c) 2011-2022 Rice University, Baylor College of Medicine, Aiden Lab
 *
 * Permission is hereby granted, free of charge, to any person obtaining a copy
 * of this software and associated documentation files (the "Software"), to deal
 * in the Software without restriction, including without limitation the rights
 * to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the Software is
 * furnished to do so, subject to the following conditions:
 *
 * The above copyright notice and this permission notice shall be included in
 * all copies or substantial portions of the Software.
 *
 *  THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 *  IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 *  FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 *  AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 *  LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 *  OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
 *  THE SOFTWARE.
 */

package mixer.utils.refinement;

import javastraw.feature1D.GenomeWide1DList;
import javastraw.reader.basics.ChromosomeHandler;
import javastraw.tools.MatrixTools;
import mixer.utils.drive.FinalMatrix;
import mixer.utils.similaritymeasures.SimilarityMetric;
import mixer.utils.tracks.Concensus2DTools;
import mixer.utils.tracks.Concensus3DTools;
import mixer.utils.tracks.SliceUtils;
import mixer.utils.tracks.SubcompartmentInterval;

import java.io.File;
import java.util.*;

public class IterativeRefinement {

    public static void resolve(FinalMatrix matrix, float[][] backupForL2, String prefix, int numClusters,
                               int z, int[] assignments1, int[] assignments2,
                               ChromosomeHandler handler, File outputDirectory) {
        Map<String, List<Integer>> superClusterToIndices = new HashMap<>();
        int[][] summary = populateSummary(numClusters, assignments1, assignments2,
                superClusterToIndices);

        MatrixTools.saveMatrixTextNumpy(new File(outputDirectory, numClusters + "_summary.npy").getAbsolutePath(),
                summary);

        List<String> finalKeys = new ArrayList<>(numClusters);
        for (int q = 0; q < numClusters; q++) {
            int[] coords = Concensus2DTools.getMaxCoordinates(summary);
            if (coords != null) {
                clearSectionSlices(summary, coords);
                finalKeys.add(makeKey(coords));
            } else {
                System.err.println("Weird error getting consensus; skipping");
                return;
            }
        }

        int[] hubAssignment = Concensus3DTools.assignHubs(matrix, numClusters, superClusterToIndices, finalKeys);
        assignRemainingEntries(hubAssignment, matrix.matrix, backupForL2, numClusters);

        GenomeWide1DList<SubcompartmentInterval> finalCompartments = matrix.getClusteringResult(hubAssignment, handler);
        SliceUtils.collapseGWList(finalCompartments);
        File outBedFile = new File(outputDirectory, prefix + "_k" + numClusters + "_consensus_clusters.bed");
        finalCompartments.simpleExport(outBedFile);
    }

    public static void assignRemainingEntries(int[] hubAssignment, float[][] matrixL1, float[][] matrixL2,
                                              int numClusters) {
        List<Integer> unassigned = getUnassignedIndices(hubAssignment);
        List<NearestNeighbors> neighborhoods = getNearestNeighbors(unassigned, matrixL1, matrixL2, hubAssignment);
        sortAndAssignRemaining(neighborhoods, hubAssignment, numClusters);
    }

    public static void assignRemainingEntries(int[] hubAssignment, float[][] matrix, SimilarityMetric metric,
                                              int numClusters) {
        List<Integer> unassigned = getUnassignedIndices(hubAssignment);
        List<NearestNeighbors> neighborhoods = getNearestNeighbors(unassigned, matrix, metric, hubAssignment);
        sortAndAssignRemaining(neighborhoods, hubAssignment, numClusters);
    }

    private static void sortAndAssignRemaining(List<NearestNeighbors> neighborhoods, int[] hubAssignment, int numClusters) {
        neighborhoods.sort(Comparator.comparingDouble(NearestNeighbors::getPercentUnassigned));
        for (NearestNeighbors neighborhood : neighborhoods) {
            hubAssignment[neighborhood.getIndex()] = neighborhood.getMajorityNeighborAssignment(hubAssignment,
                    numClusters);
        }
    }

    private static List<NearestNeighbors> getNearestNeighbors(List<Integer> unassigned, float[][] matrix,
                                                              SimilarityMetric metric, int[] hubAssignment) {
        List<NearestNeighbors> neighborhoods = new ArrayList<>(unassigned.size());
        for (int i : unassigned) {
            neighborhoods.add(new NearestNeighbors(i, matrix, metric, hubAssignment));
        }
        return neighborhoods;
    }

    private static List<NearestNeighbors> getNearestNeighbors(List<Integer> unassigned,
                                                              float[][] matrixL1, float[][] matrixL2, int[] hubAssignment) {
        List<NearestNeighbors> neighborhoods = new ArrayList<>(unassigned.size());
        for (int i : unassigned) {
            neighborhoods.add(new NearestNeighbors(i, matrixL1, matrixL2, hubAssignment));
        }
        return neighborhoods;
    }

    private static List<Integer> getUnassignedIndices(int[] hubAssignment) {
        List<Integer> indices = new LinkedList<>();
        for (int i = 0; i < hubAssignment.length; i++) {
            if (hubAssignment[i] < 0) {
                indices.add(i);
            }
        }
        return indices;
    }

    private static int[][] populateSummary(int numClusters, int[] assignments1, int[] assignments2,
                                           Map<String, List<Integer>> superClusterToIndices) {
        int[][] summary = new int[numClusters][numClusters];
        for (int x = 0; x < assignments1.length; x++) {
            int i = assignments1[x];
            int j = assignments2[x];
            if (i < 0 || j < 0) continue;
            summary[i][j]++;
            String key = makeKey(i, j);
            if (!superClusterToIndices.containsKey(key)) {
                superClusterToIndices.put(key, new LinkedList<>());
            }
            superClusterToIndices.get(key).add(x);
        }
        return summary;
    }

    private static void clearSectionSlices(int[][] matrix, int[] coords) {
        for (int j = 0; j < matrix.length; j++) {
            matrix[coords[0]][j] = -1;
        }
        for (int i = 0; i < matrix.length; i++) {
            matrix[i][coords[1]] = -1;
        }
    }

    private static String makeKey(int[] coords) {
        return coords[0] + "_" + coords[1];
    }

    private static String makeKey(int i, int j) {
        return i + "_" + j;
    }

    /*
    public static void refine(FinalMatrix fMatrix, float[][] data, SimilarityMetric metric,
                              int numClusters, int[] assignments,
                              ChromosomeHandler handler, File outputDirectory, String outputName) {

        refineEntries(assignments, data, numClusters, metric);
        GenomeWide1DList<SubcompartmentInterval> finalCompartments = fMatrix.getClusteringResult(assignments, handler);
        SliceUtils.collapseGWList(finalCompartments);
        File outBedFile = new File(outputDirectory, outputName);
        finalCompartments.simpleExport(outBedFile);
    }

    private static void refineEntries(int[] assignment, float[][] matrix, int numClusters, SimilarityMetric metric) {
        List<NearestNeighbors> neighborhoods = getAllNearestNeighbors(matrix, assignment, metric);
        int count = 0;
        for (NearestNeighbors neighborhood : neighborhoods) {
            if (++count % 1000 == 0) System.out.print(".");
            int newAssignment = neighborhood.getMajorityNeighborAssignment(assignment,
                    numClusters, 0.55f);
            if (newAssignment > -1) {
                assignment[neighborhood.getIndex()] = newAssignment;
            }
        }
        System.out.println(".");
    }

    private static List<NearestNeighbors> getAllNearestNeighbors(float[][] matrix, int[] assignment,
                                                                 SimilarityMetric metric) {
        List<NearestNeighbors> neighborhoods = new ArrayList<>(matrix.length);
        for (int i = 0; i < matrix.length; i++) {
            neighborhoods.add(new NearestNeighbors(i, matrix, metric, nullx));
        }
        return neighborhoods;
    }
    */
}
