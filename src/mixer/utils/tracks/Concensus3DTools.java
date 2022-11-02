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

package mixer.utils.tracks;

import javastraw.reader.basics.ChromosomeHandler;
import mixer.utils.drive.FinalMatrix;
import mixer.utils.kmeans.ClusteringMagic;
import mixer.utils.refinement.IterativeRefinement;
import mixer.utils.similaritymeasures.RobustEuclideanDistance;
import mixer.utils.similaritymeasures.RobustManhattanDistance;
import mixer.utils.similaritymeasures.SimilarityMetric;

import java.util.*;

public class Concensus3DTools {
    public static void resolve(int[][] results, ClusteringMagic clusteringMagic,
                               int z, String prefix, boolean useKMedians, int numClusters,
                               FinalMatrix matrix, ChromosomeHandler handler) {

        Map<String, List<Integer>> superClusterToIndices = new HashMap<>();
        int[][][] summary = populateSummary(numClusters, matrix.getNumRows(), results, superClusterToIndices);

        List<String> finalKeys = new ArrayList<>(numClusters);
        for (int q = 0; q < numClusters; q++) {
            int[] coords = getMaxCoordinates(summary);
            if (coords != null) {
                clearSectionSlices(summary, coords);
                finalKeys.add(makeKey(coords));
            } else {
                System.err.println("Could not get consensus; skipping refinement");
                return;
            }
        }

        SimilarityMetric metric = RobustEuclideanDistance.SINGLETON;
        if (useKMedians) metric = RobustManhattanDistance.SINGLETON;

        int[] hubAssignment = assignHubs(matrix, numClusters, superClusterToIndices, finalKeys);
        IterativeRefinement.assignRemainingEntries(hubAssignment, matrix.matrix, metric, numClusters);

        clusteringMagic.exportKMeansClusteringResults(z, prefix + "_refined", useKMedians,
                matrix.getClusteringResult(hubAssignment, handler), null);
    }

    public static int[] assignHubs(FinalMatrix matrix, int numClusters, Map<String, List<Integer>> superClusterToIndices, List<String> finalKeys) {
        int[] hubAssignment = new int[matrix.getNumRows()];
        Arrays.fill(hubAssignment, -1);
        for (int q = 0; q < numClusters; q++) {
            for (int i : superClusterToIndices.get(finalKeys.get(q))) {
                hubAssignment[i] = q;
            }
        }
        return hubAssignment;
    }

    private static int[][][] populateSummary(int numClusters, int numRows, int[][] results,
                                             Map<String, List<Integer>> superClusterToIndices) {
        int[][][] summary = new int[numClusters][numClusters][numClusters];

        for (int x = 0; x < numRows; x++) {
            int i = results[0][x];
            int j = results[1][x];
            int k = results[2][x];
            //if(i < 0 || j < 0 || k < 0) continue;
            summary[i][j][k]++;
            String key = makeKey(i, j, k);
            if (!superClusterToIndices.containsKey(key)) {
                superClusterToIndices.put(key, new LinkedList<>());
            }
            superClusterToIndices.get(key).add(x);
        }
        return summary;
    }

    private static void clearSectionSlices(int[][][] tensor, int[] coords) {
        for (int j = 0; j < tensor.length; j++) {
            for (int k = 0; k < tensor.length; k++) {
                tensor[coords[0]][j][k] = -1;
            }
        }
        for (int i = 0; i < tensor.length; i++) {
            for (int k = 0; k < tensor.length; k++) {
                tensor[i][coords[1]][k] = -1;
            }
        }
        for (int i = 0; i < tensor.length; i++) {
            for (int j = 0; j < tensor.length; j++) {
                tensor[i][j][coords[2]] = -1;
            }
        }
    }

    private static int[] getMaxCoordinates(int[][][] tensor) {
        int maxVal = 0;
        int[] coords = null;
        for (int i = 0; i < tensor.length; i++) {
            for (int j = 0; j < tensor.length; j++) {
                for (int k = 0; k < tensor.length; k++) {
                    if (tensor[i][j][k] > maxVal) {
                        maxVal = tensor[i][j][k];
                        coords = new int[]{i, j, k};
                    }
                }
            }
        }
        return coords;
    }

    private static String makeKey(int[] coords) {
        return coords[0] + "_" + coords[1] + "_" + coords[2];
    }

    private static String makeKey(int i, int j, int k) {
        return i + "_" + j + "_" + k;
    }
}
