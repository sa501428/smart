/*
 * The MIT License (MIT)
 *
 * Copyright (c) 2011-2020 Rice University, Baylor College of Medicine, Aiden Lab
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

package mixer.utils.shuffle;

import javastraw.featurelist.GenomeWideList;
import javastraw.reader.Dataset;
import javastraw.reader.basics.Chromosome;
import javastraw.type.NormalizationType;
import mixer.utils.common.FloatMatrixTools;
import mixer.utils.slice.structures.SliceUtils;
import mixer.utils.slice.structures.SubcompartmentInterval;

import java.io.File;
import java.io.FileWriter;
import java.util.*;

public class ShuffleMatrix {

    private final Random generator = new Random(0);
    private final Dataset ds;
    private final NormalizationType norm;
    private final int compressionFactor;
    private final int resolution;
    private final int numRounds = 50;
    private final InterOnlyMatrix.InterMapType[] mapTypes = {InterOnlyMatrix.InterMapType.ODDS_VS_EVENS,
            InterOnlyMatrix.InterMapType.SKIP_BY_TWOS, InterOnlyMatrix.InterMapType.FIRST_HALF_VS_SECOND_HALF};
    private final double[][] baseline = new double[mapTypes.length][numRounds];
    private final double[][] scores = new double[mapTypes.length][numRounds];
    final double[] ratios = new double[mapTypes.length];

    public ShuffleMatrix(Dataset ds, NormalizationType norm, int resolution, int compressionFactor) {
        this.resolution = resolution;
        this.compressionFactor = compressionFactor;
        this.ds = ds;
        this.norm = norm;
    }

    public void runAnalysis(GenomeWideList<SubcompartmentInterval> subcompartments, File outfolder) {
        SliceUtils.collapseGWList(subcompartments);
        for (int y = 0; y < mapTypes.length; y++) {
            final InterOnlyMatrix interMatrix = new InterOnlyMatrix(ds, norm, resolution, mapTypes[y]);
            final Map<Integer, List<Integer>> clusterToRowIndices = new HashMap<>();
            final Map<Integer, List<Integer>> clusterToColIndices = new HashMap<>();
            populateClusterToIndexMaps(interMatrix, subcompartments, clusterToRowIndices, clusterToColIndices);

            float[][] baselineM = determineScore(interMatrix, clusterToRowIndices, clusterToColIndices, true, baseline[y]);
            float[][] shuffleM = determineScore(interMatrix, clusterToRowIndices, clusterToColIndices, false, scores[y]);

            FloatMatrixTools.saveMatrixTextNumpy(new File(outfolder, mapTypes[y].toString() + "_baseline.npy").getAbsolutePath(), baselineM);
            FloatMatrixTools.saveMatrixTextNumpy(new File(outfolder, mapTypes[y].toString() + "_shuffle.npy").getAbsolutePath(), shuffleM);
            FloatMatrixTools.saveMatrixTextNumpy(new File(outfolder, mapTypes[y].toString() + "_matrix.npy").getAbsolutePath(), interMatrix.getMatrix());

            double sumN = 0, sumD = 9;
            for (int i = 0; i < numRounds; i++) {
                sumN += scores[y][i];
                sumD += baseline[y][i];
            }
            ratios[y] = sumN / sumD;
        }
    }

    public void savePlotsAndResults(File outfolder) {
        try {
            FileWriter myWriter = new FileWriter(new File(outfolder, "scores.txt"));
            for (int y = 0; y < mapTypes.length; y++) {
                myWriter.write(mapTypes[y].toString() + " score: " + ratios[y] + "\n");
            }
            myWriter.close();
        } catch (Exception ee) {
            System.err.println("Unable to write results to text file");
        }
        //System.out.println("Scores: " + sumN + "  baseline: " + sumD + "  ratio:" + (sumN / sumD));
        //TTest test = new TTest();
        //System.out.println("Ttest " + test.pairedTTest(baseline, scores) / 2);
        //System.out.println("Ttest " + test.pairedTTest(baseline, scores, 0.1));

    }

    private float[][] determineScore(InterOnlyMatrix interMatrix,
                                     Map<Integer, List<Integer>> clusterToRowIndices,
                                     Map<Integer, List<Integer>> clusterToColIndices,
                                     boolean randomizeGW, double[] vector) {
        float[][] aggregate = null;
        for (int k = 0; k < numRounds; k++) {
            List<Integer> allRowIndices = getShuffledByClusterIndices(clusterToRowIndices, randomizeGW);
            List<Integer> allColIndices = getShuffledByClusterIndices(clusterToColIndices, randomizeGW);
            float[][] matrix = getShuffledMatrix(interMatrix, allRowIndices, allColIndices);

            if (aggregate == null) {
                aggregate = matrix;
            } else {
                FloatMatrixTools.addBToA(aggregate, matrix);
            }
            vector[k] = scoreMatrix(matrix);
        }

        FloatMatrixTools.scaleBy(aggregate, 1f / ((float) numRounds));
        return aggregate;
    }

    private List<Integer> getShuffledByClusterIndices(Map<Integer, List<Integer>> clusterToIndices, boolean randomizeGW) {
        List<Integer> allIndices = new ArrayList<>();

        List<Integer> order = new ArrayList<>(clusterToIndices.keySet());
        Collections.sort(order);

        for (Integer index : order) {
            List<Integer> indexList = clusterToIndices.get(index);
            Collections.shuffle(indexList, generator);
            int numToUse = (indexList.size() / compressionFactor) * compressionFactor;
            for (int z = 0; z < numToUse; z++) {
                allIndices.add(indexList.get(z));
            }
        }
        if (randomizeGW) {
            Collections.shuffle(allIndices, generator);
        }
        return allIndices;
    }

    private float[][] getShuffledMatrix(InterOnlyMatrix interMatrix, List<Integer> allRowIndices, List<Integer> allColIndices) {
        int numRows = allRowIndices.size() / compressionFactor;
        int numCols = allColIndices.size() / compressionFactor;
        int numRowsKept = numRows * compressionFactor;
        int numColsKept = numCols * compressionFactor;
        float[][] original = interMatrix.getMatrix();

        float[][] result = new float[numRows][numCols];
        for (int i = 0; i < numRowsKept; i++) {
            final int i0 = allRowIndices.get(i);
            for (int j = 0; j < numColsKept; j++) {
                final int j0 = allColIndices.get(j);
                result[i / compressionFactor][j / compressionFactor] += original[i0][j0];
            }
        }

        return result;
    }

    private double scoreMatrix(float[][] matrix) {
        double diff = 0;
        for (int i = 0; i < matrix.length - 1; i++) {
            for (int j = 0; j < matrix[i].length - 1; j++) {
                diff += Math.abs(matrix[i][j] - matrix[i + 1][j]);
                diff += Math.abs(matrix[i][j] - matrix[i][j + 1]);
                diff += Math.abs(matrix[i][j] - matrix[i + 1][j + 1]);
                diff += Math.abs(matrix[i + 1][j] - matrix[i][j + 1]);
            }
        }
        return diff / ((matrix.length - 1) * (matrix[0].length - 1));
    }

    private void populateClusterToIndexMaps(InterOnlyMatrix interMatrix, GenomeWideList<SubcompartmentInterval> intraSubcompartments, Map<Integer, List<Integer>> clusterToRowIndices, Map<Integer, List<Integer>> clusterToColIndices) {
        populateCluster(interMatrix.getRowChromosomes(), interMatrix.getRowOffsets(), intraSubcompartments, clusterToRowIndices);
        populateCluster(interMatrix.getColChromosomes(), interMatrix.getColOffsets(), intraSubcompartments, clusterToColIndices);
    }

    private void populateCluster(Chromosome[] chromosomes, int[] offsets, GenomeWideList<SubcompartmentInterval> subcompartments,
                                 Map<Integer, List<Integer>> clusterToIndices) {
        for (int x = 0; x < chromosomes.length; x++) {
            Chromosome chrom = chromosomes[x];
            List<SubcompartmentInterval> intervalList = subcompartments.getFeatures("" + chrom.getIndex());
            for (SubcompartmentInterval interval : intervalList) {
                int xStart = interval.getX1() / resolution;
                int xEnd = interval.getX2() / resolution;
                int clusterID = interval.getClusterID();

                List<Integer> tempList = new ArrayList<>();
                for (int k = xStart; k < xEnd; k++) {
                    final int actualPosition = k + offsets[x];
                    tempList.add(actualPosition);
                }

                if (clusterToIndices.containsKey(clusterID)) {
                    clusterToIndices.get(clusterID).addAll(tempList);
                } else {
                    clusterToIndices.put(clusterID, tempList);
                }
            }
        }
    }
}
