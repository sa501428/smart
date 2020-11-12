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
import java.util.*;

public class ShuffleMatrix {

    private final InterOnlyMatrix interMatrix;
    private final Map<Integer, List<Integer>> clusterToRowIndices = new HashMap<>();
    private final Map<Integer, List<Integer>> clusterToColIndices = new HashMap<>();
    private final Random generator = new Random(0);
    private final int compressionFactor;
    private final int resolution;
    private final int numRounds = 10;
    private final double[] baseline = new double[numRounds];
    private final double[] scores = new double[numRounds];
    private final Dataset ds;
    private final File outfolder;

    public ShuffleMatrix(Dataset ds, NormalizationType norm, int resolution,
                         GenomeWideList<SubcompartmentInterval> subcompartments, int compressionFactor,
                         File outfolder) {
        this.ds = ds;
        this.resolution = resolution;
        this.compressionFactor = compressionFactor;
        this.outfolder = outfolder;
        interMatrix = new InterOnlyMatrix(ds, norm, resolution, InterOnlyMatrix.InterMapType.ODDS_VS_EVENS);


        SliceUtils.collapseGWList(subcompartments);
        populateClusterToIndexMaps(subcompartments);

        determineBaseline();
        determineScore();

        double sumN = 0, sumD = 0;
        for (int i = 0; i < numRounds; i++) {
            sumN += scores[i];
            sumD += baseline[i];
        }
        System.out.println("Scores: " + sumN + "  baseline: " + sumD + "  ratio:" + (sumN / sumD));
    }

    private void determineScore() {
        for (int k = 0; k < numRounds; k++) {
            List<Integer> allRowIndices = getShuffledByClusterIndices(clusterToRowIndices);
            List<Integer> allColIndices = getShuffledByClusterIndices(clusterToColIndices);
            float[][] matrix = getShuffledMatrix(allRowIndices, allColIndices);

            File temp = new File(outfolder, "scored" + k + ".npy");
            FloatMatrixTools.saveMatrixTextNumpy(temp.getAbsolutePath(), matrix);

            scores[k] = scoreMatrix(matrix);
            System.out.println("Score " + scores[k]);
        }
    }

    private List<Integer> getShuffledByClusterIndices(Map<Integer, List<Integer>> clusterToIndices) {
        List<Integer> allIndices = new ArrayList<>();
        for (List<Integer> indexList : clusterToIndices.values()) {
            Collections.shuffle(indexList, generator);
            int numToUse = (indexList.size() / compressionFactor) * compressionFactor;
            for (int z = 0; z < numToUse; z++) {
                allIndices.add(indexList.get(z));
            }
        }
        return allIndices;
    }

    private void determineBaseline() {
        List<Integer> allRowIndices = new ArrayList<>();
        for (List<Integer> rowList : clusterToRowIndices.values()) {
            allRowIndices.addAll(rowList);
        }

        List<Integer> allColIndices = new ArrayList<>();
        for (List<Integer> colList : clusterToColIndices.values()) {
            allColIndices.addAll(colList);
        }

        for (int k = 0; k < numRounds; k++) {
            // create random permutation
            Collections.shuffle(allRowIndices, generator);
            Collections.shuffle(allColIndices, generator);
            float[][] matrix = getShuffledMatrix(allRowIndices, allColIndices);
            baseline[k] = scoreMatrix(matrix);

            File temp = new File(outfolder, "baseline" + k + ".npy");
            FloatMatrixTools.saveMatrixTextNumpy(temp.getAbsolutePath(), matrix);

            System.out.println("Cost " + baseline[k]);
        }
    }

    private float[][] getShuffledMatrix(List<Integer> allRowIndices, List<Integer> allColIndices) {
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

    private void populateClusterToIndexMaps(GenomeWideList<SubcompartmentInterval> intraSubcompartments) {
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
