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

package mixer.utils.cleaning;

import javastraw.expected.ZScoreArray;
import javastraw.reader.basics.Chromosome;
import mixer.clt.ParallelizedMixerTools;
import mixer.utils.common.LogTools;
import mixer.utils.common.ParallelizedStatTools;
import mixer.utils.common.ZScoreTools;
import mixer.utils.drive.Mappings;
import mixer.utils.drive.MatrixAndWeight;
import mixer.utils.magic.FinalScale;
import mixer.utils.magic.SymmLLInterMatrix;

import java.util.concurrent.atomic.AtomicInteger;

public class MatrixPreprocessor {

    public static void clean(MatrixAndWeight matrix, Mappings mappings, Chromosome[] chromosomes,
                             boolean doScale, boolean doZscore, boolean doSecondaryCompression,
                             long seed) {
        int numRows = mappings.getNumRows();
        int numCols = mappings.getNumCols();

        int[] totalDistribution = getSumOfAllLoci(mappings, numCols, chromosomes);
        System.arraycopy(totalDistribution, 0, matrix.weights, 0, numCols);

        if (doScale) {
            //scaleMatrixColumns(matrix.matrix, totalDistribution);
            matrix.matrix = FinalScale.scaleMatrix(new SymmLLInterMatrix(matrix.matrix),
                    createTargetVector(totalDistribution, numRows, numCols, matrix.matrix));
        }

        ParallelizedStatTools.setZerosToNan(matrix.matrix);
        //ParallelizedStatTools.scaleDown(matrix.matrix, matrix.weights);
        LogTools.simpleLogWithCleanup(matrix.matrix, Float.NaN);
        removeHighGlobalThresh(matrix.matrix, 5);
        normalize(matrix.matrix, -3, 3);
        LogTools.simpleExpm1(matrix.matrix);

        // todo normalizeMatrix(matrix, mappings, chromosomes);
        //float[] coverage = new float[numRows];
        //updateCoverage(matrix, coverage);
        //scaleCoverage(coverage);
        //todo matrix = EmptyRowCleaner.cleanUpMatrix(matrix, rowIndexToIntervalMap, coverage);

        //LogTools.simpleExpm1(matrix);

        if (doZscore) {
            ZScoreTools.inPlaceZscorePositivesDownColAndSetZeroToNan(matrix.matrix);
        }

        if (doSecondaryCompression) {
            ClusteringCompressor.process(matrix, seed);
        }
    }

    private static long[] createTargetVector(int[] weights, int numRows, int numCols, float[][] matrix) {
        long[] colTarget = getSummedWeightsColumns(weights, matrix, numCols);
        long[] rowTarget = getSummedWeightsRows(weights, matrix, numRows);

        long[] target = new long[numCols + numRows];
        System.arraycopy(colTarget, 0, target, 0, numCols);
        System.arraycopy(rowTarget, 0, target, numCols, numRows);
        return target;
    }

    private static long[] getSummedWeightsColumns(int[] weights, float[][] matrix, int numCols) {
        long[] colSums = new long[numCols];
        for (int i = 0; i < matrix.length; i++) {
            for (int j = 0; j < matrix[i].length; j++) {
                if (matrix[i][j] > -1) {
                    colSums[j] += weights[j];
                }
            }
        }
        return colSums;
    }

    private static long[] getSummedWeightsRows(int[] weights, float[][] matrix, int numRows) {
        long[] rowSums = new long[numRows];
        for (int i = 0; i < matrix.length; i++) {
            for (int j = 0; j < matrix[i].length; j++) {
                if (matrix[i][j] > -1) {
                    rowSums[i] += weights[j];
                }
            }
        }
        return rowSums;
    }

    private static int[] getSumOfAllLoci(Mappings mappings, int numCols, Chromosome[] chromosomes) {
        int[] totalLoci = new int[numCols];
        for (Chromosome chromosome : chromosomes) {
            int[] row = mappings.getDistributionForChrom(chromosome);
            for (int z = 0; z < row.length; z++) {
                totalLoci[z] += row[z];
            }
        }
        return totalLoci;
    }

    private static void normalizeMatrix(float[][] matrix, Mappings mappings, Chromosome[] chromosomes) {

        for (Chromosome chromosome : chromosomes) {
            int[] globalIndices = mappings.getGlobalIndex(chromosome);
            int[] divisor = mappings.getDistributionForChrom(chromosome);
            for (int i : globalIndices) {
                if (i > -1) {
                    divide(matrix[i], divisor);
                }
            }
        }
    }

    private static void divide(float[] row, int[] totalLoci) {
        for (int k = 0; k < row.length; k++) {
            if (totalLoci[k] > 0) {
                row[k] = row[k] / totalLoci[k];
            } else if (row[k] > 0) {
                System.err.println("Impossible situation reached: row val: " + row[k] + " but expect no entries: " + totalLoci[k]);
            }
        }
    }

    private static void scaleMatrixColumns(float[][] matrix, int[] scalars) {
        for (int i = 0; i < matrix.length; i++) {
            for (int z = 0; z < matrix[i].length; z++) {
                matrix[i][z] *= scalars[z];
            }
        }
    }

    private static void removeHighGlobalThresh(float[][] data, int cutoff) {
        ZScoreArray zscores = ZScoreTools.getZscores(data);

        AtomicInteger totalNumFixed = new AtomicInteger();
        AtomicInteger index = new AtomicInteger(0);
        ParallelizedMixerTools.launchParallelizedCode(() -> {
            int i = index.getAndIncrement();
            int numFixed = 0;
            while (i < data.length) {
                for (int j = 0; j < data[i].length; j++) {
                    if (data[i][j] > 0) {
                        if (zscores.getZscore(j, data[i][j]) > cutoff) {
                            data[i][j] = Float.NaN;
                            numFixed++;
                        }
                    }
                }
                i = index.getAndIncrement();
            }
            totalNumFixed.addAndGet(numFixed);
        });
    }

    private static void normalize(float[][] data, int lowCutOff, int highCutOff) {
        ZScoreArray zscores = ZScoreTools.getZscores(data);
        fixToNormalRange(data, zscores, lowCutOff, highCutOff);
    }

    private static void fixToNormalRange(float[][] data, ZScoreArray zscores, int lowCutOff, int highCutOff) {
        AtomicInteger totalNumFixed = new AtomicInteger();
        AtomicInteger index = new AtomicInteger(0);
        ParallelizedMixerTools.launchParallelizedCode(() -> {
            int i = index.getAndIncrement();
            int numFixed = 0;
            while (i < data.length) {
                for (int j = 0; j < data[i].length; j++) {
                    if (data[i][j] > 0) {
                        double zscore = zscores.getZscore(j, data[i][j]);
                        if (zscore < lowCutOff || zscore > highCutOff) { //
                            data[i][j] = Float.NaN;
                            numFixed++;
                        }
                    }
                }
                i = index.getAndIncrement();
            }
            totalNumFixed.addAndGet(numFixed);
        });
    }
}
