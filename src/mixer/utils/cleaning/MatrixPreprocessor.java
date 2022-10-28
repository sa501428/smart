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

import javastraw.expected.Welford;
import javastraw.expected.ZScoreArray;
import javastraw.expected.Zscore;
import javastraw.reader.basics.Chromosome;
import javastraw.tools.ParallelizationTools;
import mixer.utils.common.SimpleArray2DTools;
import mixer.utils.common.ZScoreTools;
import mixer.utils.drive.MatrixAndWeight;
import mixer.utils.transform.MatrixTransform;

import java.util.concurrent.atomic.AtomicInteger;

public class MatrixPreprocessor {

    public static MatrixAndWeight clean(MatrixAndWeight matrix, Chromosome[] chromosomes,
                                        int cutoff, boolean useExp, boolean useZscore) {

        SimpleArray2DTools.setZerosToNan(matrix.matrix);

        matrix.updateWeights(chromosomes);

        SimpleArray2DTools.simpleLogWithCleanup(matrix.matrix, Float.NaN);
        removeHighGlobalThresh(matrix.matrix, 5);
        makeColumnsHaveNormalRange(matrix.matrix, -cutoff, cutoff);
        if (useExp) {
            SimpleArray2DTools.simpleExpm1(matrix.matrix);
        }
        if (useZscore) {
            ZScoreTools.inPlaceZscorePositivesDownColAndSetBelowThreshToNan(matrix.matrix, 0);
        }

        matrix.removeAllNanRows();
        return matrix;
    }

    private static void removeHighGlobalThresh(float[][] data, int cutoff) {
        ZScoreArray zscores = ZScoreTools.getZscores(data, 0);

        AtomicInteger totalNumFixed = new AtomicInteger();
        AtomicInteger index = new AtomicInteger(0);
        ParallelizationTools.launchParallelizedCode(() -> {
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

    private static void makeColumnsHaveNormalRange(float[][] data, int lowCutOff, int highCutOff) {
        ZScoreArray zscores = ZScoreTools.getZscores(data, 0);
        fixColumnsToNormalRange(data, zscores, lowCutOff, highCutOff);
    }

    private static void fixColumnsToNormalRange(float[][] data, ZScoreArray zscores, int lowCutOff, int highCutOff) {
        AtomicInteger totalNumFixed = new AtomicInteger();
        AtomicInteger index = new AtomicInteger(0);
        ParallelizationTools.launchParallelizedCode(() -> {
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

    public static MatrixAndWeight clean2(MatrixAndWeight matrix, Chromosome[] chromosomes,
                                         boolean useLog, boolean doColumnZscore, boolean doGlobalThresholding,
                                         boolean setZeroToNan, boolean doRowZscoreWithThreshold,
                                         boolean restoreEXP, boolean useCosineCompression) {
        matrix.updateWeights(chromosomes);
        matrix.divideColumnsByWeights();
        if (useLog) {
            SimpleArray2DTools.simpleLogWithCleanup(matrix.matrix, Float.NaN);
        }

        if (doGlobalThresholding) {
            // any z higher than 5 set to nan
            thresholdGlobally(matrix.matrix, 5);
        }

        if (setZeroToNan) {
            cleanUpZerosAndInfs(matrix.matrix);
        }

        if (doRowZscoreWithThreshold) {
            MatrixTransform.zscoreByRows(matrix.matrix, 3);
        }

        if (restoreEXP) {
            SimpleArray2DTools.simpleExpm1(matrix.matrix);
        }

        if (doColumnZscore) {
            ZScoreTools.inPlaceZscorePositivesDownColAndSetBelowThreshToNan(matrix.matrix, -100);
        }
        matrix.removeAllNanRows();

        if (useCosineCompression) {
            matrix.matrix = SimilarityMatrixTools.getCompressedCosineSimilarityMatrix(matrix.matrix,
                    getNumCentroids(matrix.getNumRows(), matrix.getNumCols()), 0);
        }
        return matrix;
    }

    private static int getNumCentroids(int numRows, int numCols) {
        return Math.min(numCols, numRows / 20);
    }

    private static void thresholdGlobally(float[][] matrix, int upperLimit) {
        Welford welford = new Welford();
        for (int r = 0; r < matrix.length; r++) {
            for (int c = 0; c < matrix[r].length; c++) {
                if (matrix[r][c] > 0) {
                    welford.addValue(matrix[r][c]);
                }
            }
        }
        Zscore zscore = welford.getZscore();
        for (int r = 0; r < matrix.length; r++) {
            for (int c = 0; c < matrix[r].length; c++) {
                if (zscore.getZscore(matrix[r][c]) > upperLimit) {
                    matrix[r][c] = Float.NaN;
                }
            }
        }
    }

    private static void cleanUpZerosAndInfs(float[][] matrix) {
        for (int r = 0; r < matrix.length; r++) {
            for (int c = 0; c < matrix[r].length; c++) {
                float val = matrix[r][c];
                if (val > 0) {
                    if (Float.isInfinite(val)) {
                        matrix[r][c] = Float.NaN;
                    }
                } else {
                    matrix[r][c] = Float.NaN;
                }
            }
        }
    }
}
