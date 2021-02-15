/*
 * The MIT License (MIT)
 *
 * Copyright (c) 2011-2021 Rice University, Baylor College of Medicine, Aiden Lab
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

package mixer.utils.common;

import java.util.Arrays;

@SuppressWarnings("ForLoopReplaceableByForEach")
public class ZScoreTools {

    private static final int UPPER_ZSCORE = 3;
    private static final float ZERO = 1e-10f;

    public static void inPlaceTwoStageZscoreDownCol(float[][] matrix) {
        float[] initialCutoff = new float[matrix[0].length];
        Arrays.fill(initialCutoff, Float.MAX_VALUE);

        float[] colMeans = getColNonZeroMean(matrix, initialCutoff);
        float[] colStdDevs = getColNonZeroStdDev(matrix, colMeans, initialCutoff);
        float[] zscoreThresholds = getZscoreUpperThresholds(colMeans, colStdDevs);

        postFilterZscoreDownCol(matrix, zscoreThresholds);
    }

    public static void postFilterZscoreDownCol(float[][] matrix, float[] thresholds) {
        float[] colMeans = getColNonZeroMean(matrix, thresholds);
        float[] colStdDevs = getColNonZeroStdDev(matrix, colMeans, thresholds);

        for (int i = 0; i < matrix.length; i++) {
            for (int j = 0; j < matrix[i].length; j++) {
                float val = matrix[i][j];
                if (!Float.isNaN(val)) {
                    matrix[i][j] = (val - colMeans[j]) / colStdDevs[j];
                }
            }
        }
    }

    private static float[] getZscoreUpperThresholds(float[] means, float[] std) {
        float[] cutoffs = new float[means.length];
        for (int k = 0; k < cutoffs.length; k++) {
            cutoffs[k] = means[k] + (UPPER_ZSCORE * std[k]);
        }
        return cutoffs;
    }

    public static float[] getColNonZeroMean(float[][] matrix, float[] cutoffs) {
        double[] colSums = new double[matrix[0].length];
        int[] colSize = new int[colSums.length];

        for (int i = 0; i < matrix.length; i++) {
            for (int j = 0; j < matrix[i].length; j++) {
                float val = matrix[i][j];
                if (isValid(val, cutoffs[j])) {
                    colSums[j] += val;
                    colSize[j] += 1;
                }
            }
        }

        float[] colMeans = new float[colSums.length];
        for (int k = 0; k < colSums.length; k++) {
            colMeans[k] = (float) (colSums[k] / Math.max(colSize[k], 1));
        }

        return colMeans;
    }

    public static float[] getColNonZeroStdDev(float[][] matrix, float[] means, float[] cutoffs) {

        double[] squares = new double[means.length];
        int[] colSize = new int[means.length];

        for (int i = 0; i < matrix.length; i++) {
            for (int j = 0; j < matrix[i].length; j++) {
                float val = matrix[i][j];
                if (isValid(val, cutoffs[j])) {
                    float diff = val - means[j];
                    squares[j] += diff * diff;
                    colSize[j] += 1;
                }
            }
        }

        float[] stdDev = new float[means.length];
        for (int k = 0; k < squares.length; k++) {
            stdDev[k] = (float) Math.sqrt(squares[k] / Math.max(colSize[k], 1));
        }

        return stdDev;
    }

    private static boolean isValid(float val, float threshold) {
        return !Float.isNaN(val) && val > ZERO && val < threshold;
    }
}
