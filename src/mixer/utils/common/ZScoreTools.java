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

import java.util.ArrayList;
import java.util.List;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;
import java.util.concurrent.atomic.AtomicInteger;

@SuppressWarnings("ForLoopReplaceableByForEach")
public class ZScoreTools {

    private static final float ZERO = 1e-10f;
    private static final int numCPUThreads = 30;

    public static void inPlaceZscoreDownCol(float[][] matrix, int[] weights) {
        float[] colMeans = getColMean(matrix);
        float[] colStdDevs = getColStdDev(matrix, colMeans);

        for (int i = 0; i < matrix.length; i++) {
            for (int j = 0; j < matrix[i].length; j++) {
                float val = matrix[i][j];
                if (!Float.isNaN(val)) {
                    matrix[i][j] = weights[j] * (val - colMeans[j]) / colStdDevs[j];
                }
            }
        }
    }

    public static float[] getColMean(float[][] matrix) {
        double[] colSums = new double[matrix[0].length];
        int[] colSize = new int[colSums.length];

        for (int i = 0; i < matrix.length; i++) {
            for (int j = 0; j < matrix[i].length; j++) {
                float val = matrix[i][j];
                if (isValid(val)) {
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

    public static float[] getColStdDev(float[][] matrix, float[] means) {

        double[] squares = new double[means.length];
        int[] colSize = new int[means.length];

        for (int i = 0; i < matrix.length; i++) {
            for (int j = 0; j < matrix[i].length; j++) {
                float val = matrix[i][j];
                if (isValid(val)) {
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


    /**
     * Robust Zscore Methods
     */
    public static void inPlaceRobustZscoreDownCol(float[][] matrix) {
        float[] colMedians = getParColNonZeroMedian(matrix);
        float[] colMADs = getParColNonZeroMedianAbsoluteDeviations(matrix, colMedians);

        for (int i = 0; i < matrix.length; i++) {
            for (int j = 0; j < matrix[i].length; j++) {
                float val = matrix[i][j];
                if (!Float.isNaN(val)) {
                    //matrix[i][j] = (val - colMedians[j]) / colMADs[j];
                    matrix[i][j] = val / colMADs[j]; //weights[j] *
                }
            }
        }
    }

    public static float[] getParColNonZeroMedian(float[][] matrix) {
        float[] colMedians = new float[matrix[0].length];
        AtomicInteger currColIndex = new AtomicInteger(0);
        ExecutorService executor = Executors.newFixedThreadPool(numCPUThreads);
        for (int l = 0; l < numCPUThreads; l++) {
            Runnable worker = () -> {
                int j = currColIndex.getAndIncrement();
                while (j < colMedians.length) {
                    colMedians[j] = getMedian(matrix, j);
                    j = currColIndex.getAndIncrement();
                }
            };
            executor.execute(worker);
        }
        executor.shutdown();
        while (!executor.isTerminated()) {
        }
        return colMedians;
    }

    private static float[] getParColNonZeroMedianAbsoluteDeviations(float[][] matrix, float[] medians) {
        float[] colMADs = new float[medians.length];
        AtomicInteger currColIndex = new AtomicInteger(0);
        ExecutorService executor = Executors.newFixedThreadPool(numCPUThreads);
        for (int l = 0; l < numCPUThreads; l++) {
            Runnable worker = () -> {
                int j = currColIndex.getAndIncrement();
                while (j < colMADs.length) {
                    colMADs[j] = getMedianAbsoluteDeviation(matrix, j, medians[j]);
                    j = currColIndex.getAndIncrement();
                }
            };
            executor.execute(worker);
        }
        executor.shutdown();
        while (!executor.isTerminated()) {
        }
        return colMADs;
    }

    private static float getMedian(float[][] matrix, int col) {
        List<Float> vals = new ArrayList<>();
        for (int i = 0; i < matrix.length; i++) {
            float val = matrix[i][col];
            if (isValid(val)) {
                vals.add(val);
            }
        }
        return QuickMedian.fastMedian(vals);
    }

    private static float getMedianAbsoluteDeviation(float[][] matrix, int col, float median) {
        List<Float> vals = new ArrayList<>();
        for (int i = 0; i < matrix.length; i++) {
            float val = matrix[i][col];
            if (isValid(val)) {
                vals.add(Math.abs(val - median));
            }
        }
        float mad = QuickMedian.fastMedian(vals);
        if (mad > ZERO) {
            return mad;
        }
        return 1;
    }

    private static boolean isValid(float val) {
        return !Float.isNaN(val); // && val > ZERO
    }
}
