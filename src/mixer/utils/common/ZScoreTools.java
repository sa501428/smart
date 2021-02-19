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

    public static void inPlaceRobustZscoreDownCol(float[][] matrix, int[] weights) {
        float[] colMedians = getParColNonZeroMedian(matrix);
        float[] colMADs = getParColNonZeroMedianAbsoluteDeviations(matrix, colMedians);

        for (int i = 0; i < matrix.length; i++) {
            for (int j = 0; j < matrix[i].length; j++) {
                float val = matrix[i][j];
                if (!Float.isNaN(val)) {
                    //matrix[i][j] = (val - colMedians[j]) / colMADs[j];
                    matrix[i][j] = weights[j] * val / colMADs[j];
                }
            }
        }
    }

    public static float[] getColNonZeroMedian(float[][] matrix) {
        float[] colMedians = new float[matrix[0].length];
        for (int j = 0; j < matrix[0].length; j++) {
            colMedians[j] = getMedian(matrix, j);
        }
        return colMedians;
    }

    private static float[] getColNonZeroMedianAbsoluteDeviations(float[][] matrix, float[] medians) {
        float[] colMADs = new float[medians.length];
        for (int j = 0; j < medians.length; j++) {
            colMADs[j] = getMedianAbsoluteDeviation(matrix, j, medians[j]);
        }
        return colMADs;
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
        return sortedMidpoint(vals);
    }

    private static float getMedianAbsoluteDeviation(float[][] matrix, int col, float median) {
        List<Float> vals = new ArrayList<>();
        for (int i = 0; i < matrix.length; i++) {
            float val = matrix[i][col];
            if (isValid(val)) {
                vals.add(Math.abs(val - median));
            }
        }
        return sortedMidpoint(vals);
    }

    private static float sortedMidpoint(List<Float> vals) {
        int size = vals.size();
        vals.sort(null);
        return (vals.get(size / 2) + vals.get((size - 1) / 2)) / 2;
    }

    private static boolean isValid(float val) {
        return !Float.isNaN(val) && val > ZERO;
    }
}
