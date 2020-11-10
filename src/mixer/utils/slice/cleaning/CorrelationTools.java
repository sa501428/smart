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

package mixer.utils.slice.cleaning;

import mixer.utils.shuffle.Metrics;

import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;
import java.util.concurrent.atomic.AtomicInteger;

public class CorrelationTools {

    /**
     * @param matrix
     * @param numCentroids
     * @param type         0 - cosine, 1 - corr, 2 - mse, 3 - mae, 4 - emd, 5 - KLD, 6 - JSD
     * @return
     */
    public static float[][] getMinimallySufficientNonNanSimilarityMatrix(float[][] matrix, int numCentroids, Metrics.Type type) {

        if (numCentroids == matrix.length) {
            return getNonNanDistanceMatrix(matrix, type);
        }

        float[][] centroids = QuickCentroids.generateCentroids(matrix, numCentroids, 5);
        float[][] result = new float[matrix.length][numCentroids]; // *2

        int numCPUThreads = 20;
        AtomicInteger currRowIndex = new AtomicInteger(0);
        ExecutorService executor = Executors.newFixedThreadPool(numCPUThreads);
        for (int l = 0; l < numCPUThreads; l++) {
            Runnable worker = new Runnable() {
                @Override
                public void run() {
                    int i = currRowIndex.getAndIncrement();
                    while (i < matrix.length) {
                        Metrics.fillEntries(i, 0, matrix, centroids, result, type, false);
                        i = currRowIndex.getAndIncrement();
                    }
                }
            };
            executor.execute(worker);
        }
        executor.shutdown();

        // Wait until all threads finish
        while (!executor.isTerminated()) {
        }

        return result;
        //return FloatMatrixTools.concatenate(matrix, result);
    }

    public static float[][] getNonNanDistanceMatrix(float[][] matrix, Metrics.Type type) {


        float[][] result = new float[matrix.length][matrix.length]; // *2

        int numCPUThreads = 20;
        AtomicInteger currRowIndex = new AtomicInteger(0);
        ExecutorService executor = Executors.newFixedThreadPool(numCPUThreads);
        for (int l = 0; l < numCPUThreads; l++) {
            Runnable worker = new Runnable() {
                @Override
                public void run() {
                    int i = currRowIndex.getAndIncrement();
                    while (i < matrix.length) {
                        Metrics.fillEntries(i, i, matrix, matrix, result, type, true);
                        i = currRowIndex.getAndIncrement();
                    }
                }
            };
            executor.execute(worker);
        }
        executor.shutdown();

        // Wait until all threads finish
        while (!executor.isTerminated()) {
        }

        return result;
    }
}
