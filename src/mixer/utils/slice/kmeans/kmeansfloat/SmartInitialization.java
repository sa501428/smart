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

package mixer.utils.slice.kmeans.kmeansfloat;

import mixer.utils.similaritymeasures.RobustEuclideanDistance;

import java.util.Arrays;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;
import java.util.concurrent.atomic.AtomicInteger;

public class SmartInitialization {

    private final float[][] data;
    private final int numClusters;
    private final int[] bestIndices;
    private final float[] distFromClosestPoint;

    public SmartInitialization(float[][] data, int numClusters, int initialID) {
        this.data = data;
        this.numClusters = numClusters;
        bestIndices = new int[numClusters];
        bestIndices[0] = initialID;
        distFromClosestPoint = new float[data.length];
        Arrays.fill(distFromClosestPoint, Float.MAX_VALUE);
    }

    public int[] getSmartKmeansInitialization() {

        for (int c = 0; c < numClusters - 1; c++) {
            updateDistances(bestIndices[c]);
            bestIndices[c + 1] = getIndexOfMaxVal();
        }

        return bestIndices;
    }

    private void updateDistances(Integer index) {

        int numCPUThreads = Runtime.getRuntime().availableProcessors();
        AtomicInteger currRowIndex = new AtomicInteger(0);
        ExecutorService executor = Executors.newFixedThreadPool(numCPUThreads);
        for (int l = 0; l < numCPUThreads; l++) {
            executor.execute(() -> {
                int k = currRowIndex.getAndIncrement();
                while (k < data.length) {
                    float newDist = 0;
                    if (k != index) {
                        newDist = RobustEuclideanDistance.SINGLETON.distance(data[k], data[index]);
                    }
                    distFromClosestPoint[k] = Math.min(distFromClosestPoint[k], newDist);
                    k = currRowIndex.getAndIncrement();
                }
            });
        }
        executor.shutdown();
        //noinspection StatementWithEmptyBody
        while (!executor.isTerminated()) {
        }
    }

    private int getIndexOfMaxVal() {
        float max = distFromClosestPoint[0];
        int index = 0;

        for (int i = 0; i < distFromClosestPoint.length; i++) {
            if (distFromClosestPoint[i] > max) {
                index = i;
                max = distFromClosestPoint[i];
            }
        }
        return index;
    }
}
