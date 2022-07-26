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

package mixer.utils.custom;

import mixer.MixerGlobals;
import robust.concurrent.kmeans.clustering.Cluster;
import robust.concurrent.kmeans.clustering.KMeansListener;
import robust.concurrent.kmeans.clustering.RobustConcurrentKMeans;

import java.util.Arrays;
import java.util.Random;
import java.util.concurrent.TimeUnit;
import java.util.concurrent.atomic.AtomicInteger;

public class QuickClustering {

    private final int maxIters;
    private final float[][] matrix;
    private final int numClusters;
    private final Random generator = new Random(128736);
    private final AtomicInteger numActualClusters = new AtomicInteger(0);
    private int[] assignments = null;

    public QuickClustering(float[][] matrix, int numCentroids, long seed, int maxIters) {
        this.matrix = matrix;
        this.numClusters = numCentroids;
        this.maxIters = maxIters;
        generator.setSeed(seed);
        if (matrix.length == 0 || matrix[0].length == 0) {
            System.err.println("Empty matrix provided for quick centroids");
            System.exit(5);
        }
    }

    public int[] cluster() {
        RobustConcurrentKMeans kMeans = new RobustConcurrentKMeans(matrix, numClusters, maxIters, generator.nextLong());

        KMeansListener kMeansListener = new KMeansListener() {
            @Override
            public void kmeansMessage(String s) {
                if (MixerGlobals.printVerboseComments) {
                    System.out.println(s);
                }
            }

            @Override
            public void kmeansComplete(Cluster[] clusters) {
                setAssignments(clusters);
                System.out.print(".");
                numActualClusters.set(clusters.length);
            }

            @Override
            public void kmeansError(Throwable throwable) {
                throwable.printStackTrace();
                System.err.println("Error - " + throwable.getLocalizedMessage());
                System.exit(18);
            }
        };
        kMeans.addKMeansListener(kMeansListener);
        kMeans.run();

        waitUntilDone();
        return assignments;
    }

    private void waitUntilDone() {
        while (numActualClusters.get() < 1) {
            System.out.print(".");
            try {
                TimeUnit.SECONDS.sleep(10);
            } catch (InterruptedException e) {
                e.printStackTrace();
            }
        }
    }

    private void setAssignments(Cluster[] clusters) {
        assignments = new int[matrix[0].length];
        Arrays.fill(assignments, -1);
        for (int i = 0; i < clusters.length; i++) {
            for (int k : clusters[i].getMemberIndexes()) {
                assignments[k] = i;
            }
        }
    }
}
