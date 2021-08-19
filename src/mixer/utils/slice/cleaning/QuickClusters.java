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

package mixer.utils.slice.cleaning;

import mixer.MixerGlobals;
import mixer.utils.slice.kmeans.kmeansfloat.Cluster;
import mixer.utils.slice.kmeans.kmeansfloat.ConcurrentKMeans;
import mixer.utils.slice.kmeans.kmeansfloat.KMeansListener;

import java.util.ArrayList;
import java.util.List;
import java.util.Random;
import java.util.concurrent.TimeUnit;
import java.util.concurrent.atomic.AtomicInteger;

public class QuickClusters {

    private final float[][] matrix;
    private final int initialNumClusters;
    private final Random generator = new Random(0);
    private final AtomicInteger numActualClusters = new AtomicInteger(0);
    private int maxIters = 20;
    private List<List<Integer>> indices = null;

    public QuickClusters(float[][] matrix, int numCentroids, long seed) {
        this.matrix = matrix;
        this.initialNumClusters = numCentroids;
        generator.setSeed(seed);
        if (matrix.length == 0 || matrix[0].length == 0) {
            System.err.println("Empty matrix provided for quick centroids");
            System.exit(5);
        }
    }

    public QuickClusters(float[][] matrix, int numCentroids, long seed, int numIters) {
        this(matrix, numCentroids, seed);
        this.maxIters = numIters;
    }

    public List<List<Integer>> getClusters() {
        ConcurrentKMeans kMeans = new ConcurrentKMeans(matrix, initialNumClusters, maxIters, generator.nextLong());

        KMeansListener kMeansListener = new KMeansListener() {
            @Override
            public void kmeansMessage(String s) {
                if (MixerGlobals.printVerboseComments) {
                    System.out.println(s);
                }
            }

            @Override
            public void kmeansComplete(Cluster[] clusters, long l) {
                convertClustersToIndicesList(clusters);
                System.out.print(".");
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
        return indices;
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

    private void convertClustersToIndicesList(Cluster[] initialClusters) {
        indices = new ArrayList<>();
        for (Cluster c : initialClusters) {
            List<Integer> members = new ArrayList<>();
            for (int member : c.getMemberIndexes()) {
                members.add(member);
            }
            indices.add(members);
        }

        numActualClusters.set(indices.size());
    }
}
