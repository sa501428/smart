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

package mixer.utils.impute;

import javastraw.tools.ParallelizationTools;
import mixer.SmartTools;
import robust.concurrent.kmeans.clustering.Cluster;
import robust.concurrent.kmeans.clustering.KMeansListener;
import robust.concurrent.kmeans.clustering.RobustConcurrentKMeans;
import robust.concurrent.kmeans.clustering.RobustConcurrentKMedians;

import java.util.Random;
import java.util.concurrent.TimeUnit;
import java.util.concurrent.atomic.AtomicInteger;

public class CentroidImputer {

    public static void updateBasedOnCentroids(float[][] imputed, int numClusters, Random generator) {
        AtomicInteger numActualClusters = new AtomicInteger(-1);
        RobustConcurrentKMeans kMeans = new RobustConcurrentKMedians(imputed, numClusters, 10,
                generator.nextLong(), SmartTools.NUM_ENTRIES_TO_SKIP_MEDIAN);

        KMeansListener kMeansListener = new KMeansListener() {
            @Override
            public void kmeansMessage(String s) {
                if (SmartTools.printVerboseComments) {
                    System.out.println(s);
                }
            }

            @Override
            public void kmeansComplete(Cluster[] clusters) {
                updateImputedMatrixEntries(clusters, imputed, numActualClusters);
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
        waitUntilDone(numActualClusters);
    }

    private static void updateImputedMatrixEntries(Cluster[] clusters, float[][] imputed,
                                                   AtomicInteger numActualClusters) {
        AtomicInteger clusterIndex = new AtomicInteger(0);
        ParallelizationTools.launchParallelizedCode(() -> {
            int c = clusterIndex.getAndIncrement();
            while (c < clusters.length) {
                processCluster(clusters[c], imputed);
                c = clusterIndex.getAndIncrement();
            }
        });
        numActualClusters.set(clusters.length);
    }

    private static void processCluster(Cluster cluster, float[][] imputed) {
        float[] center = cluster.getCenter();
        for (int r : cluster.getMemberIndexes()) {
            for (int j = 0; j < imputed[r].length; j++) {
                if (Float.isNaN(imputed[r][j]) && center[j] > 0) {
                    imputed[r][j] = center[j];
                }
            }
        }
    }


    private static void waitUntilDone(AtomicInteger numActualClusters) {
        while (numActualClusters.get() < 1) {
            System.out.print(".");
            try {
                TimeUnit.SECONDS.sleep(10);
            } catch (InterruptedException e) {
                e.printStackTrace();
            }
        }
    }
}
