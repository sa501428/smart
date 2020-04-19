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

package mixer.commandline.utils.drink;

import com.google.common.util.concurrent.AtomicDouble;
import mixer.MixerGlobals;
import mixer.commandline.utils.drink.kmeansfloat.Cluster;
import mixer.commandline.utils.drink.kmeansfloat.ClusterTools;
import mixer.commandline.utils.drink.kmeansfloat.ConcurrentKMeans;
import mixer.commandline.utils.drink.kmeansfloat.KMeansListener;
import mixer.data.ChromosomeHandler;
import mixer.data.feature.GenomeWideList;
import org.broad.igv.util.Pair;

import java.util.Arrays;
import java.util.concurrent.TimeUnit;
import java.util.concurrent.atomic.AtomicInteger;

public class GenomeWideKmeansRunner {

    private static Cluster[] recentClusters;
    private static int[] recentIDs;

    private final CompositeGenomeWideDensityMatrix matrix;
    private final ChromosomeHandler chromosomeHandler;
    private final AtomicInteger numActualClusters = new AtomicInteger(0);
    private final AtomicDouble withinClusterSumOfSquaresForRun = new AtomicDouble(0);
    private final int maxIters = 20000;


    private GenomeWideList<SubcompartmentInterval> finalCompartments;
    private int numClusters = 0;

    public GenomeWideKmeansRunner(ChromosomeHandler chromosomeHandler, CompositeGenomeWideDensityMatrix interMatrix) {
        matrix = interMatrix;
        this.chromosomeHandler = chromosomeHandler;
    }

    public void prepareForNewRun(int numClusters) {
        recentClusters = null;
        recentIDs = null;
        this.numClusters = numClusters;
        numActualClusters.set(0);
        withinClusterSumOfSquaresForRun.set(0);
        finalCompartments = new GenomeWideList<>(chromosomeHandler);
    }

    public void launchKmeansGWMatrix(long seed) {

        if (matrix.getLength() > 0 && matrix.getWidth() > 0) {

            if (MixerGlobals.printVerboseComments) {
                System.out.println("Using seed " + seed);
            }

            ConcurrentKMeans kMeans = new ConcurrentKMeans(matrix.getCleanedData(), numClusters, maxIters, seed);

            KMeansListener kMeansListener = new KMeansListener() {
                @Override
                public void kmeansMessage(String s) {
                    if (MixerGlobals.printVerboseComments) {
                        System.out.println(s);
                    }
                }

                @Override
                public void kmeansComplete(Cluster[] preSortedClusters, long l) {

                    Cluster[] clusters = ClusterTools.getSortedClusters(preSortedClusters);

                    System.out.print(".");
                    Pair<Double, int[]> wcssAndIds = matrix.processGWKmeansResult(clusters, finalCompartments);
                    recentClusters = ClusterTools.clone(clusters);
                    recentIDs = wcssAndIds.getSecond();
                    numActualClusters.set(clusters.length);
                    withinClusterSumOfSquaresForRun.set(wcssAndIds.getFirst());
                }

                @Override
                public void kmeansError(Throwable throwable) {
                    throwable.printStackTrace();
                    System.err.println("gw full drink - err - " + throwable.getLocalizedMessage());
                    System.exit(98);
                }
            };
            kMeans.addKMeansListener(kMeansListener);
            kMeans.run();
        }

        waitUntilDone();
    }

    private void waitUntilDone() {
        while (numActualClusters.get() < 1 && withinClusterSumOfSquaresForRun.get() == 0.0) {
            System.out.print(".");
            try {
                TimeUnit.SECONDS.sleep(10);
            } catch (InterruptedException e) {
                e.printStackTrace();
            }
        }
    }

    public int getNumActualClusters() {
        return numActualClusters.get();
    }

    public double getWithinClusterSumOfSquares() {
        return withinClusterSumOfSquaresForRun.get();
    }

    public Cluster[] getRecentClustersClone() {
        return ClusterTools.clone(recentClusters);
    }

    public int[] getRecentIDsClone() {
        return Arrays.copyOf(recentIDs, recentIDs.length);
    }

    public GenomeWideList<SubcompartmentInterval> getFinalCompartments() {
        return finalCompartments.deepClone();
    }
}
