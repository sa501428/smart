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

package mixer.utils.slice;

import com.google.common.util.concurrent.AtomicDouble;
import javastraw.featurelist.GenomeWideList;
import javastraw.reader.ChromosomeHandler;
import mixer.MixerGlobals;
import mixer.utils.common.Pair;
import mixer.utils.slice.kmeansfloat.Cluster;
import mixer.utils.slice.kmeansfloat.ClusterTools;
import mixer.utils.slice.kmeansfloat.ConcurrentKMeans;
import mixer.utils.slice.kmeansfloat.KMeansListener;
import mixer.utils.slice.matrices.CompositeGenomeWideDensityMatrix;
import mixer.utils.slice.structures.SubcompartmentInterval;

import java.util.List;
import java.util.concurrent.TimeUnit;
import java.util.concurrent.atomic.AtomicInteger;

public class GenomeWideKmeansRunner {

    private static Cluster[] recentClusters;
    private final CompositeGenomeWideDensityMatrix matrix;
    private final ChromosomeHandler chromosomeHandler;
    private final AtomicInteger numActualClusters = new AtomicInteger(0);
    private final AtomicDouble withinClusterSumOfSquaresForRun = new AtomicDouble(0);
    private final int maxIters = 1000; // changed from 20,000
    private int[][] recentIDs;
    private int[][] recentIDsForIndex;
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

            ConcurrentKMeans kMeans = new ConcurrentKMeans(matrix.getCleanedData(),
                    numClusters, maxIters, seed);

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
                    Pair<Double, List<int[][]>> wcssAndIds = matrix.processGWKmeansResult(clusters, finalCompartments);
                    recentClusters = ClusterTools.clone(clusters);
                    recentIDs = wcssAndIds.getSecond().get(0);
                    recentIDsForIndex = wcssAndIds.getSecond().get(1);
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
        int[] temp = new int[recentIDs[0].length];
        System.arraycopy(recentIDs[0], 0, temp, 0, temp.length);
        return temp;
    }

    public int[][] getRecentIDsForIndex() {
        int[][] temp = new int[recentIDsForIndex.length][recentIDsForIndex[0].length];
        for (int k = 0; k < temp.length; k++) {
            System.arraycopy(recentIDsForIndex[k], 0, temp[k], 0, temp[k].length);
        }

        return temp;
    }

    public GenomeWideList<SubcompartmentInterval> getFinalCompartments() {
        return finalCompartments.deepClone();
    }
}
