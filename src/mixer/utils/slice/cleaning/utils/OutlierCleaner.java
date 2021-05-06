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

package mixer.utils.slice.cleaning.utils;

import javastraw.tools.ParallelizedJuicerTools;
import mixer.utils.common.ArrayTools;
import mixer.utils.similaritymeasures.RobustEuclideanDistance;
import mixer.utils.similaritymeasures.RobustManhattanDistance;
import mixer.utils.similaritymeasures.SimilarityMetric;
import mixer.utils.slice.cleaning.QuickCentroids;

import java.util.Arrays;
import java.util.HashSet;
import java.util.Set;
import java.util.concurrent.atomic.AtomicInteger;

public class OutlierCleaner {
    private static final float DIST_CUTOFF = 5;
    protected final int NUM_POINTS_PER_CENTROID = 30;
    private final float[][] matrix;
    boolean useOnlyCorr;

    public OutlierCleaner(float[][] matrix, boolean useOnlyCorr) {
        this.matrix = matrix;
        this.useOnlyCorr = useOnlyCorr;
    }

    public Set<Integer> getConsistentOutliers() {
        Set<Integer> outlierIndices = getCorrOutlierIndices(matrix, useOnlyCorr);
        if (!useOnlyCorr) {
            int numInitialClusters = matrix.length / NUM_POINTS_PER_CENTROID;
            float[][] centroids = new QuickCentroids(matrix, numInitialClusters, 5L).generateCentroids(5);
            Set<Integer> distOutlierIndices = new HashSet<>();
            for (SimilarityMetric metric : new SimilarityMetric[]{RobustEuclideanDistance.SINGLETON,
                    RobustManhattanDistance.SINGLETON}) {
                distOutlierIndices.addAll(getOutlierIndices(matrix, centroids, metric));
            }
            outlierIndices.retainAll(distOutlierIndices);
        }
        System.out.println("Number of dist outliers (corr:" + useOnlyCorr + ") " + outlierIndices.size());

        return outlierIndices;
    }

    private Set<Integer> getCorrOutlierIndices(float[][] matrix, boolean onlyUseCorr) {
        Correlator correlator = new Correlator(matrix, onlyUseCorr);
        return correlator.getOutlierIndices();
    }

    private Set<Integer> getOutlierIndices(float[][] matrix, float[][] centroids, SimilarityMetric metric) {
        float[] minDist = new float[matrix.length];
        Arrays.fill(minDist, Float.MAX_VALUE);
        AtomicInteger rowIndex = new AtomicInteger(0);
        ParallelizedJuicerTools.launchParallelizedCode(() -> {
            int i = rowIndex.getAndIncrement();
            while (i < matrix.length) {
                for (int j = 0; j < centroids.length; j++) {
                    float dist = metric.distance(matrix[i], centroids[j]);
                    if (dist < minDist[i]) {
                        minDist[i] = dist;
                    }
                }
                i = rowIndex.getAndIncrement();
            }
        });

        return indexOfOutliers(minDist);
    }

    private Set<Integer> indexOfOutliers(float[] minDist) {
        Set<Integer> outliers = new HashSet<>();
        float mu = ArrayTools.getNonZeroMean(minDist);
        float std = ArrayTools.getNonZeroStd(minDist, mu);
        for (int k = 0; k < minDist.length; k++) {
            float zscore = (minDist[k] - mu) / std;
            if (zscore > DIST_CUTOFF) {
                outliers.add(k);
            }
        }
        return outliers;
    }
}
