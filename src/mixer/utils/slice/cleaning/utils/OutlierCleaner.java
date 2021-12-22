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

import javastraw.tools.MatrixTools;
import javastraw.tools.ParallelizationTools;
import mixer.MixerGlobals;
import mixer.utils.common.ArrayTools;
import mixer.utils.similaritymeasures.RobustEuclideanDistance;
import mixer.utils.similaritymeasures.RobustManhattanDistance;
import mixer.utils.similaritymeasures.SimilarityMetric;
import mixer.utils.slice.cleaning.QuickCentroids;

import java.io.File;
import java.util.Arrays;
import java.util.HashSet;
import java.util.Set;
import java.util.concurrent.atomic.AtomicInteger;

public class OutlierCleaner {
    private static final float DIST_CUTOFF = 3;
    protected final int NUM_POINTS_PER_CENTROID = 30;
    private final float[][] matrix;
    boolean useOnlyCorr;
    private final int ONE_MB = 1000000;
    int numInitialClusters = 100;

    public OutlierCleaner(float[][] matrix, boolean useOnlyCorr) {
        this.matrix = matrix;
        this.useOnlyCorr = useOnlyCorr;
    }

    public Set<Integer> getConsistentOutliers(int resolution, File outputDirectory) {
        Set<Integer> outlierIndices = getCorrOutlierIndices(matrix, useOnlyCorr);
        if (!useOnlyCorr) {
            float[][] centroids = new QuickCentroids(matrix, numInitialClusters,
                    5L).generateCentroids(ONE_MB / resolution, true);
            for (SimilarityMetric metric : new SimilarityMetric[]{RobustEuclideanDistance.SINGLETON,
                    RobustManhattanDistance.SINGLETON}) {
                outlierIndices.addAll(getOutlierIndices(matrix, centroids, metric, outputDirectory));
            }
        }
        if (MixerGlobals.printVerboseComments) {
            System.out.println("Number of dist outliers (use_only_corr:" + useOnlyCorr + ") " + outlierIndices.size());
        }
        return outlierIndices;
    }

    private Set<Integer> getCorrOutlierIndices(float[][] matrix, boolean onlyUseCorr) {
        Correlator correlator = new Correlator(matrix, onlyUseCorr);
        return correlator.getOutlierIndices();
    }

    private Set<Integer> getOutlierIndices(float[][] matrix, float[][] centroids,
                                           SimilarityMetric metric, File outputDirectory) {
        float[] minDist = new float[matrix.length];
        Arrays.fill(minDist, Float.MAX_VALUE);
        AtomicInteger rowIndex = new AtomicInteger(0);
        ParallelizationTools.launchParallelizedCode(() -> {
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

        if (MixerGlobals.printVerboseComments) {
            float[][] mtrx = new float[1][minDist.length];
            mtrx[0] = minDist;
            File outfile = new File(outputDirectory, metric.toString() + "_distances.npy");
            MatrixTools.saveMatrixTextNumpy(outfile.getAbsolutePath(), mtrx);
        }

        return indexOfOutliers(minDist);
    }

    private Set<Integer> indexOfOutliers(float[] minDist) {
        Set<Integer> outliers = new HashSet<>();
        float mu = ArrayTools.getNonZeroMean(minDist);
        float std = ArrayTools.getNonZeroStd(minDist, mu);
        for (int k = 0; k < minDist.length; k++) {
            float zScore = (minDist[k] - mu) / std;
            if (zScore > DIST_CUTOFF) {
                outliers.add(k);
            }
        }
        return outliers;
    }
}
