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

package mixer.utils.rougheval;

import mixer.clt.ParallelizedMixerTools;
import mixer.utils.common.ArrayTools;
import mixer.utils.common.FloatMatrixTools;
import mixer.utils.similaritymeasures.SimilarityMetric;

import java.util.Arrays;
import java.util.List;
import java.util.concurrent.atomic.AtomicInteger;

public class Silhouette {

    private final float[][] matrix;
    private final double[] meanIntraDistanceA, minMeanInterDistOtherClustersB, silhouette;
    private final List<List<Integer>> clusters;
    private final SimilarityMetric metric;
    private Double score = null;

    public Silhouette(float[][] matrix, List<List<Integer>> clusters, SimilarityMetric metric) {
        this.matrix = matrix;
        this.clusters = clusters;
        this.metric = metric;
        meanIntraDistanceA = new double[matrix.length];
        minMeanInterDistOtherClustersB = new double[matrix.length];
        silhouette = new double[matrix.length];
        Arrays.fill(meanIntraDistanceA, Double.NaN);
        Arrays.fill(minMeanInterDistOtherClustersB, Double.NaN);
        Arrays.fill(silhouette, Double.NaN);
    }

    public double getScore() {
        if (score == null) {
            calculateDistances();
            score = getAverageScore();
        }
        return score;
    }

    private double getAverageScore() {
        for (int z = 0; z < silhouette.length; z++) {
            silhouette[z] = calculateSilhouette(meanIntraDistanceA[z],
                    minMeanInterDistOtherClustersB[z]);
        }
        return ArrayTools.mean(silhouette);
    }

    private double calculateSilhouette(double a, double b) {
        if (Double.isNaN(a) || Double.isNaN(b)) {
            return 0;
        }
        if (a < b) {
            return 1 - (a / b);
        } else if (b < a) {
            return (b / a) - 1;
        }
        return 0;
    }

    private void calculateDistances() {
        for (int i = 0; i < clusters.size(); i++) {
            for (int j = i; j < clusters.size(); j++) {
                updateDistancesBetweenClusters(i, j);
            }
        }
    }

    private void updateDistancesBetweenClusters(int i, int j) {
        if (i == j) {
            getIntraDistancesWithinCluster(i);
        } else {
            getInterDistancesBetweenClusters(i, j);
        }
    }

    private void getInterDistancesBetweenClusters(int indexI, int indexJ) {
        List<Integer> indicesI = clusters.get(indexI);
        List<Integer> indicesJ = clusters.get(indexJ);
        int nI = indicesI.size();
        int nJ = indicesJ.size();
        final float[][] result = new float[nI][nJ];

        AtomicInteger currIndex = new AtomicInteger(0);
        ParallelizedMixerTools.launchParallelizedCode(() -> {
            int i = currIndex.getAndIncrement();
            while (i < nI) {
                for (int j = 0; j < nJ; j++) {
                    result[i][j] = getDistanceBetweenVectors(indicesI.get(i), indicesJ.get(j));
                }
                i = currIndex.getAndIncrement();
            }
        });

        populateInterClusterDistances(nI, indicesI, result);
        populateInterClusterDistances(nJ, indicesJ, FloatMatrixTools.transpose(result));
    }

    private void populateInterClusterDistances(int n, List<Integer> indices, float[][] result) {
        for (int i = 0; i < n; i++) {
            int index = indices.get(i);
            double val = ArrayTools.nanMean(result[i]);
            double originalVal = minMeanInterDistOtherClustersB[index];
            if (Double.isNaN(originalVal)) {
                minMeanInterDistOtherClustersB[index] = val;
            } else {
                minMeanInterDistOtherClustersB[index] = Math.min(originalVal, val);
            }
        }
    }

    private void getIntraDistancesWithinCluster(int index) {
        List<Integer> indices = clusters.get(index);
        int n = indices.size();
        final float[][] result = new float[n][n];

        AtomicInteger currIndex = new AtomicInteger(0);
        ParallelizedMixerTools.launchParallelizedCode(() -> {
            int i = currIndex.getAndIncrement();
            while (i < n) {
                result[i][i] = Float.NaN;
                for (int j = i + 1; j < n; j++) {
                    result[i][j] = getDistanceBetweenVectors(indices.get(i), indices.get(j));
                    result[j][i] = result[i][j];
                }
                i = currIndex.getAndIncrement();
            }
        });

        for (int i = 0; i < n; i++) {
            meanIntraDistanceA[indices.get(i)] = ArrayTools.nanMean(result[i]);
        }
    }

    private float getDistanceBetweenVectors(int i, int j) {
        return metric.distance(matrix[i], matrix[j]);
    }
}
