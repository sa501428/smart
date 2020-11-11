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

import tagbio.umap.metric.Metric;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

public class QuickCentroids {
    public static float[][] generateCentroids(float[][] matrix, int numCentroids, int numIters, Metric metric) {

        List<Integer> bestIndices = getInitialIndices(matrix, numCentroids, metric);

        int vectorLength = matrix[0].length;
        float[][] centroids = new float[bestIndices.size()][vectorLength];
        for (int i = 0; i < numCentroids; i++) {
            System.arraycopy(matrix[bestIndices.get(i)], 0, centroids[i], 0, vectorLength);
        }

        for (int z = 0; z < numIters; z++) {
            centroids = getUpdatedCentroids(centroids, matrix, metric);
        }

        return centroids;
    }

    private static float[][] getUpdatedCentroids(float[][] prevCentroids, float[][] matrix, Metric metric) {

        int[] closestCentroid = getClosestCentroidsForEachVector(prevCentroids, matrix, metric);

        int n = prevCentroids.length;
        float[][] newCentroids = new float[n][matrix[0].length];
        int[][] countsForCentroid = new int[n][matrix[0].length];
        double[][] centroidTotals = new double[n][matrix[0].length];

        for (int i = 0; i < matrix.length; i++) {
            int cID = closestCentroid[i];
            if (cID > -1) {
                for (int j = 0; j < matrix[i].length; j++) {
                    if (!Float.isNaN(matrix[i][j])) {
                        centroidTotals[cID][j] += matrix[i][j];
                        countsForCentroid[cID][j] += 1;
                    }
                }
            }
        }

        for (int i = 0; i < newCentroids.length; i++) {
            for (int j = 0; j < newCentroids[i].length; j++) {
                if (countsForCentroid[i][j] > 0) {
                    newCentroids[i][j] = (float) (centroidTotals[i][j] / countsForCentroid[i][j]);
                } else {
                    newCentroids[i][j] = Float.NaN;
                }
            }
        }
        return newCentroids;
    }

    private static int[] getClosestCentroidsForEachVector(float[][] prevCentroids, float[][] matrix, Metric metric) {
        int[] closestCentroid = new int[matrix.length];
        for (int i = 0; i < matrix.length; i++) {
            closestCentroid[i] = getNearestCentroidIndex(matrix[i], prevCentroids, metric);
        }
        return closestCentroid;
    }

    private static int getNearestCentroidIndex(float[] vector, float[][] prevCentroids, Metric metric) {
        int bestIndexSoFar = -1;
        double currDist = Double.MAX_VALUE;

        for (int k = 0; k < prevCentroids.length; k++) {
            double newDist = metric.distance(prevCentroids[k], vector);
            if (newDist < currDist) {
                currDist = newDist;
                bestIndexSoFar = k;
            }
        }
        return bestIndexSoFar;
    }

    private static List<Integer> getInitialIndices(float[][] matrix, int numCentroids, Metric metric) {
        float[] distFromClosestPoint = new float[matrix.length];
        Arrays.fill(distFromClosestPoint, Float.MAX_VALUE);

        List<Integer> bestIndices = new ArrayList<>();
        bestIndices.add(getMostDenseVectorIndex(matrix));

        for (int c = 0; c < numCentroids - 1; c++) {
            updateDistances(distFromClosestPoint, matrix, bestIndices.get(c), metric);
            bestIndices.add(getIndexOfMaxVal(distFromClosestPoint));
        }

        return bestIndices;
    }

    private static void updateDistances(float[] distFromClosestPoint, float[][] matrix, Integer index, Metric metric) {
        for (int k = 0; k < matrix.length; k++) {
            float newDist = metric.distance(matrix[k], matrix[index]);
            distFromClosestPoint[k] = Math.min(distFromClosestPoint[k], newDist);
        }
    }

    private static int getIndexOfMaxVal(float[] data) {
        float max = data[0];
        int index = 0;

        for (int i = 0; i < data.length; i++) {
            if (data[i] > max) {
                index = i;
                max = data[i];
            }
        }
        return index;
    }

    public static int getMostDenseVectorIndex(float[][] matrix) {
        float[] numNonNans = new float[matrix.length];
        for (int i = 0; i < matrix.length; i++) {
            numNonNans[i] = getNumNonNanEntries(matrix[i]);
        }

        return getIndexOfMaxVal(numNonNans);
    }

    public static int getNumNonNanEntries(float[] row) {
        int numNanEntries = 0;
        for (float val : row) {
            if (Float.isNaN(val)) {
                numNanEntries++;
            }
        }
        return row.length - numNanEntries;
    }
}
