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

package mixer.utils.refinement;

import mixer.utils.similaritymeasures.RobustEuclideanDistance;
import mixer.utils.similaritymeasures.RobustManhattanDistance;
import mixer.utils.similaritymeasures.SimilarityMetric;
import org.apache.commons.math3.stat.descriptive.DescriptiveStatistics;

import java.util.HashSet;
import java.util.Set;

public class NearestNeighbors {

    private static final int percentile = 2;
    private final int index;
    private final double percentUnassigned;
    private final Set<Integer> neighbors;

    public NearestNeighbors(int index, float[][] matrixL1, float[][] matrixL2, int[] hubs) {
        this.index = index;
        neighbors = getNearestN(index, matrixL1, RobustManhattanDistance.SINGLETON);
        Set<Integer> neighbors2 = getNearestN(index, matrixL2, RobustEuclideanDistance.SINGLETON);
        neighbors.retainAll(neighbors2);
        percentUnassigned = setPercentUnassigned(hubs);
    }

    public NearestNeighbors(int index, float[][] matrix, int[] assignment, SimilarityMetric metric) {
        this.index = index;
        neighbors = getNearestN(index, matrix, metric);
        percentUnassigned = 0;
    }

    private Set<Integer> getNearestN(int index, float[][] matrix, SimilarityMetric metric) {
        double[] distances = new double[matrix.length];
        for (int r = 0; r < matrix.length; r++) {
            distances[r] = metric.distance(matrix[r], matrix[index]);
        }
        distances[index] = Double.MAX_VALUE;
        double cutoff = getNearbyIndicesFromPercentile(distances);
        return indicesLessThan(distances, cutoff);
    }

    private Set<Integer> indicesLessThan(double[] distances, double cutoff) {
        Set<Integer> indices = new HashSet<>();
        for (int i = 0; i < distances.length; i++) {
            if (distances[i] < cutoff) {
                indices.add(i);
            }
        }
        return indices;
    }

    private double getNearbyIndicesFromPercentile(double[] values) {
        DescriptiveStatistics statistics = new DescriptiveStatistics();
        for (double v : values) {
            statistics.addValue(v);
        }
        return statistics.getPercentile(percentile);
    }

    private double setPercentUnassigned(int[] hubs) {
        double notAssigned = 0;
        for (int i : neighbors) {
            if (hubs[i] < 0) {
                notAssigned++;
            }
        }
        return notAssigned / neighbors.size();
    }

    public double getPercentUnassigned() {
        return percentUnassigned;
    }

    public int getIndex() {
        return index;
    }

    public int getMajorityNeighborAssignment(int[] assignment, int numClusters) {
        int[] counts = getCounts(assignment, numClusters);

        int bestIndex = 0;
        for (int i = 1; i < numClusters; i++) {
            if (counts[i] > counts[bestIndex]) {
                bestIndex = i;
            }
        }

        return bestIndex;
    }

    private int[] getCounts(int[] assignment, int numClusters) {
        int[] counts = new int[numClusters];
        for (int i : neighbors) {
            if (assignment[i] > -1) {
                counts[assignment[i]]++;
            }
        }
        return counts;
    }

    public int getMajorityNeighborAssignment(int[] assignment, int numClusters, float fraction) {
        int[] counts = getCounts(assignment, numClusters);

        int cutoff = (int) (fraction * assignment.length);
        int bestIndex = -1;
        for (int i = 1; i < numClusters; i++) {
            if (counts[i] > cutoff) {
                bestIndex = i;
            }
        }
        return bestIndex;
    }
}
