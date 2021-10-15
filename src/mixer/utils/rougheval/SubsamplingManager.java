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

import mixer.utils.similaritymeasures.RobustEuclideanDistance;
import mixer.utils.similaritymeasures.RobustManhattanDistance;
import mixer.utils.similaritymeasures.SimilarityMetric;
import robust.concurrent.kmeans.clustering.Cluster;

import java.util.ArrayList;
import java.util.List;
import java.util.Random;

public class SubsamplingManager {

    public static final int NUMBER_OF_TYPES = 4;
    private final Cluster[] clusters;
    private final float[][] matrix;
    private final boolean useKMedians;
    private final int type;
    private final double compressionFactor;
    private final MatrixSamplingContainer sample;
    private final SimilarityMetric metric;
    private final Random generator = new Random(0);
    private final double IDEAL_NUM_ROWS = 1000.0;


    public SubsamplingManager(Cluster[] clusters, float[][] matrix, boolean useKMedians, int type) {
        this.clusters = clusters;
        this.matrix = matrix;
        this.compressionFactor = matrix.length / IDEAL_NUM_ROWS;
        this.useKMedians = useKMedians;
        if (useKMedians) {
            this.metric = RobustManhattanDistance.SINGLETON;
        } else {
            this.metric = RobustEuclideanDistance.SINGLETON;
        }
        this.type = type;
        this.sample = doSubSampling();
    }

    private MatrixSamplingContainer doSubSampling() {
        List<float[][]> subMatrices = new ArrayList<>(clusters.length);
        List<List<Integer>> clusterIndices = new ArrayList<>(clusters.length);
        for (Cluster cluster : clusters) {
            subMatrices.add(getSubMatrix(cluster));
        }

        int n = getTotalNumberOfRows(subMatrices);
        int counter = 0;
        float[][] subSampledMatrix = new float[n][matrix[0].length];
        for (float[][] subMatrix : subMatrices) {
            List<Integer> indices = new ArrayList<>();
            for (float[] row : subMatrix) {
                indices.add(counter);
                System.arraycopy(row, 0, subSampledMatrix[counter], 0, row.length);
                counter++;
            }
            clusterIndices.add(indices);
        }
        return new MatrixSamplingContainer(subSampledMatrix, clusterIndices);
    }

    private int getTotalNumberOfRows(List<float[][]> subMatrices) {
        int counter = 0;
        for (float[][] matrix : subMatrices) {
            counter += matrix.length;
        }
        return counter;
    }

    private float[][] getSubMatrix(Cluster cluster) {
        int numValues = (int) ((cluster.getMemberIndexes().length / compressionFactor) + 1);
        if (type == 0) {
            return EfficientSubsampling.subsampleQuickClustering(matrix,
                    numValues, generator.nextLong(), useKMedians,
                    cluster.getMemberIndexes());
        } else if (type == 1) {
            return EfficientSubsampling.subsampleSpreadOut(matrix,
                    numValues, generator.nextLong(), useKMedians,
                    cluster.getMemberIndexes());
        } else {
            return EfficientSubsampling.subsampleAtRandom(matrix, numValues,
                    generator.nextLong(), cluster.getMemberIndexes());
        }
    }

    public double getScore() {
        Silhouette silhouette = new Silhouette(sample.subSampledMatrix, sample.clusterIndices, metric);
        return silhouette.getScore();
    }
}
