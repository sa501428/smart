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

package mixer.utils.kmeans;

import javastraw.feature1D.GenomeWide1DList;
import javastraw.reader.basics.ChromosomeHandler;
import mixer.SmartTools;
import mixer.utils.drive.FinalMatrix;
import mixer.utils.similaritymeasures.RobustCorrelationSimilarity;
import mixer.utils.similaritymeasures.RobustEuclideanDistance;
import mixer.utils.similaritymeasures.RobustManhattanDistance;
import mixer.utils.tracks.SubcompartmentInterval;
import robust.concurrent.kmeans.clustering.Cluster;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

public class KmeansResult {

    private static final int MIN_EXPECTED_CLUSTER_SIZE = 5;
    private final int numClustersDesired;
    private final GenomeWide1DList<SubcompartmentInterval> finalCompartments;
    private final List<List<Integer>> indicesMap = new ArrayList<>();
    private int numActualClusters = 0;
    private double wcss = 0;

    public KmeansResult(int numClusters, ChromosomeHandler chromosomeHandler) {
        numClustersDesired = numClusters;
        finalCompartments = new GenomeWide1DList<>(chromosomeHandler);
    }

    public int getNumClustersDesired() {
        return numClustersDesired;
    }

    public int getNumActualClusters() {
        return numActualClusters;
    }

    public double getWithinClusterSumOfSquares() {
        return wcss;
    }

    public GenomeWide1DList<SubcompartmentInterval> getFinalCompartmentsClone() {
        return finalCompartments.deepClone();
    }

    public void processResultAndUpdateScoringMetrics(Cluster[] clusters, FinalMatrix matrix,
                                                     boolean useKMedians, boolean useCorrMatrix) {
        populateIndicesMap(clusters);
        matrix.processKMeansClusteringResult(clusters, finalCompartments);
        wcss = getWCSS(clusters, matrix, useCorrMatrix, useKMedians);
        //clusterCorrMatrix = calculateCorrelations(clusters);
        //worstCorr = getMaxOffDiag(clusterCorrMatrix);
        numActualClusters = clusters.length;
    }

    private float getMaxOffDiag(float[][] matrix) {
        int n = matrix.length;
        float maxVal = -1;
        for (int i = 0; i < n; i++) {
            for (int j = i + 1; j < n; j++) {
                maxVal = Math.max(maxVal, matrix[i][j]);
            }
        }
        return maxVal;
    }

    public double getWCSS(Cluster[] clusters, FinalMatrix matrix,
                          boolean useCorr, boolean useKMedians) {
        double withinClusterSumOfSquares = 0;

        for (int z = 0; z < clusters.length; z++) {
            Cluster cluster = clusters[z];

            if (cluster.getMemberIndexes().length < MIN_EXPECTED_CLUSTER_SIZE) {
                withinClusterSumOfSquares += Float.MAX_VALUE;
            }

            float[][] vectors = matrix.matrix;
            for (int i : cluster.getMemberIndexes()) {
                withinClusterSumOfSquares += getDistance(cluster.getCenter(), vectors[i], useKMedians);
            }
        }

        withinClusterSumOfSquares = withinClusterSumOfSquares / clusters.length;
        if (SmartTools.printVerboseComments) {
            System.out.println("Final WCSS " + withinClusterSumOfSquares);
        }

        return withinClusterSumOfSquares;
    }

    private double getDistance(float[] center, float[] vector, boolean useKMedians) {
        if (useKMedians) {
            return RobustManhattanDistance.SINGLETON.distance(center, vector);
        }
        return RobustEuclideanDistance.getNonNanMeanSquaredError(center, vector);
    }

    private void populateIndicesMap(Cluster[] clusters) {
        indicesMap.clear();
        for (Cluster cluster : clusters) {
            List<Integer> group = new ArrayList<>();
            for (int member : cluster.getMemberIndexes()) {
                group.add(member);
            }
            indicesMap.add(group);
        }
    }

    private float[][] calculateCorrelations(Cluster[] clusters) {
        int n = clusters.length;
        float[][] matrix = new float[n][n];
        for (int i = 0; i < n; i++) {
            matrix[i][i] = 1;
            for (int j = i + 1; j < n; j++) {
                matrix[i][j] = RobustCorrelationSimilarity.SINGLETON.distance(
                        clusters[i].getCenter(), clusters[j].getCenter());
                matrix[j][i] = matrix[i][j];
            }
        }
        return matrix;
    }

    public int[] getAssignments(int numRows) {
        int[] result = new int[numRows];
        Arrays.fill(result, -1);

        for (int c = 0; c < indicesMap.size(); c++) {
            for (int i : indicesMap.get(c)) {
                result[i] = c;
            }
        }

        return result;
    }
}
