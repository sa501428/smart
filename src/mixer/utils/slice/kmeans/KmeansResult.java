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

package mixer.utils.slice.kmeans;

import javastraw.feature1D.GenomeWideList;
import javastraw.reader.basics.ChromosomeHandler;
import mixer.MixerGlobals;
import mixer.utils.rougheval.SubsamplingManager;
import mixer.utils.similaritymeasures.RobustEuclideanDistance;
import mixer.utils.similaritymeasures.RobustManhattanDistance;
import mixer.utils.slice.matrices.CompositeGenomeWideMatrix;
import mixer.utils.slice.structures.SubcompartmentInterval;
import robust.concurrent.kmeans.clustering.Cluster;

import java.util.ArrayList;
import java.util.List;

public class KmeansResult {

    private static final int MIN_EXPECTED_CLUSTER_SIZE = 5;
    private static final double NUM_ITERS = 5;
    private final int numClustersDesired;
    private final GenomeWideList<SubcompartmentInterval> finalCompartments;
    private final List<List<Integer>> indicesMap = new ArrayList<>();
    private int numActualClusters = 0;
    private double wcss = 0;
    private double silhouette = 0;

    public KmeansResult(int numClusters, ChromosomeHandler chromosomeHandler) {
        numClustersDesired = numClusters;
        finalCompartments = new GenomeWideList<>(chromosomeHandler);
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

    public double getSilhouette() {
        return silhouette;
    }

    public GenomeWideList<SubcompartmentInterval> getFinalCompartmentsClone() {
        return finalCompartments.deepClone();
    }

    public List<List<Integer>> getIndicesMapClone() {
        List<List<Integer>> output = new ArrayList<>();
        for (List<Integer> group : indicesMap) {
            output.add(getDeepCopy(group));
        }
        return output;
    }

    private List<Integer> getDeepCopy(List<Integer> listA) {
        List<Integer> output = new ArrayList<>();
        for (int a : listA) {
            output.add(a);
        }
        return output;
    }

    public void processResultAndUpdateScoringMetrics(Cluster[] clusters, CompositeGenomeWideMatrix matrix, boolean useKMedians, boolean useCorrMatrix) {
        populateIndicesMap(clusters);
        matrix.processKMeansClusteringResult(clusters, finalCompartments);
        wcss = getWCSS(clusters, matrix, useCorrMatrix, useKMedians);
        silhouette = getSilhouette(clusters, matrix, useCorrMatrix, useKMedians);
        numActualClusters = clusters.length;
    }

    public double getWCSS(Cluster[] clusters, CompositeGenomeWideMatrix matrix,
                          boolean useCorr, boolean useKMedians) {
        double withinClusterSumOfSquares = 0;

        for (int z = 0; z < clusters.length; z++) {
            Cluster cluster = clusters[z];

            if (cluster.getMemberIndexes().length < MIN_EXPECTED_CLUSTER_SIZE) {
                withinClusterSumOfSquares += Float.MAX_VALUE;
            }

            float[][] vectors = matrix.getData(useCorr);
            for (int i : cluster.getMemberIndexes()) {
                withinClusterSumOfSquares += getDistance(cluster.getCenter(), vectors[i], useKMedians);
            }
        }

        withinClusterSumOfSquares = withinClusterSumOfSquares / clusters.length;
        if (MixerGlobals.printVerboseComments) {
            System.out.println("Final WCSS " + withinClusterSumOfSquares);
        }

        return withinClusterSumOfSquares;
    }

    private double getSilhouette(Cluster[] clusters, CompositeGenomeWideMatrix matrix, boolean useCorr, boolean useKMedians) {
        double score = 0;
        SubsamplingManager manager = new SubsamplingManager(clusters, matrix.getData(useCorr), useKMedians);
        for (int iter = 0; iter < NUM_ITERS; iter++) {
            score += manager.getScore();
        }
        score = score / NUM_ITERS;
        //System.out.println("Silhouette: " + score);
        return score;
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
}
