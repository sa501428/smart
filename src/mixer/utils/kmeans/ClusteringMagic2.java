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
import mixer.utils.common.FloatMatrixTools;
import mixer.utils.drive.FinalMatrix;
import mixer.utils.tracks.SubcompartmentInterval;

import java.io.File;

public class ClusteringMagic2 extends ClusteringMagic {
    public ClusteringMagic2(FinalMatrix matrix, File outputDirectory, ChromosomeHandler handler, long seed) {
        super(matrix, outputDirectory, handler, seed);
    }


    public String[] extractFinalGWSubcompartments2(String prefix) {
        System.out.println("\nKmeans clustering");
        matrix.inPlaceScaleSqrtWeightCol();
        int[][] assignments1 = runClusteringOnMatrix(prefix, false);
        float[][] backup = FloatMatrixTools.deepClone(matrix.matrix);

        System.out.println("\nKmedians clustering");
        matrix.inPlaceScaleSqrtWeightCol();
        int[][] assignments2 = runClusteringOnMatrix(prefix, true);

        System.out.println("\nNearest Neighbors Refinement");
        for (int z = 0; z < numClusterSizeKValsUsed; z++) {
            /*
            int numClusters = z + startingClusterSizeK;
            IterativeRefinement.refine(matrix, backup, RobustEuclideanDistance.SINGLETON,
                    numClusters, assignments1[z], handler, outputDirectory,
                    prefix + "_k" + numClusters + "_kmeans_refinement_clusters.bed");
            IterativeRefinement.refine(matrix, matrix.matrix, RobustManhattanDistance.SINGLETON,
                    numClusters, assignments2[z], handler, outputDirectory,
                    prefix + "_k" + numClusters + "_kmedians_refinement_clusters.bed");
            System.out.println('*');
            */
        }
        return new String[0];
    }

    private int[][] runClusteringOnMatrix(String prefix, boolean useKMedians) {
        GenomeWideKmeansRunner kmeansRunner = new GenomeWideKmeansRunner(handler, matrix,
                false, useKMedians);
        int[][] assignments = new int[numClusterSizeKValsUsed][0];
        for (int z = 0; z < numClusterSizeKValsUsed; z++) {
            int numClusters = z + startingClusterSizeK;
            assignments[z] = getGoodResults(kmeansRunner, numClusters, z, prefix, useKMedians);
        }
        System.out.println(".");
        return assignments;
    }

    private int[] getGoodResults(GenomeWideKmeansRunner kmeansRunner, int numClusters,
                                 int z, String prefix, boolean useKMedians) {
        double wcssLimit = Float.MAX_VALUE;
        GenomeWide1DList<SubcompartmentInterval> bestClusters = null;
        int[] bestAssignments = null;
        int attemptsWhichWorked = 0;
        int attempt = 0;
        while (bestAssignments == null || attemptsWhichWorked < 20) {
            System.out.print("-");
            KmeansResult currResult = resetAndRerun(kmeansRunner, generator, numClusters);
            double wcss = currResult.getWithinClusterSumOfSquares();
            if (currResult.getNumActualClusters() == numClusters && wcss < Float.MAX_VALUE) {
                if (wcss < wcssLimit) {
                    wcssLimit = wcss;
                    bestAssignments = currResult.getAssignments(matrix.getNumRows());
                    bestClusters = currResult.getFinalCompartmentsClone();
                }
                attemptsWhichWorked++;
            }
            if (++attempt > 200 && bestAssignments != null) break;
        }
        exportKMeansClusteringResults(z, prefix, useKMedians, bestClusters, null);
        return bestAssignments;
    }
}
