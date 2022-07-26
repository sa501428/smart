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

package mixer.utils.slice;

import javastraw.feature1D.GenomeWide1DList;
import mixer.utils.slice.gmm.GenomeWideGMMRunner;
import mixer.utils.slice.kmeans.FullGenomeOEWithinClusters;
import mixer.utils.slice.kmeans.GenomeWideKmeansRunner;
import mixer.utils.slice.kmeans.KmeansEvaluator;
import mixer.utils.slice.structures.SliceUtils;
import mixer.utils.slice.structures.SubcompartmentInterval;

import java.io.File;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

public class CorrMatrixClusterer {

    public static void runClusteringOnCorrMatrix(FullGenomeOEWithinClusters parent, String prefix, boolean useKmedians) {

        Map<Integer, GenomeWide1DList<SubcompartmentInterval>> gmmClustersToResults = new HashMap<>();
        Map<Integer, GenomeWide1DList<SubcompartmentInterval>> kmeansClustersToResults = new HashMap<>();
        Map<Integer, List<List<Integer>>> kmeansIndicesMap = new HashMap<>();

        GenomeWideKmeansRunner kmeansRunner = new GenomeWideKmeansRunner(parent.getChromosomeHandler(),
                parent.getSliceMatrix(), true, useKmedians);
        KmeansEvaluator evaluator = new KmeansEvaluator(FullGenomeOEWithinClusters.numClusterSizeKValsUsed);

        for (int z = 0; z < FullGenomeOEWithinClusters.numClusterSizeKValsUsed; z++) {
            parent.runRepeatedKMeansClusteringLoop(FullGenomeOEWithinClusters.numAttemptsForKMeans, kmeansRunner, evaluator, z,
                    parent.getMaxIters(), kmeansClustersToResults, kmeansIndicesMap);
            parent.exportKMeansClusteringResults(z, kmeansClustersToResults, prefix, kmeansIndicesMap, useKmedians);

            runGMMClusteringLoop(z, 20, kmeansIndicesMap.get(z), gmmClustersToResults, parent);
            exportGMMClusteringResults(z, gmmClustersToResults, prefix, parent);
        }
        parent.exportEvaluatorInfo(evaluator, useKmedians);

        System.out.println(".");
    }

    private static void exportGMMClusteringResults(int z, Map<Integer, GenomeWide1DList<SubcompartmentInterval>> numClustersToResults,
                                                   String prefix, FullGenomeOEWithinClusters parent) {
        int k = z + FullGenomeOEWithinClusters.startingClusterSizeK;
        if (numClustersToResults.containsKey(k)) {
            GenomeWide1DList<SubcompartmentInterval> gwList = numClustersToResults.get(k);
            SliceUtils.collapseGWList(gwList);
            File outBedFile = new File(parent.getOutputDirectory(), prefix + "_" + k + "_gmm_clusters.bed");
            gwList.simpleExport(outBedFile);
        }
    }

    private static void runGMMClusteringLoop(int z, int maxIters, List<List<Integer>> startingIndices,
                                             Map<Integer, GenomeWide1DList<SubcompartmentInterval>> numClustersToResults,
                                             FullGenomeOEWithinClusters parent) {
        int numClusters = z + FullGenomeOEWithinClusters.startingClusterSizeK;
        GenomeWideGMMRunner gmmRunner = new GenomeWideGMMRunner(parent.getChromosomeHandler(), parent.getSliceMatrix());
        gmmRunner.launch(numClusters, maxIters, numClustersToResults, startingIndices, false);
        System.out.print("*");
    }

}
