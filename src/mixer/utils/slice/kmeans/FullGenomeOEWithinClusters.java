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

import javastraw.featurelist.GenomeWideList;
import javastraw.reader.ChromosomeHandler;
import javastraw.reader.Dataset;
import javastraw.tools.MatrixTools;
import javastraw.type.NormalizationType;
import mixer.MixerGlobals;
import mixer.utils.similaritymeasures.SimilarityMetric;
import mixer.utils.slice.cleaning.GWBadIndexFinder;
import mixer.utils.slice.cleaning.LeftOverClusterIdentifier;
import mixer.utils.slice.kmeans.kmeansfloat.Cluster;
import mixer.utils.slice.kmeans.kmeansfloat.ClusterTools;
import mixer.utils.slice.matrices.SliceMatrix;
import mixer.utils.slice.structures.SliceUtils;
import mixer.utils.slice.structures.SubcompartmentInterval;

import java.io.File;
import java.util.*;

public class FullGenomeOEWithinClusters {
    public static int startingClusterSizeK = 2;
    public static int numClusterSizeKValsUsed = 10;
    public static int numAttemptsForKMeans = 3;
    private final List<Dataset> datasets;
    protected final File outputDirectory;
    private final ChromosomeHandler chromosomeHandler;
    private final int resolution;
    private final NormalizationType[] norms;
    private final SliceMatrix sliceMatrix;

    private final Random generator = new Random(0);

    public FullGenomeOEWithinClusters(List<Dataset> datasets, ChromosomeHandler chromosomeHandler, int resolution, NormalizationType[] norms,
                                      File outputDirectory, long seed, String[] referenceBedFiles,
                                      SimilarityMetric metric) {
        this.chromosomeHandler = chromosomeHandler;
        this.resolution = resolution;
        this.norms = norms;
        this.outputDirectory = outputDirectory;
        this.datasets = datasets;
        generator.setSeed(seed);

        GWBadIndexFinder badIndexFinder = new GWBadIndexFinder(chromosomeHandler.getAutosomalChromosomesArray(),
                resolution, norms);
        badIndexFinder.createInternalBadList(datasets, chromosomeHandler.getAutosomalChromosomesArray());

        sliceMatrix = new SliceMatrix(
                chromosomeHandler, datasets.get(0), norms[0], resolution, outputDirectory, generator.nextLong(), referenceBedFiles,
                badIndexFinder, metric);

        for (int dI = 1; dI < datasets.size(); dI++) {
            SliceMatrix additionalData = new SliceMatrix(chromosomeHandler, datasets.get(dI),
                    norms[dI], resolution, outputDirectory, generator.nextLong(), new String[]{}, badIndexFinder, metric);
            sliceMatrix.appendDataAlongExistingRows(additionalData);
        }

        sliceMatrix.cleanUpMatricesBySparsity();
    }

    public void extractFinalGWSubcompartments(List<String> inputHicFilePaths,
                                              String prefix, int index, boolean compareMaps) {

        Map<Integer, GenomeWideList<SubcompartmentInterval>> numItersToResults = new HashMap<>();

        if (MixerGlobals.printVerboseComments) {
            sliceMatrix.exportData();
        }

        GenomeWideKmeansRunner kmeansRunner = new GenomeWideKmeansRunner(chromosomeHandler, sliceMatrix);

        double[][] iterToWcssAicBic = new double[4][numClusterSizeKValsUsed];
        Arrays.fill(iterToWcssAicBic[1], Double.MAX_VALUE);
        Arrays.fill(iterToWcssAicBic[2], Double.MAX_VALUE);
        Arrays.fill(iterToWcssAicBic[3], Double.MAX_VALUE);

        System.out.println("Genomewide clustering");
        for (int z = 0; z < numClusterSizeKValsUsed; z++) {

            int k = z + startingClusterSizeK;
            Cluster[] bestClusters = null;
            int[] bestIDs = null;
            int[][] novelIDsForIndx = null;

            for (int p = 0; p < numAttemptsForKMeans; p++) {

                kmeansRunner.prepareForNewRun(k);
                kmeansRunner.launchKmeansGWMatrix(generator.nextLong());

                int numActualClustersThisAttempt = kmeansRunner.getNumActualClusters();
                double wcss = kmeansRunner.getWithinClusterSumOfSquares();

                if (wcss < iterToWcssAicBic[1][z]) {
                    setMseAicBicValues(z, iterToWcssAicBic, numActualClustersThisAttempt, wcss);
                    numItersToResults.put(k, kmeansRunner.getFinalCompartments());
                    bestClusters = kmeansRunner.getRecentClustersClone();
                    bestIDs = kmeansRunner.getRecentIDsClone();
                    novelIDsForIndx = kmeansRunner.getRecentIDsForIndex();
                }
                System.out.print(".");
            }

            ClusterTools.performStatisticalAnalysisBetweenClusters(outputDirectory, "final_gw_" + k, bestClusters, bestIDs);
            MatrixTools.saveMatrixTextNumpy((new File(outputDirectory, "novel_ids_for_index_" + k + ".npy")).getAbsolutePath(), novelIDsForIndx);
        }
        System.out.println(".");

        if (!compareMaps && false) {
            System.out.println("Post processing");
            LeftOverClusterIdentifier identifier = new LeftOverClusterIdentifier(chromosomeHandler, datasets.get(0), norms[0], resolution);
            identifier.identify(numItersToResults, sliceMatrix.getBadIndices());
        }

        MatrixTools.saveMatrixTextNumpy(new File(outputDirectory, "clusterSize_WCSS_AIC_BIC.npy").getAbsolutePath(), iterToWcssAicBic);

        String hicFileName = SliceUtils.cleanUpPath(inputHicFilePaths.get(index));
        for (Integer key : numItersToResults.keySet()) {
            GenomeWideList<SubcompartmentInterval> gwList = numItersToResults.get(key);
            SliceUtils.collapseGWList(gwList);
            File outBedFile = new File(outputDirectory, prefix + key + "_clusters_" + hicFileName + ".subcompartment.bed");
            gwList.simpleExport(outBedFile);
        }
    }

    private void setMseAicBicValues(int z, double[][] iterToWcssAicBic, int numClusters, double sumOfSquares) {
        iterToWcssAicBic[0][z] = numClusters;
        iterToWcssAicBic[1][z] = sumOfSquares;
        // AIC
        iterToWcssAicBic[2][z] = sumOfSquares + 2 * sliceMatrix.getWidth() * numClusters;
        // BIC .5*k*d*log(n)
        iterToWcssAicBic[3][z] = sumOfSquares + 0.5 * sliceMatrix.getWidth() * numClusters * Math.log(sliceMatrix.getLength());
    }
}
