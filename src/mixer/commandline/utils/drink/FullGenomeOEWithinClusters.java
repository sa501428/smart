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

package mixer.commandline.utils.drink;

import mixer.MixerGlobals;
import mixer.commandline.utils.common.DoubleMatrixTools;
import mixer.commandline.utils.drink.kmeansfloat.Cluster;
import mixer.commandline.utils.drink.kmeansfloat.ClusterTools;
import mixer.data.ChromosomeHandler;
import mixer.data.Dataset;
import mixer.data.feature.GenomeWideList;
import mixer.windowui.NormalizationType;

import java.io.File;
import java.util.*;

public class FullGenomeOEWithinClusters {
    private final Dataset ds;
    private final ChromosomeHandler chromosomeHandler;
    private final int resolution;
    private final NormalizationType norm;
    private final GenomeWideList<SubcompartmentInterval> origIntraSubcompartments;
    private final int numRounds = 10;
    private int minIntervalSizeAllowed; // 1
    private final int numAttemptsForKMeans = 5;
    private final CompositeGenomeWideDensityMatrix interMatrix;
    private final float oeThreshold;

    public FullGenomeOEWithinClusters(Dataset ds, ChromosomeHandler chromosomeHandler, int resolution, NormalizationType norm,
                                      GenomeWideList<SubcompartmentInterval> origIntraSubcompartments, float oeThreshold, int minIntervalSizeAllowed) {
        this.ds = ds;
        this.chromosomeHandler = chromosomeHandler;
        this.resolution = resolution;
        this.norm = norm;
        this.oeThreshold = oeThreshold;
        DrinkUtils.collapseGWList(origIntraSubcompartments);
        this.origIntraSubcompartments = origIntraSubcompartments;
        this.minIntervalSizeAllowed = minIntervalSizeAllowed;

        interMatrix = new CompositeGenomeWideDensityMatrix(
                chromosomeHandler, ds, norm, resolution, origIntraSubcompartments, oeThreshold, minIntervalSizeAllowed);

        System.gc();
    }

    public void appendGWDataFromAdditionalDataset(Dataset ds2) {

        CompositeGenomeWideDensityMatrix additionalData = new CompositeGenomeWideDensityMatrix(
                chromosomeHandler, ds2, norm, resolution, origIntraSubcompartments, oeThreshold, minIntervalSizeAllowed);

        interMatrix.appendDataAlongExistingRows(additionalData);
    }

    public void extractFinalGWSubcompartments(File outputDirectory, Random generator, List<String> inputHicFilePaths,
                                              String prefix, int index, double[] convolution) {

        Map<Integer, GenomeWideList<SubcompartmentInterval>> numItersToResults = new HashMap<>();

        if (MixerGlobals.printVerboseComments) {
            interMatrix.exportData(outputDirectory);
        }

        GenomeWideKmeansRunner kmeansRunner = new GenomeWideKmeansRunner(chromosomeHandler, interMatrix);

        double[][] iterToWcssAicBic = new double[4][numRounds];
        Arrays.fill(iterToWcssAicBic[1], Double.MAX_VALUE);
        Arrays.fill(iterToWcssAicBic[2], Double.MAX_VALUE);
        Arrays.fill(iterToWcssAicBic[3], Double.MAX_VALUE);

        for (int z = 0; z < numRounds; z++) {

            int k = z + 2;
            Cluster[] bestClusters = null;
            int[] bestIDs = null;

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
                }
            }

            ClusterTools.performStatisticalAnalysisBetweenClusters(outputDirectory, "final_gw_" + k, bestClusters, bestIDs);
        }

        if (minIntervalSizeAllowed > 1) {
            LeftOverClusterIdentifier.identify(chromosomeHandler, ds, norm, resolution, numItersToResults,
                    origIntraSubcompartments, minIntervalSizeAllowed, oeThreshold, convolution);
        }

        DoubleMatrixTools.saveMatrixTextNumpy(new File(outputDirectory, "clusterSize_WCSS_AIC_BIC.npy").getAbsolutePath(), iterToWcssAicBic);

        String hicFileName = DrinkUtils.cleanUpPath(inputHicFilePaths.get(index));
        for (Integer key : numItersToResults.keySet()) {
            GenomeWideList<SubcompartmentInterval> gwList = numItersToResults.get(key);
            DrinkUtils.collapseGWList(gwList);
            File outBedFile = new File(outputDirectory, prefix + key + "_clusters_" + hicFileName + ".subcompartment.bed");
            gwList.simpleExport(outBedFile);
        }
    }

    private void setMseAicBicValues(int z, double[][] iterToWcssAicBic, int numClusters, double sumOfSquares) {
        iterToWcssAicBic[0][z] = numClusters;
        iterToWcssAicBic[1][z] = sumOfSquares;
        // AIC
        iterToWcssAicBic[2][z] = sumOfSquares + 2 * interMatrix.getWidth() * numClusters;
        // BIC .5*k*d*log(n)
        iterToWcssAicBic[3][z] = sumOfSquares + 0.5 * interMatrix.getWidth() * numClusters * Math.log(interMatrix.getLength());
    }
}
