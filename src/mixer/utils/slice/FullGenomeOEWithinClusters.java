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

package mixer.utils.slice;

import javastraw.featurelist.GenomeWideList;
import javastraw.reader.ChromosomeHandler;
import javastraw.reader.Dataset;
import javastraw.type.NormalizationType;
import mixer.MixerGlobals;
import mixer.utils.common.DoubleMatrixTools;
import mixer.utils.common.IntMatrixTools;
import mixer.utils.slice.kmeansfloat.Cluster;
import mixer.utils.slice.kmeansfloat.ClusterTools;

import java.io.File;
import java.util.*;

public class FullGenomeOEWithinClusters {
    private final Dataset ds;
    private final ChromosomeHandler chromosomeHandler;
    private final int resolution;
    private final NormalizationType norm;
    protected final File outputDirectory;
    private final int numRounds = 2;//10
    private final CompositeGenomeWideDensityMatrix interMatrix;
    private final int numAttemptsForKMeans = 10;//5 //7 //5
    private final Random generator;
    
    public FullGenomeOEWithinClusters(Dataset ds, ChromosomeHandler chromosomeHandler, int resolution, NormalizationType norm,
                                      int minIntervalSizeAllowed, File outputDirectory, Random generator, String[] referenceBedFiles) {
        this.ds = ds;
        this.chromosomeHandler = chromosomeHandler;
        this.resolution = resolution;
        this.norm = norm;
        this.outputDirectory = outputDirectory;
        this.generator = generator;

        if (minIntervalSizeAllowed > 0) {
            SliceMatrix.numColumnsToPutTogether = minIntervalSizeAllowed;
        } else {
            SliceMatrix.numColumnsToPutTogether = 200000 / resolution;
        }
        interMatrix = new SliceMatrix(
                chromosomeHandler, ds, norm, resolution, outputDirectory, generator, referenceBedFiles);
        
        System.gc();
    }

    public void appendGWDataFromAdditionalDataset(Dataset ds2) {
        CompositeGenomeWideDensityMatrix additionalData;
        additionalData = new SliceMatrix(chromosomeHandler, ds2, norm, resolution, outputDirectory, generator, new String[]{});
        interMatrix.appendDataAlongExistingRows(additionalData);
    }
    
    public void extractFinalGWSubcompartments(Random generator, List<String> inputHicFilePaths,
                                              String prefix, int index) {
        
        Map<Integer, GenomeWideList<SubcompartmentInterval>> numItersToResults = new HashMap<>();
        
        if (MixerGlobals.printVerboseComments) {
            interMatrix.exportData();
        }
        
        GenomeWideKmeansRunner kmeansRunner = new GenomeWideKmeansRunner(chromosomeHandler, interMatrix);

        double[][] iterToWcssAicBic = new double[4][numRounds];
        Arrays.fill(iterToWcssAicBic[1], Double.MAX_VALUE);
        Arrays.fill(iterToWcssAicBic[2], Double.MAX_VALUE);
        Arrays.fill(iterToWcssAicBic[3], Double.MAX_VALUE);

        System.out.println("Genomewide clustering");
        for (int z = 0; z < numRounds; z++) {
    
            int k = z + 5; //5
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
            IntMatrixTools.saveMatrixTextNumpy((new File(outputDirectory, "novel_ids_for_index_" + k + ".npy")).getAbsolutePath(), novelIDsForIndx);
        }
        System.out.println(".");

        System.out.println("Post processing");
        LeftOverClusterIdentifier.identify(chromosomeHandler, ds, norm, resolution, numItersToResults,
                interMatrix.getBadIndices());

        DoubleMatrixTools.saveMatrixTextNumpy(new File(outputDirectory, "clusterSize_WCSS_AIC_BIC.npy").getAbsolutePath(), iterToWcssAicBic);
    
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
        iterToWcssAicBic[2][z] = sumOfSquares + 2 * interMatrix.getWidth() * numClusters;
        // BIC .5*k*d*log(n)
        iterToWcssAicBic[3][z] = sumOfSquares + 0.5 * interMatrix.getWidth() * numClusters * Math.log(interMatrix.getLength());
    }
}
