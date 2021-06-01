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

package mixer.utils.slice.matrices;

import javastraw.feature1D.GenomeWideList;
import javastraw.reader.Dataset;
import javastraw.reader.basics.Chromosome;
import javastraw.reader.basics.ChromosomeHandler;
import javastraw.reader.type.NormalizationType;
import javastraw.tools.MatrixTools;
import mixer.MixerGlobals;
import mixer.utils.common.FloatMatrixTools;
import mixer.utils.similaritymeasures.RobustEuclideanDistance;
import mixer.utils.similaritymeasures.SimilarityMetric;
import mixer.utils.slice.cleaning.GWBadIndexFinder;
import mixer.utils.slice.cleaning.MatrixImputer;
import mixer.utils.slice.cleaning.SliceMatrixCleaner;
import mixer.utils.slice.kmeans.kmeansfloat.Cluster;
import mixer.utils.slice.structures.SliceUtils;
import mixer.utils.slice.structures.SubcompartmentInterval;

import java.io.File;
import java.util.*;

public abstract class CompositeGenomeWideMatrix {
    protected final NormalizationType intraNorm, interNorm;
    protected final int resolution;
    protected final Map<Integer, SubcompartmentInterval> rowIndexToIntervalMap = new HashMap<>();
    protected final Chromosome[] chromosomes;
    protected final Random generator = new Random(0);
    protected final File outputDirectory;
    private float[][] gwCleanMatrix, corrData;
    private final int[] weights;
    protected final GWBadIndexFinder badIndexLocations;
    protected final SimilarityMetric metric;
    protected final int maxClusterSizeExpected;

    public CompositeGenomeWideMatrix(ChromosomeHandler chromosomeHandler, Dataset ds,
                                     NormalizationType intraNorm, NormalizationType interNorm,
                                     int resolution,
                                     File outputDirectory, long seed,
                                     GWBadIndexFinder badIndexLocations, SimilarityMetric metric,
                                     int maxClusterSizeExpected) {
        this.maxClusterSizeExpected = maxClusterSizeExpected;
        this.intraNorm = intraNorm;
        this.interNorm = interNorm;
        System.out.println("Intra: " + intraNorm.getLabel() + " Inter: " + interNorm.getLabel());
        this.resolution = resolution;
        this.outputDirectory = outputDirectory;
        this.generator.setSeed(seed);
        this.badIndexLocations = badIndexLocations;
        this.metric = metric;

        chromosomes = chromosomeHandler.getAutosomalChromosomesArray();
        MatrixAndWeight mw = makeCleanScaledInterMatrix(ds);
        gwCleanMatrix = mw.matrix;
        this.weights = mw.weights;
    }

    abstract MatrixAndWeight makeCleanScaledInterMatrix(Dataset ds);

    public void cleanUpMatricesBySparsity() {
        SliceMatrixCleaner matrixCleanupReduction = new SliceMatrixCleaner(gwCleanMatrix,
                generator.nextLong(), outputDirectory, metric);
        gwCleanMatrix = matrixCleanupReduction.getCleanFilteredZscoredMatrix(rowIndexToIntervalMap, weights);
        corrData = MatrixImputer.imputeUntilNoNansOnlyNN(gwCleanMatrix);
        //corrData = SimilarityMatrixTools.getCosinePearsonCorrMatrix(gwCleanMatrix, 50, generator.nextLong());

        File file1 = new File(outputDirectory, "zscore_matrix.npy");
        MatrixTools.saveMatrixTextNumpy(file1.getAbsolutePath(), gwCleanMatrix);

        File file2 = new File(outputDirectory, "corr_matrix.npy");
        MatrixTools.saveMatrixTextNumpy(file2.getAbsolutePath(), corrData);
    }

    public synchronized double processKMeansClusteringResult(Cluster[] clusters,
                                                             GenomeWideList<SubcompartmentInterval> subcompartments,
                                                             boolean useCorr) {

        Set<SubcompartmentInterval> subcompartmentIntervals = new HashSet<>();
        if (MixerGlobals.printVerboseComments) {
            System.out.println("GW Composite data vs clustered into " + clusters.length + " clusters");
        }

        double withinClusterSumOfSquares = 0;
        int numGoodClusters = 0;
        int genomewideCompartmentID = 0;

        Set<Integer> badIndices = new HashSet<>();

        for (int z = 0; z < clusters.length; z++) {
            Cluster cluster = clusters[z];
            int currentClusterID = ++genomewideCompartmentID;

            if (MixerGlobals.printVerboseComments) {
                System.out.println("Size of cluster " + currentClusterID + " - " + cluster.getMemberIndexes().length);
            }

            if (cluster.getMemberIndexes().length < 5) {
                for (int badIndex : cluster.getMemberIndexes()) {
                    badIndices.add(badIndex);
                }
                withinClusterSumOfSquares += Float.MAX_VALUE;
            } else {
                numGoodClusters++;
            }

            for (int i : cluster.getMemberIndexes()) {
                if (useCorr) {
                    withinClusterSumOfSquares +=
                            RobustEuclideanDistance.getNonNanMeanSquaredError(cluster.getCenter(), corrData[i]);
                } else {
                    withinClusterSumOfSquares +=
                            RobustEuclideanDistance.getNonNanMeanSquaredError(cluster.getCenter(), gwCleanMatrix[i]);
                }

                if (rowIndexToIntervalMap.containsKey(i)) {
                    SubcompartmentInterval interv = rowIndexToIntervalMap.get(i);
                    if (interv != null) {
                        subcompartmentIntervals.add(generateNewSubcompartment(interv, currentClusterID));
                    }
                }
            }
        }

        withinClusterSumOfSquares = withinClusterSumOfSquares / numGoodClusters;
        if (MixerGlobals.printVerboseComments) {
            System.out.println("Final WCSS " + withinClusterSumOfSquares);
        }

        subcompartments.addAll(new ArrayList<>(subcompartmentIntervals));
        SliceUtils.reSort(subcompartments);

        return withinClusterSumOfSquares;
    }

    public synchronized void processGMMClusteringResult(int[] clusterID,
                                                        GenomeWideList<SubcompartmentInterval> subcompartments) {

        Set<SubcompartmentInterval> subcompartmentIntervals = new HashSet<>();

        for (int i = 0; i < corrData.length; i++) {
            if (rowIndexToIntervalMap.containsKey(i)) {
                SubcompartmentInterval interv = rowIndexToIntervalMap.get(i);
                if (interv != null) {
                    subcompartmentIntervals.add(generateNewSubcompartment(interv, clusterID[i]));
                }
            }
        }

        subcompartments.addAll(new ArrayList<>(subcompartmentIntervals));
        SliceUtils.reSort(subcompartments);
    }

    protected SubcompartmentInterval generateNewSubcompartment(SubcompartmentInterval interv, int currentClusterID) {
        int chrIndex = interv.getChrIndex();
        String chrName = interv.getChrName();
        int x1 = interv.getX1();
        int x2 = interv.getX2();
        return new SubcompartmentInterval(chrIndex, chrName, x1, x2, currentClusterID);
    }

    public float[][] getData(boolean getCorrMatrix) {
        if (getCorrMatrix) {
            return corrData;
        }
        return gwCleanMatrix;
    }

    public void appendDataAlongExistingRows(CompositeGenomeWideMatrix additionalData) {
        if (gwCleanMatrix.length != additionalData.gwCleanMatrix.length) {
            System.err.println("***************************************\n" +
                    "Dimension mismatch: " + gwCleanMatrix.length + " != " + additionalData.gwCleanMatrix.length);
        } else {
            gwCleanMatrix = FloatMatrixTools.concatenate(gwCleanMatrix, additionalData.gwCleanMatrix);
        }
    }

    public GWBadIndexFinder getBadIndices() {
        return badIndexLocations;
    }
}
