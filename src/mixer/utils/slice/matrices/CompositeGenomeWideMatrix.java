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
import mixer.algos.Slice;
import mixer.utils.common.ArrayTools;
import mixer.utils.common.FloatMatrixTools;
import mixer.utils.common.ZScoreTools;
import mixer.utils.rougheval.SubsamplingManager;
import mixer.utils.similaritymeasures.RobustEuclideanDistance;
import mixer.utils.similaritymeasures.RobustManhattanDistance;
import mixer.utils.slice.cleaning.GWBadIndexFinder;
import mixer.utils.slice.cleaning.SimilarityMatrixTools;
import mixer.utils.slice.cleaning.SliceMatrixCleaner;
import mixer.utils.slice.gmm.SimpleScatterPlot;
import mixer.utils.slice.structures.SliceUtils;
import mixer.utils.slice.structures.SubcompartmentInterval;
import robust.concurrent.kmeans.clustering.Cluster;

import java.io.File;
import java.util.*;

public abstract class CompositeGenomeWideMatrix {
    private static final int MIN_EXPECTED_CLUSTER_SIZE = 5;
    private static final double NUM_ITERS = 5;
    protected final NormalizationType[] norms;
    protected final int resolution;
    protected final Map<Integer, SubcompartmentInterval> rowIndexToIntervalMap = new HashMap<>();
    protected final Chromosome[] chromosomes;
    protected final Random generator = new Random(0);
    protected final File outputDirectory;
    private MatrixAndWeight gwCleanMatrix, projectedData = null;
    private float[][] umapProjection;
    protected final GWBadIndexFinder badIndexLocations;
    protected final int maxClusterSizeExpected;

    public CompositeGenomeWideMatrix(ChromosomeHandler chromosomeHandler, Dataset ds,
                                     NormalizationType[] norms,
                                     int resolution,
                                     File outputDirectory, long seed,
                                     GWBadIndexFinder badIndexLocations,
                                     int maxClusterSizeExpected) {
        this.maxClusterSizeExpected = maxClusterSizeExpected;
        this.norms = norms;
        if (MixerGlobals.printVerboseComments) {
            System.out.println("Norms: " + Arrays.toString(norms));
        }
        this.resolution = resolution;
        this.outputDirectory = outputDirectory;
        this.generator.setSeed(seed);
        this.badIndexLocations = badIndexLocations;

        chromosomes = chromosomeHandler.getAutosomalChromosomesArray();
        gwCleanMatrix = makeCleanScaledInterMatrix(ds, norms[Slice.INTER_SCALE_INDEX]);
        //gwCleanMatrix = makeComboNormMatrix(ds, norms);
    }

    private MatrixAndWeight makeComboNormMatrix(Dataset ds, NormalizationType[] norms) {
        return concatenate(
                makeCleanScaledInterMatrix(ds, norms[Slice.INTER_SCALE_INDEX]),
                makeCleanScaledInterMatrix(ds, norms[Slice.GW_SCALE_INDEX])
        );
    }

    abstract MatrixAndWeight makeCleanScaledInterMatrix(Dataset ds, NormalizationType interNorm);

    public void cleanUpMatricesBySparsity() {

        SliceMatrixCleaner matrixCleanupReduction = new SliceMatrixCleaner(gwCleanMatrix.matrix,
                generator.nextLong(), outputDirectory, resolution);
        gwCleanMatrix = matrixCleanupReduction.getCleanFilteredZscoredMatrix(rowIndexToIntervalMap,
                gwCleanMatrix.weights);

        inPlaceScaleSqrtWeightCol();

        File file0 = new File(outputDirectory, "genome_indices.npy");
        MatrixTools.saveMatrixTextNumpy(file0.getAbsolutePath(), getGenomeIndices());

        File file1 = new File(outputDirectory, "slice_matrix.npy");
        MatrixTools.saveMatrixTextNumpy(file1.getAbsolutePath(), gwCleanMatrix.matrix);

        if (Slice.USE_INTER_CORR_CLUSTERING || Slice.PROJECT_TO_UMAP) {
            projectedData = new MatrixAndWeight(SimilarityMatrixTools.getCosinePearsonCorrMatrix(gwCleanMatrix.matrix,
                    50, generator.nextLong()), gwCleanMatrix.weights);

            umapProjection = SimpleScatterPlot.getUmapProjection2D(projectedData.matrix);

            File file2 = new File(outputDirectory, "corr_slice_matrix.npy");
            MatrixTools.saveMatrixTextNumpy(file2.getAbsolutePath(), projectedData.matrix);

            runUmapAndSaveMatricesByChrom(outputDirectory);
        }

        /*
        projectedData = MatrixImputer.imputeUntilNoNansOnlyNN(gwCleanMatrix);
        //corrData = SimilarityMatrixTools.getCosinePearsonCorrMatrix(gwCleanMatrix, 50, generator.nextLong());
        */
    }

    protected int[][] getGenomeIndices() {
        int n = gwCleanMatrix.matrix.length;
        int[][] coordinates = new int[n][3];
        for (int i = 0; i < n; i++) {
            SubcompartmentInterval interval = rowIndexToIntervalMap.get(i);
            coordinates[i][0] = interval.getChrIndex();
            coordinates[i][1] = interval.getX1();
            coordinates[i][2] = interval.getX2();
        }
        return coordinates;
    }

    public void inPlaceScaleSqrtWeightCol() {
        ZScoreTools.inPlaceScaleSqrtWeightCol(gwCleanMatrix.matrix, gwCleanMatrix.weights);
    }


    public void runUmapAndSaveMatricesByChrom(File outputDirectory) {
        int[] indices = new int[umapProjection.length];
        for (int i = 0; i < indices.length; i++) {
            SubcompartmentInterval interval = rowIndexToIntervalMap.get(i);
            indices[i] = interval.getChrIndex();
        }

        plotUmapProjection(outputDirectory, indices, "chrom");
    }

    public void plotUmapProjection(File outputDirectory, List<List<Integer>> colorList, String stem) {
        int[] colors = new int[umapProjection.length];
        Arrays.fill(colors, -1);
        for (int index = 0; index < colorList.size(); index++) {
            for (Integer val : colorList.get(index)) {
                colors[val] = index;
            }
        }

        plotUmapProjection(outputDirectory, colors, stem);
    }

    public void plotUmapProjection(File outputDirectory, int[] colors, String stem) {
        SimpleScatterPlot plotter = new SimpleScatterPlot(umapProjection);
        File outfile = new File(outputDirectory, "umap_" + stem);
        plotter.plot(colors, outfile.getAbsolutePath());
    }

    public void processKMeansClusteringResult(Cluster[] clusters,
                                              GenomeWideList<SubcompartmentInterval> subcompartments) {

        Set<SubcompartmentInterval> subcompartmentIntervals = new HashSet<>();
        if (MixerGlobals.printVerboseComments) {
            System.out.println("GW Composite data vs clustered into " + clusters.length + " clusters");
        }

        int genomewideCompartmentID = 0;
        for (int z = 0; z < clusters.length; z++) {
            Cluster cluster = clusters[z];
            int currentClusterID = ++genomewideCompartmentID;

            if (MixerGlobals.printVerboseComments) {
                System.out.println("Size of cluster " + currentClusterID + " - " + cluster.getMemberIndexes().length);
            }

            for (int i : cluster.getMemberIndexes()) {
                if (rowIndexToIntervalMap.containsKey(i)) {
                    SubcompartmentInterval interv = rowIndexToIntervalMap.get(i);
                    if (interv != null) {
                        subcompartmentIntervals.add(generateNewSubcompartment(interv, currentClusterID));
                    }
                }
            }
        }

        subcompartments.addAll(new ArrayList<>(subcompartmentIntervals));
        SliceUtils.reSort(subcompartments);
    }

    public double getWCSS(Cluster[] clusters, boolean useCorr, boolean useKMedians) {
        double withinClusterSumOfSquares = 0;

        for (int z = 0; z < clusters.length; z++) {
            Cluster cluster = clusters[z];

            if (cluster.getMemberIndexes().length < MIN_EXPECTED_CLUSTER_SIZE) {
                withinClusterSumOfSquares += Float.MAX_VALUE;
            }

            for (int i : cluster.getMemberIndexes()) {
                if (useCorr) {
                    withinClusterSumOfSquares += getDistance(cluster.getCenter(), projectedData.matrix[i], useKMedians);
                } else {
                    withinClusterSumOfSquares += getDistance(cluster.getCenter(), gwCleanMatrix.matrix[i], useKMedians);
                }
            }
        }

        withinClusterSumOfSquares = withinClusterSumOfSquares / clusters.length;
        if (MixerGlobals.printVerboseComments) {
            System.out.println("Final WCSS " + withinClusterSumOfSquares);
        }

        return withinClusterSumOfSquares;
    }

    public double getSilhouette(Cluster[] clusters, boolean useCorr, boolean useKMedians) {
        double score = 0;
        SubsamplingManager manager;
        if (useCorr) {
            manager = new SubsamplingManager(clusters, projectedData.matrix, useKMedians);
        } else {
            manager = new SubsamplingManager(clusters, gwCleanMatrix.matrix, useKMedians);
        }
        for (int iter = 0; iter < NUM_ITERS; iter++) {
            score += manager.getScore();
        }
        score = score / NUM_ITERS;
        System.out.println("Silhouette: " + score);
        return score;
    }

    protected double getDistance(float[] center, float[] vector, boolean useKMedians) {
        if (useKMedians) {
            return RobustManhattanDistance.SINGLETON.distance(center, vector);
        }
        return RobustEuclideanDistance.getNonNanMeanSquaredError(center, vector);
    }

    public synchronized void processGMMClusteringResult(int[] clusterID,
                                                        GenomeWideList<SubcompartmentInterval> subcompartments) {

        Set<SubcompartmentInterval> subcompartmentIntervals = new HashSet<>();

        for (int i = 0; i < projectedData.matrix.length; i++) {
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
        SubcompartmentInterval newInterv = (SubcompartmentInterval) interv.deepClone();
        newInterv.setClusterID(currentClusterID);
        return newInterv;
    }

    public float[][] getData(boolean getCorrMatrix) {
        if (getCorrMatrix) {
            return projectedData.matrix;
        }
        return gwCleanMatrix.matrix;
    }

    public void appendDataAlongExistingRows(CompositeGenomeWideMatrix additionalData) {
        if (gwCleanMatrix.matrix.length != additionalData.gwCleanMatrix.matrix.length) {
            System.err.println("***************************************\n" +
                    "Dimension mismatch: " + gwCleanMatrix.matrix.length + " != " + additionalData.gwCleanMatrix.matrix.length);
        } else {
            gwCleanMatrix = concatenate(gwCleanMatrix, additionalData.gwCleanMatrix);
        }
    }

    private MatrixAndWeight concatenate(MatrixAndWeight mw1, MatrixAndWeight mw2) {
        return new MatrixAndWeight(
                FloatMatrixTools.concatenate(mw1.matrix, mw2.matrix),
                ArrayTools.concatenate(mw1.weights, mw2.weights)
        );
    }

    public GWBadIndexFinder getBadIndices() {
        return badIndexLocations;
    }

    public Map<Integer, SubcompartmentInterval> getRowIndexToIntervalMap() {
        return rowIndexToIntervalMap;
    }
}
