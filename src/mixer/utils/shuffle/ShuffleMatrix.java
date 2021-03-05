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

package mixer.utils.shuffle;

import javastraw.featurelist.GenomeWideList;
import javastraw.reader.Dataset;
import javastraw.reader.basics.Chromosome;
import javastraw.type.NormalizationType;
import mixer.utils.common.FloatMatrixTools;
import mixer.utils.shuffle.scoring.*;
import mixer.utils.similaritymeasures.SimilarityMetric;
import mixer.utils.slice.structures.SliceUtils;
import mixer.utils.slice.structures.SubcompartmentInterval;

import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.util.*;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;
import java.util.concurrent.atomic.AtomicInteger;

public class ShuffleMatrix {

    private final Random generator = new Random(0);
    private final Dataset ds;
    private final NormalizationType norm;
    private final int compressionFactor;
    private final int resolution;
    private final int numRounds = 50;
    private final InterOnlyMatrix.InterMapType[] mapTypes = {InterOnlyMatrix.InterMapType.ODDS_VS_EVENS,
            InterOnlyMatrix.InterMapType.SKIP_BY_TWOS, InterOnlyMatrix.InterMapType.FIRST_HALF_VS_SECOND_HALF};
    private final String[] scoreTypes = {"Derivative", "Variation", "Kernel (Uniform)", "Kernel (Gaussian)",
            "KL Divergence Matrix||Baseline", "KL Divergence Baseline||Matrix"};
    private final double[][] baselines = new double[mapTypes.length][scoreTypes.length]; //[numRounds]
    private final double[][] shuffled = new double[mapTypes.length][scoreTypes.length];
    private final double[][] ratios = new double[mapTypes.length][scoreTypes.length];

    private final double[][] aggBaselines = new double[mapTypes.length][scoreTypes.length]; //[numRounds]
    private final double[][] aggShuffled = new double[mapTypes.length][scoreTypes.length];
    private final double[][] aggRatios = new double[mapTypes.length][scoreTypes.length];


    private final double[][] logBaselines = new double[mapTypes.length][scoreTypes.length];
    private final double[][] logShuffled = new double[mapTypes.length][scoreTypes.length];
    private final double[][] logRatios = new double[mapTypes.length][scoreTypes.length];

    private final double[][] aggLogBaselines = new double[mapTypes.length][scoreTypes.length];
    private final double[][] aggLogShuffled = new double[mapTypes.length][scoreTypes.length];
    private final double[][] aggLogRatios = new double[mapTypes.length][scoreTypes.length];

    private final SimilarityMetric metric;

    public ShuffleMatrix(Dataset ds, NormalizationType norm, int resolution, int compressionFactor,
                         SimilarityMetric metric) {
        this.resolution = resolution;
        this.compressionFactor = compressionFactor;
        this.ds = ds;
        this.norm = norm;
        this.metric = metric;
    }

    public static void updateMatrixScores(double[][] scores, int k, float[][] matrix, Integer[] rowBounds, Integer[] colBounds,
                                          boolean isBaseline) {
        scores[0][k] = (new DerivativeScoring(matrix, rowBounds, colBounds)).score(isBaseline);
        scores[1][k] = (new VarianceScoring(matrix, rowBounds, colBounds)).score(isBaseline);
        scores[2][k] = (new UniformKernelScoring(matrix, rowBounds, colBounds)).score(isBaseline);
        scores[3][k] = (new GaussianKernelScoring(matrix, rowBounds, colBounds)).score(isBaseline);
        scores[4][k] = (new KLDivergenceScoring(matrix, rowBounds, colBounds, true)).score(isBaseline);
        scores[5][k] = (new KLDivergenceScoring(matrix, rowBounds, colBounds, false)).score(isBaseline);
    }

    public static void updateAggMatrixScores(double[] scores, float[][] matrix, Integer[] rowBounds, Integer[] colBounds,
                                             boolean isBaseline) {
        scores[0] = (new DerivativeScoring(matrix, rowBounds, colBounds)).score(isBaseline);
        scores[1] = (new VarianceScoring(matrix, rowBounds, colBounds)).score(isBaseline);
        scores[2] = (new UniformKernelScoring(matrix, rowBounds, colBounds)).score(isBaseline);
        scores[3] = (new GaussianKernelScoring(matrix, rowBounds, colBounds)).score(isBaseline);
        scores[4] = (new KLDivergenceScoring(matrix, rowBounds, colBounds, true)).score(isBaseline);
        scores[5] = (new KLDivergenceScoring(matrix, rowBounds, colBounds, false)).score(isBaseline);
    }

    public void runAnalysis(GenomeWideList<SubcompartmentInterval> subcompartments, File outfolder) {
        SliceUtils.collapseGWList(subcompartments);
        // todo, using only big size? todo sorting picture
        GenomeWideStatistics statistics = new GenomeWideStatistics(ds, resolution, norm, subcompartments);
        statistics.saveInteractionMap(outfolder);
        System.out.println("Interaction summary statistics saved");

        //GenomeWideList<SubcompartmentInterval> subcompartments = statistics.reorder(defaultSubcompartments);
        for (int y = 0; y < mapTypes.length; y++) {
            final InterOnlyMatrix interMatrix = new InterOnlyMatrix(ds, norm, resolution, mapTypes[y], metric);
            Map<Integer, List<Integer>> clusterToRowIndices = populateCluster(interMatrix.getRowChromosomes(), interMatrix.getRowOffsets(), subcompartments);
            Map<Integer, List<Integer>> clusterToColIndices = populateCluster(interMatrix.getColChromosomes(), interMatrix.getColOffsets(), subcompartments);

            shuffleMap(interMatrix, clusterToRowIndices, clusterToColIndices,
                    baselines[y], logBaselines[y], aggBaselines[y], aggLogBaselines[y],
                    shuffled[y], logShuffled[y], aggShuffled[y], aggLogShuffled[y],
                    outfolder, mapTypes[y]);

            for (int l = 0; l < scoreTypes.length; l++) {
                ratios[y][l] = shuffled[y][l] / baselines[y][l];
                logRatios[y][l] = logShuffled[y][l] / logBaselines[y][l];
                aggRatios[y][l] = aggShuffled[y][l] / aggBaselines[y][l];
                aggLogRatios[y][l] = aggLogShuffled[y][l] / aggLogBaselines[y][l];
            }
        }
    }

    public void savePlotsAndResults(File outfolder, String prefix) {
        try {
            writeToFile(outfolder, "average_scores_" + prefix + ".txt", shuffled, baselines, ratios);
            writeToFile(outfolder, "average_scores_" + prefix + "_log.txt", logShuffled, logBaselines, logRatios);
            writeToFile(outfolder, "aggregate_scores_" + prefix + ".txt", aggShuffled, aggBaselines, aggRatios);
            writeToFile(outfolder, "aggregate_scores_" + prefix + "_log.txt", aggLogShuffled, aggLogBaselines, aggLogRatios);
        } catch (Exception ee) {
            System.err.println("Unable to write results to text file");
        }
    }

    private void writeToFile(File outfolder, String filename, double[][] shuffle, double[][] baseline, double[][] ratio) throws IOException {
        FileWriter myWriter = new FileWriter(new File(outfolder, filename));
        for (int y = 0; y < mapTypes.length; y++) {
            myWriter.write(mapTypes[y].toString() + "------------------\n");
            for (int z = 0; z < scoreTypes.length; z++) {
                myWriter.write(scoreTypes[z] + " Loss\n");
                myWriter.write("Shuffled  : " + shuffle[y][z] + "\n");
                myWriter.write("Baseline  : " + baseline[y][z] + "\n");
                myWriter.write("Ratio     : " + ratio[y][z] + "\n\n");

            }
            myWriter.write("----------------------------------------------------------------\n");
        }
        myWriter.close();
    }

    private void shuffleMap(InterOnlyMatrix interMatrix, Map<Integer, List<Integer>> clusterToRowIndices,
                            Map<Integer, List<Integer>> clusterToColIndices,
                            double[] scoringsBaseline, double[] logScoringsBaseline,
                            double[] aggScoringsBaseline, double[] aggLogScoringsBaseline,
                            double[] scorings, double[] logScorings,
                            double[] aggScorings, double[] aggLogScorings,

                            File outfolder, InterOnlyMatrix.InterMapType mapType) {
        boolean isBaseline = false;

        final AggregateMatrix aggregate = new AggregateMatrix(isBaseline);
        double[][] scoresBaselineForRound = new double[scoreTypes.length][numRounds];
        double[][] logScoresBaselineForRound = new double[scoreTypes.length][numRounds];
        double[][] scoresForRound = new double[scoreTypes.length][numRounds];
        double[][] logScoresForRound = new double[scoreTypes.length][numRounds];

        int numCPUThreads = Runtime.getRuntime().availableProcessors();
        AtomicInteger currRowIndex = new AtomicInteger(0);
        ExecutorService executor = Executors.newFixedThreadPool(numCPUThreads);

        final ShuffledIndices[] globalAllIndices = new ShuffledIndices[2];
        for (int l = 0; l < numCPUThreads; l++) {
            executor.execute(() -> {
                int k = currRowIndex.getAndIncrement();
                while (k < numRounds) {

                    ShuffledIndices allRowIndices = getShuffledByClusterIndices(clusterToRowIndices, isBaseline);
                    ShuffledIndices allColIndices = getShuffledByClusterIndices(clusterToColIndices, isBaseline);
                    if (k == 0) {
                        globalAllIndices[0] = getShuffledByClusterIndices(clusterToRowIndices, isBaseline);
                        globalAllIndices[1] = getShuffledByClusterIndices(clusterToColIndices, isBaseline);
                    }

                    float[][] matrix = getShuffledMatrix(interMatrix, allRowIndices.allIndices, allColIndices.allIndices);

                    updateMatrixScores(scoresBaselineForRound, k, matrix, allRowIndices.boundaries, allColIndices.boundaries, true);
                    updateMatrixScores(scoresForRound, k, matrix, allRowIndices.boundaries, allColIndices.boundaries, false);

                    FloatMatrixTools.log(matrix, 1);
                    updateMatrixScores(logScoresBaselineForRound, k, matrix, allRowIndices.boundaries, allColIndices.boundaries, true);
                    updateMatrixScores(logScoresForRound, k, matrix, allRowIndices.boundaries, allColIndices.boundaries, false);

                    aggregate.addBToA(matrix);
                    k = currRowIndex.getAndIncrement();
                }
            });
        }
        executor.shutdown();
        //noinspection StatementWithEmptyBody
        while (!executor.isTerminated()) {
        }

        for (int k = 0; k < scoreTypes.length; k++) {
            scorings[k] = calcMean(scoresForRound[k]);
            logScorings[k] = calcMean(logScoresForRound[k]);
            scoringsBaseline[k] = calcMean(scoresBaselineForRound[k]);
            logScoringsBaseline[k] = calcMean(logScoresBaselineForRound[k]);
        }

        aggregate.scaleForNumberOfRounds(numRounds);
        aggregate.saveToPNG(outfolder, mapType.toString());
        //aggregate.balanceAndSave(outfolder, mapType);

        float[][] aggMatrix = aggregate.getMatrixCopy();

        updateAggMatrixScores(aggScoringsBaseline, aggMatrix, globalAllIndices[0].boundaries, globalAllIndices[1].boundaries, true);
        updateAggMatrixScores(aggScorings, aggMatrix, globalAllIndices[0].boundaries, globalAllIndices[1].boundaries, false);
        FloatMatrixTools.log(aggMatrix, 1);
        updateAggMatrixScores(aggLogScoringsBaseline, aggMatrix, globalAllIndices[0].boundaries, globalAllIndices[1].boundaries, true);
        updateAggMatrixScores(aggLogScorings, aggMatrix, globalAllIndices[0].boundaries, globalAllIndices[1].boundaries, false);
    }

    private ShuffledIndices getShuffledByClusterIndices(Map<Integer, List<Integer>> clusterToIndices, boolean isBaseline) {
        List<Integer> allIndices = new ArrayList<>();

        List<Integer> order = new ArrayList<>(clusterToIndices.keySet());
        Collections.sort(order);

        int count = 0;
        List<Integer> boundaries = new ArrayList<>();
        boundaries.add(count);

        for (Integer index : order) {
            List<Integer> indexList = clusterToIndices.get(index);
            Collections.shuffle(indexList, generator);
            int numToUse = (indexList.size() / compressionFactor) * compressionFactor;
            for (int z = 0; z < numToUse; z++) {
                allIndices.add(indexList.get(z));
            }
            count += (numToUse / compressionFactor);
            boundaries.add(count);
        }
        if (isBaseline) {
            Collections.shuffle(allIndices, generator);
            return new ShuffledIndices(allIndices, new Integer[]{0, count});
        }
        Integer[] output = new Integer[boundaries.size()];
        return new ShuffledIndices(allIndices, boundaries.toArray(output));
    }

    private float[][] getShuffledMatrix(InterOnlyMatrix interMatrix, List<Integer> allRowIndices, List<Integer> allColIndices) {
        int numRows = allRowIndices.size() / compressionFactor;
        int numCols = allColIndices.size() / compressionFactor;
        int numRowsKept = numRows * compressionFactor;
        int numColsKept = numCols * compressionFactor;
        float[][] original = interMatrix.getMatrix();

        float[][] result = new float[numRows][numCols];
        for (int i = 0; i < numRowsKept; i++) {
            final int i0 = allRowIndices.get(i);
            for (int j = 0; j < numColsKept; j++) {
                final int j0 = allColIndices.get(j);
                result[i / compressionFactor][j / compressionFactor] += original[i0][j0];
            }
        }

        return result;
    }

    private Map<Integer, List<Integer>> populateCluster(Chromosome[] chromosomes, int[] offsets, GenomeWideList<SubcompartmentInterval> subcompartments) {
        Map<Integer, List<Integer>> clusterToIndices = new HashMap<>();

        for (int x = 0; x < chromosomes.length; x++) {
            Chromosome chrom = chromosomes[x];
            List<SubcompartmentInterval> intervalList = subcompartments.getFeatures("" + chrom.getIndex());
            for (SubcompartmentInterval interval : intervalList) {
                int xStart = interval.getX1() / resolution;
                int xEnd = interval.getX2() / resolution;
                int clusterID = interval.getClusterID();

                List<Integer> tempList = new ArrayList<>();
                for (int k = xStart; k < xEnd; k++) {
                    final int actualPosition = k + offsets[x];
                    tempList.add(actualPosition);
                }

                if (clusterToIndices.containsKey(clusterID)) {
                    clusterToIndices.get(clusterID).addAll(tempList);
                } else {
                    clusterToIndices.put(clusterID, tempList);
                }
            }
        }
        return clusterToIndices;
    }

    private double calcMean(double[] data) {
        double sum = 0;
        for (double d : data) {
            sum += d;
        }
        return sum / data.length;
    }
}
