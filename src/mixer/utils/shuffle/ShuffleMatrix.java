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

package mixer.utils.shuffle;

import javastraw.featurelist.GenomeWideList;
import javastraw.reader.Dataset;
import javastraw.reader.basics.Chromosome;
import javastraw.type.NormalizationType;
import mixer.utils.common.FloatMatrixTools;
import mixer.utils.common.Pair;
import mixer.utils.shuffle.scoring.*;
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

    private final double[][] logBaselines = new double[mapTypes.length][scoreTypes.length];
    private final double[][] logShuffled = new double[mapTypes.length][scoreTypes.length];
    private final double[][] logRatios = new double[mapTypes.length][scoreTypes.length];

    public ShuffleMatrix(Dataset ds, NormalizationType norm, int resolution, int compressionFactor) {
        this.resolution = resolution;
        this.compressionFactor = compressionFactor;
        this.ds = ds;
        this.norm = norm;
    }

    public void runAnalysis(GenomeWideList<SubcompartmentInterval> subcompartments, File outfolder) {
        SliceUtils.collapseGWList(subcompartments);
        // todo, using only big size? todo sorting picture
        //GenomeWideStatistics statistics = new GenomeWideStatistics(ds,resolution,norm,defaultSubcompartments);
        //float[][] interactionMap = statistics.getInteractionMap();
        //GenomeWideList<SubcompartmentInterval> subcompartments = statistics.reorder(defaultSubcompartments);
        for (int y = 0; y < mapTypes.length; y++) {
            final InterOnlyMatrix interMatrix = new InterOnlyMatrix(ds, norm, resolution, mapTypes[y]);
            Map<Integer, List<Integer>> clusterToRowIndices = populateCluster(interMatrix.getRowChromosomes(), interMatrix.getRowOffsets(), subcompartments);
            Map<Integer, List<Integer>> clusterToColIndices = populateCluster(interMatrix.getColChromosomes(), interMatrix.getColOffsets(), subcompartments);

            shuffleMap(interMatrix, clusterToRowIndices, clusterToColIndices, true, baselines[y],
                    logBaselines[y], outfolder, mapTypes[y]);
            shuffleMap(interMatrix, clusterToRowIndices, clusterToColIndices, false, shuffled[y],
                    logShuffled[y], outfolder, mapTypes[y]);

            for (int l = 0; l < scoreTypes.length; l++) {
                ratios[y][l] = shuffled[y][l] / baselines[y][l];
                logRatios[y][l] = logShuffled[y][l] / logBaselines[y][l];
            }
        }
    }

    public void savePlotsAndResults(File outfolder, String prefix) {
        try {
            writeToFile(outfolder, "scores_" + prefix + ".txt", shuffled, baselines, ratios);
            writeToFile(outfolder, "scores_" + prefix + "_log.txt", logShuffled, logBaselines, logRatios);
        } catch (Exception ee) {
            System.err.println("Unable to write results to text file");
        }
    }

    private void writeToFile(File outfolder, String filename, double[][] shuffle, double[][] baseline, double[][] ratio) throws IOException {
        FileWriter myWriter = new FileWriter(new File(outfolder, filename));
        for (int y = 0; y < mapTypes.length; y++) {
            myWriter.write(mapTypes[y].toString() + "------------------\n");
            for (int z = 0; z < scoreTypes.length; z++) {
                myWriter.write(scoreTypes[z] + " Score\n");
                myWriter.write("Shuffled  : " + shuffle[y][z] + "\n");
                myWriter.write("Baseline  : " + baseline[y][z] + "\n");
                myWriter.write("Ratio     : " + ratio[y][z] + "\n\n");

            }
            myWriter.write("----------------------------------------------------------------\n");
        }
        myWriter.close();
    }

    private void shuffleMap(InterOnlyMatrix interMatrix, Map<Integer, List<Integer>> clusterToRowIndices,
                            Map<Integer, List<Integer>> clusterToColIndices, boolean isBaseline,
                            double[] scorings, double[] logScorings, File outfolder, InterOnlyMatrix.InterMapType mapType) {
        final AggregateMatrix aggregate = new AggregateMatrix(isBaseline);
        double[][] scoresForRound = new double[scoreTypes.length][numRounds];
        double[][] logScoresForRound = new double[scoreTypes.length][numRounds];

        int numCPUThreads = Runtime.getRuntime().availableProcessors();
        AtomicInteger currRowIndex = new AtomicInteger(0);
        ExecutorService executor = Executors.newFixedThreadPool(numCPUThreads);
        for (int l = 0; l < numCPUThreads; l++) {
            executor.execute(() -> {
                int k = currRowIndex.getAndIncrement();
                while (k < numRounds) {

                    Pair<List<Integer>, Integer[]> allRowIndices = getShuffledByClusterIndices(clusterToRowIndices, isBaseline);
                    Pair<List<Integer>, Integer[]> allColIndices = getShuffledByClusterIndices(clusterToColIndices, isBaseline);
                    float[][] matrix = getShuffledMatrix(interMatrix, allRowIndices.getFirst(), allColIndices.getFirst());

                    aggregate.addBToA(matrix);
                    updateMatrixScores(scoresForRound, k, matrix, allRowIndices.getSecond(), allColIndices.getSecond());
                    FloatMatrixTools.log(matrix, 1);
                    updateMatrixScores(logScoresForRound, k, matrix, allRowIndices.getSecond(), allColIndices.getSecond());

                    k = currRowIndex.getAndIncrement();
                }
            });
        }
        executor.shutdown();
        while (!executor.isTerminated()) {
        }

        for (int k = 0; k < scoreTypes.length; k++) {
            scorings[k] = calcMean(scoresForRound[k]);
            logScorings[k] = calcMean(logScoresForRound[k]);
        }

        aggregate.scaleForNumberOfRounds(numRounds);
        aggregate.saveToPNG(outfolder, mapType);
    }

    private void updateMatrixScores(double[][] scores, int k, float[][] matrix, Integer[] rowBounds, Integer[] colBounds) {
        scores[0][k] = (new DerivativeScoring(matrix, rowBounds, colBounds)).score();
        scores[1][k] = (new VarianceScoring(matrix, rowBounds, colBounds)).score();
        scores[2][k] = (new UniformKernelScoring(matrix, rowBounds, colBounds)).score();
        scores[3][k] = (new GaussianKernelScoring(matrix, rowBounds, colBounds)).score();
        scores[4][k] = (new KLDivergenceScoring(matrix, rowBounds, colBounds, true)).score();
        scores[5][k] = (new KLDivergenceScoring(matrix, rowBounds, colBounds, false)).score();
    }

    private Pair<List<Integer>, Integer[]> getShuffledByClusterIndices(Map<Integer, List<Integer>> clusterToIndices, boolean isBaseline) {
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
            return new Pair<>(allIndices, new Integer[]{0, count});
        }
        Integer[] output = new Integer[boundaries.size()];
        return new Pair<>(allIndices, boundaries.toArray(output));
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
