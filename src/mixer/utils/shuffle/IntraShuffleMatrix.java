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

public class IntraShuffleMatrix {

    private final Random generator = new Random(0);
    private final Dataset ds;
    private final NormalizationType norm;
    private final int compressionFactor;
    private final int resolution;
    private final int numRounds = 1;
    private final Chromosome[] chromosomes;
    private final String[] scoreTypes = {"Derivative", "Variation", "Kernel (Uniform)", "Kernel (Gaussian)",
            "KL Divergence Matrix||Baseline", "KL Divergence Baseline||Matrix"};
    private final double[][] baselines;
    private final double[][] shuffled;
    private final double[][] ratios;

    private final double[][] logBaselines;
    private final double[][] logShuffled;
    private final double[][] logRatios;
    private final InterOnlyMatrix.INTRA_TYPE intra_type;
    private final SimilarityMetric metric;

    public IntraShuffleMatrix(Dataset ds, NormalizationType norm, int resolution, int compressionFactor,
                              InterOnlyMatrix.INTRA_TYPE intra_type, SimilarityMetric metric) {
        this.resolution = resolution;
        this.compressionFactor = compressionFactor;
        this.ds = ds;
        this.norm = norm;
        this.intra_type = intra_type;
        this.metric = metric;
        chromosomes = ds.getChromosomeHandler().getAutosomalChromosomesArray();

        baselines = new double[chromosomes.length][scoreTypes.length]; //[numRounds]
        shuffled = new double[chromosomes.length][scoreTypes.length];
        ratios = new double[chromosomes.length][scoreTypes.length];

        logBaselines = new double[chromosomes.length][scoreTypes.length];
        logShuffled = new double[chromosomes.length][scoreTypes.length];
        logRatios = new double[chromosomes.length][scoreTypes.length];
    }

    public void runAnalysis(GenomeWideList<SubcompartmentInterval> subcompartments, File outfolder) {
        SliceUtils.collapseGWList(subcompartments);
        // todo, using only big size? todo sorting picture
        GenomeWideStatistics statistics = new GenomeWideStatistics(ds, resolution, norm, subcompartments);
        statistics.saveInteractionMap(outfolder);
        System.out.println("Interaction summary statistics saved");

        //GenomeWideList<SubcompartmentInterval> subcompartments = statistics.reorder(defaultSubcompartments);
        for (int y = 0; y < chromosomes.length; y++) {
            final InterOnlyMatrix interMatrix = new InterOnlyMatrix(ds, norm, resolution, chromosomes[y],
                    intra_type, metric);
            Map<Integer, List<Integer>> clusterToRowIndices = populateCluster(interMatrix.getRowChromosomes(), interMatrix.getRowOffsets(), subcompartments);

            shuffleMap(interMatrix, clusterToRowIndices, true, baselines[y],
                    logBaselines[y], outfolder, chromosomes[y]);
            shuffleMap(interMatrix, clusterToRowIndices, false, shuffled[y],
                    logShuffled[y], outfolder, chromosomes[y]);

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
        for (int y = 0; y < chromosomes.length; y++) {
            myWriter.write(chromosomes[y].toString() + "------------------\n");
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
                            boolean isBaseline, double[] scorings, double[] logScorings,
                            File outfolder, Chromosome mapType) {
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

                    ShuffledIndices allRowIndices = getShuffledByClusterIndices(clusterToRowIndices, isBaseline);
                    float[][] matrix = getShuffledMatrix(interMatrix, allRowIndices.allIndices);

                    aggregate.addBToA(matrix);
                    ShuffleMatrix.updateMatrixScores(scoresForRound, k, matrix, allRowIndices.boundaries, allRowIndices.boundaries, false);
                    FloatMatrixTools.log(matrix, 1);
                    ShuffleMatrix.updateMatrixScores(logScoresForRound, k, matrix, allRowIndices.boundaries, allRowIndices.boundaries, false);

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
        }

        aggregate.scaleForNumberOfRounds(numRounds);
        aggregate.saveToPNG(outfolder, mapType.getName());
        //aggregate.balanceAndSave(outfolder, mapType);
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

    private float[][] getShuffledMatrix(InterOnlyMatrix interMatrix, List<Integer> allRowIndices) {
        int numRows = allRowIndices.size() / compressionFactor;
        int numRowsKept = numRows * compressionFactor;
        float[][] original = interMatrix.getMatrix();

        float[][] result = new float[numRows][numRows];
        for (int i = 0; i < numRowsKept; i++) {
            final int i0 = allRowIndices.get(i);
            for (int j = 0; j < numRowsKept; j++) {
                final int j0 = allRowIndices.get(j);
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
