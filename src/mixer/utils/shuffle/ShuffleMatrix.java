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
import java.util.*;

public class ShuffleMatrix {

    private final Random generator = new Random(0);
    private final Dataset ds;
    private final NormalizationType norm;
    private final int compressionFactor;
    private final int resolution;
    private final int numRounds = 50;
    private final InterOnlyMatrix.InterMapType[] mapTypes = {InterOnlyMatrix.InterMapType.ODDS_VS_EVENS,
            InterOnlyMatrix.InterMapType.SKIP_BY_TWOS, InterOnlyMatrix.InterMapType.FIRST_HALF_VS_SECOND_HALF};
    private final String[] scoreTypes = {"Derivative", "Variation", "Kernel (Uniform)", "Kernel (Gaussian)", "Baseline Shift"};
    private final double[][] baselines = new double[mapTypes.length][scoreTypes.length]; //[numRounds]
    private final double[][] shuffled = new double[mapTypes.length][scoreTypes.length];
    private final double[][] ratios = new double[mapTypes.length][scoreTypes.length];

    public ShuffleMatrix(Dataset ds, NormalizationType norm, int resolution, int compressionFactor) {
        this.resolution = resolution;
        this.compressionFactor = compressionFactor;
        this.ds = ds;
        this.norm = norm;
    }

    public void runAnalysis(GenomeWideList<SubcompartmentInterval> subcompartments, File outfolder) {
        SliceUtils.collapseGWList(subcompartments);

        // todo, using only big size?
        // todo sorting picture
        //GenomeWideStatistics statistics = new GenomeWideStatistics(ds,resolution,norm,defaultSubcompartments);
        //float[][] interactionMap = statistics.getInteractionMap();
        //GenomeWideList<SubcompartmentInterval> subcompartments = statistics.reorder(defaultSubcompartments);

        for (int y = 0; y < mapTypes.length; y++) {
            final InterOnlyMatrix interMatrix = new InterOnlyMatrix(ds, norm, resolution, mapTypes[y]);
            Map<Integer, List<Integer>> clusterToRowIndices = populateCluster(interMatrix.getRowChromosomes(), interMatrix.getRowOffsets(), subcompartments);
            Map<Integer, List<Integer>> clusterToColIndices = populateCluster(interMatrix.getColChromosomes(), interMatrix.getColOffsets(), subcompartments);

            float[][] baselineM = shuffleMap(interMatrix, clusterToRowIndices, clusterToColIndices,
                    true, baselines[y]);
            float[][] shuffleM = shuffleMap(interMatrix, clusterToRowIndices, clusterToColIndices,
                    false, shuffled[y]);

            calculateScores(baselineM, shuffleM, y, outfolder);
        }
        getRatioResult();
    }

    private void calculateScores(float[][] baselineM, float[][] shuffleM, int y, File outfolder) {
        //FloatMatrixTools.saveMatrixTextNumpy(new File(outfolder, mapTypes[y].toString() + "_baseline.npy").getAbsolutePath(), baselineM);
        //FloatMatrixTools.saveMatrixTextNumpy(new File(outfolder, mapTypes[y].toString() + "_shuffle.npy").getAbsolutePath(), shuffleM);
        File baselineFile = new File(outfolder, mapTypes[y].toString() + "_baseline.png");
        File shuffleFile = new File(outfolder, mapTypes[y].toString() + "_shuffle.png");
        File baselineLogFile = new File(outfolder, mapTypes[y].toString() + "_log_baseline.png");
        File shuffleLogFile = new File(outfolder, mapTypes[y].toString() + "_log_shuffle.png");

        FloatMatrixTools.saveMatrixToPNG(baselineFile, baselineM, false);
        FloatMatrixTools.saveMatrixToPNG(shuffleFile, shuffleM, false);
        FloatMatrixTools.saveMatrixToPNG(baselineLogFile, baselineM, true);
        FloatMatrixTools.saveMatrixToPNG(shuffleLogFile, shuffleM, true);
        //FloatMatrixTools.saveMatrixTextNumpy(new File(outfolder, mapTypes[y].toString() + "_matrix.npy").getAbsolutePath(), interMatrix.getMatrix());
        //EntropyCalculations ec = new EntropyCalculations(shuffleFile, baselineFile, shuffleLogFile, baselineLogFile, shuffleM);
    }

    private void getRatioResult() {
        for (int k = 0; k < mapTypes.length; k++) {
            for (int l = 0; l < scoreTypes.length; l++) {
                ratios[k][l] = shuffled[k][l] / baselines[k][l];
            }
        }
    }

    public void savePlotsAndResults(File outfolder) {
        try {
            FileWriter myWriter = new FileWriter(new File(outfolder, "scores.txt"));
            for (int y = 0; y < mapTypes.length; y++) {
                myWriter.write(mapTypes[y].toString() + "------------------\n");
                for (int z = 0; z < scoreTypes.length; z++) {
                    myWriter.write(scoreTypes[z] + " Score\n");
                    myWriter.write("Shuffled  : " + shuffled[y][0] + "\n");
                    myWriter.write("Baseline  : " + baselines[y][1] + "\n");
                    myWriter.write("Ratio     : " + ratios[y][2] + "\n\n");

                }
                myWriter.write("----------------------------------------------------------------\n");
            }
            myWriter.close();
        } catch (Exception ee) {
            System.err.println("Unable to write results to text file");
        }
    }

    private float[][] shuffleMap(InterOnlyMatrix interMatrix,
                                 Map<Integer, List<Integer>> clusterToRowIndices,
                                 Map<Integer, List<Integer>> clusterToColIndices,
                                 boolean randomizeGW, double[] scorings) {
        float[][] aggregate = null;
        double[][] scoresForRound = new double[scoreTypes.length][numRounds];

        for (int k = 0; k < numRounds; k++) {
            Pair<List<Integer>, Integer[]> allRowIndices = getShuffledByClusterIndices(clusterToRowIndices, randomizeGW);
            Pair<List<Integer>, Integer[]> allColIndices = getShuffledByClusterIndices(clusterToColIndices, randomizeGW);
            float[][] matrix = getShuffledMatrix(interMatrix, allRowIndices.getFirst(), allColIndices.getFirst());

            if (aggregate == null) {
                aggregate = matrix;
            } else {
                FloatMatrixTools.addBToA(aggregate, matrix);
            }

            scoresForRound[0][k] = (new DerivativeScoring(matrix, allRowIndices.getSecond(), allColIndices.getSecond())).score();
            scoresForRound[1][k] = (new VarianceScoring(matrix, allRowIndices.getSecond(), allColIndices.getSecond())).score();
            scoresForRound[2][k] = (new UniformKernelScoring(matrix, allRowIndices.getSecond(), allColIndices.getSecond())).score();
            scoresForRound[3][k] = (new GaussianKernelScoring(matrix, allRowIndices.getSecond(), allColIndices.getSecond())).score();
            scoresForRound[4][k] = (new EarthShiftScoring(matrix, allRowIndices.getSecond(), allColIndices.getSecond())).score();
        }

        for (int k = 0; k < scoreTypes.length; k++) {
            scorings[k] = calcMean(scoresForRound[k]);
        }

        FloatMatrixTools.scaleBy(aggregate, 1f / ((float) numRounds));
        return aggregate;
    }

    private Pair<List<Integer>, Integer[]> getShuffledByClusterIndices(Map<Integer, List<Integer>> clusterToIndices, boolean randomizeGW) {
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
        if (randomizeGW) {
            Collections.shuffle(allIndices, generator);
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
