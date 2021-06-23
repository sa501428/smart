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

package mixer.utils.matrix;

import mixer.utils.common.FloatMatrixTools;
import mixer.utils.shuffle.scoring.KLDivergenceScoring;
import mixer.utils.shuffle.scoring.VarianceScoring;

import java.io.File;
import java.io.FileWriter;
import java.io.IOException;

public class ScoreContainer {

    private final String[] scoreTypes = {"Variation", "Relative Entropy"};

    private final double[][] baselines;
    private final double[][] shuffled;
    private final double[][] ratios;

    private final double[][] aggBaselines;
    private final double[][] aggShuffled;
    private final double[][] aggRatios;

    private final double[][] logBaselines;
    private final double[][] logShuffled;
    private final double[][] logRatios;

    private final double[][] aggLogBaselines;
    private final double[][] aggLogShuffled;
    private final double[][] aggLogRatios;

    public ScoreContainer(int numMaps, int numScores) {
        baselines = new double[numMaps][numScores];
        shuffled = new double[numMaps][numScores];
        ratios = new double[numMaps][numScores];

        aggBaselines = new double[numMaps][numScores];
        aggShuffled = new double[numMaps][numScores];
        aggRatios = new double[numMaps][numScores];

        logBaselines = new double[numMaps][numScores];
        logShuffled = new double[numMaps][numScores];
        logRatios = new double[numMaps][numScores];

        aggLogBaselines = new double[numMaps][numScores];
        aggLogShuffled = new double[numMaps][numScores];
        aggLogRatios = new double[numMaps][numScores];
    }

    public static void updateAggMatrixScores(double[] scores, float[][] matrix, Integer[] rowBounds, Integer[] colBounds,
                                             boolean isBaseline) {
        scores[0] = (new VarianceScoring(matrix, rowBounds, colBounds)).score(isBaseline);
        scores[1] = (new KLDivergenceScoring(matrix, rowBounds, colBounds, true)).score(isBaseline);
    }

    public void calculateRatios() {
        for (int y = 0; y < ratios.length; y++) {
            for (int l = 0; l < ratios[0].length; l++) {
                logRatios[y][l] = logShuffled[y][l] / logBaselines[y][l];
                ratios[y][l] = shuffled[y][l] / baselines[y][l];
                aggRatios[y][l] = aggShuffled[y][l] / aggBaselines[y][l];
                aggLogRatios[y][l] = aggLogShuffled[y][l] / aggLogBaselines[y][l];
            }
        }
    }

    public void calculateMeans(double[][] scoresForRound, double[][] logScoresForRound,
                               double[][] scoresBaselineForRound, double[][] logScoresBaselineForRound, int mapIndex) {
        for (int k = 0; k < scoresForRound.length; k++) {
            shuffled[mapIndex][k] = calcMean(scoresForRound[k]);
            baselines[mapIndex][k] = calcMean(scoresBaselineForRound[k]);
            logShuffled[mapIndex][k] = calcMean(logScoresForRound[k]);
            logBaselines[mapIndex][k] = calcMean(logScoresBaselineForRound[k]);
        }
    }

    private double calcMean(double[] data) {
        double sum = 0;
        for (double d : data) {
            sum += d;
        }
        return sum / data.length;
    }

    public void updateAggregateScores(AggregateMatrix aggregate, ShuffledIndices[] globalAllIndices, int mapIndex) {
        float[][] aggMatrix = aggregate.getMatrixCopy();
        updateAggMatrixScores(aggBaselines[mapIndex], aggMatrix, globalAllIndices[0].boundaries, globalAllIndices[1].boundaries, true);
        updateAggMatrixScores(aggShuffled[mapIndex], aggMatrix, globalAllIndices[0].boundaries, globalAllIndices[1].boundaries, false);
        FloatMatrixTools.log(aggMatrix, 1);
        updateAggMatrixScores(aggLogBaselines[mapIndex], aggMatrix, globalAllIndices[0].boundaries, globalAllIndices[1].boundaries, true);
        updateAggMatrixScores(aggLogShuffled[mapIndex], aggMatrix, globalAllIndices[0].boundaries, globalAllIndices[1].boundaries, false);
    }

    public void savePlotsAndResults(File outfolder, String prefix, String[] names) {
        try {
            writeToFile(outfolder, "average_scores_" + prefix + ".txt", shuffled, baselines, ratios, names);
            writeToFile(outfolder, "aggregate_scores_" + prefix + ".txt", aggShuffled, aggBaselines, aggRatios, names);
            writeToFile(outfolder, "average_scores_" + prefix + "_log.txt", logShuffled, logBaselines, logRatios, names);
            writeToFile(outfolder, "aggregate_scores_" + prefix + "_log.txt", aggLogShuffled, aggLogBaselines, aggLogRatios, names);
        } catch (Exception ee) {
            System.err.println("Unable to write results to text file");
        }
    }

    private void writeToFile(File outfolder, String filename, double[][] shuffle, double[][] baseline, double[][] ratio,
                             String[] names) throws IOException {
        FileWriter myWriter = new FileWriter(new File(outfolder, filename));
        for (int y = 0; y < ratio.length; y++) {
            myWriter.write(names[y] + "------------------\n");
            for (int z = 0; z < ratio[y].length; z++) {
                myWriter.write(scoreTypes[z] + " Loss\n");
                myWriter.write("Shuffled  : " + shuffle[y][z] + "\n");
                myWriter.write("Baseline  : " + baseline[y][z] + "\n");
                myWriter.write("Ratio     : " + ratio[y][z] + "\n\n");

            }
            myWriter.write("----------------------------------------------------------------\n");
        }
        myWriter.close();
    }
}
