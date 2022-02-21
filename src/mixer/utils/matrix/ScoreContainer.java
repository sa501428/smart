/*
 * The MIT License (MIT)
 *
 * Copyright (c) 2011-2022 Rice University, Baylor College of Medicine, Aiden Lab
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

import mixer.utils.shuffle.scoring.KLDivergenceScoring;
import mixer.utils.shuffle.scoring.VarianceScoring;

import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.util.Arrays;

public class ScoreContainer {

    private final String[] scoreTypes = {"Variation", "Relative Entropy"};

    private final double[][] baselines;
    private final double[][] shuffled;
    private final double[][] ratios;

    public ScoreContainer(int numMaps, int numScores) {
        baselines = new double[numMaps][numScores];
        shuffled = new double[numMaps][numScores];
        ratios = new double[numMaps][numScores];
    }

    public static double[] updateAggMatrixScores(float[][] matrix, Integer[] rowBounds, Integer[] colBounds,
                                                 boolean isBaseline) {
        double[] scores = new double[2];
        scores[0] = (new VarianceScoring(matrix, rowBounds, colBounds)).score(isBaseline);
        scores[1] = (new KLDivergenceScoring(matrix, rowBounds, colBounds, true)).score(isBaseline);
        return scores;
    }

    public void calculateRatios() {
        for (int y = 0; y < ratios.length; y++) {
            for (int l = 0; l < ratios[0].length; l++) {
                ratios[y][l] = baselines[y][l] / shuffled[y][l];
            }
        }
    }

    public void updateAggregateScores(AggregateMatrix aggregate, ShuffledIndices[] globalAllIndices, int mapIndex) {
        float[][] aggMatrix = aggregate.getMatrixCopy();
        baselines[mapIndex] = updateAggMatrixScores(aggMatrix, globalAllIndices[0].boundaries, globalAllIndices[1].boundaries, true);
        shuffled[mapIndex] = updateAggMatrixScores(aggMatrix, globalAllIndices[0].boundaries, globalAllIndices[1].boundaries, false);
    }

    public void savePlotsAndResults(File outfolder, String prefix, String[] names) {
        try {
            writeToFile(outfolder, "aggregate_scores_" + prefix + "_log.txt", shuffled, baselines, ratios, names);
        } catch (Exception ee) {
            System.err.println("Unable to write results to text file");
        }
    }

    private void writeToFile(File outfolder, String filename, double[][] shuffle, double[][] baseline, double[][] ratio,
                             String[] names) throws IOException {
        FileWriter myWriter = new FileWriter(new File(outfolder, filename));
        double[] geometricMeans = new double[scoreTypes.length];
        Arrays.fill(geometricMeans, 1);

        for (int y = 0; y < ratio.length; y++) {
            myWriter.write(names[y] + "------------------\n");
            for (int z = 0; z < ratio[y].length; z++) {
                myWriter.write("CHIC Score (" + scoreTypes[z] + ")\n");
                myWriter.write("Shuffled  : " + shuffle[y][z] + "\n");
                myWriter.write("Baseline  : " + baseline[y][z] + "\n");
                myWriter.write("Score     : " + ratio[y][z] + "\n\n");
                geometricMeans[z] *= ratio[y][z];
            }
            myWriter.write("----------------------------------------------------------------\n");
        }

        for (int z = 0; z < scoreTypes.length; z++) {
            geometricMeans[z] = Math.pow(geometricMeans[z], 1.0 / names.length);
        }

        myWriter.write("Geometric Mean of CHIC Scores------------------\n\n");
        for (int z = 0; z < scoreTypes.length; z++) {
            myWriter.write("CHIC Score (" + scoreTypes[z] + "): " + geometricMeans[z] + "\n\n");
        }

        myWriter.close();
    }
}
