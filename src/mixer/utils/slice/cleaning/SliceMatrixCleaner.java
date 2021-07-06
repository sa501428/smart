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

package mixer.utils.slice.cleaning;

import mixer.MixerGlobals;
import mixer.utils.common.LogTools;
import mixer.utils.common.ZScoreTools;
import mixer.utils.slice.cleaning.utils.RowCleaner;
import mixer.utils.slice.matrices.MatrixAndWeight;
import mixer.utils.slice.structures.SubcompartmentInterval;

import java.io.File;
import java.util.Map;
import java.util.Random;

public class SliceMatrixCleaner {

    private static final int MAX_ZSCORE = 5;
    private static final int MAX_NORMAL_ZSCORE = 2;
    public static int NUM_PER_CENTROID = 100;
    protected final File outputDirectory;
    protected float[][] data;
    protected final Random generator = new Random(0);
    protected int resolution;

    public SliceMatrixCleaner(float[][] data, long seed, File outputDirectory, int resolution) {
        this.outputDirectory = outputDirectory;
        generator.setSeed(seed);
        this.resolution = resolution;
        this.data = data;
        System.out.println("Initial matrix size " + data.length + " x " + data[0].length);
    }

    public static void setZerosToNan(float[][] matrix) {
        for (int i = 0; i < matrix.length; i++) {
            for (int j = 0; j < matrix[i].length; j++) {
                if (matrix[i][j] < 1e-20) {
                    matrix[i][j] = Float.NaN;
                }
            }
        }
    }

    /*
    public float[][] justRemoveBadRows(Set<Integer> badIndices, Map<Integer, SubcompartmentInterval> rowIndexToIntervalMap, int[] weights) {
        return filterOutColumnsAndRowsGivenBadIndices(badIndices, data, rowIndexToIntervalMap);
    }
    */

    public MatrixAndWeight getCleanFilteredZscoredMatrix(Map<Integer, SubcompartmentInterval> rowIndexToIntervalMap,
                                                         int[] weights) {
        setZerosToNan(data);
        scaleDown(data, weights);
        LogTools.simpleLogWithCleanup(data, Float.NaN);
        removeHighGlobalThresh(data, weights);
        renormalize(data, weights);

        if (MixerGlobals.printVerboseComments) {
            System.out.println("Initial matrix size " + data.length + " x " + data[0].length);
        }
        //MatrixAndWeight mw = (new ColumnCleaner(data, weights)).getCleanedData();
        //data = mw.matrix;
        //weights = mw.weights;
        //System.out.println("Matrix size after column cleanup " + mw.matrix.length + " x " + mw.matrix[0].length);

        data = (new RowCleaner(data, rowIndexToIntervalMap, weights)).getCleanedData(resolution, outputDirectory).matrix;
        if (MixerGlobals.printVerboseComments) {
            System.out.println("Matrix size after row cleanup " + data.length + " x " + data[0].length);
        }

        ZScoreTools.inPlaceZscoreDownCol(data);
        ZScoreTools.inPlaceScaleSqrtWeightCol(data, weights);

        return new MatrixAndWeight(data, weights);
    }

    private void renormalize(float[][] data, int[] weights) {
        double mu = getGlobalMean(data, weights);
        double std = getGlobalStdDev(data, weights, mu);
        if (MixerGlobals.printVerboseComments) {
            System.out.println("mu " + mu + " std" + std);
        }
        fixToNormalRange(data, mu, std);
    }

    private void fixToNormalRange(float[][] data, double mu, double std) {
        int numFixed = 0;
        for (int i = 0; i < data.length; i++) {
            for (int j = 0; j < data[i].length; j++) {
                if (!Float.isNaN(data[i][j])) {
                    double zscore = (data[i][j] - mu) / std;
                    if (zscore > SliceMatrixCleaner.MAX_NORMAL_ZSCORE ||
                            zscore < -SliceMatrixCleaner.MAX_NORMAL_ZSCORE) { //
                        data[i][j] = Float.NaN;
                        numFixed++;
                    }
                }
            }
        }
        if (MixerGlobals.printVerboseComments) {
            System.out.println("Num fixed part 2: z < -2 : " + numFixed);
        }
    }

    private void removeHighGlobalThresh(float[][] data, int[] weights) {
        double mu = getGlobalMean(data, weights);
        double std = getGlobalStdDev(data, weights, mu);
        if (MixerGlobals.printVerboseComments) {
            System.out.println("mu " + mu + " std" + std);
        }
        thresholdByMax(data, mu, std, MAX_ZSCORE);
    }

    private void thresholdByMax(float[][] data, double mu, double std, int maxZscore) {
        int numFixed = 0;
        for (int i = 0; i < data.length; i++) {
            for (int j = 0; j < data[i].length; j++) {
                if (!Float.isNaN(data[i][j])) {
                    double zscore = (data[i][j] - mu) / std;
                    if (zscore > maxZscore) {
                        data[i][j] = Float.NaN;
                        numFixed++;
                    }
                }
            }
        }
        if (MixerGlobals.printVerboseComments) {
            System.out.println("Num fixed z > 5 : " + numFixed);
        }
    }

    private double getGlobalStdDev(float[][] data, int[] weights, double mu) {
        double squares = 0;
        long count = 0;
        for (int i = 0; i < data.length; i++) {
            for (int j = 0; j < data[i].length; j++) {
                if (!Float.isNaN(data[i][j])) {
                    double diff = data[i][j] - mu;
                    squares += diff * diff * weights[j];
                    count += weights[j];
                }
            }
        }
        return Math.sqrt(squares / count);
    }

    private double getGlobalMean(float[][] data, int[] weights) {
        double total = 0;
        long count = 0;
        for (int i = 0; i < data.length; i++) {
            for (int j = 0; j < data[i].length; j++) {
                if (!Float.isNaN(data[i][j])) {
                    total += data[i][j] * weights[j];
                    count += weights[j];
                }
            }
        }
        return total / count;
    }

    private void scaleDown(float[][] data, int[] weights) {
        for (int i = 0; i < data.length; i++) {
            for (int j = 0; j < data[i].length; j++) {
                if (!Float.isNaN(data[i][j])) {
                    data[i][j] /= weights[j];
                }
            }
        }
    }
}
