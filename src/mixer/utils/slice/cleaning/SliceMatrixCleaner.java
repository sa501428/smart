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
import mixer.algos.Slice;
import mixer.utils.common.LogTools;
import mixer.utils.common.ParallelizedStatTools;
import mixer.utils.common.ZScoreTools;
import mixer.utils.slice.cleaning.utils.RowCleaner;
import mixer.utils.slice.matrices.MatrixAndWeight;
import mixer.utils.slice.structures.SubcompartmentInterval;

import java.io.File;
import java.util.Map;
import java.util.Random;

public class SliceMatrixCleaner {
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
    }

    /*
    public float[][] justRemoveBadRows(Set<Integer> badIndices, Map<Integer, SubcompartmentInterval> rowIndexToIntervalMap, int[] weights) {
        return filterOutColumnsAndRowsGivenBadIndices(badIndices, data, rowIndexToIntervalMap);
    }
    */

    public MatrixAndWeight getCleanFilteredZscoredMatrix(Map<Integer, SubcompartmentInterval> rowIndexToIntervalMap,
                                                         int[] weights) {
        if (Slice.FILTER_OUTLIERS) {
            ParallelizedStatTools.setZerosToNan(data);
            ParallelizedStatTools.scaleDown(data, weights);
            LogTools.simpleLogWithCleanup(data, Float.NaN);
            removeHighGlobalThresh(data, weights, 5, Slice.USE_WEIGHTED_MEAN);
            renormalize(data, weights, -2, 2, Slice.USE_WEIGHTED_MEAN);
            LogTools.simpleExpm1(data);
        }

        if (MixerGlobals.printVerboseComments) {
            System.out.println("Initial matrix size " + data.length + " x " + data[0].length);
        }
        //MatrixAndWeight mw = (new ColumnCleaner(data, weights)).getCleanedData();
        //data = mw.matrix;
        //weights = mw.weights;
        //System.out.println("Matrix size after column cleanup " + mw.matrix.length + " x " + mw.matrix[0].length);

        if (true || MixerGlobals.printVerboseComments) {
            System.out.println("Matrix size before row cleanup " + data.length + " x " + data[0].length);
        }
        data = (new RowCleaner(data, rowIndexToIntervalMap, weights)).getCleanedData(resolution, outputDirectory).matrix;
        if (true || MixerGlobals.printVerboseComments) {
            System.out.println("Matrix size after row cleanup " + data.length + " x " + data[0].length);
        }

        ZScoreTools.inPlaceZscoreDownCol(data);

        return new MatrixAndWeight(data, weights);
    }

    private void renormalize(float[][] data, int[] weights, int lowCutOff, int highCutOff, boolean useWeights) {
        double[] muAndStd = ParallelizedStatTools.getMeanAndStandardDev(data, weights, useWeights);
        if (MixerGlobals.printVerboseComments) {
            System.out.println("mu " + muAndStd[0] + " std" + muAndStd[1]);
        }
        fixToNormalRange(data, muAndStd[0], muAndStd[1], lowCutOff, highCutOff);
    }

    private void fixToNormalRange(float[][] data, double mu, double std, int lowCutOff, int highCutOff) {
        int numFixed = 0;
        for (int i = 0; i < data.length; i++) {
            for (int j = 0; j < data[i].length; j++) {
                if (!Float.isNaN(data[i][j])) {
                    double zscore = (data[i][j] - mu) / std;
                    if (zscore < lowCutOff || zscore > highCutOff) { //
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

    private void removeHighGlobalThresh(float[][] data, int[] weights, int cutoff, boolean useWeights) {
        double[] muAndStd = ParallelizedStatTools.getMeanAndStandardDev(data, weights, useWeights);
        if (MixerGlobals.printVerboseComments) {
            System.out.println("mu " + muAndStd[0] + " std" + muAndStd[1]);
        }
        thresholdByMax(data, muAndStd[0], muAndStd[1], cutoff);
    }

    private void thresholdByMax(float[][] data, double mu, double std, int maxZscore) {
        int numFixed = 0;
        for (int i = 0; i < data.length; i++) {
            for (int j = 0; j < data[i].length; j++) {
                if (!Float.isNaN(data[i][j]) && data[i][j] > 0) {
                    double zscore = (data[i][j] - mu) / std;
                    if (zscore > maxZscore) {
                        data[i][j] = Float.NaN;
                        numFixed++;
                    }
                }
            }
        }
        if (true || MixerGlobals.printVerboseComments) {
            System.out.println("Num fixed z > " + maxZscore + " : " + numFixed);
        }
    }
}
