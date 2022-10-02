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

package mixer.utils.cleaning;

import javastraw.expected.ZScoreArray;
import javastraw.reader.basics.Chromosome;
import javastraw.tools.ParallelizationTools;
import mixer.utils.common.SimpleArray2DTools;
import mixer.utils.common.ZScoreTools;
import mixer.utils.drive.MatrixAndWeight;

import java.util.concurrent.atomic.AtomicInteger;

public class MatrixPreprocessor {

    public static MatrixAndWeight clean(MatrixAndWeight matrix, Chromosome[] chromosomes,
                                        int cutoff, boolean useExp, boolean useZscore) {

        SimpleArray2DTools.setZerosToNan(matrix.matrix);

        matrix.updateWeights(chromosomes);

        SimpleArray2DTools.simpleLogWithCleanup(matrix.matrix, Float.NaN);
        removeHighGlobalThresh(matrix.matrix, 5);
        normalize(matrix.matrix, -cutoff, cutoff);
        if (useExp) {
            SimpleArray2DTools.simpleExpm1(matrix.matrix);
        }
        if (useZscore) {
            ZScoreTools.inPlaceZscorePositivesDownColAndSetZeroToNan(matrix.matrix);
        }

        matrix.removeAllNanRows();
        return matrix;
    }

    private static void removeHighGlobalThresh(float[][] data, int cutoff) {
        ZScoreArray zscores = ZScoreTools.getZscores(data);

        AtomicInteger totalNumFixed = new AtomicInteger();
        AtomicInteger index = new AtomicInteger(0);
        ParallelizationTools.launchParallelizedCode(() -> {
            int i = index.getAndIncrement();
            int numFixed = 0;
            while (i < data.length) {
                for (int j = 0; j < data[i].length; j++) {
                    if (data[i][j] > 0) {
                        if (zscores.getZscore(j, data[i][j]) > cutoff) {
                            data[i][j] = Float.NaN;
                            numFixed++;
                        }
                    }
                }
                i = index.getAndIncrement();
            }
            totalNumFixed.addAndGet(numFixed);
        });
    }

    private static void normalize(float[][] data, int lowCutOff, int highCutOff) {
        ZScoreArray zscores = ZScoreTools.getZscores(data);
        fixToNormalRange(data, zscores, lowCutOff, highCutOff);
    }

    private static void fixToNormalRange(float[][] data, ZScoreArray zscores, int lowCutOff, int highCutOff) {
        AtomicInteger totalNumFixed = new AtomicInteger();
        AtomicInteger index = new AtomicInteger(0);
        ParallelizationTools.launchParallelizedCode(() -> {
            int i = index.getAndIncrement();
            int numFixed = 0;
            while (i < data.length) {
                for (int j = 0; j < data[i].length; j++) {
                    if (data[i][j] > 0) {
                        double zscore = zscores.getZscore(j, data[i][j]);
                        if (zscore < lowCutOff || zscore > highCutOff) { //
                            data[i][j] = Float.NaN;
                            numFixed++;
                        }
                    }
                }
                i = index.getAndIncrement();
            }
            totalNumFixed.addAndGet(numFixed);
        });
    }
}
