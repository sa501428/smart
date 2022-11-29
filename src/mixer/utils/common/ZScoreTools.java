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

package mixer.utils.common;

import javastraw.expected.WelfordArray;
import javastraw.expected.ZScoreArray;
import javastraw.tools.ParallelizationTools;

import java.util.concurrent.atomic.AtomicInteger;

@SuppressWarnings("ForLoopReplaceableByForEach")
public class ZScoreTools {


    public static void inPlaceScaleSqrtWeightCol(float[][] matrix, int[] weights) {
        if (weights.length != matrix[0].length) {
            System.err.println("Weights mismatch error " + weights.length + " vs " + matrix[0].length);
            System.exit(54);
        }

        AtomicInteger index = new AtomicInteger(0);
        ParallelizationTools.launchParallelizedCode(() -> {
            int i = index.getAndIncrement();
            while (i < matrix.length) {
                for (int j = 0; j < matrix[i].length; j++) {
                    matrix[i][j] = (float) (Math.sqrt(weights[j]) * matrix[i][j]);
                }
                i = index.getAndIncrement();
            }
        });
    }

    public static void inPlaceZscorePositivesDownColAndSetBelowThreshToNan(float[][] matrix, int minval) {
        ZScoreArray zscores = getZscores(matrix, minval);

        AtomicInteger index = new AtomicInteger(0);
        ParallelizationTools.launchParallelizedCode(() -> {
            int i = index.getAndIncrement();
            while (i < matrix.length) {
                for (int j = 0; j < matrix[i].length; j++) {
                    float val = matrix[i][j];
                    if (val > minval) {
                        matrix[i][j] = (float) zscores.getZscore(j, val);
                    } else {
                        matrix[i][j] = Float.NaN;
                    }
                }
                i = index.getAndIncrement();
            }
        });
    }

    public static ZScoreArray getZscores(float[][] matrix, int minval) {
        WelfordArray welfords = new WelfordArray(matrix[0].length);

        for (int i = 0; i < matrix.length; i++) {
            for (int j = 0; j < matrix[i].length; j++) {
                if (matrix[i][j] > minval) {
                    welfords.addValue(j, matrix[i][j]);
                }
            }
        }

        return welfords.getZscores();
    }
}
