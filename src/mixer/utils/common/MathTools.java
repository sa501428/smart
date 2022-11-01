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

import javastraw.expected.Welford;
import javastraw.expected.Zscore;

public class MathTools {

    public static void setZerosToNan(float[][] matrix) {
        for (int r = 0; r < matrix.length; r++) {
            for (int c = 0; c < matrix[0].length; c++) {
                if (matrix[r][c] <= 0) {
                    matrix[r][c] = Float.NaN;
                }
            }
        }
    }

    public static void simpleLogWithCleanup(float[][] matrix, float badVal) {
        for (int r = 0; r < matrix.length; r++) {
            for (int c = 0; c < matrix[0].length; c++) {
                float val = (float) Math.log(matrix[r][c] + 1);
                if (val > 0 && !Float.isInfinite(val)) {
                    matrix[r][c] = val;
                } else {
                    matrix[r][c] = badVal;
                }
            }
        }
    }

    public static void simpleExpm1(float[][] matrix) {
        for (int r = 0; r < matrix.length; r++) {
            for (int c = 0; c < matrix[0].length; c++) {
                float val = (float) Math.expm1(matrix[r][c]);
                if (val > 0) {
                    matrix[r][c] = val;
                }
            }
        }
    }


    public static void removeHighGlobalThresh(float[][] matrix, int maxZscoreUpperLimit) {
        Zscore zscore = getZscore(matrix);
        for (int r = 0; r < matrix.length; r++) {
            for (int c = 0; c < matrix[0].length; c++) {
                if (matrix[r][c] > 0 && zscore.getZscore(matrix[r][c]) > maxZscoreUpperLimit) {
                    matrix[r][c] = Float.NaN;
                }
            }
        }
    }

    public static void regularize(float[][] matrix, int zscoreLimit) {
        Zscore zscore = getZscore(matrix);
        for (int r = 0; r < matrix.length; r++) {
            for (int c = 0; c < matrix[0].length; c++) {
                if (matrix[r][c] > 0 && Math.abs(zscore.getZscore(matrix[r][c])) > zscoreLimit) {
                    matrix[r][c] = Float.NaN;
                }
            }
        }
    }

    private static Zscore getZscore(float[][] matrix) {
        Welford welford = new Welford();
        for (int r = 0; r < matrix.length; r++) {
            for (int c = 0; c < matrix[0].length; c++) {
                if (matrix[r][c] > 0) {
                    welford.addValue(matrix[r][c]);
                }
            }
        }
        return welford.getZscore();
    }
}
