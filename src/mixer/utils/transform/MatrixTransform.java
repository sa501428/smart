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

package mixer.utils.transform;

import javastraw.expected.Welford;
import javastraw.expected.WelfordArray;
import javastraw.expected.ZScoreArray;
import javastraw.expected.Zscore;

public class MatrixTransform {
    public static void zscoreByRows(float[][] matrix, int limit) {
        for (int i = 0; i < matrix.length; i++) {
            normalizeRow(matrix[i], limit);
        }
    }

    private static void normalizeRow(float[] row, int limit) {
        Welford welford = new Welford();
        for (float val : row) {
            if (val > 0) {
                welford.addValue(val);
            }
        }
        if (welford.getCounts() > 2) {
            Zscore zscore = welford.getZscore();
            for (int c = 0; c < row.length; c++) {
                row[c] = thresholdZscore(zscore, row[c], limit);
            }
        }
    }

    public static void zscoreByCols(float[][] matrix, int limit) {
        WelfordArray welfords = new WelfordArray(matrix[0].length);
        for (int i = 0; i < matrix.length; i++) {
            for (int j = 0; j < matrix[i].length; j++) {
                if (matrix[i][j] > 0) {
                    welfords.addValue(j, matrix[i][j]);
                }
            }
        }
        ZScoreArray zscores = welfords.getZscores();
        for (int i = 0; i < matrix.length; i++) {
            for (int j = 0; j < matrix[i].length; j++) {
                if (matrix[i][j] > 0) {
                    matrix[i][j] = thresholdZscores(zscores, j, matrix[i][j], limit);
                }
            }
        }
    }

    private static float thresholdZscores(ZScoreArray zscores, int index, double val, int limit) {
        return threshold(zscores.getZscore(index, val), limit);
    }

    private static float thresholdZscore(Zscore zscore, double val, int limit) {
        return threshold(zscore.getZscore(val), limit);
    }

    private static float threshold(double v, int limit) {
        if (v > limit) return limit;
        if (v < -limit) return -limit;
        return (float) v;
    }
}
