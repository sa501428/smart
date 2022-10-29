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
import javastraw.expected.Zscore;

public class MatrixTransform {
    public static void zscoreByRows(float[][] matrix, int limit, boolean zscoreWithNeighbors) {
        int offset = 0;
        if (zscoreWithNeighbors) offset = 2;
        for (int i = 0; i < matrix.length; i++) {
            normalizeRegion0(matrix, i, Math.max(0, i - offset), Math.min(i + offset + 1, matrix.length), limit);
        }
    }

    private static void normalizeRegion0(float[][] matrix, int r, int r0, int rF, int limit) {
        Welford welford = new Welford();
        for (int k = r0; k < rF; k++) {
            for (float val : matrix[k]) {
                if (val > 0) {
                    welford.addValue(val);
                }
            }
        }
        if (welford.getCounts() > 2) {
            Zscore zscore = welford.getZscore();
            for (int c = 0; c < matrix[r].length; c++) {
                matrix[r][c] = zscoreRow(zscore, matrix[r][c], limit);
            }
        }
    }

    private static float zscoreRow(Zscore zscore, double val, int limit) {
        float v = (float) zscore.getZscore(val);
        if (v > limit) return limit;
        if (v < -limit) return -limit;
        return v;
    }

    public static void limitInternally(float[][] matrix, float innerLimit) {
        for (int i = 0; i < matrix.length; i++) {
            for (int j = 0; j < matrix[i].length; j++) {
                if (matrix[i][j] > -innerLimit && matrix[i][j] < innerLimit) {
                    matrix[i][j] = Float.NaN;
                }
            }
        }
    }
}
