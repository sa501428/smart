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

package mixer.utils.drive;

public class ScaleTargetVector {
    public static long[] create(int[] weights, int numRows, int numCols, float[][] matrix) {
        long[] colTarget = getSummedWeightsColumns(weights, matrix, numCols);
        long[] rowTarget = getSummedWeightsRows(weights, matrix, numRows);

        long[] target = new long[numCols + numRows];
        System.arraycopy(colTarget, 0, target, 0, numCols);
        System.arraycopy(rowTarget, 0, target, numCols, numRows);
        return target;
    }

    public static long[] getSummedWeightsColumns(int[] weights, float[][] matrix, int numCols) {
        long[] colSums = new long[numCols];
        for (int i = 0; i < matrix.length; i++) {
            for (int j = 0; j < matrix[i].length; j++) {
                if (matrix[i][j] > -1) {
                    colSums[j] += weights[j];
                }
            }
        }
        return colSums;
    }

    public static long[] getSummedWeightsRows(int[] weights, float[][] matrix, int numRows) {
        long[] rowSums = new long[numRows];
        for (int i = 0; i < matrix.length; i++) {
            for (int j = 0; j < matrix[i].length; j++) {
                if (matrix[i][j] > -1) {
                    rowSums[i] += weights[j];
                }
            }
        }
        return rowSums;
    }
}
