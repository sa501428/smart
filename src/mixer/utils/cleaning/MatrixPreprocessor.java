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

import javastraw.expected.Welford;
import javastraw.expected.Zscore;
import javastraw.reader.basics.Chromosome;
import mixer.utils.common.SimpleArray2DTools;
import mixer.utils.drive.MatrixAndWeight;
import mixer.utils.transform.MatrixTransform;

public class MatrixPreprocessor {

    private static final int ZSCORE_LIMIT = 3;

    public static MatrixAndWeight clean2(MatrixAndWeight matrix, Chromosome[] chromosomes) {
        matrix.updateWeights(chromosomes);
        matrix.divideColumnsByWeights();
        SimpleArray2DTools.simpleLogWithCleanup(matrix.matrix, Float.NaN);
        MatrixTransform.zscoreByRows(matrix.matrix, ZSCORE_LIMIT);
        matrix.removeAllNanRows();
        return matrix;
    }

    private static int getNumCentroids(int numRows, int numCols) {
        return Math.min(numCols, numRows / 20);
    }

    private static void thresholdGlobally(float[][] matrix, int upperLimit) {
        Welford welford = new Welford();
        for (int r = 0; r < matrix.length; r++) {
            for (int c = 0; c < matrix[r].length; c++) {
                if (matrix[r][c] > 0) {
                    welford.addValue(matrix[r][c]);
                }
            }
        }
        Zscore zscore = welford.getZscore();
        for (int r = 0; r < matrix.length; r++) {
            for (int c = 0; c < matrix[r].length; c++) {
                if (zscore.getZscore(matrix[r][c]) > upperLimit) {
                    matrix[r][c] = Float.NaN;
                }
            }
        }
    }

    private static void cleanUpZerosAndInfs(float[][] matrix) {
        for (int r = 0; r < matrix.length; r++) {
            for (int c = 0; c < matrix[r].length; c++) {
                float val = matrix[r][c];
                if (val > 0) {
                    if (Float.isInfinite(val)) {
                        matrix[r][c] = Float.NaN;
                    }
                } else {
                    matrix[r][c] = Float.NaN;
                }
            }
        }
    }
}
