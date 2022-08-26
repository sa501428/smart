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

package mixer.utils.intra;

import java.util.Arrays;
import java.util.Set;

public class IntraMatrixCleaner {

    public static float[][] basicClean(float[][] matrix, Set<Integer> badIndices, int pixelDist) {
        nanFillTheRowsColumnsWeDontWant(badIndices, matrix);
        NearDiagonalTrim.trimDiagonalWithinPixelDist(matrix, pixelDist);
        return matrix;
    }

    private static void nanFillTheRowsColumnsWeDontWant(Set<Integer> badIndices, float[][] matrix) {
        if (badIndices.size() < 1) return;

        for (int i = 0; i < matrix.length; i++) {
            for (int j : badIndices) {
                matrix[i][j] = Float.NaN;
            }
        }

        for (int i : badIndices) {
            Arrays.fill(matrix[i], Float.NaN);
        }
    }

    public static float[][] oeClean(float[][] matrix, Set<Integer> badIndices) {
        basicClean(matrix, badIndices, 3);
        nanFillZeroEntries(matrix);
        //nanFillExtremeOEValues(matrix); todo
        return matrix;
    }

    private static void nanFillZeroEntries(float[][] matrix) {
        for (int i = 0; i < matrix.length; i++) {
            for (int j = 0; j < matrix[i].length; j++) {
                if (Math.abs(matrix[i][j]) < 1e-10) {
                    matrix[i][j] = Float.NaN;
                }
            }
        }
    }
}
