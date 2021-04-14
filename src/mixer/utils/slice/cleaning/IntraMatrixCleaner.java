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

import javastraw.reader.basics.Chromosome;

import java.util.Arrays;
import java.util.Set;

public class IntraMatrixCleaner {

    private static void eraseTheRowsColumnsWeDontWant(Set<Integer> badIndices, float[][] matrix) {
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

    public static float[][] clean(Chromosome chrom, float[][] matrix, int resolution,
                                  int smoothingInterval, Set<Integer> badIndices) {
        NearDiagonalTrim.trim(chrom, matrix, resolution);
        eraseTheRowsColumnsWeDontWant(badIndices, matrix);
        return rollingAverage(matrix, smoothingInterval);

    }

    public static float[][] rollingAverage(float[][] matrix, int smoothingInterval) {
        int bufferWidth = smoothingInterval / 2;
        for (int i = 0; i < matrix.length; i++) {
            float[] tempRow = new float[matrix[i].length];
            System.arraycopy(matrix[i], 0, tempRow, 0, tempRow.length);
            for (int j = 0; j < tempRow.length; j++) {
                if (Float.isNaN(tempRow[j])) continue;
                int startK = Math.max(j - bufferWidth, 0);
                int endK = Math.min(j + 1 + bufferWidth, tempRow.length);
                float total = 0;
                int numVals = 0;
                for (int k = startK; k < endK; k++) {
                    if (!Float.isNaN(tempRow[k])) {
                        total += tempRow[k];
                        numVals++;
                    }
                }
                matrix[i][j] = total / Math.max(numVals, 1);
            }
        }
        return matrix;
    }

    public static float[][] compress(float[][] interMatrix, int compressionFactor) {
        int width = interMatrix[0].length / compressionFactor + 1;
        float[][] result = new float[interMatrix.length][width];
        for (int i = 0; i < interMatrix.length; i++) {
            for (int j = 0; j < interMatrix[i].length; j++) {
                result[i][j / 3] += interMatrix[i][j];
            }
        }
        return result;
    }
}
