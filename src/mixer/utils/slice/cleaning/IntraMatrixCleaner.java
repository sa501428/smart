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

    private static void nanFillTheRowsColumnsWeDontWant(Set<Integer> badIndices, float[][] matrix, int resFactor) {
        if (badIndices.size() < 1) return;

        for (int i = 0; i < matrix.length; i++) {
            for (int j : badIndices) {
                matrix[i][j / resFactor] = Float.NaN;
            }
        }

        for (int i : badIndices) {
            Arrays.fill(matrix[i / resFactor], Float.NaN);
        }
    }

    public static float[][] cleanAndCompress(Chromosome chrom, float[][] matrix, int resolution,
                                             Set<Integer> badIndices, int resFactor) {
        //NearDiagonalTrim.nanFill(chrom, matrix, resolution);
        nanFillTheRowsColumnsWeDontWant(badIndices, matrix, resFactor);
        // float[][] compressedMatrix = compress(matrix, smoothingInterval);
        //nanFillZeroEntries(matrix);
        //subtractOEBy1(matrix);
        // ZScoreTools.inPlaceZscoreDownCol(matrix);
        // return rollingAverage(matrix, smoothingInterval);
        return matrix;
    }

    private static void subtractOEBy1(float[][] matrix) {
        for (int i = 0; i < matrix.length; i++) {
            for (int j = 0; j < matrix[i].length; j++) {
                if (!Float.isNaN(matrix[i][j])) {
                    matrix[i][j] -= 1; // because 1 means O = E -> 0 for corr calcs
                }
            }
        }
    }

    private static void nanFillZeroEntries(float[][] matrix) {
        for (int i = 0; i < matrix.length; i++) {
            for (int j = 0; j < matrix[i].length; j++) {
                if (matrix[i][j] < 1e-10) {
                    matrix[i][j] = Float.NaN;
                }
            }
        }
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
        int width = (int) Math.ceil(interMatrix[0].length / ((float) compressionFactor));
        float[][] result = new float[interMatrix.length][width];
        for (int i = 0; i < interMatrix.length; i++) {
            for (int j = 0; j < interMatrix[i].length; j++) {
                result[i][j / compressionFactor] += interMatrix[i][j];
            }
        }
        return result;
    }
}
