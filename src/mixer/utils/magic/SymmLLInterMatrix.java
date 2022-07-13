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

package mixer.utils.magic;

public class SymmLLInterMatrix {

    private final float[][] matrix;
    private final int offsetR, n;

    // assume we are given lower left
    SymmLLInterMatrix(float[][] matrix) {
        this.matrix = matrix;
        this.offsetR = matrix[0].length;
        n = matrix.length + matrix[0].length;
    }

    public int[] getNumNonZerosInRow() {
        int[] numNonZero = new int[n];
        for (int i = 0; i < matrix.length; i++) {
            int r = i + offsetR;
            for (int j = 0; j < matrix[0].length; j++) {
                int c = j;
                if (matrix[i][j] > 0) {
                    numNonZero[r]++;
                    numNonZero[c]++;
                }
            }
        }
        return numNonZero;
    }

    public double[] sparseMultiply(float[] vector) {
        double[] sumVector = new double[n];

        for (int i = 0; i < matrix.length; i++) {
            int r = i + offsetR;
            for (int j = 0; j < matrix[0].length; j++) {
                int c = j;
                if (matrix[i][j] > 0) {
                    sumVector[r] += matrix[i][j] * vector[r];
                    sumVector[c] += matrix[i][j] * vector[c];
                }
            }
        }

        return sumVector;
    }


}
