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

package mixer.utils.shuffle;

public class TensorTools {

    public static double[][][] log(double[][][] a) {
        double[][][] result = new double[a.length][a[0].length][a[0][0].length];
        for (int i = 0; i < a.length; i++) {
            for (int j = 0; j < a[i].length; j++) {
                for (int k = 0; k < a[i][j].length; k++) {
                    if (a[i][j][k] > 0) {
                        result[i][j][k] = Math.log(1 + a[i][j][k]);
                    }
                }
            }
        }
        return result;
    }


    public static void addBtoA(double[][][] a, double[][][] b) {
        for (int i = 0; i < a.length; i++) {
            for (int j = 0; j < a[i].length; j++) {
                for (int k = 0; k < a[i][j].length; k++) {
                    a[i][j][k] += b[i][j][k];
                }
            }
        }
    }

    public static void addBtoA(long[][][] a, long[][][] b) {
        for (int i = 0; i < a.length; i++) {
            for (int j = 0; j < a[i].length; j++) {
                for (int k = 0; k < a[i][j].length; k++) {
                    a[i][j][k] += b[i][j][k];
                }
            }
        }
    }

    public static double[][] makeSymmetric(double[][] input) {
        int n = input.length;
        double[][] results = new double[n][n];
        for (int i = 0; i < n; i++) {
            for (int j = 0; j < n; j++) {
                results[i][j] += input[i][j];
                results[j][i] += input[i][j];
            }
        }
        return results;
    }

    public static long[][] makeSymmetric(long[][] input) {
        int n = input.length;
        long[][] results = new long[n][n];
        for (int i = 0; i < n; i++) {
            for (int j = 0; j < n; j++) {
                results[i][j] += input[i][j];
                results[j][i] += input[i][j];
            }
        }
        return results;
    }

    public static double[][] divideBy(double[][] totals, long[][] areas) {
        for (int i = 0; i < totals.length; i++) {
            for (int j = 0; j < totals[i].length; j++) {
                if (areas[i][j] > 0) {
                    totals[i][j] /= areas[i][j];
                }
            }
        }
        return totals;
    }

    public static double[][][] divide(double[][][] a, long[][][] b) {
        double[][][] result = new double[a.length][a[0].length][a[0][0].length];
        for (int i = 0; i < a.length; i++) {
            for (int j = 0; j < a[i].length; j++) {
                for (int k = 0; k < a[i][j].length; k++) {
                    if (b[i][j][k] > 0) {
                        result[i][j][k] = a[i][j][k] / b[i][j][k];
                    }
                }
            }
        }
        return result;
    }

    public static float[][] concatenate(double[][][] input) {
        int n = input[0].length;
        float[][] matrix = new float[2 * n][2 * n];
        for (int i = 0; i < n; i++) {
            for (int j = 0; j < n; j++) {
                matrix[2 * i][2 * j] = (float) input[0][i][j];
                matrix[2 * i + 1][2 * j] = (float) input[1][i][j];
                matrix[2 * i][2 * j + 1] = (float) input[2][i][j];
                matrix[2 * i + 1][2 * j + 1] = (float) input[3][i][j];
            }
        }
        return matrix;
    }
}
