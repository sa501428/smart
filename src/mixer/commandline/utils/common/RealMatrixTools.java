/*
 * The MIT License (MIT)
 *
 * Copyright (c) 2011-2020 Rice University, Baylor College of Medicine, Aiden Lab
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

package mixer.commandline.utils.common;

import org.apache.commons.math.linear.Array2DRowRealMatrix;
import org.apache.commons.math.linear.RealMatrix;

import java.util.Random;

public class RealMatrixTools {

    /**
     * @return matrix initialized with 0s
     */
    public static RealMatrix cleanArray2DMatrix(int rows, int cols) {
        return presetValueMatrix(rows, cols, 0);
    }

    /**
     * @return matrix of size m x n initialized with a specified value
     */
    private static RealMatrix presetValueMatrix(int numRows, int numCols, int val) {
        RealMatrix matrix = new Array2DRowRealMatrix(numRows, numCols);
        for (int r = 0; r < numRows; r++)
            for (int c = 0; c < numCols; c++)
                matrix.setEntry(r, c, val);
        return matrix;
    }

    /**
     * Generate a matrix with randomly initialized 1s and 0s
     *
     * @param rows number of rows
     * @param cols number of columns
     * @return randomized binary matrix
     */
    private static RealMatrix randomUnitMatrix(int rows, int cols) {
        Random generator = new Random();
        RealMatrix matrix = cleanArray2DMatrix(rows, cols);
        for (int r = 0; r < rows; r++)
            for (int c = 0; c < cols; c++)
                if (generator.nextBoolean())
                    matrix.setEntry(r, c, 1);
        return matrix;
    }

    /**
     * @return matrix flipped across the antidiagonal
     */
    public static RealMatrix flipAcrossAntiDiagonal(RealMatrix matrix) {
        int n = Math.min(matrix.getColumnDimension(), matrix.getRowDimension());
        RealMatrix antiDiagFlippedMatrix = cleanArray2DMatrix(n, n);
        int maxIndex = n - 1;
        for (int i = 0; i < n; i++) {
            for (int j = 0; j < n; j++) {
                antiDiagFlippedMatrix.setEntry(maxIndex - j, maxIndex - i, matrix.getEntry(i, j));
            }
        }
        return antiDiagFlippedMatrix;
    }

    /**
     * @return matrix flipped Left-Right
     */
    public static RealMatrix flipLeftRight(RealMatrix matrix) {
        int r = matrix.getRowDimension(), c = matrix.getColumnDimension();
        RealMatrix leftRightFlippedMatrix = cleanArray2DMatrix(r, c);
        for (int i = 0; i < r; i++) {
            for (int j = 0; j < c; j++) {
                leftRightFlippedMatrix.setEntry(i, c - 1 - j, matrix.getEntry(i, j));
            }
        }
        return leftRightFlippedMatrix;
    }

    /**
     * @return matrix flipped Top-Bottom
     */
    public static RealMatrix flipTopBottom(RealMatrix matrix) {
        int r = matrix.getRowDimension(), c = matrix.getColumnDimension();
        RealMatrix topBottomFlippedMatrix = cleanArray2DMatrix(r, c);
        for (int i = 0; i < r; i++) {
            for (int j = 0; j < c; j++) {
                topBottomFlippedMatrix.setEntry(r - 1 - i, j, matrix.getEntry(i, j));
            }
        }
        return topBottomFlippedMatrix;
    }

    /**
     * @return Element-wise multiplication of matrices i.e. M.*N in Matlab
     */
    public static RealMatrix elementBasedMultiplication(RealMatrix matrix1, RealMatrix matrix2) {
        // chooses minimal intersection of dimensions
        int r = Math.min(matrix1.getRowDimension(), matrix2.getRowDimension());
        int c = Math.min(matrix1.getColumnDimension(), matrix2.getColumnDimension());

        RealMatrix elementwiseMultipliedMatrix = cleanArray2DMatrix(r, c);
        for (int i = 0; i < r; i++) {
            for (int j = 0; j < c; j++) {
                elementwiseMultipliedMatrix.setEntry(i, j, matrix1.getEntry(i, j) * matrix2.getEntry(i, j));
            }
        }
        return elementwiseMultipliedMatrix;
    }


    /**
     * @return Element-wise division of matrices i.e. M./N in Matlab
     */
    public static RealMatrix elementBasedDivision(RealMatrix matrix1, RealMatrix matrix2) {
        // chooses minimal intersection of dimensions
        int r = Math.min(matrix1.getRowDimension(), matrix2.getRowDimension());
        int c = Math.min(matrix1.getColumnDimension(), matrix2.getColumnDimension());

        RealMatrix elementwiseDividedMatrix = cleanArray2DMatrix(r, c);
        for (int i = 0; i < r; i++) {
            for (int j = 0; j < c; j++) {
                elementwiseDividedMatrix.setEntry(i, j, matrix1.getEntry(i, j) / matrix2.getEntry(i, j));
            }
        }
        return elementwiseDividedMatrix;
    }


    /**
     * Replace NaNs in given matrix with given value
     */
    public static void setNaNs(RealMatrix matrix, int val) {
        for (int i = 0; i < matrix.getRowDimension(); i++) {
            for (int j = 0; j < matrix.getColumnDimension(); j++) {
                if (Double.isNaN(matrix.getEntry(i, j))) {
                    matrix.setEntry(i, j, val);
                }
            }
        }
    }

    /**
     * Return sign of values in matrix:
     * val > 0 : 1
     * val = 0 : 0
     * val < 0 : -1
     */
    public static RealMatrix sign(RealMatrix matrix) {
        int r = matrix.getRowDimension();
        int c = matrix.getColumnDimension();
        RealMatrix signMatrix = cleanArray2DMatrix(r, c);
        for (int i = 0; i < r; i++) {
            for (int j = 0; j < c; j++) {
                double val = matrix.getEntry(i, j);
                if (val > 0) {
                    signMatrix.setEntry(i, j, 1);
                } else if (val < 0) {
                    signMatrix.setEntry(i, j, -1);
                }
            }
        }
        return signMatrix;
    }

    /**
     * Flatten a 2D double matrix into a double array
     *
     * @param matrix
     * @return 1D double array in row major order
     */
    public static double[] flattenedRowMajorOrderMatrix(RealMatrix matrix) {
        int m = matrix.getRowDimension();
        int n = matrix.getColumnDimension();
        int numElements = m * n;
        double[] flattenedMatrix = new double[numElements];

        int index = 0;
        for (int i = 0; i < m; i++) {
            System.arraycopy(matrix.getRow(i), 0, flattenedMatrix, index, n);
            index += n;
        }
        return flattenedMatrix;
    }

    /**
     * Replace all of a given value in a matrix with a new value
     */
    public static void replaceValue(RealMatrix matrix, int initialVal, int newVal) {
        for (int i = 0; i < matrix.getRowDimension(); i++) {
            for (int j = 0; j < matrix.getColumnDimension(); j++) {
                if (matrix.getEntry(i, j) == initialVal) {
                    matrix.setEntry(i, j, newVal);
                }
            }
        }
    }

    /**
     * Normalize matrix by dividing by max element
     *
     * @return matrix * (1/max_element)
     */
    public static RealMatrix normalizeByMax(RealMatrix matrix) {
        double max = calculateMax(matrix);
        return matrix.scalarMultiply(1 / max);
    }

    /**
     * @return max element in matrix
     */
    public static double calculateMax(RealMatrix matrix) {
        double max = matrix.getEntry(0, 0);
        for (int i = 0; i < matrix.getRowDimension(); i++) {
            for (int j = 0; j < matrix.getColumnDimension(); j++) {
                double val = matrix.getEntry(i, j);
                if (max < val) {
                    max = val;
                }
            }
        }
        return max;
    }

    /**
     * @return min element in matrix
     */
    public static double calculateMin(RealMatrix matrix) {
        double min = matrix.getEntry(0, 0);
        for (int i = 0; i < matrix.getRowDimension(); i++) {
            for (int j = 0; j < matrix.getColumnDimension(); j++) {
                double val = matrix.getEntry(i, j);
                if (min > val) {
                    min = val;
                }
            }
        }
        return min;
    }

    /**
     * print for matrix
     */
    public static void print(RealMatrix matrix) {
        DoubleMatrixTools.print(matrix.getData());
    }

    /**
     * @return region within matrix specified by indices
     */
    public static RealMatrix getSubMatrix(RealMatrix matrix, int[] indices) {
        return matrix.getSubMatrix(indices[0], indices[1], indices[2], indices[3]);
    }

    /**
     * Fill lower left triangle with values from upper right triangle
     *
     * @param matrix
     * @return
     */
    public static RealMatrix fillLowerLeftTriangle(RealMatrix matrix) {
        for (int r = 0; r < matrix.getRowDimension(); r++)
            for (int c = 0; c < matrix.getColumnDimension(); c++)
                matrix.setEntry(c, r, matrix.getEntry(r, c));
        return matrix;
    }

    public static void thresholdValues(RealMatrix matrix, int val) {
        for (int i = 0; i < matrix.getRowDimension(); i++) {
            for (int j = 0; j < matrix.getColumnDimension(); j++) {
                if (matrix.getEntry(i, j) > val) {
                    matrix.setEntry(i, j, val);
                }
            }
        }
    }

    public static void thresholdValuesDouble(RealMatrix matrix, double lowVal, double highVal) {
        for (int i = 0; i < matrix.getRowDimension(); i++) {
            for (int j = 0; j < matrix.getColumnDimension(); j++) {
                if (matrix.getEntry(i, j) > highVal) {
                    matrix.setEntry(i, j, highVal);
                }
                if (matrix.getEntry(i, j) < lowVal) {
                    matrix.setEntry(i, j, lowVal);
                }
            }
        }
    }


    public static void cleanUpNaNs(RealMatrix matrix) {
        for (int r = 0; r < matrix.getRowDimension(); r++) {
            for (int c = 0; c < matrix.getColumnDimension(); c++) {
                if (Double.isNaN(matrix.getEntry(r, c))) {
                    matrix.setEntry(r, c, 0);
                }
            }
        }
    }
}
