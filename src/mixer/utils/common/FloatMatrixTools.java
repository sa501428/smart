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

package mixer.utils.common;

import javastraw.tools.MatrixTools;
import org.apache.commons.math3.stat.correlation.PearsonsCorrelation;

import java.util.Arrays;
import java.util.List;


/**
 * Helper methods to handle matrix operations
 */
public class FloatMatrixTools {

    private static final PearsonsCorrelation pearsonsCorrelation = new PearsonsCorrelation();

    public static void thresholdNonZerosByZscoreToNanDownColumn(float[][] matrix, float threshold, int batchSize) {
        float[] colMeans = getColNonZeroMeansNonNan(matrix, batchSize);
        float[] colStdDevs = getColNonZeroStdDevNonNans(matrix, colMeans, batchSize);

        for (int i = 0; i < matrix.length; i++) {
            for (int j = 0; j < matrix[i].length; j++) {
                float val = matrix[i][j];
                if (!Float.isNaN(val) && val > 1e-10) {
                    int newJ = j / batchSize;
                    float newVal = (val - colMeans[newJ]) / colStdDevs[newJ];
                    if (newVal > threshold) { // || newVal < -threshold || val < 1e-10
                        matrix[i][j] = Float.NaN;
                    }
                }
            }
        }
    }

    public static void inPlaceZscoreDownColsNoNan(float[][] matrix, int batchSize) {
        float[] colMeans = getColNonZeroMeansNonNan(matrix, batchSize);
        float[] colStdDevs = getColNonZeroStdDevNonNans(matrix, colMeans, batchSize);

        for (int i = 0; i < matrix.length; i++) {
            for (int j = 0; j < matrix[i].length; j++) {
                float val = matrix[i][j];
                if (!Float.isNaN(val)) {
                    int newJ = j / batchSize;
                    matrix[i][j] = (val - colMeans[newJ]) / colStdDevs[newJ];
                }
            }
        }
    }

    public static float[] getColNonZeroStdDevNonNans(float[][] matrix, float[] means, int batchSize) {

        float[] stdDevs = new float[means.length];
        int[] colNonNans = new int[means.length];

        for (int i = 0; i < matrix.length; i++) {
            for (int j = 0; j < matrix[i].length; j++) {
                float val = matrix[i][j];
                if (!Float.isNaN(val) && val > 1e-10) {
                    int newJ = j / batchSize;
                    float diff = val - means[newJ];
                    stdDevs[newJ] += diff * diff;
                    colNonNans[newJ] += 1;
                }
            }
        }

        for (int k = 0; k < stdDevs.length; k++) {
            stdDevs[k] = (float) Math.sqrt(stdDevs[k] / Math.max(colNonNans[k], 1));
        }

        return stdDevs;
    }

    public static float[] getColNonZeroMeansNonNan(float[][] matrix, int batchSize) {
        float[] colMeans = new float[matrix[0].length / batchSize + 1];
        int[] colNonNans = new int[matrix[0].length / batchSize + 1];
        for (int i = 0; i < matrix.length; i++) {
            for (int j = 0; j < matrix[i].length; j++) {
                float val = matrix[i][j];
                if (!Float.isNaN(val) && val > 1e-10) {
                    int newJ = j / batchSize;
                    colMeans[newJ] += val;
                    colNonNans[newJ] += 1;
                }
            }
        }

        for (int k = 0; k < colMeans.length; k++) {
            colMeans[k] = colMeans[k] / Math.max(colNonNans[k], 1);
        }

        return colMeans;
    }

    public static float[][] fill(float[][] allDataForRegion, float val) {
        for (int i = 0; i < allDataForRegion.length; i++) {
            Arrays.fill(allDataForRegion[i], val);
        }
        return allDataForRegion;
    }

    public static void saveMatrixTextNumpy(String filename, float[][] matrix) {
        MatrixTools.saveMatrixTextNumpy(filename, matrix);
    }

    public static float[] flattenedRowMajorOrderMatrix(float[][] matrix) {
        int m = matrix.length;
        int n = matrix[0].length;

        int numElements = m * n;
        float[] flattenedMatrix = new float[numElements];

        int index = 0;
        for (int i = 0; i < m; i++) {
            System.arraycopy(matrix[i], 0, flattenedMatrix, index, n);
            index += n;
        }
        return flattenedMatrix;
    }

    public static float[] getRowMajorOrderFlattendedSectionFromMatrix(float[][] matrix, int numCols) {
        int numRows = matrix.length - numCols;

        int numElements = numRows * numCols;
        float[] flattenedMatrix = new float[numElements];

        int index = 0;
        for (int i = numCols; i < numRows; i++) {
            System.arraycopy(matrix[i], 0, flattenedMatrix, index, numCols);
            index += numCols;
        }
        return flattenedMatrix;
    }

    public static float[][] concatenate(float[][] matrix1, float[][] matrix2) {
        float[][] combo = new float[matrix1.length][matrix1[0].length + matrix2[0].length];
        for (int i = 0; i < matrix1.length; i++) {
            System.arraycopy(matrix1[i], 0, combo[i], 0, matrix1[i].length);
            System.arraycopy(matrix2[i], 0, combo[i], matrix1[i].length, matrix2[i].length);
        }
        return combo;
    }

    public static float[][] transpose(float[][] matrix) {
        float[][] result = new float[matrix[0].length][matrix.length];
        for (int i = 0; i < matrix.length; i++) {
            for (int j = 0; j < matrix[i].length; j++) {
                result[j][i] = matrix[i][j];
            }
        }
        return result;
    }

    // column length assumed identical and kept the same
    public static float[][] stitchMultipleMatricesTogetherByRowDim(List<float[][]> data) {
        if (data.size() == 1) return data.get(0);

        int colNums = data.get(0)[0].length;
        int rowNums = 0;
        for (float[][] mtrx : data) {
            rowNums += mtrx.length;
        }

        float[][] aggregate = new float[rowNums][colNums];

        int rowOffSet = 0;
        for (float[][] region : data) {
            copyFromAToBRegion(region, aggregate, rowOffSet, 0);
            rowOffSet += region.length;
        }

        return aggregate;
    }

    public static void copyFromAToBRegion(float[][] source, float[][] destination, int rowOffSet, int colOffSet) {
        for (int i = 0; i < source.length; i++) {
            System.arraycopy(source[i], 0, destination[i + rowOffSet], colOffSet, source[0].length);
        }
    }

    // column length assumed identical and kept the same
    public static float[][] stitchMultipleMatricesTogetherByColDim(List<float[][]> data) {
        if (data.size() == 1) return data.get(0);

        int rowNums = data.get(0).length;
        int colNums = 0;
        for (float[][] mtrx : data) {
            colNums += mtrx[0].length;
        }

        float[][] aggregate = new float[rowNums][colNums];

        int colOffSet = 0;
        for (float[][] region : data) {
            copyFromAToBRegion(region, aggregate, 0, colOffSet);
            colOffSet += region[0].length;
        }

        return aggregate;
    }

    public static float[][] cleanUpMatrix(float[][] matrix) {
        for (int r = 0; r < matrix.length; r++) {
            for (int c = 0; c < matrix[r].length; c++) {
                if (Float.isNaN(matrix[r][c]) || Float.isInfinite(matrix[r][c]) || Math.abs(matrix[r][c]) < 1E-10) {
                    matrix[r][c] = 0;
                }
            }
        }
        return matrix;
    }

    public static void addBToA(float[][] a, float[][] b) {
        if (a.length == b.length && a[0].length == b[0].length) {
            for (int i = 0; i < a.length; i++) {
                for (int j = 0; j < a[i].length; j++) {
                    a[i][j] += b[i][j];
                }
            }
        } else {
            System.err.println("dimensions incorrect " + a.length + "==" + b.length
                    + "; " + a[0].length + "==" + b[0].length);
        }
    }

    public static void scaleBy(float[][] matrix, float scalar) {
        for (int i = 0; i < matrix.length; i++) {
            for (int j = 0; j < matrix[i].length; j++) {
                matrix[i][j] *= scalar;
            }
        }
    }
}
