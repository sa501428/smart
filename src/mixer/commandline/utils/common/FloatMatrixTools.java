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

import org.jetbrains.bio.npy.NpyFile;

import java.nio.file.Paths;
import java.util.Arrays;


/**
 * Helper methods to handle matrix operations
 */
public class FloatMatrixTools {
    
    public static void thresholdInPlaceByZscoreDownCols(float[][] matrix, float threshold, int batchSize) {
        float[] colMeans = getColMeansNonNan(matrix, batchSize);
        float[] colStdDevs = getColStdDevNonNans(matrix, colMeans, batchSize);
        
        for (int i = 0; i < matrix.length; i++) {
            for (int j = 0; j < matrix[i].length; j++) {
                float val = matrix[i][j];
                if (!Float.isNaN(val)) {
                    int newJ = j / batchSize;
                    float newVal = (val - colMeans[newJ]) / colStdDevs[newJ];
                    if (newVal > threshold || newVal < -threshold) {
                        matrix[i][j] = Float.NaN;
                    }
                }
            }
        }
    }
    
    public static void inPlaceZscoreDownColsNoNan(float[][] matrix, int batchSize) {
        float[] colMeans = getColMeansNonNan(matrix, batchSize);
        float[] colStdDevs = getColStdDevNonNans(matrix, colMeans, batchSize);
        
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
    
    public static float[] getColStdDevNonNans(float[][] matrix, float[] means, int batchSize) {
        
        float[] stdDevs = new float[means.length];
        int[] colNonNans = new int[means.length];
        
        for (int i = 0; i < matrix.length; i++) {
            for (int j = 0; j < matrix[i].length; j++) {
                float val = matrix[i][j];
                if (!Float.isNaN(val)) {
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
    
    public static float[] getColMeansNonNan(float[][] matrix, int batchSize) {
        float[] colMeans = new float[matrix[0].length / batchSize + 1];
        int[] colNonNans = new int[matrix[0].length / batchSize + 1];
        for (int i = 0; i < matrix.length; i++) {
            for (int j = 0; j < matrix[i].length; j++) {
                float val = matrix[i][j];
                if (!Float.isNaN(val)) {
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
    
    public static void saveMatrixTextNumpy(String filename, float[][] matrix, int[][] breakpoints) {
        saveMatrixTextNumpy(filename, matrix, breakpoints[0]);
    }

    public static void saveMatrixTextNumpy(String filename, float[][] matrix, int[] breakpoints) {
        long totalNum = ((long) matrix.length) * ((long) matrix[0].length);
        if (totalNum >= Integer.MAX_VALUE / 2) {
            for (int k = 2; k < breakpoints.length - 2; k *= 2) {
                System.gc();
                int numCols = breakpoints[k];
                float[] flattenedArray = FloatMatrixTools.getRowMajorOrderFlattendedSectionFromMatrix(matrix, numCols);
                NpyFile.write(Paths.get(filename + "_" + numCols + ".npy"), flattenedArray, new int[]{flattenedArray.length / numCols, numCols});
            }
        } else {
            saveMatrixTextNumpy(filename, matrix);
        }
        System.gc();
    }
    
    public static void saveMatrixTextNumpy(String filename, float[][] matrix) {
        int numRows = matrix.length;
        int numCols = matrix[0].length;
        float[] flattenedArray = FloatMatrixTools.flattenedRowMajorOrderMatrix(matrix);
        NpyFile.write(Paths.get(filename), flattenedArray, new int[]{numRows, numCols});
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
}
