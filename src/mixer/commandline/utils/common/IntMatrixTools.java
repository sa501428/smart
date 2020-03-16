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

import java.io.*;
import java.nio.charset.StandardCharsets;
import java.nio.file.Paths;
import java.util.Arrays;

public class IntMatrixTools {

    public static int[] flattenedRowMajorOrderMatrix(int[][] matrix) {
        int m = matrix.length;
        int n = matrix[0].length;

        int numElements = m * n;
        int[] flattenedMatrix = new int[numElements];

        int index = 0;
        for (int i = 0; i < m; i++) {
            System.arraycopy(matrix[i], 0, flattenedMatrix, index, n);
            index += n;
        }
        return flattenedMatrix;
    }

    public static int[][] normalizeMatrixUsingColumnSum(int[][] matrix) {
        int[][] newMatrix = new int[matrix.length][matrix[0].length];
        int[] columnSum = new int[matrix[0].length];
        for (int[] row : matrix) {
            for (int i = 0; i < row.length; i++) {
                columnSum[i] += row[i];
            }
        }

        for (int i = 0; i < matrix.length; i++) {
            for (int j = 0; j < matrix[i].length; j++) {
                newMatrix[i][j] = matrix[i][j] / columnSum[j];
            }
        }

        return newMatrix;
    }

    public static int[][] normalizeMatrixUsingRowSum(int[][] matrix) {
        int[][] newMatrix = new int[matrix.length][matrix[0].length];
        int[] rowSum = getRowSums(matrix);

        for (int i = 0; i < matrix.length; i++) {
            for (int j = 0; j < matrix[i].length; j++) {
                newMatrix[i][j] = matrix[i][j] / rowSum[i];
            }
        }

        return newMatrix;
    }

    public static int[] getRowSums(int[][] matrix) {
        int[] rowSum = new int[matrix.length];
        for (int i = 0; i < matrix.length; i++) {
            for (int val : matrix[i]) {
                rowSum[i] += val;
            }
        }
        return rowSum;
    }

    public static int[] getAbsValColSums(int[][] matrix) {
        int[] colSum = new int[matrix[0].length];
        for (int i = 0; i < matrix.length; i++) {
            for (int j = 0; j < matrix[i].length; j++) {
                colSum[j] += Math.abs(matrix[i][j]);
            }
        }
        return colSum;
    }

    public static void saveMatrixTextV2(String filename, int[][] matrix) {
        Writer writer = null;
        try {
            writer = new BufferedWriter(new OutputStreamWriter(new FileOutputStream(filename), StandardCharsets.UTF_8));
            for (int[] row : matrix) {
                String s = Arrays.toString(row);//.replaceAll().replaceAll("]","").trim();
                s = s.replaceAll("\\[", "").replaceAll("\\]", "").trim();
                writer.write(s + "\n");
            }
        } catch (IOException ex) {
            ex.printStackTrace();
        } finally {
            try {
                if (writer != null)
                    writer.close();
            } catch (Exception ex) {
                ex.printStackTrace();
            }
        }
    }

    public static void saveMatrixTextNumpy(String filename, int[][] matrix) {
        int numRows = matrix.length;
        int numCols = matrix[0].length;
        int[] flattenedArray = IntMatrixTools.flattenedRowMajorOrderMatrix(matrix);

        NpyFile.write(Paths.get(filename), flattenedArray, new int[]{numRows, numCols});
    }

    public static void saveMatrixTextNumpy(String filename, int[] matrix) {
        NpyFile.write(Paths.get(filename), matrix, new int[]{1, matrix.length});
    }

    public static void labelRegionWithOnes(int[][] labelsMatrix, int rowLength, int numRows, int colLength, int numCols, int startRowOf1, int startColOf1) {
        for (int i = 0; i < Math.min(rowLength, numRows); i++) {
            for (int j = 0; j < Math.min(colLength, numCols); j++) {
                labelsMatrix[startRowOf1 + i][startColOf1 + j] = 1;
            }
        }
    }
}
