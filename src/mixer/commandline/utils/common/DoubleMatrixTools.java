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

import mixer.data.ContactRecord;
import org.jetbrains.bio.npy.NpyFile;

import java.io.*;
import java.nio.charset.StandardCharsets;
import java.nio.file.Paths;
import java.text.DecimalFormat;
import java.util.Arrays;
import java.util.List;

public class DoubleMatrixTools {


    /**
     * @return minimal positive entry in the matrix greater than 0
     */
    private static double minimumPositive(double[][] data) {
        double minVal = Double.MAX_VALUE;
        for (double[] row : data) {
            for (double val : row) {
                if (val > 0 && val < minVal)
                    minVal = val;
            }
        }
        if (minVal == Double.MAX_VALUE)
            minVal = 0;
        return minVal;
    }

    private static double mean(double[][] data) {
        double average = 0;
        if (data.length > 0) {
            double total = 0;
            for (double[] vals : data) {
                for (double val : vals) {
                    total += val;
                }
            }
            average = (total / data.length) / data[0].length;
        }
        return average;
    }


    public static double[] flattenedRowMajorOrderMatrix(double[][] matrix) {
        int m = matrix.length;
        int n = matrix[0].length;

        int numElements = m * n;
        double[] flattenedMatrix = new double[numElements];

        int index = 0;
        for (int i = 0; i < m; i++) {
            System.arraycopy(matrix[i], 0, flattenedMatrix, index, n);
            index += n;
        }
        return flattenedMatrix;
    }


    /**
     * print for 2D array
     */
    public static void print(double[][] data) {
        for (double[] row : data) {
            System.out.println(Arrays.toString(row));
        }
    }


    public static double[] getRowSums(double[][] matrix) {
        double[] rowSum = new double[matrix.length];
        for (int i = 0; i < matrix.length; i++) {
            for (double val : matrix[i]) {
                rowSum[i] += val;
            }
        }
        return rowSum;
    }


    public static double[] getRowSums(List<ContactRecord> unNormedRecordList, double scalar, double[] normVector) {
        double[] rowSum = new double[normVector.length];
        for (ContactRecord record : unNormedRecordList) {
            int x = record.getBinX();
            int y = record.getBinY();
            float counts = record.getCounts();

            double normVal = counts * scalar / (normVector[x] * normVector[y]);
            rowSum[x] += normVal;
            if (x != y) {
                rowSum[y] += normVal;
            }

        }
        return rowSum;
    }


    public static void cleanUpNaNs(double[][] matrix) {
        for (int r = 0; r < matrix.length; r++) {
            for (int c = 0; c < matrix[r].length; c++) {
                if (Double.isNaN(matrix[r][c])) {
                    matrix[r][c] = 0;
                }
            }
        }
    }


    public static double sum(double[][] data) {
        double sum = 0;
        for (double[] row : data) {
            for (double val : row) {
                sum += val;
            }
        }
        return sum;
    }

    public static void exportData(double[][] data, File file) {
        try {
            DecimalFormat df = new DecimalFormat("##.###");

            final FileWriter fw = new FileWriter(file);
            for (double[] row : data) {
                for (double val : row) {
                    if (Double.isNaN(val)) {
                        fw.write("NaN, ");
                    } else {
                        fw.write(Double.valueOf(df.format(val)) + ", ");
                    }
                }
                fw.write("0\n");
            }
            fw.close();
        } catch (Exception e) {
            System.err.println("Error exporting matrix");
            e.printStackTrace();
            System.exit(86);
        }
    }

    public static double[][] transpose(double[][] matrix) {
        int h0 = matrix.length;
        int w0 = matrix[0].length;
        double[][] transposedMatrix = new double[w0][h0];

        for (int i = 0; i < h0; i++) {
            for (int j = 0; j < w0; j++) {
                transposedMatrix[j][i] = matrix[i][j];
            }
        }
        return transposedMatrix;
    }

    public static double[][] convertToDoubleMatrix(boolean[][] adjacencyMatrix) {
        double[][] matrix = new double[adjacencyMatrix.length][adjacencyMatrix[0].length];
        for (int i = 0; i < adjacencyMatrix.length; i++) {
            for (int j = 0; j < adjacencyMatrix[0].length; j++) {
                if (adjacencyMatrix[i][j]) {
                    matrix[i][j] = 1;
                }
            }
        }
        return matrix;
    }

    public static double[][] convertToDoubleMatrix(int[][] adjacencyMatrix) {
        double[][] matrix = new double[adjacencyMatrix.length][adjacencyMatrix[0].length];
        for (int i = 0; i < adjacencyMatrix.length; i++) {
            for (int j = 0; j < adjacencyMatrix[0].length; j++) {
                matrix[i][j] = adjacencyMatrix[i][j];
            }
        }
        return matrix;
    }

    public static float[][] convertToFloatMatrix(double[][] dataMatrix) {
        float[][] matrix = new float[dataMatrix.length][dataMatrix[0].length];
        for (int i = 0; i < dataMatrix.length; i++) {
            for (int j = 0; j < dataMatrix[0].length; j++) {
                matrix[i][j] = (float) dataMatrix[i][j];
            }
        }
        return matrix;
    }


    public static void saveMatrixTextV2(String filename, double[][] matrix) {
        Writer writer = null;
        try {
            writer = new BufferedWriter(new OutputStreamWriter(new FileOutputStream(filename), StandardCharsets.UTF_8));
            for (double[] row : matrix) {
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

    public static void saveMatrixTextNumpy(String filename, double[][] matrix) {
        int numRows = matrix.length;
        int numCols = matrix[0].length;
        double[] flattenedArray = flattenedRowMajorOrderMatrix(matrix);

        NpyFile.write(Paths.get(filename), flattenedArray, new int[]{numRows, numCols});
    }

    public static double[][] deepClone(double[][] data) {
        double[][] copy = new double[data.length][data[0].length];
        for (int i = 0; i < data.length; i++) {
            System.arraycopy(data[i], 0, copy[i], 0, data[i].length);
        }
        return copy;
    }


    // column length assumed identical and kept the same
    public static double[][] stitchMultipleMatricesTogetherByRowDim(List<double[][]> data) {
        // todo currently assuming each one identical...

        int colNums = data.get(0)[0].length;
        int rowNums = 0;
        for (double[][] mtrx : data) {
            rowNums += mtrx.length;
        }

        double[][] aggregate = new double[rowNums][colNums];

        int rowOffSet = 0;
        for (double[][] region : data) {
            copyFromAToBRegion(region, aggregate, rowOffSet, 0);
            rowOffSet += region.length;
        }

        return aggregate;
    }

    // column length assumed identical and kept the same
    public static double[][] stitchMultipleMatricesTogetherByColDim(List<double[][]> data) {
        // todo currently assuming each one identical...

        int rowNums = data.get(0).length;
        int colNums = 0;
        for (double[][] mtrx : data) {
            colNums += mtrx[0].length;
        }

        double[][] aggregate = new double[rowNums][colNums];

        int colOffSet = 0;
        for (double[][] region : data) {
            copyFromAToBRegion(region, aggregate, 0, colOffSet);
            colOffSet += region[0].length;
        }

        return aggregate;
    }


    public static double[][] smoothAndAppendDerivativeDownColumn(double[][] data, double[] convolution) {

        int numColumns = data[0].length;
        if (convolution != null && convolution.length > 1) {
            numColumns -= (convolution.length - 1);
        }

        double[][] appendedDerivative = new double[data.length][2 * numColumns - 1];

        if (convolution != null && convolution.length > 1) {
            for (int i = 0; i < data.length; i++) {
                for (int j = 0; j < numColumns; j++) {
                    for (int k = 0; k < convolution.length; k++) {
                        appendedDerivative[i][j] += convolution[k] * data[i][j + k];
                    }
                }
            }
        } else {
            for (int i = 0; i < data.length; i++) {
                System.arraycopy(data[i], 0, appendedDerivative[i], 0, numColumns);
            }
        }

        for (int i = 0; i < data.length; i++) {
            for (int j = 0; j < numColumns - 1; j++) {
                appendedDerivative[i][numColumns + j] = appendedDerivative[i][j] - appendedDerivative[i][j + 1];
            }
        }

        return appendedDerivative;
    }

    public static void copyFromAToBRegion(double[][] source, double[][] destination, int rowOffSet, int colOffSet) {
        for (int i = 0; i < source.length; i++) {
            System.arraycopy(source[i], 0, destination[i + rowOffSet], colOffSet, source[0].length);
        }
    }

    public static void labelEnrichedRegionWithOnes(int[][] labelsMatrix, double[][] data, int rowLength, int numRows, int colLength, int numCols, int startRowOf1, int startColOf1) {
        double total = 0;
        int numVals = 0;

        for (int i = 0; i < Math.min(rowLength, numRows); i++) {
            for (int j = 0; j < Math.min(colLength, numCols); j++) {
                total += data[startRowOf1 + i][startColOf1 + j];
                numVals++;
            }
        }
        double average = total / numVals;

        for (int i = 0; i < Math.min(rowLength, numRows); i++) {
            for (int j = 0; j < Math.min(colLength, numCols); j++) {
                if (data[startRowOf1 + i][startColOf1 + j] > average) {
                    labelsMatrix[startRowOf1 + i][startColOf1 + j] = 1;
                }
            }
        }
    }


    public static double[][] takeDerivativeDownColumn(double[][] data) {
        double[][] derivative = new double[data.length][data[0].length - 1];

        for (int i = 0; i < data.length; i++) {
            System.arraycopy(data[i], 0, derivative[i], 0, derivative[i].length);
        }
        for (int i = 0; i < derivative.length; i++) {
            for (int j = 0; j < derivative[i].length; j++) {
                derivative[i][j] -= data[i][j + 1];
            }
        }

        return derivative;
    }

}
