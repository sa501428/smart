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

import org.apache.commons.math3.stat.descriptive.rank.Median;
import org.jetbrains.bio.npy.NpyFile;

import java.io.*;
import java.nio.charset.StandardCharsets;
import java.nio.file.Paths;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;


/**
 * Helper methods to handle matrix operations
 */
public class FloatMatrixTools {


    private static float mean(float[][] data) {
        double average = 0;
        if (data.length > 0) {
            double total = 0;
            for (float[] vals : data) {
                for (float val : vals) {
                    total += val;
                }
            }
            average = (total / data.length) / data[0].length;
        }
        return (float) average;
    }

    public static float[][] add(float[][] first, float[][] second, float scaleFirst, float scaleSecond) {
        float[][] answer = new float[first.length][first[0].length];
        for (int i = 0; i < answer.length; i++) {
            for (int j = 0; j < answer[i].length; j++) {
                answer[i][j] = scaleFirst * first[i][j] + scaleSecond * second[i][j];
            }
        }
        return answer;
    }

    public static float[][] max(float[][] first, float[][] second) {
        float[][] answer = new float[first.length][first[0].length];
        for (int i = 0; i < answer.length; i++) {
            for (int j = 0; j < answer[i].length; j++) {
                answer[i][j] = Math.max(first[i][j], second[i][j]);
            }
        }
        return answer;
    }

    public static void inPlaceZscoreDownCols(float[][] matrix) {
        float[] colMeans = getColMeansNonNan(matrix);
        float[] colStdDevs = getColStdDevNonNans(matrix, colMeans);

        for (int i = 0; i < matrix.length; i++) {
            for (int j = 0; j < matrix[i].length; j++) {
                float val = matrix[i][j];
                if (!Float.isNaN(val)) {
                    float newVal = (val - colMeans[j]) / colStdDevs[j];
                    newVal = Math.min(5, newVal);
                    newVal = Math.max(-5, newVal);
                    matrix[i][j] = newVal;
                }
            }
        }
    }

    private static float[] getColStdDevNonNans(float[][] matrix, float[] means) {

        float[] stdDevs = new float[means.length];
        int[] colNonNans = new int[means.length];

        for (int i = 0; i < matrix.length; i++) {
            for (int j = 0; j < matrix[i].length; j++) {
                float val = matrix[i][j];
                if (!Float.isNaN(val)) {
                    float diff = val - means[j];
                    stdDevs[j] += diff * diff;
                    colNonNans[j] += 1;
                }
            }
        }

        for (int k = 0; k < stdDevs.length; k++) {
            stdDevs[k] = (float) Math.sqrt(stdDevs[k] / Math.max(colNonNans[k], 1));
        }

        return stdDevs;
    }

    private static float[] getColMeansNonNan(float[][] matrix) {
        float[] colMeans = new float[matrix[0].length];
        int[] colNonNans = new int[matrix[0].length];
        for (int i = 0; i < matrix.length; i++) {
            for (int j = 0; j < matrix[i].length; j++) {
                float val = matrix[i][j];
                if (!Float.isNaN(val)) {
                    colMeans[j] += val;
                    colNonNans[j] += 1;
                }
            }
        }

        for (int k = 0; k < colMeans.length; k++) {
            colMeans[k] = colMeans[k] / Math.max(colNonNans[k], 1);
        }

        return colMeans;
    }

    public static float[][] inPlaceZscoreDownRows(float[][] matrix, float threshold) {
        float[] rowMeans = getRowSums(matrix);
        for (int k = 0; k < rowMeans.length; k++) {
            rowMeans[k] = rowMeans[k] / matrix[k].length;
        }

        float[] stdDevs = new float[matrix.length];
        for (int k = 0; k < matrix.length; k++) {
            stdDevs[k] = getStdDev(matrix[k], rowMeans[k]);
        }

        for (int i = 0; i < matrix.length; i++) {
            for (int j = 0; j < matrix[i].length; j++) {
                matrix[i][j] = (matrix[i][j] - rowMeans[i]) / stdDevs[i];
            }
        }

        for (int i = 0; i < matrix.length; i++) {
            for (int j = 0; j < matrix[i].length; j++) {
                matrix[i][j] = Math.min(Math.max(matrix[i][j], -threshold), threshold);
            }
        }

        return matrix;
    }

    private static float getStdDev(float[] data, float mean, int[] counts, int totalNumEntries) {
        double stddev = 0;

        for (int i = 0; i < data.length; i++) {
            float val = data[i];
            float diff = val - mean;
            stddev += counts[i] * diff * diff;
        }
        stddev = (stddev / totalNumEntries);

        return (float) Math.sqrt(stddev);
    }

    private static float getStdDev(float[] data, float mean) {
        double stddev = 0;

        for (float val : data) {
            float diff = val - mean;
            stddev += diff * diff;
        }
        stddev = (stddev / data.length);

        return (float) Math.sqrt(stddev);
    }

    public static float[][] log(float[][] matrix) {
        for (int i = 0; i < matrix.length; i++) {
            for (int j = 0; j < matrix[i].length; j++) {
                matrix[i][j] = (float) Math.log(matrix[i][j]);
            }
        }

        return matrix;
    }

    public static float[][] getMatrixModIndicesOfColumns(float[][] originalData, int modIndx, int base) {
        int numOrigColumns = originalData[0].length;
        int numModColumns = numOrigColumns / base + (numOrigColumns % base);
        float[][] modData = new float[originalData.length][numModColumns];

        for (int i = 0; i < originalData.length; i++) {
            int counter = 0;
            for (int j = modIndx; j < numOrigColumns; j += base) {
                modData[i][counter] = originalData[i][j];
                counter++;
            }
        }
        return modData;
    }

    public static float[][] getHalfOfMatrix(float[][] originalData, boolean getFirstHalf) {
        int numNewColumns = originalData[0].length / 2;
        float[][] modData = new float[originalData.length][numNewColumns];

        int offset = numNewColumns;
        if (getFirstHalf) {
            offset = 0;
        }

        for (int i = 0; i < originalData.length; i++) {
            for (int j = 0; j < numNewColumns; j++) {
                modData[i][j] = originalData[i][j + offset];
            }
        }
        return modData;
    }

    public static float[][] getRoundedLog(float[][] matrix) {
        float[][] matrix2 = new float[matrix.length][matrix[0].length];

        for (int i = 0; i < matrix.length; i++) {
            for (int j = 0; j < matrix[i].length; j++) {
                float val = (float) Math.round((float) 10 * Math.log(matrix[i][j])) / 10f;
                if (Float.isInfinite(val) || Float.isNaN(val)) {
                    val = 0;
                }
                matrix2[i][j] = val;
            }
        }
        return matrix2;
    }

    public static void scaleValuesByCount(float[][] matrix, int[] counts) {
        double[] countsSqrt = DoubleMatrixTools.sqrt(counts);

        for (int i = 0; i < matrix.length; i++) {
            for (int j = 0; j < countsSqrt.length; j++) {
                matrix[i][j] = (float) (matrix[i][j] * countsSqrt[j]);
            }
        }
    }

    public static void scaleValuesInPlaceByCountAndZscore(float[][] matrix, int[] counts) {

        int numTotalEntries = 0;
        for (int val : counts) {
            numTotalEntries += val;
        }

        float[] rowMeans = getRowMeans(matrix, counts, numTotalEntries);
        float[] rowStdDevs = getRowStandardDeviations(matrix, rowMeans, counts, numTotalEntries);
        double[] countsSqrt = DoubleMatrixTools.sqrt(counts);

        // zscore   (x-mu)/std
        for (int i = 0; i < matrix.length; i++) {
            for (int j = 0; j < matrix[i].length; j++) {
                matrix[i][j] = (float) (countsSqrt[j] * ((matrix[i][j] - rowMeans[i]) / rowStdDevs[i]));
            }
        }
    }

    private static float[] getRowStandardDeviations(float[][] matrix, float[] rowMeans, int[] counts, int numTotalEntries) {
        float[] rowStdDevs = new float[matrix.length];
        for (int k = 0; k < matrix.length; k++) {
            rowStdDevs[k] = getStdDev(matrix[k], rowMeans[k], counts, numTotalEntries);
        }
        return rowStdDevs;
    }

    public static float[][] fill(float[][] allDataForRegion, float val) {
        for (int i = 0; i < allDataForRegion.length; i++) {
            Arrays.fill(allDataForRegion[i], val);
        }
        return allDataForRegion;
    }

    public static void cleanUpNansInfinitesNegatives(float[][] matrix) {
        for (int i = 0; i < matrix.length; i++) {
            for (int j = 0; j < matrix[0].length; j++) {
                if (Float.isNaN(matrix[i][j]) || Float.isInfinite(matrix[i][j]) || matrix[i][j] < 1E-10) {
                    matrix[i][j] = 0;
                }
            }
        }
    }


    public float standardDeviation(float[][] data, float mean) {
        double stddev = 0;

        for (float[] vals : data) {
            for (float val : vals) {
                stddev += (val - mean) * (val - mean);
            }
        }
        stddev = (stddev / data.length) / data[0].length;

        return (float) Math.sqrt(stddev);
    }

    public void inPlaceZscore(float[][] data) {
        float mean = mean(data);
        float stddev = standardDeviation(data, mean);

        for (int i = 0; i < data.length; i++) {
            for (int j = 0; j < data[i].length; j++) {
                data[i][j] = (data[i][j] - mean) / stddev;
            }
        }
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

    /**
     * Reshape array into a matrix
     *
     * @param flatMatrix
     * @param numRows
     * @param numCols
     * @return properly dimensioned matrix
     */
    private static float[][] reshapeFlatMatrix(float[] flatMatrix, int numRows, int numCols) {
        float[][] matrix = new float[numRows][numCols];

        for (int i = 0; i < numRows; i++) {
            System.arraycopy(flatMatrix, i * numCols, matrix[i], 0, numCols);
        }
        return matrix;
    }

    /**
     * From Matrix M, extract out M[r1:r2,c1:c2]
     * r2, c2 not inclusive (~python numpy)
     *
     * @return extracted matrix region M[r1:r2,c1:c2]
     */
    public static float[][] extractLocalMatrixRegion(float[][] matrix, int r1, int r2, int c1, int c2) {

        int numRows = r2 - r1;
        int numColumns = c2 - c1;
        float[][] extractedRegion = new float[numRows][numColumns];

        for (int i = 0; i < numRows; i++) {
            System.arraycopy(matrix[r1 + i], c1, extractedRegion[i], 0, numColumns);
        }

        return extractedRegion;
    }

    /**
     * print for 2D array
     */
    private static void print(float[][] data) {
        for (float[] row : data) {
            System.out.println(Arrays.toString(row));
        }
    }





    public static float[] getAbsValColSums(float[][] matrix) {
        float[] colSum = new float[matrix[0].length];
        for (int i = 0; i < matrix.length; i++) {
            for (int j = 0; j < matrix[i].length; j++) {
                colSum[j] += Math.abs(matrix[i][j]);
            }
        }
        return colSum;
    }


    public static float[] getRowSums(float[][] matrix) {
        float[] rowSum = new float[matrix.length];
        for (int i = 0; i < matrix.length; i++) {
            for (float val : matrix[i]) {
                rowSum[i] += val;
            }
        }
        return rowSum;
    }

    public static float[] getRowMeans(float[][] matrix, int[] colCounts, int numTotalCounts) {
        float[] rowSum = new float[matrix.length];
        for (int i = 0; i < matrix.length; i++) {
            for (int j = 0; j < colCounts.length; j++) {
                rowSum[i] += (matrix[i][j] * colCounts[j]);
            }
        }
        for (int k = 0; k < rowSum.length; k++) {
            rowSum[k] = rowSum[k] / numTotalCounts;
        }
        return rowSum;
    }


    public static void cleanUpNaNs(float[][] matrix) {
        for (int r = 0; r < matrix.length; r++) {
            for (int c = 0; c < matrix[r].length; c++) {
                if (Float.isNaN(matrix[r][c])) {
                    matrix[r][c] = 0;
                }
            }
        }
    }

    public static float[][] transpose(float[][] matrix) {
        int h0 = matrix.length;
        int w0 = matrix[0].length;
        float[][] transposedMatrix = new float[w0][h0];

        for (int i = 0; i < h0; i++) {
            for (int j = 0; j < w0; j++) {
                transposedMatrix[j][i] = matrix[i][j];
            }
        }
        return transposedMatrix;
    }

    public static void copyFromAToBRegion(float[][] source, float[][] destination, int rowOffSet, int colOffSet) {
        for (int i = 0; i < source.length; i++) {
            System.arraycopy(source[i], 0, destination[i + rowOffSet], colOffSet, source[0].length);
        }
    }




    public static void saveMatrixTextV2(String filename, float[][] matrix) {
        Writer writer = null;
        try {
            writer = new BufferedWriter(new OutputStreamWriter(new FileOutputStream(filename), StandardCharsets.UTF_8));
            for (float[] row : matrix) {
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

    public static void saveMatrixTextNumpy(String filename, float[][] matrix) {
        int numRows = matrix.length;
        int numCols = matrix[0].length;
        float[] flattenedArray = FloatMatrixTools.flattenedRowMajorOrderMatrix(matrix);

        NpyFile.write(Paths.get(filename), flattenedArray, new int[]{numRows, numCols});
    }

    public static float[][] generateCompositeMatrixWithNansCleaned(float[][] matrixDiag1, float[][] matrixDiag2, float[][] matrix1vs2) {
        int newLength = matrixDiag1.length + matrixDiag2.length;
        float[][] compositeMatrix = new float[newLength][newLength];

        copyFromAToBRegion(matrixDiag1, compositeMatrix, 0, 0);
        copyFromAToBRegion(matrixDiag2, compositeMatrix, matrixDiag1.length, matrixDiag1.length);

        for (int i = 0; i < matrix1vs2.length; i++) {
            for (int j = 0; j < matrix1vs2[0].length; j++) {
                compositeMatrix[i][matrixDiag1.length + j] = matrix1vs2[i][j];
                compositeMatrix[matrixDiag1.length + j][i] = matrix1vs2[i][j];
            }
        }

        FloatMatrixTools.cleanUpNaNs(compositeMatrix);
        return compositeMatrix;
    }


    public static float[][] deepClone(float[][] data) {
        float[][] copy = new float[data.length][data[0].length];
        for (int i = 0; i < data.length; i++) {
            System.arraycopy(data[i], 0, copy[i], 0, data[i].length);
        }
        return copy;
    }

    public static float[][] getWithAppendedDerivative(float[][] data) {

        int numColumns = data[0].length;
        float[][] appendedDerivative = new float[data.length][2 * numColumns - 1];

        for (int i = 0; i < data.length; i++) {
            System.arraycopy(data[i], 0, appendedDerivative[i], 0, numColumns);
        }

        for (int i = 0; i < data.length; i++) {
            for (int j = 0; j < numColumns - 1; j++) {
                appendedDerivative[i][numColumns + j] = appendedDerivative[i][j] - appendedDerivative[i][j + 1];
            }
        }

        return appendedDerivative;
    }

    public static float[][] getNormalizedThresholdedAndAppendedDerivativeDownColumn(float[][] data, float maxVal, float scaleDerivFactor, float derivativeThreshold) {

        double[] averageVal = new double[data.length];
        for (int i = 0; i < data.length; i++) {
            for (float val : data[i]) {
                averageVal[i] += val;
            }
        }

        for (int i = 0; i < data.length; i++) {
            averageVal[i] = averageVal[i] / data[i].length;
        }

        float[][] thresholdedData = new float[data.length][data[0].length];
        for (int i = 0; i < data.length; i++) {
            for (int j = 0; j < data[i].length; j++) {
                thresholdedData[i][j] = (float) Math.min(maxVal, data[i][j] / averageVal[i]);
            }
        }

        return getMainAppendedDerivativeScaledPosDownColumn(thresholdedData, scaleDerivFactor, derivativeThreshold);
    }

    public static float[][] getNormalizedThresholdedByMedian(float[][] data, float maxVal) {

        double[] medianVal = new double[data.length];
        for (int i = 0; i < data.length; i++) {
            medianVal[i] = getMedian(data[i]);
        }

        float[][] thresholdedData = new float[data.length][data[0].length];
        for (int i = 0; i < data.length; i++) {
            for (int j = 0; j < data[i].length; j++) {
                thresholdedData[i][j] = (float) Math.min(maxVal, data[i][j] / medianVal[i]);
            }
        }

        return thresholdedData;
    }

    public static double getMedian(float[] values) {
        double[] array = new double[values.length];
        for (int k = 0; k < values.length; k++) {
            array[k] = values[k];
        }
        Median median = new Median();
        return median.evaluate(array);
    }


    public static float[][] getMainAppendedDerivativeScaledPosDownColumn(float[][] data, float scaleDerivFactor, float threshold) {

        int numColumns = data[0].length;
        float[][] derivative = getRelevantDerivativeScaledPositive(data, scaleDerivFactor, threshold);
        float[][] appendedDerivative = new float[data.length][numColumns + derivative[0].length];
        for (int i = 0; i < data.length; i++) {
            System.arraycopy(data[i], 0, appendedDerivative[i], 0, numColumns);
        }

        for (int i = 0; i < data.length; i++) {
            System.arraycopy(derivative[i], 0, appendedDerivative[i], numColumns, derivative[i].length);
        }

        return appendedDerivative;
    }

    public static float[][] getMainAppendedDerivativeDownColumnV2(float[][] data, float scaleDerivFactor, float threshold) {

        int numColumns = data[0].length;
        float[][] derivative = onlyGetRelevantDerivative(data, 1, threshold);
        float[][] appendedDerivative = new float[data.length][numColumns + derivative[0].length];
        for (int i = 0; i < data.length; i++) {
            System.arraycopy(data[i], 0, appendedDerivative[i], 0, numColumns);
        }

        for (int i = 0; i < data.length; i++) {
            for (int j = 0; j < data[i].length; j++) {
                appendedDerivative[i][j] = Math.min(.5f, Math.max(-.5f, appendedDerivative[i][j]));
            }
        }

        for (int i = 0; i < data.length; i++) {
            System.arraycopy(derivative[i], 0, appendedDerivative[i], numColumns, derivative[i].length);
        }

        return appendedDerivative;
    }

    public static float[][] getTrimmedMatrix(float[][] data) {
        List<Integer> indicesToUse = getImportantIndices(data);
        int n = indicesToUse.size();

        float[][] mainData = new float[data.length][n];

        for (int i = 0; i < data.length; i++) {
            for (int k = 0; k < n; k++) {
                int indexToUse = indicesToUse.get(k);
                mainData[i][k] = data[i][indexToUse];
            }
        }

        return mainData;
    }

    public static float[][] getTrimmedMatrixWithAppendedDerivativeDownColumn(float[][] data, float threshold) {
        List<Integer> indicesToUse = getImportantIndices(data);
        int n = indicesToUse.size();

        float[][] importantDataAndDerivative = new float[data.length][2 * n];

        for (int i = 0; i < data.length; i++) {
            for (int k = 0; k < n; k++) {
                int indexToUse = indicesToUse.get(k);
                importantDataAndDerivative[i][k] = Math.min(threshold, Math.max(-threshold, data[i][indexToUse]));
            }
        }

        for (int i = 0; i < data.length; i++) {
            for (int k = 0; k < n; k++) {
                int indexToUse = indicesToUse.get(k);
                importantDataAndDerivative[i][n + k] = Math.min(threshold, Math.max(-threshold, data[i][indexToUse] - data[i][indexToUse + 1]));
            }
        }

        return importantDataAndDerivative;
    }

    private static List<Integer> getImportantIndices(float[][] data) {
        int n = Math.max(data.length / 5, 1);
        float[][] derivative = new float[n][data[0].length - 1];
        for (int i = 0; i < n; i++) {
            for (int j = 0; j < data[0].length - 1; j++) {
                derivative[i][j] = data[i][j] - data[i][j + 1];
            }
        }

        float[] columnSums = getAbsValColSums(derivative);
        List<Integer> indicesToUse = new ArrayList<>();
        for (int k = 0; k < columnSums.length; k++) {
            if (columnSums[k] > 0) {
                indicesToUse.add(k);
            }
        }
        return indicesToUse;
    }

    public static float[][] getFullMatrixWithAppendedDerivative(float[][] data, float scaleDerivFactor, float threshold) {

        int numColumns = data[0].length;
        float[][] derivative = onlyGetRelevantDerivative(data, 1, threshold);
        float[][] appendedDerivative = new float[data.length][numColumns + derivative[0].length];
        for (int i = 0; i < data.length; i++) {
            System.arraycopy(data[i], 0, appendedDerivative[i], 0, numColumns);
        }

        for (int i = 0; i < data.length; i++) {
            System.arraycopy(derivative[i], 0, appendedDerivative[i], numColumns, derivative[i].length);
        }

        return appendedDerivative;
    }

    public static float[][] getRelevantDerivativeScaledPositive(float[][] data, float scaleDerivFactor, float threshold) {

        float[][] derivative = new float[data.length][data[0].length - 1];
        for (int i = 0; i < data.length; i++) {
            for (int j = 0; j < data[0].length - 1; j++) {
                derivative[i][j] = data[i][j] - data[i][j + 1];
            }
        }

        float[] columnSums = getAbsValColSums(derivative);
        List<Integer> indicesToUse = new ArrayList<>();
        for (int k = 0; k < columnSums.length; k++) {
            if (columnSums[k] > 0) {
                indicesToUse.add(k);
            }
        }

        float[][] importantDerivative = new float[data.length][indicesToUse.size()];

        for (int i = 0; i < data.length; i++) {
            for (int k = 0; k < indicesToUse.size(); k++) {
                int indexToUse = indicesToUse.get(k);
                importantDerivative[i][k] = Math.min(threshold, Math.max(-threshold, derivative[i][indexToUse] * scaleDerivFactor)) + threshold;
            }
        }

        return importantDerivative;
    }

    public static float[][] onlyGetFullDerivative(float[][] data, float threshold) {

        float[][] derivative = new float[data.length][data[0].length - 1];
        for (int i = 0; i < data.length; i++) {
            for (int j = 0; j < data[0].length - 1; j++) {
                derivative[i][j] = Math.min(threshold, Math.max(-threshold, data[i][j] - data[i][j + 1]));
            }
        }

        return derivative;
    }

    public static float[][] onlyGetRelevantDerivative(float[][] data, float scaleDerivFactor, float threshold) {

        float[][] derivative = new float[data.length][data[0].length - 1];
        for (int i = 0; i < data.length; i++) {
            for (int j = 0; j < data[0].length - 1; j++) {
                derivative[i][j] = data[i][j] - data[i][j + 1];
            }
        }

        float[] columnSums = getAbsValColSums(derivative);
        List<Integer> indicesToUse = new ArrayList<>();
        for (int k = 0; k < columnSums.length; k++) {
            if (columnSums[k] > 0) {
                indicesToUse.add(k);
            }
        }

        float[][] importantDerivative = new float[data.length][indicesToUse.size()];

        for (int i = 0; i < data.length; i++) {
            for (int k = 0; k < indicesToUse.size(); k++) {
                int indexToUse = indicesToUse.get(k);
                importantDerivative[i][k] = Math.min(threshold, Math.max(-threshold, derivative[i][indexToUse] * scaleDerivFactor));
            }
        }

        return importantDerivative;
    }

    public static float[] runSlidingAverageOnArray(int radius, float[] values) {

        float[] newValues = new float[values.length];
        for (int i = 0; i < values.length; i++) {
            float sum = 0;
            int numVals = 0;
            for (int j = Math.max(i - radius, 0); j < Math.min(i + radius, values.length); j++) {
                if (values[j] > 0) {
                    sum += values[j];
                    numVals++;
                }
            }
            if (numVals == 0) {
                newValues[i] = 0;
            } else {
                newValues[i] = sum / numVals;
            }

        }
        return newValues;
    }
}
