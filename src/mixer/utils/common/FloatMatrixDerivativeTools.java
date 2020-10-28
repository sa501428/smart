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

import org.apache.commons.math3.stat.descriptive.rank.Median;

import java.util.ArrayList;
import java.util.List;

public class FloatMatrixDerivativeTools {

    public static float[][] getWithAppendedNonNanDerivative(float[][] data) {

        int numColumns = data[0].length;
        float[][] appendedDerivative = new float[data.length][2 * numColumns - 1];

        for (int i = 0; i < data.length; i++) {
            System.arraycopy(data[i], 0, appendedDerivative[i], 0, numColumns);
        }

        for (int i = 0; i < data.length; i++) {
            for (int j = 0; j < numColumns - 1; j++) {
                appendedDerivative[i][numColumns + j] = normalizedNonNanDerivative(appendedDerivative[i][j], appendedDerivative[i][j + 1]);
            }
        }

        return appendedDerivative;
    }

    private static float normalizedNonNanDerivative(float v0, float v1) {
        if (Float.isNaN(v0) || Float.isNaN(v1)) {
            return Float.NaN;
        }
        return 2 * (v1 - v0) / (v1 + v0 + 1e-5f);
    }

    public static void clearAndFillInSmoothDeriv(float[][] destination, float[][] input) {
        int[] factors = new int[]{1, 1, -1, -1};
        FloatMatrixTools.fill(destination, 0);
        for (int i = 0; i < input.length; i++) {
            for (int j = 0; j < input[i].length - factors.length + 1; j++) {
                for (int k = 0; k < factors.length; k++) {
                    destination[i][j] += factors[k] * input[i][j + k];
                }
            }
        }
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


    public static double getMedian(float[] values) {
        double[] array = new double[values.length];
        for (int k = 0; k < values.length; k++) {
            array[k] = values[k];
        }
        Median median = new Median();
        return median.evaluate(array);
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

    public static float[] getAbsValColSums(float[][] matrix) {
        float[] colSum = new float[matrix[0].length];
        for (int i = 0; i < matrix.length; i++) {
            for (int j = 0; j < matrix[i].length; j++) {
                colSum[j] += Math.abs(matrix[i][j]);
            }
        }
        return colSum;
    }
}
