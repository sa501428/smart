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

package mixer.commandline.utils.drink.obsolete;

import mixer.commandline.utils.common.DoubleMatrixTools;
import mixer.commandline.utils.common.FloatMatrixTools;
import mixer.commandline.utils.grind.GrindFloatMatrixTools;

import java.util.Arrays;

public class OldTools {
	
	
	public static float[][] inPlaceDerivAndThresholdDownCols(float[][] matrix, float threshold) {
		float[][] concatenatedMatrix = getFullMatrixWithAppendedSmoothDerivative(matrix);
		FloatMatrixTools.thresholdInPlaceByZscoreDownCols(concatenatedMatrix, threshold, 1);
		return concatenatedMatrix;
	}
	
	public static float[][] getFullMatrixWithAppendedSmoothDerivative(float[][] data) {
		
		int numColumns = data[0].length;
		float[][] appendedDerivative = new float[data.length][2 * numColumns - 3];
		for (int i = 0; i < data.length; i++) {
			System.arraycopy(data[i], 0, appendedDerivative[i], 0, numColumns);
		}
		
		for (int i = 0; i < data.length; i++) {
			for (int j = 0; j < numColumns - 3; j++) {
				appendedDerivative[i][numColumns + j] = data[i][j] + data[i][j + 1] - data[i][j + 2] - data[i][j + 3];
			}
		}
		
		return appendedDerivative;
	}
	
	public static void inPlaceZscoreDownRowsNoNan(float[][] matrix) {
		float[] rowMeans = getRowMeansNonNan(matrix);
		float[] rowStdDevs = getRowStdDevNonNans(matrix, rowMeans);
		
		for (int i = 0; i < matrix.length; i++) {
			for (int j = 0; j < matrix[i].length; j++) {
				float val = matrix[i][j];
				if (!Float.isNaN(val)) {
					matrix[i][j] = (val - rowMeans[i]) / rowStdDevs[i];
				}
			}
		}
	}
	
	
	public static float[][] inPlaceZscoreDownRows(float[][] matrix, float threshold) {
		float[] rowMeans = GrindFloatMatrixTools.getRowSums(matrix);
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
	
	public static void cleanUpNansInfinitesNegatives(float[][] matrix) {
		for (int i = 0; i < matrix.length; i++) {
			for (int j = 0; j < matrix[0].length; j++) {
				if (Float.isNaN(matrix[i][j]) || Float.isInfinite(matrix[i][j]) || matrix[i][j] < 1E-10) {
					matrix[i][j] = 0;
				}
			}
		}
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
	
	private static float[] getRowStandardDeviations(float[][] matrix, float[] rowMeans, int[] counts, int numTotalEntries) {
		float[] rowStdDevs = new float[matrix.length];
		for (int k = 0; k < matrix.length; k++) {
			rowStdDevs[k] = getStdDev(matrix[k], rowMeans[k], counts, numTotalEntries);
		}
		return rowStdDevs;
	}
	
	public static void thresholdByZscoreToNanDownRow(float[][] matrix, float threshold) {
		float[] rowMeans = getRowMeansNonNan(matrix);
		float[] rowStdDevs = getRowStdDevNonNans(matrix, rowMeans);
		
		for (int i = 0; i < matrix.length; i++) {
			for (int j = 0; j < matrix[i].length; j++) {
				float val = matrix[i][j];
				if (!Float.isNaN(val)) {
					float newVal = (val - rowMeans[i]) / rowStdDevs[i];
					if (newVal > threshold || newVal < -threshold) {
						matrix[i][j] = Float.NaN;
					}
				}
			}
		}
	}
	
	public static float[] getRowMeansNonNan(float[][] matrix) {
		float[] rowMeans = new float[matrix.length];
		int[] rowNonNans = new int[matrix.length];
		for (int i = 0; i < matrix.length; i++) {
			for (int j = 0; j < matrix[i].length; j++) {
				float val = matrix[i][j];
				if (!Float.isNaN(val)) {
					rowMeans[i] += val;
					rowNonNans[i] += 1;
				}
			}
		}
		
		for (int k = 0; k < rowMeans.length; k++) {
			rowMeans[k] = rowMeans[k] / Math.max(rowNonNans[k], 1);
		}
		
		return rowMeans;
	}
	
	public static float[] getRowStdDevNonNans(float[][] matrix, float[] means) {
		
		float[] stdDevs = new float[means.length];
		int[] rowNonNans = new int[means.length];
		
		for (int i = 0; i < matrix.length; i++) {
			for (int j = 0; j < matrix[i].length; j++) {
				float val = matrix[i][j];
				if (!Float.isNaN(val)) {
					float diff = val - means[i];
					stdDevs[i] += diff * diff;
					rowNonNans[i] += 1;
				}
			}
		}
		
		for (int k = 0; k < stdDevs.length; k++) {
			stdDevs[k] = (float) Math.sqrt(stdDevs[k] / Math.max(rowNonNans[k], 1));
		}
		
		return stdDevs;
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
	
	public void inPlaceZscore(float[][] data) {
		float mean = mean(data);
		float stddev = standardDeviation(data, mean);
		
		for (int i = 0; i < data.length; i++) {
			for (int j = 0; j < data[i].length; j++) {
				data[i][j] = (data[i][j] - mean) / stddev;
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
	
	
}
