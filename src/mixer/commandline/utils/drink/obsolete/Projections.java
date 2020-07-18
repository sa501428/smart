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

import mixer.MixerGlobals;
import mixer.commandline.utils.common.FloatMatrixTools;
import mixer.commandline.utils.drink.MatrixCleanup;
import mixer.commandline.utils.drink.SubcompartmentInterval;
import org.broad.igv.util.Pair;

import java.io.File;
import java.util.*;

public class Projections extends MatrixCleanup {
	
	private final static float sqrt3 = (float) Math.sqrt(3);
	private final static float cutoff_2_3 = 2f / 3f;
	private final static float cutoff_5_6 = 5f / 6f;
	private final static int[] factors = new int[]{1, 1, -1, -1};
	public static boolean USE_RANDOM3_TYPE = false;
	public static boolean USE_DERIV = true, DERIV_ON_ZSCORE = false, DERIV_THEN_ZSCORE = false;
	public static int reductionScalar = 4;
	
	public Projections(float[][] interMatrix, long seed, File outputDirectory) {
		super(interMatrix, seed, outputDirectory);
	}
	
	public static float[][] filterOutColumnsAndRows(float[][] interMatrix, Set<Integer> badIndices, Map<Integer, SubcompartmentInterval> original) {
		int[] newIndexToOrigIndex = new int[interMatrix.length - badIndices.size()];
		float[][] newMatrix = new float[newIndexToOrigIndex.length][newIndexToOrigIndex.length];
		
		int counter = 0;
		for (int i = 0; i < interMatrix.length; i++) {
			if (!badIndices.contains(i)) {
				newIndexToOrigIndex[counter++] = i;
			}
		}
		
		for (int i = 0; i < newMatrix.length; i++) {
			for (int j = 0; j < newMatrix.length; j++) {
				newMatrix[i][j] = interMatrix[newIndexToOrigIndex[i]][newIndexToOrigIndex[j]];
			}
		}
		
		Map<Integer, SubcompartmentInterval> newRowIndexToIntervalMap = new HashMap<>();
		for (int i = 0; i < newIndexToOrigIndex.length; i++) {
			newRowIndexToIntervalMap.put(i, (SubcompartmentInterval) original.get(newIndexToOrigIndex[i]).deepClone());
		}
		
		original.clear();
		original.putAll(newRowIndexToIntervalMap);
		
		return newMatrix;
	}
	
	private static void inPlaceMatrixMultiplication(float[][] inputMatrix, float[][] randomMatrix, float[][] resultMatrix,
													int resultColOffset, int inputColOffset) {
		for (int i = 0; i < inputMatrix.length; i++) {
			for (int z = 0; z < randomMatrix.length; z++) {
				for (int j = 0; j < randomMatrix[0].length; j++) {
					resultMatrix[i][resultColOffset + z] += inputMatrix[i][inputColOffset + j] * randomMatrix[z][j];
				}
			}
		}
	}
	
	private static float[][] createRandomMatrix(Random generator, int numRows, int numCols, boolean isDerivMatrix) {
		float[][] randMatrix = new float[numRows][numCols];
		
		for (int i = 0; i < numRows; i++) {
			for (int j = 0; j < numCols; j++) {
				randMatrix[i][j] = generator.nextFloat();
			}
		}
		
		if (USE_RANDOM3_TYPE) {
			for (int i = 0; i < numRows; i++) {
				for (int j = 0; j < numCols; j++) {
					float val = randMatrix[i][j];
					if (val < cutoff_2_3) {
						randMatrix[i][j] = 0;
					} else if (val < cutoff_5_6) {
						randMatrix[i][j] = sqrt3;
					} else {
						randMatrix[i][j] = -sqrt3;
					}
				}
			}
		} else {
			for (int i = 0; i < numRows; i++) {
				for (int j = 0; j < numCols; j++) {
					float val = randMatrix[i][j];
					if (val < .5) {
						randMatrix[i][j] = -1;
					} else {
						randMatrix[i][j] = 1;
					}
				}
			}
		}
		
		if (isDerivMatrix) {
			return createRandomDerivativeMatrix(randMatrix);
		}
		
		return randMatrix;
	}
	
	private static float[][] createRandomDerivativeMatrix(float[][] randMatrix) {
		float[][] randDerivMatrix = new float[randMatrix.length][randMatrix[0].length];
		
		for (int i = 0; i < randMatrix.length; i++) {
			for (int j = 0; j < randMatrix[0].length - 3; j++) {
				for (int k = 0; k < factors.length; k++) {
					randDerivMatrix[i][j + k] += randMatrix[i][j] * factors[k];
				}
			}
		}
		
		return randDerivMatrix;
	}
	
	private static Pair<int[], int[]> fixDimensions(int[][] dimensions, Set<Integer> badIndices) {
		
		
		int[] offsets = dimensions[0];
		int[] widthOfRegion = dimensions[1];
		
		for (int badIndex : badIndices) {
			for (int i = 0; i < offsets.length; i++) {
				if (i == offsets.length - 1 && badIndex >= offsets[i]) {
					widthOfRegion[i]--;
					break;
				} else if (badIndex >= offsets[i] && badIndex < offsets[i + 1]) {
					widthOfRegion[i]--;
					break;
				}
			}
		}
		
		for (int k = 1; k < offsets.length; k++) {
			offsets[k] = offsets[k - 1] + widthOfRegion[k - 1];
		}
		
		return new Pair<>(offsets, widthOfRegion);
	}
	
	private static void movAvgOnMatrix(float[][] data) {
		for (int i = 0; i < data.length; i++) {
			float[] rowCopy = new float[data[i].length];
			System.arraycopy(data[i], 0, rowCopy, 0, data[i].length);
			Arrays.fill(data[i], 0);
			
			for (int j = 0; j < rowCopy.length; j++) {
				
				if (Float.isNaN(rowCopy[j])) {
					data[i][j] = Float.NaN;
				} else {
					int kStart = -1, kEnd = 2;
					if (j == 0) {
						kStart = 0;
					} else if (j == rowCopy.length - 1) {
						kEnd = 1;
					}
					
					int divisor = 0;
					for (int k = kStart; k < kEnd; k++) {
						float val = rowCopy[j + k];
						if (!Float.isNaN(val)) {
							divisor++;
							data[i][j] += val;
						}
					}
					
					data[i][j] /= divisor;
				}
			}
		}
	}
	
	public float[][] getCleanedMatrix(Map<Integer, SubcompartmentInterval> rowIndexToIntervalMap, int[][] origDimensions) {
		
		thresholdByZscoreToNanDownColumn(data, zScoreThreshold, 1);
		Set<Integer> badIndices = getBadIndices(data);
		data = filterOutColumnsAndRows(data, badIndices, rowIndexToIntervalMap);
		System.out.println("matrix size " + data.length + " x " + data[0].length);
		
		Pair<int[], int[]> dimensionParams = fixDimensions(origDimensions, badIndices);
		if (MixerGlobals.printVerboseComments) {
			FloatMatrixTools.saveMatrixTextNumpy(new File(outputDirectory, "post_filter_data_matrix.npy").getAbsolutePath(), data, dimensionParams.getFirst());
		}
		
		return randomProjection(data, dimensionParams, USE_DERIV, generator); // dimensions.getSecond()
	}
	
	private float[][] randomProjection(float[][] interMatrix, Pair<int[], int[]> dimensionParams, boolean includeSmoothDerivative, Random generator) {
		
		int[] offsets = dimensionParams.getFirst();
		int[] widthOfRegion = dimensionParams.getSecond();
		
		// create holder for new matrix
		int newWidth = 0;
		for (int width : widthOfRegion) {
			int reducedWidth = width / reductionScalar;
			if (includeSmoothDerivative) {
				reducedWidth *= 2;
			}
			newWidth += reducedWidth;
		}
		float[][] reducedMatrix = new float[interMatrix.length][newWidth];
		
		System.out.println("matrix size " + interMatrix.length + " x " + interMatrix[0].length);
		System.out.println("reduced matrix size " + reducedMatrix.length + " x " + reducedMatrix[0].length);
		
		// clone original matrix
		//float[][] clone = cloneAndMovAvgOnMatrix(interMatrix);
		//interMatrix = null;
		FloatMatrixTools.thresholdInPlaceByZscoreDownCols(interMatrix, zScoreThreshold, 1);
		movAvgOnMatrix(interMatrix);
		
		if (MixerGlobals.printVerboseComments) {
			//FloatMatrixTools.saveMatrixTextNumpy(new File(outputDirectory, "smooth_matrix.npy").getAbsolutePath(), interMatrix);
			System.gc();
		}
		
		if (includeSmoothDerivative) {
			if (DERIV_THEN_ZSCORE) {
				//clone = cloneAndMovAvgOnMatrix(interMatrix);
				//FloatMatrixTools.clearAndFillInSmoothDeriv(interMatrix, clone);
				//FloatMatrixTools.thresholdInPlaceByZscoreDownCols(interMatrix, zScoreThreshold);
				//FloatMatrixTools.inPlaceZscoreDownColsNoNan(interMatrix);
			}
			
			int reducedMatrixOffset = 0;
			for (int k = 0; k < offsets.length; k++) {
				int reducedWidth = widthOfRegion[k] / reductionScalar;
				reducedMatrixOffset += reducedWidth;
				
				if (DERIV_THEN_ZSCORE) {
					//float[][] randomMatrix = createRandomMatrix(generator, reducedWidth, widthOfRegion[k], false);
					//inPlaceMatrixMultiplication(interMatrix, randomMatrix, reducedMatrix, reducedMatrixOffset, offsets[k]);
				} else if (DERIV_ON_ZSCORE) {
					//float[][] randomMatrixDeriv = createRandomMatrix(generator, reducedWidth, widthOfRegion[k], true);
					//inPlaceMatrixMultiplication(clone, randomMatrixDeriv, reducedMatrix, reducedMatrixOffset, offsets[k]);
				} else {
					float[][] randomMatrixDeriv = createRandomMatrix(generator, reducedWidth, widthOfRegion[k], true);
					inPlaceMatrixMultiplication(interMatrix, randomMatrixDeriv, reducedMatrix, reducedMatrixOffset, offsets[k]);
				}
				reducedMatrixOffset += reducedWidth;
			}
		}
		
		//float[][] clone = FloatMatrixTools.deepClone(interMatrix);
		
		FloatMatrixTools.inPlaceZscoreDownColsNoNan(interMatrix, 1);
		
		if (MixerGlobals.printVerboseComments) {
			//FloatMatrixTools.saveMatrixTextNumpy(new File(outputDirectory, "smooth_zscore_matrix.npy").getAbsolutePath(), interMatrix);
			System.gc();
		}
		
		int reducedMatrixOffset = 0;
		for (int k = 0; k < offsets.length; k++) {
			int reducedWidth = widthOfRegion[k] / reductionScalar;
			float[][] randomMatrix = createRandomMatrix(generator, reducedWidth, widthOfRegion[k], false);
			inPlaceMatrixMultiplication(interMatrix, randomMatrix, reducedMatrix, reducedMatrixOffset, offsets[k]);
			reducedMatrixOffset += reducedWidth;
			
			if (includeSmoothDerivative) {
				reducedMatrixOffset += reducedWidth;
			}
		}
		
		if (MixerGlobals.printVerboseComments) {
			FloatMatrixTools.saveMatrixTextNumpy(new File(outputDirectory, "projected_matrix.npy").getAbsolutePath(), reducedMatrix);
			System.gc();
		}
		
		return reducedMatrix;
	}
}
