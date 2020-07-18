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

package mixer.commandline.utils.grind;

import java.io.*;
import java.nio.charset.StandardCharsets;
import java.util.Arrays;

public class GrindFloatMatrixTools {
	public static float[][] deepClone(float[][] data) {
		float[][] copy = new float[data.length][data[0].length];
		for (int i = 0; i < data.length; i++) {
			System.arraycopy(data[i], 0, copy[i], 0, data[i].length);
		}
		return copy;
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
		
		cleanUpNaNs(compositeMatrix);
		return compositeMatrix;
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
	
	public static void cleanUpNaNs(float[][] matrix) {
		for (int r = 0; r < matrix.length; r++) {
			for (int c = 0; c < matrix[r].length; c++) {
				if (Float.isNaN(matrix[r][c])) {
					matrix[r][c] = 0;
				}
			}
		}
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
}
