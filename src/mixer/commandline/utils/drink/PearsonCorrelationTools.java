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

package mixer.commandline.utils.drink;

import org.apache.commons.math3.stat.regression.SimpleRegression;

import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;
import java.util.concurrent.atomic.AtomicInteger;

public class PearsonCorrelationTools {
	
	/**
	 * @param matrix
	 * @param rowNumDivisor i.e. if 1, get every row, if two skip every other, etc
	 * @return
	 */
	public static float[][] getNonNanPearsonCorrelationMatrix(float[][] matrix, int rowNumDivisor) {
		int numRows = matrix.length;
		float[][] result = new float[numRows][(numRows - 1) / rowNumDivisor + 1];
		
		int numCPUThreads = 20;
		AtomicInteger index = new AtomicInteger(0);
		ExecutorService executor = Executors.newFixedThreadPool(numCPUThreads);
		for (int l = 0; l < numCPUThreads; l++) {
			Runnable worker = new Runnable() {
				@Override
				public void run() {
					int i = index.getAndIncrement();
					while (i < numRows) {
						for (int j = i; j < numRows; j++) {
							if (j % rowNumDivisor == 0) {
								if (j == i) {
									result[i][j / rowNumDivisor] = 1;
								} else {
									float val = getNonNanPearsonCorrelation(matrix[i], matrix[j]);
									//System.out.println(i+" - "+j+" "+rowNumDivisor+" "+numRows+" "+matrix[0].length);
									result[i][j / rowNumDivisor] = val;
									if (i % rowNumDivisor == 0) {
										result[j][i / rowNumDivisor] = val;
									}
								}
							} else if (i % rowNumDivisor == 0) {
								float val = getNonNanPearsonCorrelation(matrix[i], matrix[j]);
								result[j][i / rowNumDivisor] = val;
							}
						}
						i = index.getAndIncrement();
						if (i % 100 == 0) {
							System.out.println("Completed: " + i * 100f / numRows + "%");
						}
					}
				}
			};
			executor.execute(worker);
		}
		executor.shutdown();
		
		// Wait until all threads finish
		while (!executor.isTerminated()) {
		}
		
		return result;
	}
	
	private static float getNonNanPearsonCorrelation(float[] array1, float[] array2) {
		SimpleRegression regression = new SimpleRegression();
		for (int i = 0; i < array1.length; i++) {
			if (!(Float.isNaN(array1[i]) || Float.isNaN(array2[i]))) {
				regression.addData(array1[i], array2[i]);
			}
		}
		
		return (float) regression.getR();
	}
	
	public static int[] getReSortedIndexOrder(float[][] matrix) {
		CorrelationBlockBuilder builder = new CorrelationBlockBuilder(matrix);
		return builder.getSortedIndexOrder();
	}
}
