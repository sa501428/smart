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

package mixer.utils.slice;

import org.apache.commons.math.stat.regression.SimpleRegression;

import java.io.File;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;
import java.util.concurrent.atomic.AtomicInteger;

public class CorrelationTools {

	public static float getNonNanPearsonCorrelation(float[] array1, float[] array2) {
		SimpleRegression regression = new SimpleRegression();
		for (int i = 0; i < array1.length; i++) {
			boolean entryIsBad = Float.isNaN(array1[i]) || Float.isNaN(array2[i]);
			if (!entryIsBad) {
				regression.addData(array1[i], array2[i]);
			}
		}

		return (float) regression.getR();
	}

	public static float[][] getMinimallySufficientNonNanPearsonCorrelationMatrix(float[][] matrix, int numCentroids, File outputDirectory) {

		float[][] centroids = QuickCentroids.generateCentroids(matrix, numCentroids);
		float[][] result = new float[matrix.length][numCentroids]; // *2

		int numCPUThreads = 20;
		AtomicInteger currRowIndex = new AtomicInteger(0);
		ExecutorService executor = Executors.newFixedThreadPool(numCPUThreads);
		for (int l = 0; l < numCPUThreads; l++) {
			Runnable worker = new Runnable() {
				@Override
				public void run() {
					int i = currRowIndex.getAndIncrement();
					while (i < matrix.length) {
						
						for (int j = 0; j < centroids.length; j++) {
							float val = getNonNanPearsonCorrelation(matrix[i], centroids[j]);
							result[i][j] = arctanh(val);
							
							//result[i][centroids.length+j] = cosineSimilarity(matrix[i], centroids[j])*weight;
							//result[i][j] = cosineSimilarity(matrix[i], centroids[j]);
						}
						
						i = currRowIndex.getAndIncrement();
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
	
	private static float arctanh(float x) {
		float val = Math.max(x, -.99f);
		val = Math.min(val, .99f);
		val = (float) (Math.log(1 + val) - Math.log(1 - val)) / 2;
		if (Float.isInfinite(val)) {
			val = Float.NaN;
		}
		return val;
	}
	
	private static float cosineSimilarity(float[] vectorA, float[] vectorB) {
		double dotProduct = 0.0;
		double normA = 0.0;
		double normB = 0.0;
		for (int i = 0; i < vectorA.length; i++) {
			boolean entryIsBad = Float.isNaN(vectorA[i]) || Float.isNaN(vectorB[i]);
			if (!entryIsBad) {
				dotProduct += vectorA[i] * vectorB[i];
				normA += vectorA[i] * vectorA[i];
				normB += vectorB[i] * vectorB[i];
			}
		}
		return (float) (dotProduct / (Math.sqrt(normA) * Math.sqrt(normB)));
	}
}
