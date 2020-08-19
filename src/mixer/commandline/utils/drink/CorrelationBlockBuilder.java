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

import java.util.*;

public class CorrelationBlockBuilder {
	private final float[][] matrix;
	private final float threshold = .5f;
	private final Map<Integer, Set<Integer>> nearbyNeightbors = new HashMap<>();
	private final boolean[] indexStatusNotDone;
	private final List<Set<Integer>> ordering = new ArrayList<>();
	private int latestUnseenIndex = 0;
	
	public CorrelationBlockBuilder(float[][] matrix) {
		this.matrix = matrix;
		indexStatusNotDone = new boolean[matrix.length];
		Arrays.fill(indexStatusNotDone, true);
		
		for (int j = 0; j < matrix[0].length; j++) {
			Set<Integer> batch = new HashSet<>();
			for (int i = 0; i < matrix.length; i++) {
				if (matrix[i][j] > threshold) {
					batch.add(i);
				}
			}
			
			if (batch.size() > 0) {
				System.out.println("batch size " + batch.size());
			}
			for (int indx : batch) {
				nearbyNeightbors.getOrDefault(indx, new HashSet<>()).addAll(batch);
				//System.out.println("neighbors :"+nearbyNeightbors.getOrDefault(indx,new HashSet<>()).size());
			}
		}
	}
	
	public int[] getSortedIndexOrder() {
		int currIndex = getLatestIndex();
		Set<Integer> currentBatch = new HashSet<>();
		Queue<Integer> queue = new LinkedList<>();
		
		while (currIndex > -1) {
			currentBatch.add(currIndex);
			for (Integer nextIndx : nearbyNeightbors.getOrDefault(currIndex, new HashSet<>())) {
				currentBatch.add(nextIndx);
				if (indexStatusNotDone[nextIndx]) {
					queue.add(nextIndx);
				}
			}
			indexStatusNotDone[currIndex] = false;
			if (queue.size() > 0) {
				currIndex = queue.poll();
			} else {
				// start new set
				if (currentBatch.size() > 0) {
					System.out.println("Size " + currentBatch.size());
					ordering.add(currentBatch);
					currentBatch = new HashSet<>();
				}
				currIndex = getLatestIndex();
			}
		}
		if (currentBatch.size() > 0) {
			ordering.add(currentBatch);
		}
		
		System.out.println("Estimating " + ordering.size() + " clusters may be reasonable ¯\\_(ツ)_/¯");
		
		int[] newOrder = new int[matrix.length];
		int counter = 0;
		for (Set<Integer> batch : ordering) {
			for (int indx : batch) {
				newOrder[counter] = indx;
				counter++;
			}
		}
		
		if (counter != matrix.length) {
			System.err.println("Count doesn't add up Matrix: " + matrix.length + " Count:" + counter);
		}
		
		return newOrder;
	}
	
	private int getLatestIndex() {
		for (int k = latestUnseenIndex; k < matrix.length; k++) {
			if (indexStatusNotDone[k]) {
				latestUnseenIndex = k;
				return k;
			}
		}
		return -1;
	}
	
	
}
