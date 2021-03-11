/*
 * The MIT License (MIT)
 *
 * Copyright (c) 2011-2021 Rice University, Baylor College of Medicine, Aiden Lab
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

package mixer.utils.umap;

import mixer.utils.common.FloatMatrixTools;

import java.util.Arrays;

public class ZeroRowRemover {

    private final int[] newPosition;
    private final float[][] cleanData;
    private final int[][] cleanIndices;

    public ZeroRowRemover(float[][] data, int[][] initialIndexToIDs) {
        newPosition = new int[data.length];

        int numToKeep = getNewIndices(FloatMatrixTools.getNumZerosInRow(data),
                data[0].length / 2);

        cleanData = cleanUpEmptyRows(data, numToKeep);
        cleanIndices = cleanUpIndices(initialIndexToIDs, numToKeep);
    }

    private float[][] cleanUpEmptyRows(float[][] data, int numToKeep) {
        float[][] newMatrix = new float[numToKeep][data[0].length];
        for (int i = 0; i < data.length; i++) {
            if (newPosition[i] > -1) {
                System.arraycopy(data[i], 0, newMatrix[newPosition[i]], 0, data[i].length);
            }
        }
        return newMatrix;
    }

    private int[][] cleanUpIndices(int[][] indices, int numToKeep) {
        int[][] newIndices = new int[numToKeep][indices[0].length];
        for (int i = 0; i < indices.length; i++) {
            if (newPosition[i] > -1) {
                System.arraycopy(indices[i], 0, newIndices[newPosition[i]], 0, indices[i].length);
            }
        }
        return newIndices;
    }

    private int getNewIndices(int[] numZerosInRow, int threshold) {
        Arrays.fill(newPosition, -1);

        int counter = 0;
        for (int k = 0; k < numZerosInRow.length; k++) {
            if (numZerosInRow[k] < threshold) {
                // good row
                newPosition[k] = counter;
                counter++;
            }
        }

        return counter;
    }

    public float[][] getCleanData() {
        return cleanData;
    }

    public int[][] getCleanIndices() {
        return cleanIndices;
    }
}
