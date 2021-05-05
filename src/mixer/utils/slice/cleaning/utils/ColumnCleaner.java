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

package mixer.utils.slice.cleaning.utils;

import mixer.MixerGlobals;
import mixer.utils.common.FloatMatrixTools;
import mixer.utils.slice.matrices.MatrixAndWeight;

import java.util.Set;

public class ColumnCleaner extends DimensionCleaner {
    public ColumnCleaner(float[][] data, int[] weights) {
        super(data, weights);
    }

    @Override
    protected MatrixAndWeight filterOutBadIndices(Set<Integer> badIndices, float[][] matrix, int[] weights) {
        if (MixerGlobals.printVerboseComments) {
            System.out.println("interMatrix.length " + matrix.length + " badIndices.size() " + badIndices.size());
        }

        int counter = 0;
        int[] newIndexToOrigIndex = new int[matrix[0].length - badIndices.size()];
        for (int i = 0; i < matrix[0].length; i++) {
            if (!badIndices.contains(i)) {
                newIndexToOrigIndex[counter++] = i;
            }
        }

        float[][] newMatrix = new float[matrix.length][newIndexToOrigIndex.length];
        for (int i = 0; i < newMatrix.length; i++) {
            for (int j = 0; j < newIndexToOrigIndex.length; j++) {
                int origJ = newIndexToOrigIndex[j];
                newMatrix[i][j] = matrix[i][origJ];
            }
        }

        int[] newWeights = new int[newIndexToOrigIndex.length];
        for (int j = 0; j < newIndexToOrigIndex.length; j++) {
            int origJ = newIndexToOrigIndex[j];
            newWeights[j] = weights[origJ];
        }

        return new MatrixAndWeight(newMatrix, newWeights);
    }

    @Override
    protected int[] getNumberOfNansInDimension(float[][] matrix) {
        // invalid columns
        int[] numInvalids = new int[matrix[0].length];
        for (int i = 0; i < matrix.length; i++) {
            for (int j = 0; j < matrix[i].length; j++) {
                if (Float.isNaN(matrix[i][j]) || matrix[i][j] < ZERO) {
                    numInvalids[j]++;
                }
            }
        }
        return numInvalids;
    }

    @Override
    protected int getLimit() {
        return (int) Math.ceil(PERCENT_NAN_ALLOWED * data.length);
    }

    @Override
    protected boolean useOnlyCorr() {
        return true;
    }

    @Override
    protected float[][] getAppropriatelyFlippedMatrix() {
        return FloatMatrixTools.transpose(data);
    }
}
