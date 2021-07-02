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

import mixer.utils.slice.matrices.MatrixAndWeight;

import java.io.File;
import java.util.HashSet;
import java.util.Set;

public abstract class DimensionCleaner {

    protected static final float ZERO = 1e-10f;
    protected final float[][] data;
    protected final int[] weights;

    public DimensionCleaner(float[][] data, int[] weights) {
        this.data = data;
        this.weights = weights;
    }

    public MatrixAndWeight getCleanedData(int resolution, File outputDirectory) {
        return filterOutDimension(data, resolution, outputDirectory);
    }

    private MatrixAndWeight filterOutDimension(float[][] matrix, int resolution, File outputDirectory) {
        Set<Integer> badIndices = getSparseIndices(matrix);
        Set<Integer> outlierIndices = (new OutlierCleaner(getAppropriatelyFlippedMatrix(),
                useOnlyCorr())).getConsistentOutliers(resolution, outputDirectory);
        badIndices.addAll(outlierIndices);

        if (badIndices.size() == 0) {
            return new MatrixAndWeight(matrix, weights);
        }

        return filterOutBadIndices(badIndices, matrix, weights);
    }

    protected abstract boolean useOnlyCorr();

    protected abstract float[][] getAppropriatelyFlippedMatrix();

    protected Set<Integer> getSparseIndices(float[][] matrix) {
        // sparse rows
        Set<Integer> badIndices = new HashSet<>();
        int[] numNans = getNumberOfNansInDimension(matrix);
        int maxBadEntriesAllowed = getLimit();
        for (int i = 0; i < numNans.length; i++) {
            if (numNans[i] > maxBadEntriesAllowed) {
                badIndices.add(i);
            }
        }
        return badIndices;
    }

    protected abstract int getLimit();

    protected abstract int[] getNumberOfNansInDimension(float[][] matrix);

    protected abstract MatrixAndWeight filterOutBadIndices(Set<Integer> badIndices, float[][] matrix, int[] weights);

}
