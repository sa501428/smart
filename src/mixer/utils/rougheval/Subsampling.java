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

package mixer.utils.rougheval;

import mixer.utils.slice.cleaning.QuickCentroids;

import java.util.*;

public class Subsampling {
    public float[][] subsampleAtRandom(float[][] matrix, int numValues, long seed) {
        int numCols = matrix[0].length;
        float[][] result = new float[numValues][numCols];
        List<Integer> randomIndices = getRandomIndices(numValues, seed, matrix.length);
        for (int i = 0; i < randomIndices.size(); i++) {
            System.arraycopy(matrix[randomIndices.get(i)], 0, result[i], 0, numCols);
        }
        return result;
    }

    private List<Integer> getRandomIndices(int length, long seed, int bound) {
        Random generator = new Random(seed);
        Set<Integer> indices = new HashSet<>(length);
        while (indices.size() < length) {
            indices.add(generator.nextInt(bound));
        }
        return new ArrayList<>(indices);
    }

    public float[][] subsampleSpreadOut(float[][] matrix, int numValues, long seed, boolean useKMedians) {
        return new QuickCentroids(matrix, numValues, seed, 0).generateCentroids(0);
    }

    public float[][] subsampleQuickClustering(float[][] matrix, int numValues, long seed, boolean useKMedians) {
        return new QuickCentroids(matrix, numValues, seed).generateCentroids(0);
    }


}
