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

public class EfficientSubsampling {
    public static float[][] subsampleAtRandom(float[][] matrix, int numValues, long seed,
                                              int[] memberIndexes) {
        int numCols = matrix[0].length;
        float[][] result = new float[numValues][numCols];
        List<Integer> randomIndices = getRandomIndices(numValues, seed, memberIndexes);
        for (int i = 0; i < randomIndices.size(); i++) {
            System.arraycopy(matrix[randomIndices.get(i)], 0, result[i], 0, numCols);
        }
        return result;
    }

    private static List<Integer> getRandomIndices(int length, long seed, int[] memberIndexes) {
        Random generator = new Random(seed);
        Set<Integer> indices = new HashSet<>(length);
        while (indices.size() < length) {
            int actualIndex = memberIndexes[generator.nextInt(memberIndexes.length)];
            indices.add(actualIndex);
        }
        return new ArrayList<>(indices);
    }

    public static float[][] subsampleSpreadOut(float[][] matrix, int numValues, long seed,
                                               boolean useKMedians, int[] memberIndexes) {
        System.err.println("ALPHA " + numValues + " , " + memberIndexes.length);
        return new QuickCentroids(getTempSubMatrix(matrix, memberIndexes), numValues, seed,
                0).generateCentroids(0, useKMedians);
    }

    public static float[][] subsampleQuickClustering(float[][] matrix, int numValues, long seed,
                                                     boolean useKMedians, int[] memberIndexes) {
        System.err.println("BETA " + numValues + " , " + memberIndexes.length);
        return new QuickCentroids(getTempSubMatrix(matrix, memberIndexes), numValues,
                seed).generateCentroids(0, useKMedians);
    }

    private static float[][] getTempSubMatrix(float[][] matrix, int[] memberIndexes) {
        float[][] result = new float[memberIndexes.length][matrix[0].length];
        for (int i = 0; i < memberIndexes.length; i++) {
            int k0 = memberIndexes[i];
            System.arraycopy(matrix[k0], 0, result[i], 0, result[i].length);
        }
        return result;
    }


}
