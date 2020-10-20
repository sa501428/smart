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

import mixer.utils.common.FloatMatrixTools;

import java.io.File;
import java.util.*;

public class MatrixCleanup {
    private final static float PERCENT_ZERO_ALLOWED = .5f;
    public static final int BATCHED_NUM_ROWS = 1; // 10
    protected static final float zScoreThreshold = 3f;
    protected final File outputDirectory;
    protected final Random generator;

    private final static float PERCENT_NAN_ALLOWED = .5f;
    protected float[][] data;

    public MatrixCleanup(float[][] interMatrix, long seed, File outputDirectory) {
        data = interMatrix;
        System.out.println("matrix size " + data.length + " x " + data[0].length);
        generator = new Random(seed);
        this.outputDirectory = outputDirectory;
    }

    public static float[][] filterOutColumnsAndRowsNonSymmetricMatrix(float[][] interMatrix, Map<Integer, SubcompartmentInterval> original) {
        Set<Integer> badIndices = getBadIndices(interMatrix);
        if (badIndices.size() == 0) {
            return interMatrix;
        }

        System.out.println("interMatrix.length " + interMatrix.length + " badIndices.size() " + badIndices.size());

        int counter = 0;
        int[] newIndexToOrigIndex = new int[interMatrix.length - badIndices.size()];
        for (int i = 0; i < interMatrix.length; i++) {
            if (!badIndices.contains(i)) {
                newIndexToOrigIndex[counter++] = i;
            }
        }

        float[][] newMatrix = new float[newIndexToOrigIndex.length][interMatrix[0].length];
        Map<Integer, SubcompartmentInterval> newRowIndexToIntervalMap = new HashMap<>();
        for (int i = 0; i < newMatrix.length; i++) {
            int tempI = newIndexToOrigIndex[i];
            System.arraycopy(interMatrix[tempI], 0, newMatrix[i], 0, newMatrix[0].length);
            newRowIndexToIntervalMap.put(i, (SubcompartmentInterval) original.get(newIndexToOrigIndex[i]).deepClone());
        }

        original.clear();
        original.putAll(newRowIndexToIntervalMap);

        return newMatrix;
    }
    
    public static Set<Integer> getBadIndices(float[][] matrix) {
        
        Set<Integer> badIndices = new HashSet<>();
        
        int[] numNans = getNumberOfNansInRow(matrix);
        int[] numZerosNotNans = getNumberZerosNotNansInRow(matrix);
        
        int[] numNonNan = new int[numNans.length];
        int n = matrix[0].length;

        for (int i = 0; i < numNans.length; i++) {
            numNonNan[i] = n - numNans[i];
            if (((float) numNans[i]) / ((float) n) > PERCENT_NAN_ALLOWED) {
                badIndices.add(i);
            }
    
            if (((float) numZerosNotNans[i]) / ((float) numNonNan[i]) > PERCENT_ZERO_ALLOWED) {
                badIndices.add(i);
            }
        }

        return badIndices;
    }

    private static int[] getNumberZerosNotNansInRow(float[][] matrix) {
        int[] numZeros = new int[matrix.length];
        for (int i = 0; i < matrix.length; i++) {
            for (int j = 0; j < matrix[i].length; j++) {
                if (!Float.isNaN(matrix[i][j]) && matrix[i][j] < 1e-5) {
                    numZeros[i]++;
                }
            }
        }
        return numZeros;
    }

    private static int[] getNumberOfNansInRow(float[][] matrix) {
        int[] numNans = new int[matrix.length];
        for (int i = 0; i < matrix.length; i++) {
            for (int j = 0; j < matrix[i].length; j++) {
                if (Float.isNaN(matrix[i][j])) {
                    numNans[i]++;
                }
            }
        }
        return numNans;
    }

    public float[][] getSimpleCleaningOfMatrixAppendCorr(Map<Integer, SubcompartmentInterval> rowIndexToIntervalMap, boolean useCorrelation) {
        FloatMatrixTools.thresholdByZscoreToNanDownColumn(data, zScoreThreshold, BATCHED_NUM_ROWS);

        data = filterOutColumnsAndRowsNonSymmetricMatrix(data, rowIndexToIntervalMap);
        System.out.println("matrix size " + data.length + " x " + data[0].length);

        FloatMatrixTools.inPlaceZscoreDownColsNoNan(data, BATCHED_NUM_ROWS);

        if (useCorrelation) {
            data = CorrelationTools.getMinimallySufficientNonNanPearsonCorrelationMatrix(data, data[0].length / 20, outputDirectory);
        }
        return data;
    }
}
