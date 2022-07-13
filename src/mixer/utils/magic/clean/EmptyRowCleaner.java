/*
 * The MIT License (MIT)
 *
 * Copyright (c) 2011-2022 Rice University, Baylor College of Medicine, Aiden Lab
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

package mixer.utils.magic.clean;

import mixer.utils.common.ArrayTools;
import mixer.utils.slice.cleaning.utils.MatrixRowCleaner;
import mixer.utils.slice.structures.SubcompartmentInterval;

import java.util.HashSet;
import java.util.Map;
import java.util.Set;

public class EmptyRowCleaner {

    private static final int MIN_ZSCORE_CUTOFF = -2;

    public static float[][] cleanUpMatrix(float[][] data, Map<Integer, SubcompartmentInterval> rowIndexToIntervalMap,
                                          float[] coverage) {
        Set<Integer> badIndices = getBadIndices(data);
        System.out.println("initial magic matrix num rows: " + data.length + " badIndices: " + badIndices.size());

        addRowsWithBadCoverage(badIndices, coverage);
        System.out.println("initial magic matrix num rows: " + data.length + " badIndices+badCoverage: " + badIndices.size());

        return MatrixRowCleaner.makeNewMatrixAndUpdateIndices(data, rowIndexToIntervalMap, badIndices);
    }

    private static void addRowsWithBadCoverage(Set<Integer> badIndices, float[] coverage) {
        float mean = ArrayTools.getNonZeroMean(coverage);
        float stdDev = ArrayTools.getNonZeroStd(coverage, mean);
        addRowsWithBadCoverage(coverage, mean, stdDev, badIndices);
    }

    private static Set<Integer> getBadIndices(float[][] matrix) {
        int minGoodColsRequired = Math.max(1, matrix[0].length / 2);
        Set<Integer> badIndices = new HashSet<>();
        for (int i = 0; i < matrix.length; i++) {
            if (isBadRow(matrix[i], minGoodColsRequired)) {
                badIndices.add(i);
            }
        }
        return badIndices;
    }

    private static boolean isBadRow(float[] row, int limit) {
        int numGoodEntries = 0;
        for (float val : row) {
            if (val > 0) {
                numGoodEntries++;
            }
        }
        return numGoodEntries < limit;
    }

    private static void addRowsWithBadCoverage(float[] coverage, double mean, double stdDev,
                                               Set<Integer> badIndices) {
        for (int k = 0; k < coverage.length; k++) {
            if (coverage[k] > 0) {
                double zval = (coverage[k] - mean) / stdDev;
                if (zval < MIN_ZSCORE_CUTOFF) {
                    badIndices.add(k);
                }
            } else {
                badIndices.add(k);
            }
        }
    }

    private void removeExtremeCoverages(float[] coverage) {

    }
}
