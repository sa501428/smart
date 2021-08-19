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

package mixer.utils.slice.cleaning;

import javastraw.tools.MatrixTools;
import javastraw.tools.ParallelizedJuicerTools;
import mixer.utils.common.ArrayTools;
import mixer.utils.similaritymeasures.RobustCorrelationSimilarity;
import mixer.utils.similaritymeasures.RobustEuclideanDistance;
import mixer.utils.similaritymeasures.SimilarityMetric;
import org.apache.commons.math3.stat.regression.SimpleRegression;

import java.util.*;
import java.util.concurrent.atomic.AtomicInteger;

public class MatrixImputer {
    public static float[][] getImputedMatrix(float[][] initialData) {
        float[][] imputed = MatrixTools.deepClone(initialData);
        fillInImputedMatrix(imputed, initialData);
        return imputed;
    }

    private static void fillInImputedMatrix(float[][] imputed, float[][] initialData) {
        System.out.println("Imputing...");

        AtomicInteger index = new AtomicInteger(0);
        ParallelizedJuicerTools.launchParallelizedCode(() -> {
            int i = index.getAndIncrement();
            while (i < imputed.length) {
                BitSet yIsNan = getIsNan(imputed[i]);
                if (yIsNan.cardinality() > 0) {
                    int[] rowsToUse = getTwoNearestRows(yIsNan, imputed[i], initialData);
                    updateImputedMatrixEntries(yIsNan, imputed[i], initialData, rowsToUse);
                }
                i = index.getAndIncrement();
            }
        });
        System.out.println(".");
        System.out.println("Done imputing; Nans: " + numNans(initialData) + " -> " + numNans(imputed));
    }

    private static void updateImputedMatrixEntries(BitSet yIsNan, float[] imputed, float[][] data, int[] rowsToUse) {
        for (int j = 0; j < imputed.length; j++) {
            if (yIsNan.get(j)) {
                imputed[j] = getPredictedValue(rowsToUse, j, data);
            }
        }
    }

    private static float getPredictedValue(int[] rowIndices, int col, float[][] data) {
        float val = 0;
        int count = 0;
        for (int row : rowIndices) {
            if (!Float.isNaN(data[row][col])) {
                val += data[row][col];
                count++;
            }
        }
        return val / count;
    }

    private static Map<Integer, SimpleRegression> generateAllRegressions(float[] imputed, float[][] data,
                                                                         BitSet yIsNan, List<Integer> colsToUse) {
        Map<Integer, SimpleRegression> regressions = new HashMap<>();
        for (int k : colsToUse) {
            regressions.put(k, generateSingleRegression(data[k], imputed, yIsNan));
        }
        return regressions;
    }

    private static SimpleRegression generateSingleRegression(float[] x, float[] y, BitSet yIsNan) {
        SimpleRegression regression = new SimpleRegression();
        for (int j = 0; j < y.length; j++) {
            if (!yIsNan.get(j) && !Float.isNaN(x[j])) {
                regression.addData(x[j], y[j]);
            }
        }
        return regression;
    }

    private static BitSet getIsNan(float[] vector) {
        BitSet indexIsNan = new BitSet(vector.length);
        for (int i = 0; i < vector.length; i++) {
            if (Float.isNaN(vector[i])) {
                indexIsNan.set(i);
            }
        }
        return indexIsNan;
    }

    private static int[] getTwoNearestRows(BitSet yIsNan, float[] vector, float[][] initialData) {

        List<Integer> potentialIndicesToCheck = new ArrayList<>();
        int colToCheck = yIsNan.nextSetBit(0);
        if (colToCheck > -1) {
            for (int i = 0; i < initialData.length; i++) {
                if (Float.isNaN(initialData[i][colToCheck])) {
                    continue;
                }
                potentialIndicesToCheck.add(i);
            }
        }

        return getTwoNearestRowsFromList(potentialIndicesToCheck, vector, initialData);
    }

    private static int[] getTwoNearestRowsFromList(List<Integer> potentialIndicesToCheck, float[] vector, float[][] initialData) {
        int[] closestIndices = new int[2];
        float[] closestVals = new float[2];
        Arrays.fill(closestIndices, -1);
        Arrays.fill(closestVals, Float.MAX_VALUE);

        for (int i : potentialIndicesToCheck) {
            float dist = RobustEuclideanDistance.SINGLETON.distance(vector, initialData[i]);
            if (dist < closestVals[1]) {
                if (dist < closestVals[0]) {
                    closestVals[1] = closestVals[0];
                    closestIndices[1] = closestIndices[0];
                    closestVals[0] = dist;
                    closestIndices[0] = i;
                } else {
                    closestVals[1] = dist;
                    closestIndices[1] = i;
                }
            }
        }

        return closestIndices;
    }

    private static float[][] getR2Matrix(float[][] initialData) {
        float[][] r2 = new float[initialData.length][initialData.length];

        AtomicInteger index = new AtomicInteger(0);
        ParallelizedJuicerTools.launchParallelizedCode(() -> {
            SimilarityMetric metric = RobustCorrelationSimilarity.SINGLETON;
            int i = index.getAndIncrement();
            while (i < r2.length) {
                for (int j = i + 1; j < r2.length; j++) {
                    float r = metric.distance(initialData[i], initialData[j]);
                    r2[i][j] = r * r;
                    r2[j][i] = r2[i][j];
                }
                i = index.getAndIncrement();
            }
        });

        return r2;
    }

    private static String numNans(float[][] data) {
        int counter = 0;
        int[] numNanRows = new int[data.length];
        int[] numNanCols = new int[data[0].length];

        for (int i = 0; i < data.length; i++) {
            for (int j = 0; j < data[i].length; j++) {
                if (Float.isNaN(data[i][j])) {
                    counter++;
                    numNanRows[i]++;
                    numNanCols[j]++;
                }
            }
        }
        int maxNanInRows = ArrayTools.max(numNanRows);
        int maxNanInCols = ArrayTools.max(numNanCols);
        int avgNanInRows = ArrayTools.mean(numNanRows);
        int avgNanInCols = ArrayTools.mean(numNanCols);

        return "" + counter + "/(" + maxNanInRows + " : " + avgNanInRows + ")/(" + maxNanInCols + " : " + avgNanInCols + ")";
    }

    private static int getNumNans(float[][] data) {
        int counter = 0;

        for (int i = 0; i < data.length; i++) {
            for (int j = 0; j < data[i].length; j++) {
                if (Float.isNaN(data[i][j])) {
                    counter++;
                }
            }
        }

        return counter;
    }

    public static float[][] imputeUntilNoNansOnlyNN(float[][] data) {
        float[][] newData = data;
        do {
            newData = getImputedMatrix(newData);
        } while (checkIfHasNan(newData));
        return newData;
    }

    private static boolean checkIfHasNan(float[][] data) {
        for (int i = 0; i < data.length; i++) {
            for (int j = 0; j < data[i].length; j++) {
                if (Float.isNaN(data[i][j])) {
                    return true;
                }
            }
        }
        return false;
    }
}
