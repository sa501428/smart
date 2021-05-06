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
import mixer.utils.similaritymeasures.SimilarityMetric;
import org.apache.commons.math3.stat.regression.SimpleRegression;

import java.util.*;
import java.util.concurrent.atomic.AtomicInteger;

public class MatrixImputer {
    private static float R2_CUTOFF = 0.5f;
    private static boolean DO_TRANSPOSE = false;

    public static float[][] getImputedMatrix(float[][] data, boolean doTranspose) {
        DO_TRANSPOSE = doTranspose;
        float[][] initialData = data;
        if (doTranspose) {
            initialData = MatrixTools.transpose(data);
        }
        float[][] imputed = MatrixTools.deepClone(initialData);
        fillInImputedMatrix(imputed, initialData);
        if (doTranspose) {
            return MatrixTools.transpose(imputed);
        }
        return imputed;
    }

    private static void fillInImputedMatrix(float[][] imputed, float[][] initialData) {
        float[][] r2Matrix = getR2Matrix(initialData);
        System.out.println("Imputing...");

        AtomicInteger index = new AtomicInteger(0);
        ParallelizedJuicerTools.launchParallelizedCode(() -> {
            int i = index.getAndIncrement();
            while (i < imputed.length) {
                BitSet yIsNan = getIsNan(imputed[i]);
                if (yIsNan.cardinality() > 0) {
                    List<Integer> colsToUse = getColumnsWithStrongCorr(r2Matrix[i]);
                    if (colsToUse.size() > 0) {
                        updateImputedMatrixEntries(yIsNan, imputed[i], initialData, colsToUse);
                        //System.out.print(".");
                    }
                }
                i = index.getAndIncrement();
            }
        });
        System.out.println(".");
        System.out.println("Done imputing; Nans: " + numNans(initialData) + " -> " + numNans(imputed));
    }

    private static void updateImputedMatrixEntries(BitSet yIsNan, float[] imputed, float[][] data, List<Integer> colsToUse) {
        Map<Integer, SimpleRegression> regressions = generateAllRegressions(imputed, data, yIsNan, colsToUse);

        for (int j = 0; j < imputed.length; j++) {
            if (yIsNan.get(j)) {
                imputed[j] = getPredictedValue(regressions, j, data);
            }
        }
    }

    private static float getPredictedValue(Map<Integer, SimpleRegression> regressions, int col, float[][] data) {
        double accum = 0;
        double totalWeights = 0;
        for (Integer row : regressions.keySet()) {
            if (!Float.isNaN(data[row][col])) {
                SimpleRegression regression = regressions.get(row);
                double newY = regression.predict(data[row][col]);
                double weight = regression.getRSquare();
                accum += (newY * weight);
                totalWeights += weight;
            }
        }

        if (totalWeights > .1) {
            return (float) (accum / totalWeights);
        }
        return Float.NaN;
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

    private static List<Integer> getColumnsWithStrongCorr(float[] r2Vals) {
        List<Integer> strongCorrs = new ArrayList<>();
        for (int i = 0; i < r2Vals.length; i++) {
            if (r2Vals[i] >= R2_CUTOFF) {
                strongCorrs.add(i);
            }
        }
        return strongCorrs;
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

        if (DO_TRANSPOSE) {
            return "" + counter + "/(" + maxNanInCols + " : " + avgNanInCols + ")/(" + maxNanInRows + " : " + avgNanInRows + ")";
        }
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

    public static float[][] imputeUntilNoNans(float[][] data) {
        R2_CUTOFF = 0.5f;
        int numNans = getNumNans(data);
        float[][] newData = data;
        do {
            int priorNum = numNans;
            newData = getImputedMatrix(newData, true);
            newData = getImputedMatrix(newData, false);
            numNans = getNumNans(newData);

            if (priorNum == numNans) {
                R2_CUTOFF /= 1.1f;
                System.out.println("New R^2 Cutoff " + R2_CUTOFF);
            }

        } while (numNans > 0);
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
