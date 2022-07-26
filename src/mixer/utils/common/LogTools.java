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

package mixer.utils.common;

import javastraw.tools.ParallelizationTools;

import java.util.concurrent.atomic.AtomicInteger;

public class LogTools {
    public static void simpleLogWithCleanup(float[][] matrix, float badVal) {
        AtomicInteger index = new AtomicInteger(0);
        ParallelizationTools.launchParallelizedCode(() -> {
            int i = index.getAndIncrement();
            while (i < matrix.length) {
                for (int j = 0; j < matrix[i].length; j++) {
                    float val = matrix[i][j];
                    if (!Float.isNaN(val)) {
                        val = (float) Math.log(val + 1);
                        if (Float.isInfinite(val)) {
                            matrix[i][j] = badVal;
                        } else {
                            matrix[i][j] = val;
                        }
                    }
                }
                i = index.getAndIncrement();
            }
        });
    }

    public static float getMaxAbsLogVal(float[][] matrix) {
        double maxVal = Math.abs(Math.log(matrix[0][0]));
        for (int i = 0; i < matrix.length; i++) {
            for (int j = 0; j < matrix[i].length; j++) {
                double temp = Math.abs(Math.log(matrix[i][j]));
                if (temp > maxVal) {
                    maxVal = temp;
                }
            }
        }
        return (float) maxVal;
    }

    public static void applySimpleLog(float[][] interMatrix) {
        for (int i = 0; i < interMatrix.length; i++) {
            for (int j = 0; j < interMatrix[i].length; j++) {
                float val = interMatrix[i][j];
                if (Float.isNaN(val)) {
                    interMatrix[i][j] = 0;
                } else {
                    val = (float) Math.log(val + 1);
                    if (Float.isInfinite(val)) {
                        interMatrix[i][j] = 0;
                    } else {
                        interMatrix[i][j] = val;
                    }
                }
            }
        }
    }

    public static void scaleDownThenLogThenScaleUp(float[][] matrix, int[] weights) {
        for (int i = 0; i < matrix.length; i++) {
            for (int j = 0; j < matrix[i].length; j++) {
                float val = matrix[i][j] / weights[j];
                if (!Float.isNaN(val)) {
                    val = (float) Math.log(val + 1);
                    if (Float.isInfinite(val)) {
                        matrix[i][j] = Float.NaN;
                    } else {
                        matrix[i][j] = val * weights[j];
                    }
                }
            }
        }
    }

    public static void expInPlace(float[][] interMatrix) {
        for (int i = 0; i < interMatrix.length; i++) {
            for (int j = 0; j < interMatrix[i].length; j++) {
                float val = interMatrix[i][j];
                if (!Float.isNaN(val)) {
                    val = expTanhZ(val);
                    if (Float.isInfinite(val)) {
                        val = Float.NaN;
                    }
                    interMatrix[i][j] = val;
                }
            }
        }
    }

    private static float expTanhZ(float val) {
        double val1 = (val + 1) / 10; //+1
        double val2 = 3 * Math.tanh(val1);
        //return (float) val2;
        return (float) Math.expm1(val2); // Math.expm1(val)
    }

    public static void expInPlaceType2(float[][] interMatrix) {
        for (int i = 0; i < interMatrix.length; i++) {
            for (int j = 0; j < interMatrix[i].length; j++) {
                float val = interMatrix[i][j];
                if (!Float.isNaN(val)) {
                    val = (float) (Math.expm1(val - 1)); // Math.tanh
                    if (Float.isInfinite(val)) {
                        val = Float.NaN;
                    }
                    interMatrix[i][j] = val;
                }
            }
        }
    }

    public static void expInPlaceType3(float[][] interMatrix) {
        for (int i = 0; i < interMatrix.length; i++) {
            for (int j = 0; j < interMatrix[i].length; j++) {
                float val = interMatrix[i][j];
                if (!Float.isNaN(val)) {
                    val = (float) Math.tanh(Math.expm1(val - 1)); //
                    if (Float.isInfinite(val)) {
                        val = Float.NaN;
                    }
                    interMatrix[i][j] = val;
                }
            }
        }
    }

    public static void expInPlaceType4(float[][] interMatrix) {
        for (int i = 0; i < interMatrix.length; i++) {
            for (int j = 0; j < interMatrix[i].length; j++) {
                float val = interMatrix[i][j];
                if (!Float.isNaN(val)) {
                    val = (float) Math.tanh(val - 1); //
                    if (Float.isInfinite(val)) {
                        val = Float.NaN;
                    }
                    interMatrix[i][j] = val;
                }
            }
        }
    }

    public static void expInPlaceType5(float[][] interMatrix) {
        for (int i = 0; i < interMatrix.length; i++) {
            for (int j = 0; j < interMatrix[i].length; j++) {
                float val = interMatrix[i][j];
                if (!Float.isNaN(val)) {
                    val = (float) Math.tanh(val); //
                    if (Float.isInfinite(val)) {
                        val = Float.NaN;
                    }
                    interMatrix[i][j] = val;
                }
            }
        }
    }

    public static void eluInPlaceType6(float[][] interMatrix) {
        for (int i = 0; i < interMatrix.length; i++) {
            for (int j = 0; j < interMatrix[i].length; j++) {
                float val = interMatrix[i][j];
                if (!Float.isNaN(val)) {
                    val = elu(val - 1); //
                    if (Float.isInfinite(val)) {
                        val = Float.NaN;
                    }
                    interMatrix[i][j] = val;
                }
            }
        }
    }

    private static float elu(float v) {
        if (v < 0) {
            return (float) (Math.expm1(v));
        }
        return v;
    }

    public static void reluInPlaceType7(float[][] interMatrix) {
        for (int i = 0; i < interMatrix.length; i++) {
            for (int j = 0; j < interMatrix[i].length; j++) {
                float val = interMatrix[i][j];
                if (!Float.isNaN(val)) {
                    val = relu(val - 0.5); //
                    if (Float.isInfinite(val)) {
                        val = Float.NaN;
                    }
                    interMatrix[i][j] = val;
                }
            }
        }
    }

    private static float relu(double v) {
        if (v < 0) {
            return 0f;
        }
        return (float) v;
    }

    public static void eluInPlaceType8(float[][] interMatrix) {
        for (int i = 0; i < interMatrix.length; i++) {
            for (int j = 0; j < interMatrix[i].length; j++) {
                float val = interMatrix[i][j];
                if (!Float.isNaN(val)) {
                    val = elu(val - 0.5f); //
                    if (Float.isInfinite(val)) {
                        val = Float.NaN;
                    }
                    interMatrix[i][j] = val;
                }
            }
        }
    }

    public static void eluInPlaceType9(float[][] interMatrix) {
        for (int i = 0; i < interMatrix.length; i++) {
            for (int j = 0; j < interMatrix[i].length; j++) {
                float val = interMatrix[i][j];
                if (!Float.isNaN(val)) {
                    val = elu(val); //
                    if (Float.isInfinite(val)) {
                        val = Float.NaN;
                    }
                    interMatrix[i][j] = val;
                }
            }
        }
    }

    public static void simpleExpm1(float[][] matrix) {
        AtomicInteger index = new AtomicInteger(0);
        ParallelizationTools.launchParallelizedCode(() -> {
            int i = index.getAndIncrement();
            while (i < matrix.length) {
                for (int j = 0; j < matrix[i].length; j++) {
                    float val = matrix[i][j];
                    if (!Float.isNaN(val)) {
                        matrix[i][j] = (float) Math.expm1(val);
                    }
                }
                i = index.getAndIncrement();
            }
        });
    }
}
