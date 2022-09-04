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

package mixer.utils.common;

import javastraw.tools.ParallelizationTools;

import java.util.concurrent.atomic.AtomicInteger;

public class SimpleArray2DTools {
    public static void simpleLogWithCleanup(float[][] matrix, float badVal) {
        AtomicInteger index = new AtomicInteger(0);
        ParallelizationTools.launchParallelizedCode(() -> {
            int i = index.getAndIncrement();
            while (i < matrix.length) {
                for (int j = 0; j < matrix[i].length; j++) {
                    float val = matrix[i][j];
                    if (val > 0) {
                        val = (float) Math.log(val + 1);
                        if (Float.isInfinite(val)) {
                            matrix[i][j] = badVal;
                        } else {
                            matrix[i][j] = val;
                        }
                    } else {
                        matrix[i][j] = badVal;
                    }
                }
                i = index.getAndIncrement();
            }
        });
    }

    public static void simpleLogWithCleanup(double[][] matrix, double badVal) {
        AtomicInteger index = new AtomicInteger(0);
        ParallelizationTools.launchParallelizedCode(() -> {
            int i = index.getAndIncrement();
            while (i < matrix.length) {
                for (int j = 0; j < matrix[i].length; j++) {
                    double val = matrix[i][j];
                    if (val > 0) {
                        val = (float) Math.log(val + 1);
                        if (Double.isInfinite(val)) {
                            matrix[i][j] = badVal;
                        } else {
                            matrix[i][j] = val;
                        }
                    } else {
                        matrix[i][j] = badVal;
                    }
                }
                i = index.getAndIncrement();
            }
        });
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

    public static void setZerosToNan(float[][] matrix) {
        AtomicInteger index = new AtomicInteger(0);
        ParallelizationTools.launchParallelizedCode(() -> {
            int i = index.getAndIncrement();
            while (i < matrix.length) {
                for (int j = 0; j < matrix[i].length; j++) {
                    if (matrix[i][j] < 1e-20) {
                        matrix[i][j] = Float.NaN;
                    }
                }
                i = index.getAndIncrement();
            }
        });
    }
}
