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

import java.util.Collection;

public class ArrayTools {

    private static final float ZERO = 1e-10f;

    public static float getNonZeroStd(float[] numNonZeros, float mean) {
        int count = 0;
        float total = 0;
        for (float val : numNonZeros) {
            if (val > 0) {
                float diff = val - mean;
                total += (diff * diff);
                count++;
            }
        }
        return (float) Math.sqrt(total / count);
    }

    /**
     * @param numNonZeros
     * @return
     */
    public static float getNonZeroMean(float[] numNonZeros) {
        int count = 0;
        float total = 0;
        for (float val : numNonZeros) {
            if (val > ZERO) {
                total += val;
                count++;
            }
        }
        return total / count;
    }

    public static float getNonZeroMean(Collection<float[]> allArrays) {
        double total = 0;
        long count = 0;

        for (float[] array : allArrays) {
            for (float val : array) {
                if (val > ZERO) {
                    total += val;
                    count++;
                }
            }
        }
        return (float) (total / count);
    }

    public static float getNonZeroStd(Collection<float[]> allArrays, float mean) {
        long count = 0;
        double total = 0;
        for (float[] array : allArrays) {
            for (float val : array) {
                if (val > 0) {
                    float diff = val - mean;
                    total += (diff * diff);
                    count++;
                }
            }
        }
        return (float) Math.sqrt(total / count);
    }

    public static float percentNaN(float[] array) {
        float counter = 0;
        for (float val : array) {
            if (Float.isNaN(val)) {
                counter++;
            }
        }
        return counter / array.length;
    }

    public static int max(int[] array) {
        int maxVal = array[0];
        for (int val : array) {
            maxVal = Math.max(maxVal, val);
        }
        return maxVal;
    }

    public static int mean(int[] array) {
        int sum = 0;
        for (int val : array) {
            sum += val;
        }
        return sum / array.length;
    }
}
