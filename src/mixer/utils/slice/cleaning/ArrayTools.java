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

package mixer.utils.slice.cleaning;

public class ArrayTools {
    public static float getNonZeroStdIntArray(int[] numNonZeros, float mean) {
        int count = 0;
        double total = 0;
        for (int val : numNonZeros) {
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
    public static float getNonZeroMeanIntArray(int[] numNonZeros) {
        int count = 0;
        float total = 0;
        for (int val : numNonZeros) {
            if (val > 0) {
                total += val;
                count++;
            }
        }
        return total / count;
    }

    public static double getNonZeroStd(double[] numNonZeros, double mean) {
        int count = 0;
        double total = 0;
        for (double val : numNonZeros) {
            if (val > 0) {
                double diff = val - mean;
                total += (diff * diff);
                count++;
            }
        }
        return Math.sqrt(total / count);
    }

    /**
     * @param numNonZeros
     * @return
     */
    public static double getNonZeroMean(double[] numNonZeros) {
        int count = 0;
        double total = 0;
        for (double val : numNonZeros) {
            if (val > 0) {
                total += val;
                count++;
            }
        }
        return total / count;
    }
}
