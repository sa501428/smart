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
}
