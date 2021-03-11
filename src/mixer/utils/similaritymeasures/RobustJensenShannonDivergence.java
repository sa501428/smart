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

package mixer.utils.similaritymeasures;

public class RobustJensenShannonDivergence extends SimilarityMetric {
    public static final RobustJensenShannonDivergence SINGLETON = new RobustJensenShannonDivergence();

    private RobustJensenShannonDivergence() {
        super(true);
    }

    private static float nonNanJSDistance(float[] a, float[] b) {
        double sumA = 0, sumB = 0;
        for (int k = 0; k < a.length; k++) {
            boolean isBad = Float.isNaN(a[k] + b[k]);
            if (!isBad) {
                sumA += a[k];
                sumB += b[k];
            }
        }
        final double sumC = sumA + sumB;

        double distP = 0, distQ = 0;
        for (int i = 0; i < a.length; i++) {
            boolean isBad = Float.isNaN(a[i] + b[i]);
            if (!isBad) {
                final double p = a[i] / sumA;
                final double q = b[i] / sumB;
                final double m = (a[i] + b[i]) / sumC;
                distP += (p * Math.log(p / m));
                distQ += (q * Math.log(q / m));
            }
        }
        return (float) ((distP + distQ) / 2f);
    }

    @Override
    public float distance(final float[] x, final float[] y) {
        return nonNanJSDistance(x, y);
    }
}
