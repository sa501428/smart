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
package mixer.utils.custommetrics;

import tagbio.umap.metric.Metric;

/**
 * Correlation distance.
 */
public final class RobustCorrelationMetric extends Metric {

    /**
     * Correlation distance.
     */
    public static final RobustCorrelationMetric SINGLETON = new RobustCorrelationMetric();

    private RobustCorrelationMetric() {
        super(true);
    }

    @Override
    public float distance(final float[] x, final float[] y) {
        int counter = 0;
        double sumX = 0;
        double sumY = 0;
        for (int i = 0; i < x.length; i++) {
            boolean entryIsBad = Float.isNaN(x[i]) || Float.isNaN(y[i]);
            if (!entryIsBad) {
                sumX += x[i];
                sumY += y[i];
                counter++;
            }
        }
        double muX = sumX / counter;
        double muY = sumY / counter;

        double dotProduct = 0.0;
        double normX = 0.0;
        double normY = 0.0;
        for (int i = 0; i < x.length; i++) {
            boolean entryIsBad = Float.isNaN(x[i]) || Float.isNaN(y[i]);
            if (!entryIsBad) {
                double nX = x[i] - muX;
                double nY = y[i] - muY;

                dotProduct += nX * nY;
                normX += nX * nX;
                normY += nY * nY;
            }
        }

        if (normX == 0.0 && normY == 0.0) {
            return 0;
        } else if (dotProduct == 0.0) {
            return 1;
        } else {
            return (float) (1 - (dotProduct / Math.sqrt(normX * normY)));
        }
    }
}
