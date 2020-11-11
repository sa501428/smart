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
 * Euclidean distance.
 */
public final class RobustEuclideanMetric extends Metric {

  /**
   * Euclidean metric.
   */
  public static final RobustEuclideanMetric SINGLETON = new RobustEuclideanMetric();

  private RobustEuclideanMetric() {
    super(false);
  }

  @Override
  public float distance(final float[] x, final float[] y) {
    //  D(x, y) = \sqrt{\sum_i (x_i - y_i)^2}
    float result = getNonNanMeanSquaredError(x, y);
    return (float) Math.sqrt(result * x.length);
  }

  private float getNonNanMeanSquaredError(float[] x, float[] y) {
    double sumOfSquares = 0;
    int numDiffs = 0;
    for (int i = 0; i < x.length; i++) {
      if (!Float.isNaN(x[i]) && !Float.isNaN(y[i])) {
        final float v = x[i] - y[i];
        sumOfSquares += (v * v);
        numDiffs++;
      }
    }
    numDiffs = Math.max(numDiffs, 1);
    return (float) (sumOfSquares / numDiffs);
  }
}
