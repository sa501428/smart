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
 * Bray Curtis distance.
 */
public final class RobustEarthMoversMetric extends Metric {

  /**
   * Bray Curtis distance.
   */
  public static final RobustEarthMoversMetric SINGLETON = new RobustEarthMoversMetric();

  private RobustEarthMoversMetric() {
    super(false);
  }

  public static float nonNanEarthMoversDistance(float[] a, float[] b) {
    double sumA = 0, sumB = 0;
    boolean[] isValid = new boolean[a.length];
    for (int k = 0; k < a.length; k++) {
      isValid[k] = !(Float.isNaN(a[k]) || Float.isNaN(b[k]));
      if (isValid[k]) {
        sumA += a[k];
        sumB += b[k];
      }
    }
    sumA = sumA / 100;
    sumB = sumB / 100;

    double lastDistance = 0;
    double totalDistance = 0;
    int numValid = 0;
    for (int i = 0; i < a.length; i++) {
      if (isValid[i]) {
        double tempA = a[i] / sumA;
        double tempB = b[i] / sumB;
        final double currentDistance = (tempA + lastDistance) - tempB;
        totalDistance += Math.abs(currentDistance);
        lastDistance = currentDistance;
        numValid++;
      }
    }
    numValid = Math.max(numValid, 1);
    return (float) (a.length * totalDistance / numValid);
  }

  @Override
  public float distance(final float[] x, final float[] y) {
    return nonNanEarthMoversDistance(x, y);
  }
}
