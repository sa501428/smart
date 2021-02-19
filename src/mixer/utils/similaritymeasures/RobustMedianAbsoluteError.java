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

import mixer.utils.common.QuickMedian;

import java.util.ArrayList;
import java.util.List;

/**
 * Euclidean distance.
 */
public final class RobustMedianAbsoluteError extends SimilarityMetric {

  /**
   * Euclidean metric.
   */
  public static final RobustMedianAbsoluteError SINGLETON = new RobustMedianAbsoluteError();
  private static final float ZERO = 1e-10f;

  private RobustMedianAbsoluteError() {
    super(true);
  }

  private static float sortedMidpoint(List<Float> vals) {
    int size = vals.size();
    vals.sort(null);
    return (vals.get(size / 2) + vals.get((size - 1) / 2)) / 2;
  }

  @Override
  public float distance(final float[] x, final float[] y) {
    return getNonNanMedianAbsDeviation(x, y);
  }

  private float getNonNanMedianAbsDeviation(float[] x, float[] y) {
    List<Float> vals = new ArrayList<>();
    for (int i = 0; i < x.length; i++) {
      final float v = Math.abs(x[i] - y[i]);
      if (!Float.isNaN(v) && v > ZERO) {
        vals.add(v);
      }
    }
    //return sortedMidpoint(vals);
    if (vals.size() > 0) {
      return QuickMedian.fastMedian(vals);
    } else {
      return 0f;
    }
  }
}
