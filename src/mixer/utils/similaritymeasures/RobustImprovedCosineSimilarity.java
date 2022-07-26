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

public final class RobustImprovedCosineSimilarity extends SimilarityMetric {

  public static final RobustImprovedCosineSimilarity SINGLETON = new RobustImprovedCosineSimilarity();

  private RobustImprovedCosineSimilarity() {
    super(true);
  }

  @Override
  public float distance(final float[] x, final float[] y) {
    double dotProduct = 0.0;
    double norm1X = 0.0;
    double norm1Y = 0.0;
    for (int i = 0; i < x.length; i++) {
      double product = Math.sqrt(x[i] * y[i]);
      if (!Double.isNaN(product)) {
        dotProduct += product;
        norm1X += x[i];
        norm1Y += y[i];
      }
    }

    return (float) (dotProduct / Math.sqrt(norm1X * norm1Y));
  }
}
