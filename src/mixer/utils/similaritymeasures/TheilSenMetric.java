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

public class TheilSenMetric extends SimilarityMetric {

    public static final TheilSenMetric SINGLETON = new TheilSenMetric();
    private static final float ZERO = 1e-10f;

    private TheilSenMetric() {
        super(false);
    }

    private static float getTheilSenSlope(float[] xx, float[] yy) {
        List<Float> slopesList = new ArrayList<>();
        for (int i = 0; i < xx.length; i++) {
            float x = xx[i];
            float y = yy[i];
            if (Float.isNaN(x) || Float.isNaN(y)) continue;
            for (int j = i + 1; j < xx.length; j++) {
                if (x != xx[j]) { // x must be different, otherwise slope becomes infinite
                    float slope = (yy[j] - y) / (xx[j] - x);
                    if (!Float.isNaN(slope)) {
                        slopesList.add(slope);
                    }
                }
            }
        }

        return QuickMedian.fastMedian(slopesList);
    }

    @Override
    public float distance(float[] x, float[] y, int index, int skip) {
        return getTheilSenSlope(x, y);
    }
}
