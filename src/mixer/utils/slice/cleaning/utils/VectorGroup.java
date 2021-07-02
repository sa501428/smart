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

package mixer.utils.slice.cleaning.utils;

import mixer.utils.similaritymeasures.RobustCorrelationSimilarity;

import java.util.ArrayList;
import java.util.HashSet;
import java.util.List;
import java.util.Set;

public class VectorGroup {
    private final float threshold = 0.5f;
    private final Set<Integer> indices = new HashSet<>();
    private final List<float[]> vectors = new ArrayList<>();

    public VectorGroup(int index, float[] vector) {
        append(index, vector);
    }

    public void append(int index, float[] vector) {
        indices.add(index);
        vectors.add(vector);
    }

    public boolean shouldInclude(float[] vector) {
        for (float[] vec : vectors) {
            float corr = RobustCorrelationSimilarity.SINGLETON.distance(vec, vector);
            if (Math.abs(corr) > threshold) {
                return true;
            }
        }

        return false;
    }

    public int size() {
        return indices.size();
    }

    public Set<Integer> getMembers() {
        return indices;
    }
}
