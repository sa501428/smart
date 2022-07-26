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

package mixer.utils.remedian;

import mixer.utils.common.QuickMedian;

import java.util.ArrayList;
import java.util.List;

public class Remedian {
    private final int numValsPerSet;
    private final List<Float> values;
    private final List<Float> medians = new ArrayList<>();
    private int total = 0;
    private Float theMainMedian = null;

    public Remedian(int numValsPerSet) {
        this.numValsPerSet = numValsPerSet;
        values = new ArrayList<>(numValsPerSet);
    }

    public void addVal(float value) {
        values.add(value);
        total++;
        if (total % numValsPerSet == 0) {
            getCurrentMedianAndAppend();
        }
    }

    public float getMedian() {
        if (values.size() > 0) {
            getCurrentMedianAndAppend();
            theMainMedian = null;
        }
        if (theMainMedian == null) {
            theMainMedian = QuickMedian.fastMedian(medians);
        }
        return theMainMedian;
    }

    private void getCurrentMedianAndAppend() {
        float median = QuickMedian.fastMedian(values);
        medians.add(median);
        values.clear();
    }
}
