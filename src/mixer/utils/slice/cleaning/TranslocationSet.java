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

package mixer.utils.slice.cleaning;

import javastraw.reader.basics.Chromosome;
import mixer.utils.common.ArrayTools;
import mixer.utils.common.QuickMedian;

import java.util.HashMap;
import java.util.HashSet;
import java.util.Map;
import java.util.Set;

public class TranslocationSet {
    private static final float TRANSLOCATION_CUTOFF = 3;
    private final Set<String> translocations = new HashSet<>();
    Map<String, Float> maxInRegion = new HashMap<>();

    public void put(Chromosome chr1, Chromosome chr2, float val) {
        maxInRegion.put(key(chr1, chr2), val);
    }

    public String key(Chromosome chr1, Chromosome chr2) {
        if (chr1.getIndex() < chr2.getIndex()) {
            return chr1.getIndex() + "-" + chr2.getIndex();
        } else {
            return chr2.getIndex() + "-" + chr1.getIndex();
        }
    }

    public boolean getHasTranslocation(Chromosome chr1, Chromosome chr2) {
        return translocations.contains(key(chr1, chr2));
    }

    public void determineTranslocations() {
        float[] values = ArrayTools.toArray(maxInRegion.values());
        float mean = ArrayTools.getNonZeroMean(values);
        float median = QuickMedian.fastMedian(values);
        float stdDev = ArrayTools.getNonZeroStd(values, mean);
        float mad = 1.4826f * QuickMedian.fastMedianAbsDeviation(values, median);

        for (String key : maxInRegion.keySet()) {
            Float val = maxInRegion.get(key);
            float zscore = (val - mean) / stdDev;
            float robustZscore = (val - median) / mad;
            //System.out.println("Translocation: "+key +" z: "+zscore +"  rz: "+robustZscore);
            if (zscore > TRANSLOCATION_CUTOFF) {
                translocations.add(key);
                System.out.println("*** Translocation: " + key + " z: " + zscore + "  rz: " + robustZscore + "");
            }
        }
    }
}
