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

import java.awt.*;
import java.util.List;
import java.util.*;

public class TranslocationSet {
    private final Set<String> translocations = new HashSet<>();
    Map<String, List<Rectangle>> maxInRegion = new HashMap<>();

    public void put(Chromosome chr1, Chromosome chr2, List<Rectangle> vals) {
        maxInRegion.put(key(chr1, chr2), vals);
        translocations.add(key(chr1, chr2));
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
}
