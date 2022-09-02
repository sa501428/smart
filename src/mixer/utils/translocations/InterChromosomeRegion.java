/*
 * The MIT License (MIT)
 *
 * Copyright (c) 2011-2022 Rice University, Baylor College of Medicine, Aiden Lab
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

package mixer.utils.translocations;

import javastraw.reader.basics.Chromosome;

import java.util.Objects;

public class InterChromosomeRegion {
    Chromosome c1, c2;

    public InterChromosomeRegion(Chromosome c1, Chromosome c2) {
        if (c1.getIndex() <= c2.getIndex()) {
            this.c1 = c1;
            this.c2 = c2;
        } else {
            this.c1 = c2;
            this.c2 = c1;
        }
    }

    public boolean is(Chromosome a, Chromosome b) {
        if (c1.getIndex() == a.getIndex() && c2.getIndex() == b.getIndex()) return true;
        return c2.getIndex() == a.getIndex() && c1.getIndex() == b.getIndex();
    }

    @Override
    public boolean equals(Object obj) {
        if (this == obj) {
            return true;
        }
        if (obj instanceof InterChromosomeRegion) {
            InterChromosomeRegion o = (InterChromosomeRegion) obj;
            return o.is(c1, c2);
        }
        return false;
    }

    @Override
    public int hashCode() {
        return Objects.hash(c1, c2);
    }
}
