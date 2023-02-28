/*
 * The MIT License (MIT)
 *
 * Copyright (c) 2011-2023 Rice University, Baylor College of Medicine, Aiden Lab
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

package mixer.utils.tracks;

import javastraw.feature1D.Feature1D;
import javastraw.reader.basics.Chromosome;

public class EigenvectorInterval extends SimpleInterval {

    private final float value;

    public EigenvectorInterval(int chrIndex, String chrName,
                               int x1, int x2,
                               float value, int maxX) {
        super(chrIndex, chrName, x1, x2, maxX);
        this.value = value;
    }

    public EigenvectorInterval(Chromosome chromosome, int x1, int x2, float val, int maxX) {
        this(chromosome.getIndex(), chromosome.getName(), x1, x2, val, maxX);
    }

    public EigenvectorInterval(SubcompartmentInterval interv, float val) {
        this(interv.getChrIndex(), interv.getChrName(), interv.getX1(), interv.getX2(), val, interv.getX2());
    }

    @Override
    public String toString() { // bedgraph format
        String initialName = getChrName().replace("chr", "");
        return "chr" + initialName + "\t" + getX1() + "\t" + getX2() + "\t" + value;
    }

    @Override
    public Feature1D deepClone() {
        return new EigenvectorInterval(getChrIndex(), getChrName(), getX1(), getX2(), value, getX2());
    }
}