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

package mixer.utils.intra;

import javastraw.reader.Dataset;
import javastraw.reader.basics.Chromosome;
import javastraw.reader.block.ContactRecord;
import javastraw.reader.mzd.MatrixZoomData;
import javastraw.reader.type.HiCZoom;
import javastraw.reader.type.NormalizationType;

import java.util.Iterator;

public class OETools {
    public static float[][] getOEMatrix(Dataset ds, MatrixZoomData zd, Chromosome chrom, int resolution, NormalizationType norm) {
        int length = (int) (chrom.getLength() / resolution + 1);
        float[][] matrix = new float[length][length];
        double[] ev = ds.getExpectedValues(new HiCZoom(resolution), norm,
                false).getExpectedValuesWithNormalization(chrom.getIndex()).getValues().get(0);
        Iterator<ContactRecord> recordIterator = zd.getNormalizedIterator(norm);
        while (recordIterator.hasNext()) {
            ContactRecord record = recordIterator.next();
            int dist = Math.abs(record.getBinX() - record.getBinY());
            float oe = (float) ((record.getCounts() + 1) / (ev[dist] + 1));
            matrix[record.getBinX()][record.getBinY()] = oe;
            matrix[record.getBinY()][record.getBinX()] = oe;
        }
        return matrix;
    }
}
