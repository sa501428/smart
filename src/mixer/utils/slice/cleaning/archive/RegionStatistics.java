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

package mixer.utils.slice.cleaning.archive;

import javastraw.reader.block.Block;
import javastraw.reader.block.ContactRecord;
import mixer.utils.slice.cleaning.CRecordUtils;

import java.util.List;

class RegionStatistics {
    final float[] rowSums;
    final float[] colSums;
    final float[] rowNonZeros;
    final float[] colNonZeros;

    RegionStatistics(int numRows, int numCols, List<Block> blocks) {
        rowSums = new float[numRows];
        colSums = new float[numCols];
        rowNonZeros = new float[numRows];
        colNonZeros = new float[numCols];

        for (Block b : blocks) {
            if (b != null) {
                for (ContactRecord cr : b.getContactRecords()) {
                    CRecordUtils.add(cr, rowSums, colSums, rowNonZeros, colNonZeros);
                }
            }
        }
    }
}
