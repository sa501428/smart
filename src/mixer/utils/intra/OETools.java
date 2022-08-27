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

import javastraw.reader.basics.Chromosome;
import javastraw.reader.block.ContactRecord;
import javastraw.reader.mzd.MatrixZoomData;
import javastraw.reader.type.NormalizationType;
import mixer.utils.common.LogExpectedSpline;

import java.util.Arrays;
import java.util.Iterator;
import java.util.Set;

public class OETools {

    private static final int FIVE_MB = 5000000, FIFTY_MB = 50000000;

    public static float[][] getCleanOEMatrix(MatrixZoomData zd, Chromosome chrom, int resolution,
                                             NormalizationType norm, Set<Integer> badIndices) {

        LogExpectedSpline spline = new LogExpectedSpline(zd, norm, (int) (chrom.getLength() / resolution) + 1);

        int length = (int) (chrom.getLength() / resolution + 1);
        int minDist = FIVE_MB / resolution;

        float[][] matrix = new float[length][length];
        for (float[] row : matrix) {
            Arrays.fill(row, Float.NaN);
        }
        Iterator<ContactRecord> recordIterator = zd.getNormalizedIterator(norm);
        while (recordIterator.hasNext()) {
            ContactRecord record = recordIterator.next();
            int dist = Math.abs(record.getBinX() - record.getBinY());
            if (dist > minDist) {
                float oe = (float) Math.log(record.getCounts() / spline.getExpectedFromUncompressedBin(dist));
                matrix[record.getBinX()][record.getBinY()] = oe;
                matrix[record.getBinY()][record.getBinX()] = oe;
            }
        }

        IntraMatrixCleaner.nanFillNearDiagonal(matrix, FIVE_MB / resolution);
        IntraMatrixCleaner.nanFillBadRowsColumns(badIndices, matrix);
        //IntraMatrixCleaner.nanFillZeroEntries(matrix);
        return matrix;
    }
}
