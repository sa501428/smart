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

import javastraw.expected.ExpectedModel;
import javastraw.expected.ExpectedUtils;
import javastraw.expected.LogExpectedSpline;
import javastraw.reader.basics.Chromosome;
import javastraw.reader.block.ContactRecord;
import javastraw.reader.mzd.MatrixZoomData;
import javastraw.reader.type.NormalizationType;
import mixer.utils.common.FloatMatrixTools;
import mixer.utils.drive.LogExpectedSubset;

import java.util.*;

public class OETools {

    private static final int FIVE_MB = 5000000, FIFTY_MB = 50000000;

    public static float[][] getCleanOEMatrix(MatrixZoomData zd, Chromosome chrom, int resolution,
                                             NormalizationType norm, Set<Integer> badIndices, int resFactor,
                                             boolean takeLog, boolean skipNearDiagonal, boolean useExpandedIntraOE) {

        LogExpectedSpline spline = new LogExpectedSpline(zd, norm, chrom, resolution);

        int length = (int) (chrom.getLength() / resolution + 1);

        float[][] matrix = newNanMatrix(length);
        float[][] matrix2 = newNanMatrix(length);

        List<ContactRecord> filteredContacts = filter(resolution, zd.getNormalizedIterator(norm));
        for (ContactRecord record : filteredContacts) {
            int dist = ExpectedModel.getDist(record);
            float oe = (float) (record.getCounts() / spline.getExpectedFromUncompressedBin(dist));
            if (takeLog) {
                oe = (float) Math.log(oe);
            }
            matrix[record.getBinX()][record.getBinY()] = oe;
            matrix[record.getBinY()][record.getBinX()] = oe;
        }

        if (useExpandedIntraOE) {
            LogExpectedSubset expected = new LogExpectedSubset(filteredContacts, chrom, resolution);
            for (ContactRecord record : filteredContacts) {
                if (expected.isInInterval(record)) {
                    float z = expected.getZscoreForObservedUncompressedBin(record);
                    if (Math.abs(z) < 5) {
                        matrix2[record.getBinX()][record.getBinY()] = z;
                        matrix2[record.getBinY()][record.getBinX()] = z;
                    }
                }
            }
        }

        if (skipNearDiagonal) {
            IntraMatrixCleaner.nanFillNearDiagonal(matrix, FIVE_MB / resolution);
            if (useExpandedIntraOE) {
                IntraMatrixCleaner.nanFillNearDiagonal(matrix2, FIVE_MB / resolution);
            }
        }

        IntraMatrixCleaner.nanFillBadRowsColumns(badIndices, matrix, resFactor);
        if (useExpandedIntraOE) {
            IntraMatrixCleaner.nanFillBadRowsColumns(badIndices, matrix2, resFactor);
            return FloatMatrixTools.concatenate(matrix, matrix2);
        }
        return matrix;
    }

    private static float[][] newNanMatrix(int length) {
        float[][] matrix = new float[length][length];
        for (float[] row : matrix) {
            Arrays.fill(row, Float.NaN);
        }
        return matrix;
    }

    public static List<ContactRecord> filter(int resolution, Iterator<ContactRecord> iterator) {
        int minDist = FIVE_MB / resolution;
        List<ContactRecord> records = new LinkedList<>();
        while (iterator.hasNext()) {
            ContactRecord cr = iterator.next();
            if (cr.getCounts() > 0) {
                if (ExpectedUtils.getDist(cr) > minDist) {
                    records.add(cr);
                }
            }
        }
        return records;
    }
}
