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

package mixer.utils.matrix;

import javastraw.reader.Dataset;
import javastraw.reader.basics.Chromosome;
import javastraw.reader.basics.ChromosomeHandler;
import javastraw.reader.block.ContactRecord;
import javastraw.reader.mzd.MatrixZoomData;
import javastraw.reader.type.NormalizationType;

import java.util.Iterator;

public class InterOnlyMatrix extends HiCMatrix {

    public InterOnlyMatrix(Dataset ds, NormalizationType norm, int resolution,
                           Chromosome[] rowsChromosomes, Chromosome[] colsChromosomes) {
        super(ds, norm, resolution, rowsChromosomes, colsChromosomes);
    }

    @Override
    protected void fillInChromosomeRegion(float[][] matrix, MatrixZoomData zd, Chromosome chr1, int offsetIndex1,
                                          Chromosome chr2, int offsetIndex2, boolean needToFlip) {
        Iterator<ContactRecord> iterator = zd.getNormalizedIterator(norm);
        while (iterator.hasNext()) {
            ContactRecord record = iterator.next();
            if (record.getCounts() > 0) {
                if (needToFlip) {
                    matrix[offsetIndex1 + record.getBinY()][offsetIndex2 + record.getBinX()] = record.getCounts();
                } else {
                    matrix[offsetIndex1 + record.getBinX()][offsetIndex2 + record.getBinY()] = record.getCounts();
                }
            }
        }
    }

    public static InterOnlyMatrix getMatrix(Dataset ds, NormalizationType norm, int resolution, InterMapType mapType) {
        ChromosomeHandler chromosomeHandler = ds.getChromosomeHandler();
        Chromosome[] rowsChromosomes, colsChromosomes;
        switch (mapType) {
            case SKIP_BY_TWOS: // but start with CHR 1 separate
                rowsChromosomes = chromosomeHandler.splitAutosomesAndSkipByTwos().a;
                colsChromosomes = chromosomeHandler.splitAutosomesAndSkipByTwos().b;
                break;
            case FIRST_HALF_VS_SECOND_HALF:
                rowsChromosomes = chromosomeHandler.splitAutosomesIntoHalves().a;
                colsChromosomes = chromosomeHandler.splitAutosomesIntoHalves().b;
                break;
            case ODDS_VS_EVENS:
            default:
                rowsChromosomes = chromosomeHandler.extractOddOrEvenAutosomes(true);
                colsChromosomes = chromosomeHandler.extractOddOrEvenAutosomes(false);
                break;
        }

        return new InterOnlyMatrix(ds, norm, resolution, rowsChromosomes, colsChromosomes);
    }

    public enum InterMapType {ODDS_VS_EVENS, FIRST_HALF_VS_SECOND_HALF, SKIP_BY_TWOS}
}
