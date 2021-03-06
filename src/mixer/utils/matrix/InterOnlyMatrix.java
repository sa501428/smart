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

package mixer.utils.matrix;

import javastraw.reader.ChromosomeHandler;
import javastraw.reader.Dataset;
import javastraw.reader.HiCFileTools;
import javastraw.reader.MatrixZoomData;
import javastraw.reader.basics.Chromosome;
import javastraw.type.NormalizationType;
import mixer.utils.common.FloatMatrixTools;
import mixer.utils.similaritymeasures.SimilarityMetric;

public class InterOnlyMatrix extends HiCMatrix {

    public InterOnlyMatrix(Dataset ds, NormalizationType norm, int resolution,
                           Chromosome[] rowsChromosomes, Chromosome[] colsChromosomes,
                           SimilarityMetric metric) {
        super(ds, norm, resolution, rowsChromosomes, colsChromosomes, metric, false, INTRA_TYPE.DEFAULT);
    }

    public static InterOnlyMatrix getMatrix(Dataset ds, NormalizationType norm, int resolution, InterMapType mapType,
                                            SimilarityMetric metric) {
        ChromosomeHandler chromosomeHandler = ds.getChromosomeHandler();
        Chromosome[] rowsChromosomes, colsChromosomes;
        switch (mapType) {
            case SKIP_BY_TWOS: // but start with CHR 1 separate
                rowsChromosomes = chromosomeHandler.splitAutosomesAndSkipByTwos().getFirst();
                colsChromosomes = chromosomeHandler.splitAutosomesAndSkipByTwos().getSecond();
                break;
            case FIRST_HALF_VS_SECOND_HALF:
                rowsChromosomes = chromosomeHandler.splitAutosomesIntoHalves().getFirst();
                colsChromosomes = chromosomeHandler.splitAutosomesIntoHalves().getSecond();
                break;
            case ODDS_VS_EVENS:
            default:
                rowsChromosomes = chromosomeHandler.extractOddOrEvenAutosomes(true);
                colsChromosomes = chromosomeHandler.extractOddOrEvenAutosomes(false);
                break;
        }

        return new InterOnlyMatrix(ds, norm, resolution, rowsChromosomes, colsChromosomes, metric);
    }

    protected void fillInInterChromosomeRegion(Dataset ds, float[][] matrix, MatrixZoomData zd, Chromosome chr1, int offsetIndex1,
                                               Chromosome chr2, int offsetIndex2, boolean needToFlip) {

        int lengthChr1 = (int) Math.ceil((float) chr1.getLength() / resolution);
        int lengthChr2 = (int) Math.ceil((float) chr2.getLength() / resolution);

        float[][] allDataForRegion = null;
        try {
            if (needToFlip) {
                float[][] allDataForRegionMatrix = HiCFileTools.extractLocalBoundedRegionFloatMatrix(zd, 0, lengthChr2,
                        0, lengthChr1, lengthChr2, lengthChr1, norm, false);
                allDataForRegion = FloatMatrixTools.transpose(allDataForRegionMatrix);
            } else {
                allDataForRegion = HiCFileTools.extractLocalBoundedRegionFloatMatrix(zd, 0, lengthChr1,
                        0, lengthChr2, lengthChr1, lengthChr2, norm, false);
            }
        } catch (Exception e) {
            e.printStackTrace();
            System.exit(99);
        }

        for (int i = 0; i < allDataForRegion.length; i++) {
            System.arraycopy(allDataForRegion[i], 0,
                    matrix[offsetIndex1 + i], offsetIndex2, allDataForRegion[i].length);
        }
    }


    public enum InterMapType {ODDS_VS_EVENS, FIRST_HALF_VS_SECOND_HALF, SKIP_BY_TWOS}
}
