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

package mixer.utils.shuffle;

import javastraw.reader.*;
import javastraw.reader.basics.Chromosome;
import javastraw.type.NormalizationType;
import mixer.utils.common.FloatMatrixTools;
import mixer.utils.similaritymeasures.SimilarityMetric;
import mixer.utils.slice.cleaning.MatrixCleanerAndProjector;
import mixer.utils.slice.cleaning.SimilarityMatrixTools;
import mixer.utils.slice.matrices.Dimension;

public class InterOnlyMatrix {

    private final INTRA_TYPE intra_type;
    private final boolean isIntra;
    private final SimilarityMetric metric;
    private final NormalizationType norm;
    private final int resolution;
    private final float[][] interMatrix;
    private final Chromosome[] rowsChromosomes;
    private final Chromosome[] colsChromosomes;
    private Dimension rowsDimension, colsDimension;

    public InterOnlyMatrix(Dataset ds, NormalizationType norm, int resolution, InterMapType mapType,
                           SimilarityMetric metric) {
        this.norm = norm;
        this.resolution = resolution;
        isIntra = false;
        intra_type = INTRA_TYPE.DEFAULT; // unused
        this.metric = metric;
        ChromosomeHandler chromosomeHandler = ds.getChromosomeHandler();

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

        interMatrix = makeCleanScaledInterMatrix(ds);
    }

    public InterOnlyMatrix(Dataset ds, NormalizationType norm, int resolution, Chromosome chromosome,
                           INTRA_TYPE intra_type, SimilarityMetric metric) {
        this.norm = norm;
        this.resolution = resolution;
        isIntra = true;
        this.intra_type = intra_type;
        this.metric = metric;

        rowsChromosomes = new Chromosome[1];
        colsChromosomes = new Chromosome[1];
        rowsChromosomes[0] = chromosome;
        colsChromosomes[0] = chromosome;

        interMatrix = makeCleanScaledInterMatrix(ds);
    }

    public static INTRA_TYPE getIntraType(int mapTypeOption) {
        switch (mapTypeOption) {
            case 3:
                return INTRA_TYPE.LOG_BASE_EXPECTED;
            case 2:
                return INTRA_TYPE.TRUE_OE;
            case 1:
                return INTRA_TYPE.JUST_LOG;
            case 0:
            default:
                return INTRA_TYPE.DEFAULT;
        }
    }

    private float[][] makeCleanScaledInterMatrix(Dataset ds) {
        rowsDimension = new Dimension(rowsChromosomes, resolution);
        colsDimension = new Dimension(colsChromosomes, resolution);
        float[][] interMatrix = new float[rowsDimension.length][colsDimension.length];

        for (int i = 0; i < rowsChromosomes.length; i++) {
            Chromosome chr1 = rowsChromosomes[i];
            for (int j = 0; j < colsChromosomes.length; j++) {
                Chromosome chr2 = colsChromosomes[j];
                final MatrixZoomData zd = HiCFileTools.getMatrixZoomData(ds, chr1, chr2, resolution);
                if (zd == null) continue;

                // will need to flip across diagonal
                boolean needToFlip = chr2.getIndex() < chr1.getIndex();
                fillInInterChromosomeRegion(ds, interMatrix, zd, chr1, rowsDimension.offset[i], chr2, colsDimension.offset[j], needToFlip);
            }
            System.out.print(".");
        }
        System.out.println(".");


        return FloatMatrixTools.cleanUpMatrix(interMatrix);
    }

    private void fillInInterChromosomeRegion(Dataset ds, float[][] matrix, MatrixZoomData zd, Chromosome chr1, int offsetIndex1,
                                             Chromosome chr2, int offsetIndex2, boolean needToFlip) {

        int lengthChr1 = (int) Math.ceil((float) chr1.getLength() / resolution);
        int lengthChr2 = (int) Math.ceil((float) chr2.getLength() / resolution);

        float[][] allDataForRegion = null;
        try {
            if (chr1.getIndex() == chr2.getIndex()) {
                if (intra_type == INTRA_TYPE.DEFAULT) {
                    allDataForRegion = HiCFileTools.extractLocalBoundedRegionFloatMatrix(zd, 0, lengthChr1,
                            0, lengthChr2, lengthChr1, lengthChr2, norm, true);
                } else if (intra_type == INTRA_TYPE.LOG_BASE_EXPECTED) {
                    allDataForRegion = HiCFileTools.getOEMatrixForChromosome(ds, zd, chr1, resolution,
                            norm, 10, ExtractingOEDataUtils.ThresholdType.LOG_BASE_EXP_OF_OBS,
                            true, 1, 0);
                } else if (intra_type == INTRA_TYPE.TRUE_OE) {
                    allDataForRegion = HiCFileTools.getOEMatrixForChromosome(ds, zd, chr1, resolution,
                            norm, 10, ExtractingOEDataUtils.ThresholdType.TRUE_OE,
                            true, 1, 0);
                } else if (intra_type == INTRA_TYPE.JUST_LOG) {
                    allDataForRegion = HiCFileTools.extractLocalBoundedRegionFloatMatrix(zd, 0, lengthChr1,
                            0, lengthChr2, lengthChr1, lengthChr2, norm, true);
                    MatrixCleanerAndProjector.simpleLogWithCleanup(allDataForRegion, 0);
                } else {
                    System.err.println("Invalid Matrix type " + intra_type);
                    System.exit(9);
                }
                if (metric != null) {
                    allDataForRegion = SimilarityMatrixTools.getNonNanSimilarityMatrix(allDataForRegion,
                            metric, 1, 283746L);
                }

            } else if (needToFlip) {
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

    public enum INTRA_TYPE {DEFAULT, JUST_LOG, TRUE_OE, LOG_BASE_EXPECTED}

    public Chromosome[] getRowChromosomes() {
        return rowsChromosomes;
    }

    public Chromosome[] getColChromosomes() {
        return colsChromosomes;
    }

    public float[][] getMatrix() {
        return interMatrix;
    }

    public int[] getRowOffsets() {
        return rowsDimension.offset;
    }

    public int[] getColOffsets() {
        return colsDimension.offset;
    }

    public void applySimpleLog() {
        for (int i = 0; i < interMatrix.length; i++) {
            for (int j = 0; j < interMatrix[i].length; j++) {
                float val = interMatrix[i][j];
                if (Float.isNaN(val)) {
                    interMatrix[i][j] = 0;
                } else {
                    val = (float) Math.log(val + 1);
                    if (Float.isInfinite(val)) {
                        interMatrix[i][j] = 0;
                    } else {
                        interMatrix[i][j] = val;
                    }
                }
            }
        }
    }

    public enum InterMapType {ODDS_VS_EVENS, FIRST_HALF_VS_SECOND_HALF, SKIP_BY_TWOS}
}
