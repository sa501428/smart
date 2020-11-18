/*
 * The MIT License (MIT)
 *
 * Copyright (c) 2011-2020 Rice University, Baylor College of Medicine, Aiden Lab
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

import javastraw.reader.ChromosomeHandler;
import javastraw.reader.Dataset;
import javastraw.reader.HiCFileTools;
import javastraw.reader.MatrixZoomData;
import javastraw.reader.basics.Chromosome;
import javastraw.type.NormalizationType;
import mixer.utils.common.FloatMatrixTools;
import mixer.utils.common.Pair;

public class InterOnlyMatrix {

    private final NormalizationType norm;
    private final int resolution;
    private final float[][] interMatrix;
    private final Chromosome[] rowsChromosomes;
    private final Chromosome[] colsChromosomes;
    private Pair<Integer, int[]> rowsDimension, colsDimension;

    public InterOnlyMatrix(Dataset ds, NormalizationType norm, int resolution, InterMapType mapType) {
        this.norm = norm;
        this.resolution = resolution;
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

    /**
     * methods for constructing initial matrix
     */

    private Pair<Integer, int[]> calculateDimensionInterMatrix(Chromosome[] chromosomes) {
        int total = 0;
        int[] indices = new int[chromosomes.length];

        for (int i = 0; i < chromosomes.length; i++) {
            total += (int) (chromosomes[i].getLength() / resolution + 1);
            if (i < chromosomes.length - 1) {
                indices[i + 1] = total;
            }
        }

        return new Pair<>(total, indices);
    }

    private float[][] makeCleanScaledInterMatrix(Dataset ds) {
        rowsDimension = calculateDimensionInterMatrix(rowsChromosomes);
        colsDimension = calculateDimensionInterMatrix(colsChromosomes);

        float[][] interMatrix = new float[rowsDimension.getFirst()][colsDimension.getFirst()];

        for (int i = 0; i < rowsChromosomes.length; i++) {
            Chromosome chr1 = rowsChromosomes[i];
            for (int j = 0; j < colsChromosomes.length; j++) {
                Chromosome chr2 = colsChromosomes[j];

                if (chr1.getIndex() == chr2.getIndex()) continue;
                final MatrixZoomData zd = HiCFileTools.getMatrixZoomData(ds, chr1, chr2, resolution);
                if (zd == null) continue;

                // will need to flip across diagonal
                boolean needToFlip = chr2.getIndex() < chr1.getIndex();
                fillInInterChromosomeRegion(interMatrix, zd, chr1, rowsDimension.getSecond()[i], chr2, colsDimension.getSecond()[j], needToFlip);
            }
            System.out.print(".");
        }
        System.out.println(".");


        return FloatMatrixTools.cleanUpMatrix(interMatrix);
    }

    private void fillInInterChromosomeRegion(float[][] matrix, MatrixZoomData zd, Chromosome chr1, int offsetIndex1,
                                             Chromosome chr2, int offsetIndex2, boolean needToFlip) {

        int lengthChr1 = (int) (chr1.getLength() / resolution);
        int lengthChr2 = (int) (chr2.getLength() / resolution);

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
                    matrix[offsetIndex1 + i], offsetIndex2 + 0, allDataForRegion[i].length);
        }
    }

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
        return rowsDimension.getSecond();
    }

    public int[] getColOffsets() {
        return colsDimension.getSecond();
    }

    public enum InterMapType {ODDS_VS_EVENS, FIRST_HALF_VS_SECOND_HALF, SKIP_BY_TWOS}
}
