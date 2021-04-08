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

import javastraw.reader.Dataset;
import javastraw.reader.basics.Chromosome;
import javastraw.reader.mzd.MatrixZoomData;
import javastraw.reader.type.NormalizationType;
import javastraw.tools.HiCFileTools;
import mixer.utils.common.FloatMatrixTools;
import mixer.utils.common.LogTools;
import mixer.utils.common.ZScoreTools;
import mixer.utils.similaritymeasures.SimilarityMetric;
import mixer.utils.slice.cleaning.IntraMatrixCleaner;
import mixer.utils.slice.cleaning.SimilarityMatrixTools;
import mixer.utils.slice.matrices.Dimension;

import java.util.Arrays;

abstract public class HiCMatrix {
    public static boolean USE_ZSCORE = false;
    public static int NUM_CENTROIDS = 1;
    protected final NormalizationType norm;
    protected final int resolution;
    protected final float[][] interMatrix;
    protected final Chromosome[] rowsChromosomes;
    protected final Chromosome[] colsChromosomes;
    protected final Dimension rowsDimension, colsDimension;
    protected final SimilarityMetric metric;
    protected final boolean isIntra;
    protected final INTRA_TYPE intraType;

    public HiCMatrix(Dataset ds, NormalizationType norm, int resolution,
                     Chromosome[] rowsChromosomes, Chromosome[] colsChromosomes,
                     SimilarityMetric metric, boolean isIntra, INTRA_TYPE intraType,
                     boolean shouldZeroOutNans, int compressionFactor) {
        this.norm = norm;
        this.resolution = resolution;
        this.metric = metric;
        this.rowsChromosomes = rowsChromosomes;
        this.colsChromosomes = colsChromosomes;
        this.isIntra = isIntra;
        this.intraType = intraType;
        rowsDimension = new Dimension(rowsChromosomes, resolution);
        colsDimension = new Dimension(colsChromosomes, resolution);
        interMatrix = makeCleanScaledInterMatrix(ds, shouldZeroOutNans, compressionFactor);
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

    private float[][] makeCleanScaledInterMatrix(Dataset ds, boolean shouldZeroOutNans,
                                                 int compressionFactor) {

        float[][] interMatrix = new float[rowsDimension.length][colsDimension.length];
        for (int i = 0; i < rowsChromosomes.length; i++) {
            Chromosome chr1 = rowsChromosomes[i];
            for (int j = 0; j < colsChromosomes.length; j++) {
                Chromosome chr2 = colsChromosomes[j];
                final MatrixZoomData zd = HiCFileTools.getMatrixZoomData(ds, chr1, chr2, resolution);
                if (zd == null) continue;

                // will need to flip across diagonal
                boolean needToFlip = chr2.getIndex() < chr1.getIndex();
                fillInChromosomeRegion(ds, interMatrix, zd, chr1, rowsDimension.offset[i],
                        chr2, colsDimension.offset[j], needToFlip);
            }
            System.out.print(".");
        }
        System.out.println(".");

        FloatMatrixTools.cleanUpMatrix(interMatrix, shouldZeroOutNans);

        if (compressionFactor > 0) {
            interMatrix = IntraMatrixCleaner.compress(interMatrix, compressionFactor);
            FloatMatrixTools.log(interMatrix, 1);
        }

        if (USE_ZSCORE) {
            int[] weights = new int[interMatrix[0].length];
            Arrays.fill(weights, 1);
            ZScoreTools.inPlaceZscoreDownCol(interMatrix, weights);
        }

        if (metric != null) {
            return SimilarityMatrixTools.getNonNanSimilarityMatrix(interMatrix,
                    metric, NUM_CENTROIDS, 283746L);
        }

        return interMatrix;
    }

    abstract protected void fillInChromosomeRegion(Dataset ds, float[][] matrix, MatrixZoomData zd,
                                                   Chromosome chr1, int offsetIndex1,
                                                   Chromosome chr2, int offsetIndex2, boolean needToFlip);

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
        LogTools.applySimpleLog(interMatrix);
    }

    public enum INTRA_TYPE {DEFAULT, JUST_LOG, TRUE_OE, LOG_BASE_EXPECTED}
}
