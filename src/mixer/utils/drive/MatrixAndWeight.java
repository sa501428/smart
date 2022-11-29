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

package mixer.utils.drive;

import javastraw.reader.basics.Chromosome;
import mixer.utils.common.FloatMatrixTools;
import mixer.utils.transform.MatrixTransform;

public class MatrixAndWeight {
    public float[][] matrix, intra;
    public int[] weights;
    private final Mappings mappings;

    public MatrixAndWeight(float[][] interMatrix1, float[][] intraMatrix1, int[] weights, Mappings mappings) {
        this.matrix = interMatrix1;
        this.intra = intraMatrix1;
        this.weights = weights;
        this.mappings = mappings;
    }

    public void divideColumnsByWeights() {
        FloatMatrixTools.divideColumnsByWeights(matrix, weights);
    }

    public int[] getSumOfAllLoci(Chromosome[] chromosomes) {
        int[] totalLoci = new int[mappings.getNumCols()];
        for (Chromosome chromosome : chromosomes) {
            int[] row = mappings.getDistributionForChrom(chromosome);
            for (int z = 0; z < row.length; z++) {
                totalLoci[z] += row[z];
            }
        }
        return totalLoci;
    }

    public void updateWeights(Chromosome[] chromosomes) {
        int[] totalDistribution = getSumOfAllLoci(chromosomes);
        System.arraycopy(totalDistribution, 0, weights, 0, mappings.getNumCols());
    }

    public MatrixAndWeight deepCopy() {
        if (intra != null && intra.length > 0) {
            return new MatrixAndWeight(FloatMatrixTools.deepClone(matrix), FloatMatrixTools.deepClone(intra),
                    FloatMatrixTools.deepClone(weights), mappings.deepCopy());
        }
        return new MatrixAndWeight(FloatMatrixTools.deepClone(matrix), null,
                FloatMatrixTools.deepClone(weights), mappings.deepCopy());
    }

    public void zscoreByCols(int zscoreLimit) {
        MatrixTransform.zscoreByCols(matrix, zscoreLimit);
    }

    public FinalMatrix getFinalMatrix(boolean includeIntra) {
        if (includeIntra && intra != null) {
            return new FinalMatrix(FloatMatrixTools.concatenate(matrix, intra),
                    FloatMatrixTools.concatenate(weights, weights),
                    mappings);
        }
        return new FinalMatrix(matrix, weights, mappings);
    }

    public void doSimpleVCNorm() {
        float[] rowSums = new float[matrix.length];
        float[] colSums = new float[matrix[0].length];
        double sum1 = 0;
        for (int i = 0; i < matrix.length; i++) {
            for (int j = 0; j < matrix[i].length; j++) {
                if (matrix[i][j] > 0) {
                    rowSums[i] += matrix[i][j];
                    colSums[j] += matrix[i][j];
                    sum1 += matrix[i][j];
                }
            }
        }

        double sum2 = 0;
        for (int i = 0; i < matrix.length; i++) {
            for (int j = 0; j < matrix[i].length; j++) {
                double denom = rowSums[i] * colSums[j];
                if (denom > 0 && matrix[i][j] > 0) {
                    matrix[i][j] = (float) (matrix[i][j] / denom);
                    sum2 += matrix[i][j];
                } else {
                    matrix[i][j] = Float.NaN;
                }
            }
        }

        double scale = sum1 / sum2;

        for (int i = 0; i < matrix.length; i++) {
            for (int j = 0; j < matrix[i].length; j++) {
                matrix[i][j] *= scale;
            }
        }
    }
}

