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
import mixer.utils.common.SimpleArray2DTools;
import mixer.utils.tracks.SubcompartmentInterval;
import mixer.utils.transform.MatrixTransform;

import java.util.HashMap;
import java.util.Map;

public class MatrixAndWeight {
    public float[][] matrix, intra, matrix2;
    public int[] weights;
    private final Map<Integer, SubcompartmentInterval> map = new HashMap<>();
    private final Mappings mappings;

    public MatrixAndWeight(float[][] interMatrix1, float[][] intraMatrix1,
                           float[][] interMatrix2, int[] weights, Mappings mappings) {
        this.matrix = interMatrix1;
        this.intra = intraMatrix1;
        this.matrix2 = interMatrix2;
        this.weights = weights;
        this.mappings = mappings;
        if (mappings != null) populateRowIndexToIntervalMap(mappings);
    }

    private void populateRowIndexToIntervalMap(Mappings mappings) {
        int resolution = mappings.getResolution();
        Chromosome[] chromosomes = mappings.getChromosomes();
        for (Chromosome chromosome : chromosomes) {
            int maxGenomeLen = (int) chromosome.getLength();
            int[] globalIndices = mappings.getGlobalIndex(chromosome);
            for (int i = 0; i < globalIndices.length; i++) {
                if (globalIndices[i] > -1) {
                    int coord = globalIndices[i];
                    int x1 = i * resolution;
                    int x2 = Math.min(x1 + resolution, maxGenomeLen);
                    SubcompartmentInterval newRInterval = new SubcompartmentInterval(chromosome, x1, x2, -1);
                    map.put(coord, newRInterval);
                }
            }
        }
    }

    public void divideColumnsByWeights(boolean useBothNorms) {
        FloatMatrixTools.divideColumnsByWeights(matrix, weights);
        if (useBothNorms) {
            FloatMatrixTools.divideColumnsByWeights(matrix2, weights);
        }
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
        return new MatrixAndWeight(FloatMatrixTools.deepClone(matrix), FloatMatrixTools.deepClone(intra),
                FloatMatrixTools.deepClone(matrix2), FloatMatrixTools.deepClone(weights), mappings.deepCopy());
    }

    public void putIntraIntoMainMatrix() {
        for (int i = 0; i < intra.length; i++) {
            for (int j = 0; j < intra[i].length; j++) {
                if (intra[i][j] > -10) {
                    matrix[i][j] = intra[i][j];
                    matrix2[i][j] = intra[i][j];
                }
            }
        }
    }

    public void zscoreByRows(int zscoreLimit, boolean useBothNorms) {
        MatrixTransform.zscoreByRows(matrix, zscoreLimit);
        if (useBothNorms) {
            MatrixTransform.zscoreByRows(matrix2, zscoreLimit);
        }
    }

    public void applyLog(boolean useBothNorms) {
        SimpleArray2DTools.simpleLogWithCleanup(matrix, Float.NaN);
        if (useBothNorms) {
            SimpleArray2DTools.simpleLogWithCleanup(matrix2, Float.NaN);
        }
    }

    public FinalMatrix getFinalMatrix(boolean useBothNorms, boolean appendIntra) {
        if (useBothNorms) {
            if (appendIntra) {
                return new FinalMatrix(FloatMatrixTools.concatenate(matrix, matrix2, intra),
                        FloatMatrixTools.concatenate(weights, weights, mutiply(weights, 2)),
                        mappings, map);
            } else {
                return new FinalMatrix(FloatMatrixTools.concatenate(matrix, matrix2),
                        FloatMatrixTools.concatenate(weights, weights),
                        mappings, map);
            }
        } else if (appendIntra) {
            return new FinalMatrix(FloatMatrixTools.concatenate(matrix, intra),
                    FloatMatrixTools.concatenate(weights, weights),
                    mappings, map);
        } else {
            return new FinalMatrix(matrix, weights, mappings, map);
        }
    }

    private int[] mutiply(int[] input, int scalar) {
        int[] copy = new int[input.length];
        for (int i = 0; i < copy.length; i++) {
            copy[i] = input[i] * scalar;
        }
        return copy;
    }
}

