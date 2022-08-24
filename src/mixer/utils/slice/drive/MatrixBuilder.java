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

package mixer.utils.slice.drive;

import javastraw.reader.Dataset;
import javastraw.reader.Matrix;
import javastraw.reader.basics.Chromosome;
import javastraw.reader.basics.ChromosomeHandler;
import javastraw.reader.block.ContactRecord;
import javastraw.reader.mzd.MatrixZoomData;
import javastraw.reader.type.HiCZoom;
import javastraw.reader.type.NormalizationType;
import mixer.utils.common.ArrayTools;
import mixer.utils.common.LogTools;
import mixer.utils.magic.FinalScale;
import mixer.utils.magic.SymmLLInterMatrix;
import mixer.utils.slice.matrices.MatrixAndWeight;

import java.util.Arrays;
import java.util.HashMap;
import java.util.Iterator;
import java.util.Map;

public class MatrixBuilder {
    public static MatrixAndWeight populateMatrix(Dataset ds, ChromosomeHandler handler, int resolution,
                                                 NormalizationType norm, Mappings mappings,
                                                 boolean doScale) {
        int numRows = mappings.getNumRows();
        int numCols = mappings.getNumCols();
        int[] weights = new int[numCols];
        float[][] matrix = new float[numRows][numCols];

        Chromosome[] chromosomes = handler.getAutosomalChromosomesArray();
        Map<Integer, int[]> genomewideDistributionForChrom = new HashMap<>();
        for (Chromosome chromosome : chromosomes) {
            genomewideDistributionForChrom.put(chromosome.getIndex(), new int[numCols]);
        }

        for (int i = 0; i < chromosomes.length; i++) {
            for (int j = i + 1; j < chromosomes.length; j++) {
                //if (shouldSkipRegion(chromosomes[i], chromosomes[j])) continue;

                Matrix m1 = ds.getMatrix(chromosomes[i], chromosomes[j]);
                if (m1 == null) continue;
                MatrixZoomData zd = m1.getZoomData(new HiCZoom(resolution));
                if (zd == null) continue;

                if (norm.getLabel().equalsIgnoreCase("none")) {
                    populateMatrixFromIterator(matrix, zd.getDirectIterator(), mappings, chromosomes[i], chromosomes[j]);
                } else {
                    populateMatrixFromIterator(matrix, zd.getNormalizedIterator(norm), mappings, chromosomes[i], chromosomes[j]);
                }

                updateNumberOfLoci(genomewideDistributionForChrom.get(chromosomes[i].getIndex()),
                        mappings.getDistributionForChrom(chromosomes[j]));
                updateNumberOfLoci(genomewideDistributionForChrom.get(chromosomes[j].getIndex()),
                        mappings.getDistributionForChrom(chromosomes[i]));

                System.out.print(".");
            }
            System.out.println(".");
        }

        //ParallelizedStatTools.setZerosToNan(data);
        //ParallelizedStatTools.scaleDown(data, weights);
        LogTools.simpleLogWithCleanup(matrix, Float.NaN);
        //removeHighGlobalThresh(data, weights, 5, Slice.USE_WEIGHTED_MEAN);
        //renormalize(data, weights, -2, 2, Slice.USE_WEIGHTED_MEAN);
        //LogTools.simpleExpm1(data);

        normalizeMatrix(matrix, mappings, chromosomes);

        float[] coverage = new float[numRows];
        updateCoverage(matrix, coverage);
        scaleCoverage(coverage);

        int[] totalDistribution = getSumOfAllLoci(mappings, numCols, chromosomes);
        System.arraycopy(totalDistribution, 0, weights, 0, Math.max(weights.length, totalDistribution.length));
        scaleMatrixColumns(matrix, totalDistribution);

        //todo matrix = EmptyRowCleaner.cleanUpMatrix(matrix, rowIndexToIntervalMap, coverage);

        LogTools.simpleExpm1(matrix);

        if (doScale) {
            return new MatrixAndWeight(
                    FinalScale.scaleMatrix(new SymmLLInterMatrix(matrix),
                            createTargetVector(totalDistribution, numRows, numCols)),
                    weights);
        } else {
            return new MatrixAndWeight(matrix, weights);
        }
    }


    // todo move to mapping?
    private static void populateMatrixFromIterator(float[][] matrix, Iterator<ContactRecord> iterator,
                                                   Mappings mappings, Chromosome c1, Chromosome c2) {

        int[] binToClusterID1 = mappings.getProtocluster(c1);
        int[] binToClusterID2 = mappings.getProtocluster(c2);
        int[] binToGlobalIndex1 = mappings.getGlobalIndex(c1);
        int[] binToGlobalIndex2 = mappings.getGlobalIndex(c2);

        while (iterator.hasNext()) {
            ContactRecord cr = iterator.next();
            if (cr.getCounts() > 0) {
                int r = cr.getBinX();
                int c = cr.getBinY();
                if (binToClusterID1[r] > -1 && binToClusterID2[c] > -1) {
                    matrix[binToGlobalIndex1[r]][binToClusterID2[c]] += cr.getCounts();
                    matrix[binToGlobalIndex2[c]][binToClusterID1[r]] += cr.getCounts();
                }
            }
        }
    }

    /*
    private static boolean shouldSkipRegion(Chromosome c1, Chromosome c2) {
        for (InterChromosomeRegion region : regionsToIgnore) {
            if (region.is(c1, c2)) {
                return true;
            }
        }
        return false;
    }
    */


    private static int[] createTargetVector(int[] colSums, int numRows, int numCols) {
        int[] target = new int[numRows + numCols];
        Arrays.fill(target, 1);
        for (int i = 0; i < numCols; i++) {
            target[i] = colSums[i];
        }
        return target;
    }

    private static void scaleCoverage(float[] coverage) {
        float mean = ArrayTools.getNonZeroMean(coverage);
        for (int k = 0; k < coverage.length; k++) {
            coverage[k] /= mean;
        }
    }

    private static void updateCoverage(float[][] matrix, float[] coverage) {
        for (int i = 0; i < matrix.length; i++) {
            for (int j = 0; j < matrix[i].length; j++) {
                coverage[i] += matrix[i][j];
            }
        }
    }


    private static void updateNumberOfLoci(int[] totalLoci, int[] lociForRegion) {
        for (int i = 0; i < totalLoci.length; i++) {
            totalLoci[i] += lociForRegion[i];
        }
    }

    private static int[] getSumOfAllLoci(Mappings mappings, int numCols, Chromosome[] chromosomes) {


        int[] totalLoci = new int[numCols];
        for (Chromosome chromosome : chromosomes) {
            int[] row = mappings.getDistributionForChrom(chromosome);
            for (int z = 0; z < row.length; z++) {
                totalLoci[z] += row[z];
            }
        }
        return totalLoci;
    }

    private static void normalizeMatrix(float[][] matrix, Mappings mappings, Chromosome[] chromosomes) {

        for (Chromosome chromosome : chromosomes) {
            int[] globalIndices = mappings.getGlobalIndex(chromosome);
            int[] divisor = mappings.getDistributionForChrom(chromosome);
            for (int i : globalIndices) {
                if (i > -1) {
                    divide(matrix[i], divisor);
                }
            }
        }
    }

    private static void divide(float[] row, int[] totalLoci) {
        for (int k = 0; k < row.length; k++) {
            if (totalLoci[k] > 0) {
                row[k] = row[k] / totalLoci[k];
            } else if (row[k] > 0) {
                System.err.println("Impossible situation reached: row val: " + row[k] + " but expect no entries: " + totalLoci[k]);
            }
        }
    }

    private static void scaleMatrixColumns(float[][] matrix, int[] scalars) {
        for (int i = 0; i < matrix.length; i++) {
            for (int z = 0; z < matrix[i].length; z++) {
                matrix[i][z] *= scalars[z];
            }
        }
    }
}
