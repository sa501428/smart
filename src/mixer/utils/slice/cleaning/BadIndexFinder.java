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

package mixer.utils.slice.cleaning;

import javastraw.reader.Dataset;
import javastraw.reader.HiCFileTools;
import javastraw.reader.MatrixZoomData;
import javastraw.reader.basics.Block;
import javastraw.reader.basics.Chromosome;
import javastraw.reader.basics.ContactRecord;
import javastraw.type.NormalizationType;

import java.io.IOException;
import java.util.*;

public class BadIndexFinder {
    // ridiculously high coverage filter
    private static final float ZSCORE_COVERAGE_MAX_ALLOWED_INTER = 10;
    private final int resolution;
    private final Map<Integer, Set<Integer>> badIndices = new HashMap<>();
    private final NormalizationType[] norms;
    //private final double[] globalSum, globalSumOfSquares, GLOBAL_MAX_RAW_VALUE;
    //private final long[] globalCounts;

    public BadIndexFinder(List<Dataset> datasets, Chromosome[] chromosomes, int resolution,
                          NormalizationType[] norms) {
        this.resolution = resolution;
        this.norms = norms;
        for (Chromosome chrom : chromosomes) {
            badIndices.put(chrom.getIndex(), new HashSet<>());
        }

        /*
        globalSum = new double[datasets.size()];
        globalSumOfSquares = new double[datasets.size()];
        GLOBAL_MAX_RAW_VALUE = new double[datasets.size()];
        Arrays.fill(GLOBAL_MAX_RAW_VALUE, Double.MAX_VALUE);
        globalCounts = new long[datasets.size()];
         */

        createInternalBadList(datasets, chromosomes);

        /*
        float[] globalInterMean = new float[datasets.size()];
        float[] globalInterStdDev = new float[datasets.size()];
        for (int k = 0; k < datasets.size(); k++) {
            globalInterMean[k] = (float) (globalSum[k] / globalCounts[k]);
            double variance = (globalSumOfSquares[k] / globalCounts[k]) - globalInterMean[k] * globalInterMean[k];
            globalInterStdDev[k] = (float) Math.sqrt(variance);
            GLOBAL_MAX_RAW_VALUE[k] = Math.exp(ZSCORE_MAX_ALLOWED_INTER * globalInterStdDev[k] + globalInterMean[k]) - 1;
        }
         */
    }

    private void createInternalBadList(List<Dataset> datasets, Chromosome[] chromosomes) {
        for (int z = 0; z < datasets.size(); z++) {
            for (int i = 0; i < chromosomes.length; i++) {
                Chromosome chr1 = chromosomes[i];
                for (int j = i; j < chromosomes.length; j++) {
                    Chromosome chr2 = chromosomes[j];
                    final MatrixZoomData zd = HiCFileTools.getMatrixZoomData(datasets.get(z), chr1, chr2, resolution);
                    try {
                        determineBadIndicesForRegion(chr1, chr2, zd, i == j, z);
                    } catch (IOException e) {
                        e.printStackTrace();
                    }
                }
            }
        }
    }

    private RegionStatistics getNumberOfNonZeros(List<Block> blocks, int numRows, int numCols,
                                                 boolean isIntra, int dIndex) {
        RegionStatistics stats = new RegionStatistics(numRows, numCols);

        for (Block b : blocks) {
            if (b != null) {
                for (ContactRecord cr : b.getContactRecords()) {
                    float val = (float) Math.log(cr.getCounts() + 1);
                    if (Float.isNaN(val) || val < 1e-10 || Float.isInfinite(val)) {
                        continue;
                    }

                    stats.numNonZerosRows[cr.getBinX()]++;
                    stats.numNonZerosCols[cr.getBinY()]++;

                    if (isIntra) {
                        if (cr.getBinX() != cr.getBinY()) {
                            stats.numNonZerosRows[cr.getBinY()]++;
                            stats.numNonZerosCols[cr.getBinX()]++;
                        }
                    } else {
                        stats.rowSums[cr.getBinX()] += val;
                        stats.colSums[cr.getBinY()] += val;

                        //globalSum[dIndex] += val;
                        //globalSumOfSquares[dIndex] += val * val;
                        //globalCounts[dIndex] += 1;
                    }
                }
            }
        }

        return stats;
    }

    private void determineBadIndicesForRegion(Chromosome chr1, Chromosome chr2, MatrixZoomData zd,
                                              boolean isIntra, int dIndex) throws IOException {

        int lengthChr1 = (int) (chr1.getLength() / resolution + 1);
        int lengthChr2 = (int) (chr2.getLength() / resolution + 1);

        determineBadIndicesForRegion(chr1, chr2,
                HiCFileTools.getAllRegionBlocks(zd, 0, lengthChr1, 0, lengthChr2, norms[dIndex], isIntra),
                lengthChr1, lengthChr2, isIntra, dIndex);
    }

    private void determineBadIndicesForRegion(Chromosome chr1, Chromosome chr2, List<Block> blocks,
                                              int numRows, int numCols, boolean isIntra, int dIndex) {
        RegionStatistics stats = getNumberOfNonZeros(blocks, numRows, numCols, isIntra, dIndex);
        removeSparserIndicesZscoredCount(chr1, stats.numNonZerosRows, isIntra);
        if (!isIntra) {
            removeExtremeCoverage(stats.rowSums, chr1);
            removeExtremeCoverage(stats.colSums, chr2);
            removeSparserIndicesZscoredCount(chr2, stats.numNonZerosCols, false);
        }
    }

    private void removeExtremeCoverage(double[] sums, Chromosome chrom) {
        double mean = ArrayTools.getNonZeroMean(sums);
        double stdDev = ArrayTools.getNonZeroStd(sums, mean);
        getBadCoverageIndicesByZscore(chrom, sums, mean, stdDev);
    }

    private void removeSparserIndicesZscoredCount(Chromosome chr1, int[] numNonZeros, boolean isIntra) {
        float mean = ArrayTools.getNonZeroMeanIntArray(numNonZeros);
        float stdDev = ArrayTools.getNonZeroStdIntArray(numNonZeros, mean);
        getBadIndicesByZscore(chr1, numNonZeros, mean, stdDev, isIntra);
    }

    private void getBadIndicesByZscore(Chromosome chr1, int[] numNonZeros, float mean, float stdDev, boolean isIntra) {
        for (int k = 0; k < numNonZeros.length; k++) {
            if (numNonZeros[k] < 1) {
                badIndices.get(chr1.getIndex()).add(k);
            }
        }
    }

    private void getBadCoverageIndicesByZscore(Chromosome chr1, double[] sums, double mean, double stdDev) {
        // todo ideally needs to be genome wide check
        for (int k = 0; k < sums.length; k++) {
            if (sums[k] > 0) {
                double zval = (sums[k] - mean) / stdDev;
                if (zval > ZSCORE_COVERAGE_MAX_ALLOWED_INTER) {
                    badIndices.get(chr1.getIndex()).add(k);
                }
            } else {
                badIndices.get(chr1.getIndex()).add(k);
            }
        }
    }

    public Set<Integer> getBadIndices(Chromosome chrom) {
        return badIndices.get(chrom.getIndex());
    }

    private class RegionStatistics {
        int[] numNonZerosRows;
        int[] numNonZerosCols;
        double[] rowSums;
        double[] colSums;

        RegionStatistics(int numRows, int numCols) {
            numNonZerosRows = new int[numRows];
            numNonZerosCols = new int[numCols];
            rowSums = new double[numRows];
            colSums = new double[numCols];
        }
    }
}
