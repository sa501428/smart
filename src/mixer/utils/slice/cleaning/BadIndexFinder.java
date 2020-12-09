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
import mixer.utils.common.ArrayTools;

import java.io.IOException;
import java.util.*;

public class BadIndexFinder {
    // ridiculously high coverage filter
    private static final float ZSCORE_COVERAGE_MAX_ALLOWED_INTER = 10;
    private final int resolution;
    private final Map<Integer, Set<Integer>> badIndices = new HashMap<>();
    private final NormalizationType[] norms;

    public BadIndexFinder(List<Dataset> datasets, Chromosome[] chromosomes, int resolution,
                          NormalizationType[] norms) {
        this.resolution = resolution;
        this.norms = norms;
        for (Chromosome chrom : chromosomes) {
            badIndices.put(chrom.getIndex(), new HashSet<>());
        }

        createInternalBadList(datasets, chromosomes);
    }

    private void createInternalBadList(List<Dataset> datasets, Chromosome[] chromosomes) {
        for (int z = 0; z < datasets.size(); z++) {
            for (int i = 0; i < chromosomes.length; i++) {
                Chromosome chr1 = chromosomes[i];
                for (int j = i + 1; j < chromosomes.length; j++) {
                    Chromosome chr2 = chromosomes[j];
                    final MatrixZoomData zd = HiCFileTools.getMatrixZoomData(datasets.get(z), chr1, chr2, resolution);
                    try {
                        determineBadIndicesForRegion(chr1, chr2, zd, z);
                    } catch (IOException e) {
                        e.printStackTrace();
                    }
                }
            }
        }
    }

    private void determineBadIndicesForRegion(Chromosome chr1, Chromosome chr2, MatrixZoomData zd,
                                              int dIndex) throws IOException {
        int lengthChr1 = (int) (chr1.getLength() / resolution + 1);
        int lengthChr2 = (int) (chr2.getLength() / resolution + 1);

        List<Block> blocks = HiCFileTools.getAllRegionBlocks(zd, 0, lengthChr1, 0, lengthChr2,
                norms[dIndex], false);
        RegionStatistics stats = getSums(blocks, lengthChr1, lengthChr2);
        removeExtremeCoverage(stats.rowSums, chr1);
        removeExtremeCoverage(stats.colSums, chr2);
    }

    private RegionStatistics getSums(List<Block> blocks, int numRows, int numCols) {
        RegionStatistics stats = new RegionStatistics(numRows, numCols);

        for (Block b : blocks) {
            if (b != null) {
                for (ContactRecord cr : b.getContactRecords()) {
                    // todo check log after sum
                    float val = (float) Math.log(cr.getCounts() + 1);
                    if (Float.isNaN(val) || val < 1e-10 || Float.isInfinite(val)) {
                        continue;
                    }
                    stats.rowSums[cr.getBinX()] += val;
                    stats.colSums[cr.getBinY()] += val;
                }
            }
        }

        return stats;
    }

    private void removeExtremeCoverage(double[] sums, Chromosome chrom) {
        double mean = ArrayTools.getNonZeroMean(sums);
        double stdDev = ArrayTools.getNonZeroStd(sums, mean);
        getBadCoverageIndicesByZscore(chrom, sums, mean, stdDev);
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
        double[] rowSums;
        double[] colSums;

        RegionStatistics(int numRows, int numCols) {
            rowSums = new double[numRows];
            colSums = new double[numCols];
        }
    }
}
