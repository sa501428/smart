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

package mixer.utils.slice.cleaning;

import javastraw.reader.Dataset;
import javastraw.reader.HiCFileTools;
import javastraw.reader.MatrixZoomData;
import javastraw.reader.basics.Block;
import javastraw.reader.basics.Chromosome;
import javastraw.type.NormalizationType;
import mixer.utils.common.ArrayTools;

import java.io.IOException;
import java.util.*;

public class BadIndexFinder {
    // ridiculously high coverage filter
    protected static final float ZSCORE_COVERAGE_MAX_ALLOWED_INTER = 5;
    protected static final float ZSCORE_MIN_NONZERO_COVERAGE_NEEDED_INTER = -3;
    protected final int resolution;
    protected final Map<Integer, Set<Integer>> badIndices = new HashMap<>();
    protected final NormalizationType[] norms;

    public BadIndexFinder(Chromosome[] chromosomes, int resolution,
                          NormalizationType[] norms) {
        this.resolution = resolution;
        this.norms = norms;
        for (Chromosome chrom : chromosomes) {
            badIndices.put(chrom.getIndex(), new HashSet<>());
        }
    }

    public void createInternalBadList(List<Dataset> datasets, Chromosome[] chromosomes) {
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

    protected void determineBadIndicesForRegion(Chromosome chr1, Chromosome chr2, MatrixZoomData zd,
                                                int dIndex) throws IOException {
        int lengthChr1 = (int) (chr1.getLength() / resolution + 1);
        int lengthChr2 = (int) (chr2.getLength() / resolution + 1);

        List<Block> blocks = HiCFileTools.getAllRegionBlocks(zd, 0, lengthChr1, 0, lengthChr2,
                norms[dIndex], false);
        RegionStatistics stats = new RegionStatistics(lengthChr1, lengthChr2, blocks);

        removeExtremeCoverages(stats.rowSums, chr1, false);
        removeExtremeCoverages(stats.colSums, chr2, false);
        removeExtremeCoverages(stats.rowNonZeros, chr1, true);
        removeExtremeCoverages(stats.colNonZeros, chr2, true);
    }

    private void removeExtremeCoverages(float[] sums, Chromosome chrom, boolean removeLowCoverage) {
        float mean = ArrayTools.getNonZeroMean(sums);
        float stdDev = ArrayTools.getNonZeroStd(sums, mean);
        getBadCoverageIndices(chrom, sums, mean, stdDev, removeLowCoverage);
    }

    protected void getBadCoverageIndices(Chromosome chr1, float[] sums, double mean, double stdDev,
                                         boolean removeLowCoverage) {
        if (sums == null) {
            System.err.println("Skipping 0 " + chr1.getName() + " " + removeLowCoverage);
            return;
        }
        for (int k = 0; k < sums.length; k++) {
            if (sums[k] > 0) {
                double zval = (sums[k] - mean) / stdDev;
                if (removeLowCoverage) {
                    if (zval < ZSCORE_MIN_NONZERO_COVERAGE_NEEDED_INTER) {
                        badIndices.get(chr1.getIndex()).add(k);
                    }
                } else if (zval > ZSCORE_COVERAGE_MAX_ALLOWED_INTER) {
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
}
