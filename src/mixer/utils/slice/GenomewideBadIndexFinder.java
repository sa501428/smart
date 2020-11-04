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

package mixer.utils.slice;

import javastraw.reader.Dataset;
import javastraw.reader.HiCFileTools;
import javastraw.reader.MatrixZoomData;
import javastraw.reader.basics.Block;
import javastraw.reader.basics.Chromosome;
import javastraw.reader.basics.ContactRecord;
import javastraw.type.NormalizationType;
import mixer.utils.common.Pair;

import java.io.IOException;
import java.util.*;

public class GenomewideBadIndexFinder {

    private static final float ZSCORE_MIN_SPARSE_THRESHOLD_LOWER_INTRA = -3;//-5; // 3 bottom 0.15% dropped
    private static final float ZSCORE_MAX_ALLOWED_INTER = 3;
    public static float ZSCORE_MIN_SPARSE_THRESHOLD_INTER_HIGHER = -1;//-2; //-1.28f; // bottom 10% dropped
    private final double GLOBAL_MAX_RAW_VALUE;
    private final int resolution;
    private final Map<Integer, Set<Integer>> badIndices = new HashMap<>();
    private final Map<Integer, Set<Integer>> worstIndices = new HashMap<>();
    private final NormalizationType norm;
    private final float globalInterMean, globalInterStdDev;
    private double globalSum = 0, globalSumOfSquares = 0;
    private long globalCounts = 0;

    public GenomewideBadIndexFinder(Dataset ds, Chromosome[] chromosomes, int resolution, NormalizationType norm) {
        this.resolution = resolution;
        this.norm = norm;
        createInternalBadList(ds, chromosomes);

        globalInterMean = (float) (globalSum / globalCounts);
        double variance = (globalSumOfSquares / globalCounts) - globalInterMean * globalInterMean;
        globalInterStdDev = (float) Math.sqrt(variance);
        GLOBAL_MAX_RAW_VALUE = Math.exp(ZSCORE_MAX_ALLOWED_INTER * globalInterStdDev + globalInterMean) - 1;
    }

    private void createInternalBadList(Dataset ds, Chromosome[] chromosomes) {

        for (Chromosome chr1 : chromosomes) {
            badIndices.put(chr1.getIndex(), new HashSet<Integer>());
            worstIndices.put(chr1.getIndex(), new HashSet<Integer>());
        }

        for (int i = 0; i < chromosomes.length; i++) {
            Chromosome chr1 = chromosomes[i];
            for (int j = i; j < chromosomes.length; j++) {
                Chromosome chr2 = chromosomes[j];
                final MatrixZoomData zd = HiCFileTools.getMatrixZoomData(ds, chr1, chr2, resolution);
                try {
                    determineBadIndicesForRegion(chr1, chr2, zd, i == j);
                } catch (IOException e) {
                    e.printStackTrace();
                }
            }
        }
    }

    /**
     * @return
     */
    private Pair<int[], int[]> getNumberOfNonZeros(List<Block> blocks, int numRows, int numCols, boolean isIntra) {
        int[] numNonZerosRows = new int[numRows];
        int[] numNonZerosCols = new int[numCols];

        for (Block b : blocks) {
            if (b != null) {
                for (ContactRecord cr : b.getContactRecords()) {
                    float val = (float) Math.log(cr.getCounts() + 1);
                    if (Float.isNaN(val) || val < 1e-10 || Float.isInfinite(val)) {
                        continue;
                    }

                    numNonZerosRows[cr.getBinX()]++;
                    numNonZerosCols[cr.getBinY()]++;


                    if (isIntra) {
                        if (cr.getBinX() != cr.getBinY()) {
                            numNonZerosRows[cr.getBinY()]++;
                            numNonZerosCols[cr.getBinX()]++;
                        }
                    } else {
                        globalSum += val;
                        globalSumOfSquares += val * val;
                        globalCounts += 1;
                    }
                }
            }
        }

        return new Pair<>(numNonZerosRows, numNonZerosCols);
    }

    private void determineBadIndicesForRegion(Chromosome chr1, Chromosome chr2, MatrixZoomData zd, boolean isIntra) throws IOException {

        int lengthChr1 = (int) (chr1.getLength() / resolution + 1);
        int lengthChr2 = (int) (chr2.getLength() / resolution + 1);

        determineBadIndicesForRegion(chr1, chr2,
                HiCFileTools.getAllRegionBlocks(zd, 0, lengthChr1, 0, lengthChr2, norm, isIntra),
                lengthChr1, lengthChr2, isIntra);
    }

    private void determineBadIndicesForRegion(Chromosome chr1, Chromosome chr2, List<Block> blocks, int numRows, int numCols, boolean isIntra) {
        Pair<int[], int[]> results = getNumberOfNonZeros(blocks, numRows, numCols, isIntra);
        removeSparserIndicesZscoredCount(chr1, results.getFirst(), isIntra);
        if (!isIntra) {
            removeSparserIndicesZscoredCount(chr2, results.getSecond(), isIntra);
        }
    }

    private void removeSparserIndicesZscoredCount(Chromosome chr1, int[] numNonZeros, boolean isIntra) {
        float mean = getNonZeroMean(numNonZeros);
        float stdDev = getNonZeroStd(numNonZeros, mean);
        getBadIndicesByZscore(chr1, numNonZeros, mean, stdDev, isIntra);
    }

    private void getBadIndicesByZscore(Chromosome chr1, int[] numNonZeros, float mean, float stdDev, boolean isIntra) {
        for (int k = 0; k < numNonZeros.length; k++) {
            if (numNonZeros[k] > 0) {
                float zval = (numNonZeros[k] - mean) / stdDev;
                if (isIntra) {
                    if (zval < ZSCORE_MIN_SPARSE_THRESHOLD_LOWER_INTRA) {
                        worstIndices.get(chr1.getIndex()).add(k);
                        badIndices.get(chr1.getIndex()).add(k);
                    }
                } else if (zval < ZSCORE_MIN_SPARSE_THRESHOLD_INTER_HIGHER) {
                    badIndices.get(chr1.getIndex()).add(k);
                }
            } else {
                badIndices.get(chr1.getIndex()).add(k);
                worstIndices.get(chr1.getIndex()).add(k);
            }
        }
    }

    private float getNonZeroStd(int[] numNonZeros, float mean) {
        int count = 0;
        double total = 0;
        for (int val : numNonZeros) {
            if (val > 0) {
                float diff = val - mean;
                total += (diff * diff);
                count++;
            }
        }
        return (float) Math.sqrt(total / count);
    }

    /**
     * @param numNonZeros
     * @return
     */
    private float getNonZeroMean(int[] numNonZeros) {
        int count = 0;
        float total = 0;
        for (int val : numNonZeros) {
            if (val > 0) {
                total += val;
                count++;
            }
        }
        return total / count;
    }

    public Set<Integer> getBadIndices(Chromosome chrom) {
        return badIndices.get(chrom.getIndex());
    }

    public Set<Integer> getWorstIndices(Chromosome chrom) {
        return worstIndices.get(chrom.getIndex());
    }

    public boolean getExceedsAllowedGlobalValue(float val) {
        return val > GLOBAL_MAX_RAW_VALUE;
    }
}
