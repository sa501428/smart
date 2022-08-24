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

import javastraw.reader.basics.Chromosome;
import mixer.utils.slice.structures.SubcompartmentInterval;

import java.util.Arrays;
import java.util.HashMap;
import java.util.Map;

public class BinMappings implements Mappings {
    // chromosome to (bin to proto-cluster index)
    // chromosome to (bin to global index)
    private final Map<Integer, int[]> chromToBinToProtocluster = new HashMap<>();
    private final Map<Integer, int[]> chromToBinToGlobalIndex = new HashMap<>();
    private final Map<Integer, int[]> chromToDistributionForChromosome = new HashMap<>();
    private int numRows = 0;
    private int numCols = 0;

    public void putBinToProtoCluster(Chromosome chrom, int[] binToProtocluster) {
        chromToBinToProtocluster.put(chrom.getIndex(), binToProtocluster);
    }

    public void calculateGlobalIndices(Chromosome[] chromosomes) {
        int counter = 0;
        int maxCol = 0;
        for (Chromosome chromosome : chromosomes) {
            int[] protoCluster = chromToBinToProtocluster.get(chromosome.getIndex());
            int[] globalIndex = new int[protoCluster.length];
            Arrays.fill(globalIndex, -1);
            for (int i = 0; i < protoCluster.length; i++) {
                int val = protoCluster[i];
                if (val > -1) {
                    globalIndex[i] = counter;
                    counter++;
                    maxCol = Math.max(maxCol, val);
                }
            }
            chromToBinToGlobalIndex.put(chromosome.getIndex(), globalIndex);
        }
        numRows = counter;
        numCols = maxCol + 1;

        for (Chromosome chromosome : chromosomes) {
            int[] counts = new int[numCols];
            int[] protoCluster = chromToBinToProtocluster.get(chromosome.getIndex());
            for (int val : protoCluster) {
                if (val > -1) {
                    counts[val]++;
                }
            }
            chromToDistributionForChromosome.put(chromosome.getIndex(), counts);
        }
    }

    public int[] getProtoclusterAssignments(Chromosome chrom) {
        return chromToBinToProtocluster.get(chrom.getIndex());
    }

    @Override
    public int getNumRows() {
        return numRows;
    }

    @Override
    public int getNumCols() {
        return numCols;
    }

    @Override
    public int[] getDistributionForChrom(Chromosome chromosome) {
        return chromToDistributionForChromosome.get(chromosome.getIndex());
    }

    @Override
    public int[] getProtocluster(Chromosome c1) {
        return chromToBinToProtocluster.get(c1);
    }

    @Override
    public int[] getGlobalIndex(Chromosome c1) {
        return chromToBinToGlobalIndex.get(c1);
    }

    @Override
    protected int[][] getGenomeIndices() {
        int[][] coordinates = new int[numRows][3];
        for (int i = 0; i < n; i++) {
            SubcompartmentInterval interval = rowIndexToIntervalMap.get(i);
            coordinates[i][0] = interval.getChrIndex();
            coordinates[i][1] = interval.getX1();
            coordinates[i][2] = interval.getX2();
        }
        return coordinates;
    }
}