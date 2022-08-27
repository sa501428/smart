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
    private static final int IGNORE = -1;
    private final int resolution;

    public BinMappings(int resolution) {
        this.resolution = resolution;
    }

    public void putBinToProtoCluster(Chromosome chrom, int[] binToProtocluster) {
        chromToBinToProtocluster.put(chrom.getIndex(), binToProtocluster);
    }

    public void calculateGlobalIndices(Chromosome[] chromosomes) {
        int counter = 0;
        int maxCol = 0;
        for (Chromosome chromosome : chromosomes) {
            int[] protoCluster = chromToBinToProtocluster.get(chromosome.getIndex());
            int[] globalIndex = new int[protoCluster.length];
            Arrays.fill(globalIndex, IGNORE);
            for (int i = 0; i < protoCluster.length; i++) {
                int val = protoCluster[i];
                if (val > IGNORE) {
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
        return chromToBinToProtocluster.get(c1.getIndex());
    }

    @Override
    public int[] getGlobalIndex(Chromosome c1) {
        return chromToBinToGlobalIndex.get(c1.getIndex());
    }

    @Override
    public boolean contains(Chromosome c1) {
        return chromToBinToProtocluster.containsKey(c1.getIndex()) && chromToBinToGlobalIndex.containsKey(c1.getIndex());
    }

    @Override
    public void printStatus() {
        for (Integer i : chromToBinToProtocluster.keySet()) {
            System.out.println("key1 " + i + " " + Arrays.toString(chromToBinToProtocluster.get(i)));
        }
        for (Integer i : chromToBinToGlobalIndex.keySet()) {
            System.out.println("key2 " + i + " " + Arrays.toString(chromToBinToGlobalIndex.get(i)));
        }
    }

    protected int[][] getGenomeIndices() {
        int[][] coordinates = new int[numRows][3];
        for (Integer chrIndex : chromToBinToGlobalIndex.keySet()) {
            int[] globalIndices = chromToBinToGlobalIndex.get(chrIndex);
            for (int i = 0; i < globalIndices.length; i++) {
                if (globalIndices[i] > IGNORE) {
                    int coord = globalIndices[i];
                    coordinates[coord][0] = chrIndex;
                    coordinates[coord][1] = i * resolution;
                    coordinates[coord][2] = coordinates[coord][1] + resolution;
                }
            }
        }
        return coordinates;
    }
}
