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
import mixer.utils.tracks.SubcompartmentInterval;

import java.util.*;

public class BinMappings implements Mappings {
    // chromosome to (bin to proto-cluster index)
    // chromosome to (bin to global index)
    private final Map<Integer, int[]> chromToBinToProtocluster = new HashMap<>();
    private final Map<Integer, int[]> chromToBinToGlobalIndex = new HashMap<>();
    private final Map<Integer, int[]> chromToDistributionForChromosome = new HashMap<>();
    private int numRows = 0;
    private int numCols = 0;
    protected static final int IGNORE = -1;
    private final int resolution;
    private final Chromosome[] chromosomes;

    public BinMappings(int resolution, Chromosome[] chromosomes) {
        this.resolution = resolution;
        this.chromosomes = chromosomes;
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

    public static void copyFromAtoB(Map<Integer, int[]> a, Map<Integer, int[]> b) {
        for (Integer key : a.keySet()) {
            b.put(key, FloatMatrixTools.deepClone(a.get(key)));
        }
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
    public int getResolution() {
        return resolution;
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

    @Override
    public Chromosome[] getChromosomes() {
        return chromosomes;
    }

    @Override
    public void updateInternalDataStructures(Set<Integer> badIndices) {

        List<Integer> badIndx = new ArrayList<>(badIndices);
        Collections.sort(badIndx);
        Collections.reverse(badIndx);

        for (int bi : badIndx) {
            for (int[] array : chromToBinToGlobalIndex.values()) {
                for (int i = 0; i < array.length; i++) {

                    if (array[i] > bi) {
                        array[i]--;
                    } else if (array[i] == bi) {
                        array[i] = IGNORE;
                    }
                    // else is < bi, and isn't affected
                }
            }
        }

        numRows -= badIndices.size();
    }

    @Override
    public Mappings deepCopy() {
        BinMappings newMapping = new BinMappings(resolution, chromosomes);
        newMapping.numRows = numRows;
        newMapping.numCols = numCols;
        copyFromAtoB(chromToBinToProtocluster, newMapping.chromToBinToProtocluster);
        copyFromAtoB(chromToBinToGlobalIndex, newMapping.chromToBinToGlobalIndex);
        copyFromAtoB(chromToDistributionForChromosome, newMapping.chromToDistributionForChromosome);
        return newMapping;
    }

    @Override
    public int[] getProtoclusterAssignments(Chromosome chrom) {
        return chromToBinToProtocluster.get(chrom.getIndex());
    }

    @Override
    public Map<Integer, SubcompartmentInterval> populateRowIndexToIntervalMap() {
        Map<Integer, SubcompartmentInterval> map = new HashMap<>();
        int resolution = getResolution();
        Chromosome[] chromosomes = getChromosomes();
        for (Chromosome chromosome : chromosomes) {
            int maxGenomeLen = (int) chromosome.getLength();
            int[] globalIndices = getGlobalIndex(chromosome);
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
        return map;
    }
}
