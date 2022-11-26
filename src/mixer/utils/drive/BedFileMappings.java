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

import javastraw.feature1D.GenomeWide1DList;
import javastraw.reader.basics.Chromosome;
import mixer.utils.tracks.Concensus2DTools;
import mixer.utils.tracks.SubcompartmentInterval;

import java.util.*;

public class BedFileMappings implements Mappings {
    private static final int IGNORE = -1;
    // chromosome to (bin to proto-cluster index)
    // chromosome to (bin to global index)
    private final Map<Integer, int[]> chromToBinToProtocluster = new HashMap<>();
    private final Map<Integer, int[]> chromToBinToGlobalIndex = new HashMap<>();
    private final Map<Integer, int[]> chromToDistributionForChromosome = new HashMap<>();
    private final int resolution;
    private final Chromosome[] chromosomes;
    private int numRows = 0;
    private int numCols = 0;

    private BedFileMappings(int resolution, Chromosome[] chromosomes) {
        this.resolution = resolution;
        this.chromosomes = chromosomes;
    }

    public BedFileMappings(Chromosome[] chromosomes, int resolution,
                           GenomeWide1DList<SubcompartmentInterval> clusters) {
        this(resolution, chromosomes);

        Map<String, Map<Integer, Integer>> map1 = Concensus2DTools.summarize(clusters);
        int numClusters = Concensus2DTools.getMaxId(map1) + 1;

        int counter = 0;
        for (Chromosome chrom : chromosomes) {
            int length = (int) (chrom.getLength() / resolution) + 1;
            int[] binToProtocluster = new int[length];
            Arrays.fill(binToProtocluster, -1);
            fillInClusterAssigments(binToProtocluster, map1.get("" + chrom.getIndex()), counter, resolution);
            putBinToProtoCluster(chrom, binToProtocluster);
            counter += numClusters;
        }

        calculateGlobalIndices(chromosomes);
    }

    private void fillInClusterAssigments(int[] binToProtocluster, Map<Integer, Integer> posToID,
                                         int counter, int resolution) {
        for (Integer pos : posToID.keySet()) {
            int id = posToID.get(pos);
            binToProtocluster[pos / resolution] = id + counter;
        }
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

    @Override
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
                }
            }
        }

        numRows -= badIndices.size();
    }

    @Override
    public Mappings deepCopy() {
        BedFileMappings newMapping = new BedFileMappings(resolution, chromosomes);
        newMapping.numRows = numRows;
        newMapping.numCols = numCols;
        BinMappings.copyFromAtoB(chromToBinToProtocluster, newMapping.chromToBinToProtocluster);
        BinMappings.copyFromAtoB(chromToBinToGlobalIndex, newMapping.chromToBinToGlobalIndex);
        BinMappings.copyFromAtoB(chromToDistributionForChromosome, newMapping.chromToDistributionForChromosome);
        return newMapping;
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
