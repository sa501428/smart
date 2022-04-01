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

package mixer.utils.bed;

import javastraw.reader.basics.Chromosome;
import javastraw.reader.basics.ChromosomeHandler;
import mixer.utils.slice.structures.SubcompartmentInterval;

import java.util.HashMap;
import java.util.List;
import java.util.Map;

public class BedFileMappings {

    private final int[] offsets;
    private final Map<Integer, Integer> binToColumnNumber;
    private final Map<Integer, Integer> clusterNumToColIndex;
    private int numCols = 0;
    private int numRows = 0;


    public BedFileMappings(String bedpath, ChromosomeHandler handler, int resolution) {

        offsets = generateFromChromosomes(handler, resolution);

        SubcompartmentBedFile bedFile = new SubcompartmentBedFile(bedpath, handler);
        if (bedFile.getNumTotalFeatures() < 1) {
            System.err.println("bed file is empty or incorrect path provided.");
            System.exit(3);
        }

        binToColumnNumber = new HashMap<>();
        clusterNumToColIndex = new HashMap<>();

        Chromosome[] chromosomes = handler.getAutosomalChromosomesArray();
        for (int z = 0; z < chromosomes.length; z++) {
            Chromosome chromosome = chromosomes[z];
            List<SubcompartmentInterval> intervals = bedFile.get(chromosome.getIndex());
            int offset = offsets[z];
            for (SubcompartmentInterval si : intervals) {
                for (int x = si.getX1() / resolution; x < si.getX2() / resolution; x++) {
                    int currIndex = offset + x;
                    int cluster = si.getClusterID();
                    if (!clusterNumToColIndex.containsKey(cluster)) {
                        clusterNumToColIndex.put(cluster, numCols);
                        numCols++;
                    }
                    int col = clusterNumToColIndex.get(cluster);
                    binToColumnNumber.put(currIndex, col);
                }
            }
        }
        System.out.println("Bed file loaded");
    }

    private int[] generateFromChromosomes(ChromosomeHandler handler, int resolution) {
        Chromosome[] chromosomes = handler.getAutosomalChromosomesArray();
        int[] offsets = new int[chromosomes.length];
        int offset = 0;
        for (int z = 0; z < chromosomes.length; z++) {
            Chromosome chromosome = chromosomes[z];
            offsets[z] = offset;
            offset += (int) (chromosome.getLength() / resolution) + 1;
        }
        numRows = offset;
        return offsets;
    }

    public int getNumRows() {
        return numRows;
    }

    public int getNumCols() {
        return numCols;
    }

    public int[] getOffsets() {
        return offsets;
    }

    public Map<Integer, Integer> getBinToClusterID() {
        return binToColumnNumber;
    }
}
