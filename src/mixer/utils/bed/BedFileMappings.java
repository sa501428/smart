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

import javastraw.reader.Dataset;
import javastraw.reader.basics.Chromosome;
import javastraw.reader.basics.ChromosomeHandler;
import javastraw.reader.norm.NormalizationVector;
import javastraw.reader.type.HiCZoom;
import javastraw.reader.type.NormalizationType;
import mixer.utils.tracks.SubcompartmentInterval;

import java.util.Arrays;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

public class BedFileMappings {

    private final int[] offsets;
    private final int[] binIndexToChromIndex;
    private final int[] indexToClusterID;
    private final Map<Integer, Integer> bedFileCIDToColumnID = new HashMap<>();
    private final Map<Integer, int[]> chromIndexToDistributionForChromosome = new HashMap<>();
    private int numCols = 0;
    private int numRows = 0;


    public BedFileMappings(String bedpath, ChromosomeHandler handler, int resolution, Dataset ds,
                           NormalizationType norm) {

        offsets = generateFromChromosomes(handler, resolution);
        numRows = offsets[offsets.length - 1];
        binIndexToChromIndex = new int[numRows];
        indexToClusterID = new int[numRows];
        Arrays.fill(binIndexToChromIndex, -1);
        Arrays.fill(indexToClusterID, -1);

        SubcompartmentBedFile bedFile = new SubcompartmentBedFile(bedpath, handler);
        if (bedFile.getNumTotalFeatures() < 1) {
            System.err.println("bed file is empty or incorrect path provided.");
            System.exit(3);
        }

        Chromosome[] chromosomes = handler.getAutosomalChromosomesArray();
        for (int z = 0; z < chromosomes.length; z++) {
            Chromosome chromosome = chromosomes[z];

            NormalizationVector nv = ds.getNormalizationVector(chromosome.getIndex(), new HiCZoom(resolution), norm);
            double[] vec = nv.getData().getValues().get(0);

            List<SubcompartmentInterval> intervals = bedFile.get(chromosome.getIndex());
            int offset = offsets[z];
            for (SubcompartmentInterval si : intervals) {
                for (int x = si.getX1() / resolution; x < si.getX2() / resolution; x++) {
                    if (vec[x] > 0) {
                        int currIndex = offset + x;
                        binIndexToChromIndex[currIndex] = chromosome.getIndex();
                        int tempCID = si.getClusterID();
                        if (!bedFileCIDToColumnID.containsKey(tempCID)) {
                            bedFileCIDToColumnID.put(tempCID, numCols);
                            numCols++;
                        }
                        int colID = bedFileCIDToColumnID.get(tempCID);
                        indexToClusterID[currIndex] = colID;
                    }
                }
            }
        }
        System.out.println("Bed file loaded");

        for (Chromosome chromosome : chromosomes) {
            NormalizationVector nv = ds.getNormalizationVector(chromosome.getIndex(), new HiCZoom(resolution), norm);
            double[] vec = nv.getData().getValues().get(0);

            int[] counts = new int[numCols];
            for (SubcompartmentInterval si : bedFile.get(chromosome.getIndex())) {
                int colID = bedFileCIDToColumnID.get(si.getClusterID());
                for (int x = si.getX1() / resolution; x < si.getX2() / resolution; x++) {
                    if (vec[x] > 0) {
                        counts[colID]++;
                    }
                }
            }
            chromIndexToDistributionForChromosome.put(chromosome.getIndex(), counts);
        }
    }

    private int[] generateFromChromosomes(ChromosomeHandler handler, int resolution) {
        Chromosome[] chromosomes = handler.getAutosomalChromosomesArray();
        int[] offsets = new int[chromosomes.length + 1];
        int offset = 0;
        for (int z = 0; z < chromosomes.length; z++) {
            Chromosome chromosome = chromosomes[z];
            offsets[z] = offset;
            offset += (int) (chromosome.getLength() / resolution) + 1;
        }
        offsets[chromosomes.length] = offset;
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

    public int[] getIndexToClusterID() {
        return indexToClusterID;
    }

    public Map<Integer, int[]> getChromIndexToDistributionForChromosome() {
        return chromIndexToDistributionForChromosome;
    }

    public int[] getBinIndexToChromIndex() {
        return binIndexToChromIndex;
    }
}
