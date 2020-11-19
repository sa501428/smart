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

package mixer.utils.slice.matrices;

import javastraw.reader.ChromosomeHandler;
import javastraw.reader.Dataset;
import javastraw.reader.HiCFileTools;
import javastraw.reader.MatrixZoomData;
import javastraw.reader.basics.Block;
import javastraw.reader.basics.Chromosome;
import javastraw.reader.basics.ContactRecord;
import javastraw.type.NormalizationType;
import mixer.MixerGlobals;
import mixer.utils.common.Pair;
import mixer.utils.similaritymeasures.SimilarityMetric;
import mixer.utils.slice.cleaning.BadIndexFinder;
import mixer.utils.slice.cleaning.IndexOrderer;
import mixer.utils.slice.structures.SubcompartmentInterval;

import java.io.File;
import java.util.*;

public class SliceMatrix extends CompositeGenomeWideDensityMatrix {

    public static int numColumnsToPutTogether = 2;

    /**
     * for SLICE, minIntervalSize become how many bins to collapse together (e.g. 5 means 5 bins together)
     *
     * @param chromosomeHandler
     * @param ds
     * @param norm
     * @param resolution
     * @param outputDirectory
     * @param generator
     * @param referenceBedFiles
     */
    public SliceMatrix(ChromosomeHandler chromosomeHandler, Dataset ds, NormalizationType norm, int resolution, File outputDirectory,
                       Random generator, String[] referenceBedFiles, BadIndexFinder badIndexLocations,
                       int datasetIndex, SimilarityMetric metric) {
        super(chromosomeHandler, ds, norm, resolution, outputDirectory, generator, referenceBedFiles,
                badIndexLocations, datasetIndex, metric);
    }

    float[][] makeCleanScaledInterMatrix(Dataset ds) {

        // height/weight chromosomes
        Map<Integer, Integer> indexToFilteredLength = calculateActualLengthForChromosomes(chromosomes);
        Map<Integer, Integer> indexToCompressedLength = calculateCompressedLengthForChromosomes(indexToFilteredLength);

        IndexOrderer orderer = null;
        if (numColumnsToPutTogether > 1) {
            orderer = new IndexOrderer(ds, chromosomes, resolution, norm, numColumnsToPutTogether,
                    badIndexLocations);
            indexToCompressedLength = calculateCompressedLengthForChromosomes(orderer.getIndexToRearrangedLength());
        }

        Pair<Integer, int[][]> dimensions = calculateDimensionInterMatrix(chromosomes, indexToFilteredLength);
        Pair<Integer, int[][]> compressedDimensions = calculateDimensionInterMatrix(chromosomes, indexToCompressedLength);

        if (MixerGlobals.printVerboseComments) {
            System.out.println(dimensions.getFirst() + " by " + compressedDimensions.getFirst());
        }

        float[][] interMatrix = new float[dimensions.getFirst()][compressedDimensions.getFirst()];
        try {
            if (interMatrix.length == 0 || interMatrix[0].length == 0) {
                System.err.println("No matrix created; map is likely too sparse. " +
                        "Try a lower resolution or higher compression window.");
                System.exit(9);
            }
        } catch (Exception e) {
            System.err.println("No matrix created; map is likely too sparse. Try a lower resolution.");
            System.exit(9);
        }

        System.out.println("Indexing complete");

        for (int i = 0; i < chromosomes.length; i++) {
            Chromosome chr1 = chromosomes[i];

            for (int j = i; j < chromosomes.length; j++) {
                Chromosome chr2 = chromosomes[j];

                final MatrixZoomData zd = HiCFileTools.getMatrixZoomData(ds, chr1, chr2, resolution);

                fillInChromosomeRegion(interMatrix, badIndexLocations, zd, chr1, dimensions.getSecond()[0][i], compressedDimensions.getSecond()[0][i],
                        chr2, dimensions.getSecond()[0][j], compressedDimensions.getSecond()[0][j], i == j, orderer);
                System.out.print(".");
            }
        }
        System.out.println(".");

        return interMatrix;
    }

    protected Map<Integer, Integer> calculateActualLengthForChromosomes(Chromosome[] chromosomes) {
        Map<Integer, Integer> indexToFilteredLength = new HashMap<>();
        for (Chromosome chrom : chromosomes) {
            indexToFilteredLength.put(chrom.getIndex(), (int) Math.ceil((float) chrom.getLength() / resolution) - badIndexLocations.getBadIndices(chrom).size());
        }

        return indexToFilteredLength;
    }

    /**
     * @param initialMap
     * @return
     */
    private Map<Integer, Integer> calculateCompressedLengthForChromosomes(Map<Integer, Integer> initialMap) {
        Map<Integer, Integer> indexToCompressedLength = new HashMap<>();
        for (Integer key : initialMap.keySet()) {
            int val = (int) Math.ceil((double) initialMap.get(key) / numColumnsToPutTogether);
            //System.out.println("size of " + key + " " + val + " was (" + initialMap.get(key) + ") num cols " + numColumnsToPutTogether);
            indexToCompressedLength.put(key, val);
        }

        return indexToCompressedLength;
    }

    /**
     * @param matrix
     * @param zd
     * @param chr1
     * @param offsetIndex1
     * @param compressedOffsetIndex1
     * @param chr2
     * @param offsetIndex2
     * @param compressedOffsetIndex2
     * @param isIntra
     */
    private void fillInChromosomeRegion(float[][] matrix, BadIndexFinder badIndices, MatrixZoomData zd, Chromosome chr1, int offsetIndex1, int compressedOffsetIndex1,
                                        Chromosome chr2, int offsetIndex2, int compressedOffsetIndex2, boolean isIntra,
                                        IndexOrderer orderer) {

        int lengthChr1 = (int) (chr1.getLength() / resolution + 1);
        int lengthChr2 = (int) (chr2.getLength() / resolution + 1);

        List<Block> blocks = null;
        try {
            if (!isIntra) {
                blocks = HiCFileTools.getAllRegionBlocks(zd, 0, lengthChr1, 0, lengthChr2, norm, false);

                if (blocks.size() < 1) {
                    System.err.println("Missing Interchromosomal Data " + zd.getKey());
                    System.exit(98);
                }
            }
        } catch (Exception e) {
            e.printStackTrace();
            System.exit(99);
        }

        Map<Integer, Integer> rowPosChrom1 = makeLocalIndexMap(chr1, badIndices.getBadIndices(chr1), offsetIndex1, 1);
        Map<Integer, Integer> rowPosChrom2 = makeLocalIndexMap(chr2, badIndices.getBadIndices(chr2), offsetIndex2, 1);

        Map<Integer, Integer> colPosChrom1, colPosChrom2;
        if (orderer != null) {
            colPosChrom1 = makeLocalReorderedIndexMap(chr1, badIndices.getBadIndices(chr1), compressedOffsetIndex1,
                    numColumnsToPutTogether, orderer.get(chr1));
            colPosChrom2 = makeLocalReorderedIndexMap(chr2, badIndices.getBadIndices(chr2), compressedOffsetIndex2,
                    numColumnsToPutTogether, orderer.get(chr2));
        } else {
            colPosChrom1 = makeLocalIndexMap(chr1, badIndices.getBadIndices(chr1), compressedOffsetIndex1, numColumnsToPutTogether);
            colPosChrom2 = makeLocalIndexMap(chr2, badIndices.getBadIndices(chr2), compressedOffsetIndex2, numColumnsToPutTogether);
        }

        if (isIntra) {
            updateSubcompartmentMap(chr1, badIndices.getBadIndices(chr1), offsetIndex1, rowIndexToIntervalMap);
        }

        copyValuesToArea(matrix, blocks, badIndices,
                rowPosChrom1, colPosChrom1, rowPosChrom2, colPosChrom2, isIntra);
    }


    private void copyValuesToArea(float[][] matrix, List<Block> blocks, BadIndexFinder badIndices,
                                  Map<Integer, Integer> rowPosChrom1, Map<Integer, Integer> colPosChrom1,
                                  Map<Integer, Integer> rowPosChrom2, Map<Integer, Integer> colPosChrom2, boolean isIntra) {
        if (isIntra) {
            for (int binX : rowPosChrom1.keySet()) {
                for (int binY : rowPosChrom2.keySet()) {
                    matrix[rowPosChrom1.get(binX)][colPosChrom2.get(binY)] = Float.NaN;
                }
            }
        } else {
            for (Block b : blocks) {
                if (b != null) {
                    for (ContactRecord cr : b.getContactRecords()) {
                        float val = cr.getCounts();
                        if (badIndices.getExceedsAllowedGlobalValue(val, datasetIndex)) {
                            val = Float.NaN;
                        }

                        int binX = cr.getBinX();
                        int binY = cr.getBinY();

                        if (rowPosChrom1.containsKey(binX) && rowPosChrom2.containsKey(binY)) {
                            matrix[rowPosChrom1.get(binX)][colPosChrom2.get(binY)] += val;
                            matrix[rowPosChrom2.get(binY)][colPosChrom1.get(binX)] += val;
                        }
                    }
                }
            }
        }
    }

    private void updateSubcompartmentMap(Chromosome chromosome, Set<Integer> badIndices, int offsetIndex1, Map<Integer, SubcompartmentInterval> rowIndexToIntervalMap) {
        int counter = 0;
        int chrLength = (int) (chromosome.getLength() / resolution + 1);
        for (int i = 0; i < chrLength; i++) {
            if (badIndices.contains(i)) {
                continue;
            }
            int newX = i * resolution;
            SubcompartmentInterval newRInterval = new SubcompartmentInterval(chromosome.getIndex(), chromosome.getName(), newX, newX + resolution, counter);
            rowIndexToIntervalMap.put(offsetIndex1 + (counter), newRInterval);
            counter++;
        }
    }

    private Map<Integer, Integer> makeLocalIndexMap(Chromosome chrom, Set<Integer> badIndices, int offsetIndex, int divisor) {
        Map<Integer, Integer> binToLocalMap = new HashMap<>();
        int counter = 0;

        int chrLength = (int) (chrom.getLength() / resolution + 1);
        for (int i = 0; i < chrLength; i++) {
            if (badIndices.contains(i)) {
                continue;
            }

            binToLocalMap.put(i, offsetIndex + (counter / divisor));
            counter++;
        }

        return binToLocalMap;
    }

    private Map<Integer, Integer> makeLocalReorderedIndexMap(Chromosome chrom, Set<Integer> badIndices, int offsetIndex, int divisor, int[] newOrder) {
        Map<Integer, Integer> binToLocalMap = new HashMap<>();

        int chrLength = (int) (chrom.getLength() / resolution + 1);
        for (int i = 0; i < chrLength; i++) {
            if (badIndices.contains(i)) {
                continue;
            }

            binToLocalMap.put(i, offsetIndex + (newOrder[i] / divisor));
        }

        return binToLocalMap;
    }
}
