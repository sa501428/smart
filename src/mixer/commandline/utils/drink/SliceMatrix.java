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

package mixer.commandline.utils.drink;

import mixer.MixerGlobals;
import mixer.commandline.utils.common.FloatMatrixTools;
import mixer.data.*;
import mixer.windowui.NormalizationType;
import org.broad.igv.feature.Chromosome;
import org.broad.igv.util.Pair;

import java.io.File;
import java.util.*;

public class SliceMatrix extends CompositeGenomeWideDensityMatrix {
    
    public static int numColumnsToPutTogether = 2;
    public static boolean USE_DERIV = false;
    
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
    public SliceMatrix(ChromosomeHandler chromosomeHandler, Dataset ds, NormalizationType norm, int resolution, File outputDirectory, Random generator, String[] referenceBedFiles) {
        super(chromosomeHandler, ds, norm, resolution, outputDirectory, generator, referenceBedFiles);
    }
    
    float[][] makeCleanScaledInterMatrix(Dataset ds) {
        
        // height/weight chromosomes
        Map<Integer, Integer> indexToFilteredLength = calculateActualLengthForChromosomes(chromosomes);
        Map<Integer, Integer> indexToCompressedLength = calculateCompressedLengthForChromosomes(indexToFilteredLength);
        
        Pair<Integer, int[][]> dimensions = calculateDimensionInterMatrix(chromosomes, indexToFilteredLength);
        Pair<Integer, int[][]> compressedDimensions = calculateDimensionInterMatrix(chromosomes, indexToCompressedLength);
        
        System.out.println(dimensions.getFirst() + " by " + compressedDimensions.getFirst());
        
        float[][] interMatrix = new float[dimensions.getFirst()][compressedDimensions.getFirst()];
        for (int i = 0; i < chromosomes.length; i++) {
            Chromosome chr1 = chromosomes[i];
            
            for (int j = i; j < chromosomes.length; j++) {
                Chromosome chr2 = chromosomes[j];
                
                final MatrixZoomData zd = HiCFileTools.getMatrixZoomData(ds, chr1, chr2, resolution);
                
                fillInChromosomeRegion(interMatrix, badIndexLocations, zd, chr1, dimensions.getSecond()[0][i], compressedDimensions.getSecond()[0][i],
                        chr2, dimensions.getSecond()[0][j], compressedDimensions.getSecond()[0][j], i == j);
                System.out.print(".");
            }
        }
        System.out.println(".");
        
        interMatrix = ExtractingOEDataUtils.simpleLogWithCleanup(interMatrix);
        System.out.println(".");
        
        if (MixerGlobals.printVerboseComments) {
            FloatMatrixTools.saveMatrixTextNumpy(new File(outputDirectory, "pre_data_matrix.npy").getAbsolutePath(), interMatrix, dimensions.getSecond());
        }
        
        MatrixCleanup matrixCleanupReduction = new MatrixCleanup(interMatrix, generator.nextLong(), outputDirectory);
        
        return matrixCleanupReduction.getSimpleCleaningOfMatrixAppendDeriv(rowIndexToIntervalMap, dimensions.getSecond(), USE_DERIV);
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
            System.out.println("size of " + key + " " + val + " was (" + initialMap.get(key) + ") num cols " + numColumnsToPutTogether);
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
    private void fillInChromosomeRegion(float[][] matrix, GenomewideBadIndexFinder badIndices, MatrixZoomData zd, Chromosome chr1, int offsetIndex1, int compressedOffsetIndex1,
                                        Chromosome chr2, int offsetIndex2, int compressedOffsetIndex2, boolean isIntra) {
    
        int lengthChr1 = chr1.getLength() / resolution + 1;
        int lengthChr2 = chr2.getLength() / resolution + 1;
    
        List<Block> blocks = null;
        try {
            if (!isIntra) {
                blocks = HiCFileTools.getAllRegionBlocks(zd, 0, lengthChr1, 0, lengthChr2, norm, isIntra);
            
                if (blocks == null || blocks.size() < 1) {
                    System.err.println("Missing Interchromosomal Data " + zd.getKey());
                    System.exit(98);
                }
            }
        } catch (Exception e) {
            e.printStackTrace();
            System.exit(99);
        }
    
        Map<Integer, Integer> genomePosToLocal1 = makeLocalIndexMap(chr1, badIndices.getBadIndices(chr1), offsetIndex1, 1);
        Map<Integer, Integer> genomePosToLocal2 = makeLocalIndexMap(chr2, badIndices.getBadIndices(chr2), offsetIndex2, 1);
        Map<Integer, Integer> genomePosToCompressed1 = makeLocalIndexMap(chr1, badIndices.getBadIndices(chr1), compressedOffsetIndex1, numColumnsToPutTogether);
        Map<Integer, Integer> genomePosToCompressed2 = makeLocalIndexMap(chr2, badIndices.getBadIndices(chr2), compressedOffsetIndex2, numColumnsToPutTogether);
    
        if (isIntra) {
            updateSubcompartmentMap(chr1, badIndices.getBadIndices(chr1), offsetIndex1, rowIndexToIntervalMap);
        }
    
        copyValuesToArea(matrix, blocks,
                genomePosToLocal1, genomePosToCompressed1, genomePosToLocal2, genomePosToCompressed2, isIntra);
    }
    
    
    private void copyValuesToArea(float[][] matrix, List<Block> blocks,
                                  Map<Integer, Integer> genomePosToLocal1, Map<Integer, Integer> genomePosToCompressed1,
                                  Map<Integer, Integer> genomePosToLocal2, Map<Integer, Integer> genomePosToCompressed2, boolean isIntra) {
        if (isIntra) {
            for (int binX : genomePosToLocal1.keySet()) {
                for (int binY : genomePosToLocal2.keySet()) {
                    matrix[genomePosToLocal1.get(binX)][genomePosToCompressed2.get(binY)] = Float.NaN;
                }
            }
        } else {
            for (Block b : blocks) {
                if (b != null) {
                    for (ContactRecord cr : b.getContactRecords()) {
                        float val = cr.getCounts();
                        if (Float.isNaN(val) || val < 1e-10 || Float.isInfinite(val)) {
                            continue;
                        }
                        
                        int binX = cr.getBinX();
                        int binY = cr.getBinY();
                        
                        if (genomePosToLocal1.containsKey(binX) && genomePosToLocal2.containsKey(binY)) {
                            matrix[genomePosToLocal1.get(binX)][genomePosToCompressed2.get(binY)] += val;
                            matrix[genomePosToLocal2.get(binY)][genomePosToCompressed1.get(binX)] += val;
                        }
                    }
                }
            }
        }
    }
    
    private void updateSubcompartmentMap(Chromosome chromosome, Set<Integer> badIndices, int offsetIndex1, Map<Integer, SubcompartmentInterval> rowIndexToIntervalMap) {
        int counter = 0;
        int chrLength = chromosome.getLength() / resolution + 1;
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
        
        int chrLength = chrom.getLength() / resolution + 1;
        for (int i = 0; i < chrLength; i++) {
            if (badIndices.contains(i)) {
                continue;
            }
            
            binToLocalMap.put(i, offsetIndex + (counter / divisor));
            counter++;
        }
        
        return binToLocalMap;
    }
}
