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
import mixer.commandline.utils.drink.kmeansfloat.Cluster;
import mixer.commandline.utils.drink.kmeansfloat.ClusterTools;
import mixer.data.ChromosomeHandler;
import mixer.data.Dataset;
import mixer.data.HiCFileTools;
import mixer.data.MatrixZoomData;
import mixer.data.feature.GenomeWideList;
import mixer.windowui.NormalizationType;
import org.broad.igv.feature.Chromosome;
import org.broad.igv.util.Pair;

import java.io.File;
import java.util.*;

public class CompositeGenomeWideDensityMatrix {
    private final NormalizationType norm;
    private final int resolution;
    private final GenomeWideList<SubcompartmentInterval> intraSubcompartments;
    private final float[][] gwCleanMatrix;
    private final Map<Integer, SubcompartmentInterval> rowIndexToIntervalMap = new HashMap<>();
    private final Map<Integer, Map<Integer, Integer>> chrIndxTorowIndexToGoldIDMap, chrIndxTorowIndexToSilverIDMap;

    private final Chromosome[] chromosomes;
    private final float threshold;
    private final int minIntervalSizeAllowed;
    private final File outputDirectory;
    private final Random generator;

    public CompositeGenomeWideDensityMatrix(ChromosomeHandler chromosomeHandler, Dataset ds, NormalizationType norm, int resolution,
                                            GenomeWideList<SubcompartmentInterval> intraSubcompartments, float oeThreshold,
                                            int minIntervalSizeAllowed, File outputDirectory, Random generator) {
        this.minIntervalSizeAllowed = minIntervalSizeAllowed;
        this.norm = norm;
        this.resolution = resolution;
        this.intraSubcompartments = intraSubcompartments;
        this.outputDirectory = outputDirectory;
        this.generator = generator;
        threshold = oeThreshold;

        chrIndxTorowIndexToGoldIDMap = DrinkUtils.createGoldStandardLookup("/Users/muhammad/Desktop/drinks/existingmethods/GSE63525_GM12878_subcompartments.bed");
        chrIndxTorowIndexToSilverIDMap = DrinkUtils.createGoldStandardLookup("/Users/muhammad/Desktop/drinks/existingmethods/SCI_GM12878_SCI_sub_compartments.bed");

        chromosomes = chromosomeHandler.getAutosomalChromosomesArray();
        gwCleanMatrix = makeCleanScaledInterMatrix(ds);
    }

    private float[][] makeCleanScaledInterMatrix(Dataset ds) {

        // height/weight chromosomes
        Map<Integer, Integer> indexToFilteredLength = calculateActualLengthForChromosomes(chromosomes, intraSubcompartments);
        Map<Integer, Integer> indexToCompressedLength = calculateCompressedLengthForChromosomes(chromosomes, intraSubcompartments);

        Pair<Integer, int[]> dimensions = calculateDimensionInterMatrix(chromosomes, indexToFilteredLength);
        Pair<Integer, int[]> compressedDimensions = calculateDimensionInterMatrix(chromosomes, indexToCompressedLength);

        //float[][] interMatrix = new float[dimensions.getFirst()][compressedDimensions.getFirst()];
        float[][] interMatrix = new float[dimensions.getFirst()][dimensions.getFirst()];

        int[] numCountsForCol = new int[compressedDimensions.getFirst()];
        for (int i = 0; i < chromosomes.length; i++) {
            Chromosome chr1 = chromosomes[i];

            for (int j = i; j < chromosomes.length; j++) {
                Chromosome chr2 = chromosomes[j];

                final MatrixZoomData zd = HiCFileTools.getMatrixZoomData(ds, chr1, chr2, resolution);

                fillInChromosomeRegion(interMatrix, numCountsForCol, ds, zd, chr1, dimensions.getSecond()[i], compressedDimensions.getSecond()[i], chr2, dimensions.getSecond()[j], compressedDimensions.getSecond()[j], i == j);
            }
        }

        interMatrix = ExtractingOEDataUtils.simpleLog(interMatrix);

        for (int i = 0; i < interMatrix.length; i++) {
            for (int j = 0; j < interMatrix[0].length; j++) {
                if (Float.isInfinite(interMatrix[i][j]) || Math.abs(interMatrix[i][j]) < 1E-10) {
                    interMatrix[i][j] = 0;
                }
            }
        }

        FloatMatrixTools.inPlaceZscoreDownCols(interMatrix);

        //FloatMatrixTools.scaleValuesByCount(interMatrix, numCountsForCol);
        //FloatMatrixTools.scaleValuesInPlaceByCountAndZscore(interMatrix, numCountsForCol);

        return interMatrix;
    }

    private Map<Integer, Integer> calculateActualLengthForChromosomes(Chromosome[] chromosomes, GenomeWideList<SubcompartmentInterval> intraSubcompartments) {
        Map<Integer, Integer> indexToFilteredLength = new HashMap<>();
        for (Chromosome chrom : chromosomes) {
            int val = 0;
            List<SubcompartmentInterval> intervals = intraSubcompartments.getFeatures("" + chrom.getIndex());
            for (SubcompartmentInterval interval : intervals) {
                int numCols = interval.getWidthForResolution(resolution);
                if (numCols >= minIntervalSizeAllowed) {
                    val += numCols;
                }
            }
            indexToFilteredLength.put(chrom.getIndex(), val);
        }

        return indexToFilteredLength;
    }

    private Map<Integer, Integer> calculateCompressedLengthForChromosomes(Chromosome[] chromosomes, GenomeWideList<SubcompartmentInterval> intraSubcompartments) {
        Map<Integer, Integer> indexToCompressedLength = new HashMap<>();
        for (Chromosome chrom : chromosomes) {
            int val = 0;
            List<SubcompartmentInterval> intervals = intraSubcompartments.getFeatures("" + chrom.getIndex());
            for (SubcompartmentInterval interval : intervals) {
                int numCols = interval.getWidthForResolution(resolution);
                if (numCols >= minIntervalSizeAllowed) {
                    val += 1;
                }
            }
            indexToCompressedLength.put(chrom.getIndex(), val);
        }

        return indexToCompressedLength;
    }

    private Pair<Integer, int[]> calculateDimensionInterMatrix(Chromosome[] chromosomes, Map<Integer, Integer> indexToFilteredLength) {
        int total = 0;
        int[] indices = new int[chromosomes.length];

        for (int i = 0; i < chromosomes.length; i++) {
            total += indexToFilteredLength.get(chromosomes[i].getIndex());
            if (i < chromosomes.length - 1) {
                indices[i + 1] = total;
            }
        }

        return new Pair<>(total, indices);
    }

    private void fillInChromosomeRegion(float[][] matrix, int[] numCountsForCol, Dataset ds, MatrixZoomData zd, Chromosome chr1, int offsetIndex1, int compressedOffsetIndex1,
                                        Chromosome chr2, int offsetIndex2, int compressedOffsetIndex2, boolean isIntra) {

        int lengthChr1 = chr1.getLength() / resolution + 1;
        int lengthChr2 = chr2.getLength() / resolution + 1;
        List<SubcompartmentInterval> intervals1 = intraSubcompartments.getFeatures("" + chr1.getIndex());
        List<SubcompartmentInterval> intervals2 = intraSubcompartments.getFeatures("" + chr2.getIndex());

        if (intervals1.size() == 0 || intervals2.size() == 0) {
            System.err.println("Missing interval data " + zd.getKey());
            System.exit(97);
        }
        float[][] allDataForRegion = null;
        try {

            if (!isIntra) {
                allDataForRegion = HiCFileTools.extractLocalBoundedRegioFloatMatrix(zd, 0,
                        lengthChr1, 0, lengthChr2, lengthChr1, lengthChr2, norm, isIntra);

                if (allDataForRegion == null) {
                    System.err.println("Missing Interchromosomal Data " + zd.getKey());
                    System.exit(98);
                }

                for (int i = 0; i < allDataForRegion.length; i++) {
                    for (int j = 0; j < allDataForRegion[0].length; j++) {
                        if (Float.isNaN(allDataForRegion[i][j]) || Float.isInfinite(allDataForRegion[i][j]) || allDataForRegion[i][j] < 1E-10) {
                            allDataForRegion[i][j] = 0;
                        }
                    }
                }

            }
        } catch (Exception e) {
            e.printStackTrace();
            System.exit(99);
        }

        int internalOffset1 = offsetIndex1;
        for (SubcompartmentInterval interv1 : intervals1) {
            int numRows = interv1.getWidthForResolution(resolution);
            if (numRows >= minIntervalSizeAllowed) {
                int internalOffset2 = offsetIndex2;
                for (SubcompartmentInterval interv2 : intervals2) {
                    int numCols = interv2.getWidthForResolution(resolution);
                    if (numCols >= minIntervalSizeAllowed) {
                        copyValuesToArea(matrix, interv1, internalOffset1, numRows, interv2, internalOffset2, numCols, isIntra, allDataForRegion);
                        //updateMasterMatrixWithRegionalDensities(matrix, interv1, internalOffset1, numRows, interv2, internalOffset2, numCols, isIntra, allDataForRegion);
                        internalOffset2 += numCols;
                    }
                }
                internalOffset1 += numRows;
            }
        }
    }

    private void copyValuesToArea(float[][] matrix, SubcompartmentInterval interv1, int offsetIndex1, int numRows,
                                  SubcompartmentInterval interv2, int offsetIndex2, int numCols, boolean isIntra, float[][] allDataForRegion) {

        int binXStart = interv1.getX1() / resolution;
        int binXEnd = interv1.getX2() / resolution;

        int binYStart = interv2.getX1() / resolution;
        int binYEnd = interv2.getX2() / resolution;

        if (isIntra) {
            for (int i = 0; i < numRows; i++) {
                int newX1 = (binXStart + i) * resolution;
                SubcompartmentInterval newRInterval = new SubcompartmentInterval(interv1.getChrIndex(), interv1.getChrName(), newX1, newX1 + resolution, interv1.getClusterID());
                rowIndexToIntervalMap.put(offsetIndex1 + i, newRInterval);
                for (int j = 0; j < numCols; j++) {
                    int newY1 = (binYStart + i) * resolution;
                    SubcompartmentInterval newCInterval = new SubcompartmentInterval(interv2.getChrIndex(), interv2.getChrName(), newY1, newY1 + resolution, interv1.getClusterID());
                    rowIndexToIntervalMap.put(offsetIndex2 + j, newCInterval);

                    matrix[offsetIndex1 + i][offsetIndex2 + j] = Float.NaN;
                }
            }
        } else {
            float[][] tempCopy = new float[numRows][numCols];
            for (int i = binXStart; i < binXEnd; i++) {
                for (int j = binYStart; j < binYEnd; j++) {
                    tempCopy[i - binXStart][j - binYStart] = allDataForRegion[i][j];
                }
            }

            //tempCopy = augmentMatrixRegion(tempCopy);

            for (int i = 0; i < numRows; i++) {
                int newX1 = (binXStart + i) * resolution;
                SubcompartmentInterval newRInterval = new SubcompartmentInterval(interv1.getChrIndex(), interv1.getChrName(), newX1, newX1 + resolution, interv1.getClusterID());
                rowIndexToIntervalMap.put(offsetIndex1 + i, newRInterval);
                for (int j = 0; j < numCols; j++) {
                    int newY1 = (binYStart + i) * resolution;
                    SubcompartmentInterval newCInterval = new SubcompartmentInterval(interv2.getChrIndex(), interv2.getChrName(), newY1, newY1 + resolution, interv1.getClusterID());
                    rowIndexToIntervalMap.put(offsetIndex2 + j, newCInterval);

                    matrix[offsetIndex1 + i][offsetIndex2 + j] = tempCopy[i][j];
                }
            }

            for (int i = 0; i < numRows; i++) {
                for (int j = 0; j < numCols; j++) {
                    matrix[offsetIndex2 + j][offsetIndex1 + i] = tempCopy[i][j];
                }
            }
        }
    }

    private float[][] augmentMatrixRegion(float[][] input) {

        List<Float> values = new ArrayList<>();
        for (int i = 0; i < input.length; i++) {
            for (int j = 0; j < input[i].length; j++) {
                float val = input[i][j];
                if (val > 0) {
                    values.add(val);
                }
            }
        }

        int numNonZero = values.size();
        float percentFilled = numNonZero / ((float) (input.length * input[0].length));

        if (numNonZero > 3 && percentFilled > .3) {
            Collections.sort(values);
            if (values.get(values.size() - 1) - values.get(0) < .01) {
                System.err.println("Weird - not enhancing");
            } else {
                float[][] output = new float[input.length][input[0].length];
                float mean = getMean(values);
                float stdDev = getSampleStdDev(values, mean);

                for (int i = 0; i < input.length; i++) {
                    for (int j = 0; j < input[i].length; j++) {
                        float val = input[i][j];
                        if (val < 1e-3) {
                            input[i][j] = getNewRandVal(mean, stdDev);
                        }
                    }
                }
            }
        }

        return input;
    }

    private float getNewRandVal(float mean, float stdDev) {
        float zVal = 2 * generator.nextFloat() - 1;
        return mean + stdDev * zVal;
    }

    private float getMean(List<Float> values) {
        float total = 0;
        for (Float val : values) {
            total += val;
        }
        return total / values.size();
    }

    private float getSampleStdDev(List<Float> values, float mean) {
        double total = 0;
        for (Float val : values) {
            float diff = val - mean;
            total += (diff * diff);
        }
        return (float) (total / (values.size() - 1));
    }


    public synchronized Pair<Double, List<int[][]>> processGWKmeansResult(Cluster[] clusters, GenomeWideList<SubcompartmentInterval> subcompartments) {

        Set<SubcompartmentInterval> subcompartmentIntervals = new HashSet<>();
        if (MixerGlobals.printVerboseComments) {
            System.out.println("GW Composite data vs clustered into " + clusters.length + " clusters");
        }

        double withinClusterSumOfSquares = 0;
        int genomewideCompartmentID = 0;

        int[][] ids = new int[1][clusters.length];
        int[][] idsForIndex = new int[3][gwCleanMatrix.length];
        Arrays.fill(idsForIndex[0], -1);
        Arrays.fill(idsForIndex[1], -1);
        Arrays.fill(idsForIndex[2], -1);

        for (int z = 0; z < clusters.length; z++) {
            Cluster cluster = clusters[z];
            int currentClusterID = ++genomewideCompartmentID;
            ids[0][z] = currentClusterID;

            if (MixerGlobals.printVerboseComments) {
                System.out.println("Size of cluster " + currentClusterID + " - " + cluster.getMemberIndexes().length);
            }

            for (int i : cluster.getMemberIndexes()) {
                withinClusterSumOfSquares += ClusterTools.getPositiveVectorMSEDifference(cluster.getCenter(), gwCleanMatrix[i]);

                try {
                    SubcompartmentInterval interv;

                    if (rowIndexToIntervalMap.containsKey(i)) {
                        interv = rowIndexToIntervalMap.get(i);
                        if (interv == null) continue; // probably a zero row

                        int chrIndex = interv.getChrIndex();
                        String chrName = interv.getChrName();
                        int x1 = interv.getX1();
                        int x2 = interv.getX2();

                        subcompartmentIntervals.add(
                                new SubcompartmentInterval(chrIndex, chrName, x1, x2, currentClusterID));

                        idsForIndex[2][i] = currentClusterID;

                        if (chrIndxTorowIndexToGoldIDMap.containsKey(chrIndex)) {
                            if (chrIndxTorowIndexToGoldIDMap.get(chrIndex).containsKey(x1)) {
                                idsForIndex[0][i] = chrIndxTorowIndexToGoldIDMap.get(chrIndex).get(x1);
                            }
                        }

                        if (chrIndxTorowIndexToSilverIDMap.containsKey(chrIndex)) {
                            if (chrIndxTorowIndexToSilverIDMap.get(chrIndex).containsKey(x1)) {
                                idsForIndex[1][i] = chrIndxTorowIndexToSilverIDMap.get(chrIndex).get(x1);
                            }
                        }

                    } else {
                        System.err.println("********* is weird error?");
                    }

                } catch (Exception e) {
                    e.printStackTrace();
                    System.exit(87);
                }
            }
        }

        if (MixerGlobals.printVerboseComments) {
            System.out.println("Final WCSS " + withinClusterSumOfSquares);
        }

        subcompartments.addAll(new ArrayList<>(subcompartmentIntervals));
        DrinkUtils.reSort(subcompartments);

        List<int[][]> outputs = new ArrayList<>();
        outputs.add(ids);
        outputs.add(idsForIndex);

        return new Pair<>(withinClusterSumOfSquares, outputs);
    }

    public float[][] getCleanedData() {
        return gwCleanMatrix;
    }

    public int getLength() {
        return gwCleanMatrix.length;
    }

    public int getWidth() {
        return gwCleanMatrix[0].length;
    }

    public void exportData() {
        System.out.println(getLength() + " -v- " + getWidth());
        FloatMatrixTools.saveMatrixTextNumpy(new File(outputDirectory, "data_matrix.npy").getAbsolutePath(), getCleanedData());
    }

    public void appendDataAlongExistingRows(CompositeGenomeWideDensityMatrix additionalData) {
        if (getLength() != additionalData.getLength()) {
            System.err.println("***************************************\n" +
                    "Dimension mismatch: " + getLength() + " != " + additionalData.getLength());
        } else {

        }
    }
}
