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
    private final Chromosome[] chromosomes;
    private final float threshold;
    private final int minIntervalSizeAllowed;
    private final File outputDirectory;

    public CompositeGenomeWideDensityMatrix(ChromosomeHandler chromosomeHandler, Dataset ds, NormalizationType norm, int resolution,
                                            GenomeWideList<SubcompartmentInterval> intraSubcompartments, float oeThreshold,
                                            int minIntervalSizeAllowed, File outputDirectory) {
        this.minIntervalSizeAllowed = minIntervalSizeAllowed;
        this.norm = norm;
        this.resolution = resolution;
        this.intraSubcompartments = intraSubcompartments;
        this.outputDirectory = outputDirectory;
        threshold = oeThreshold;
        chromosomes = chromosomeHandler.getAutosomalChromosomesArray();
        gwCleanMatrix = makeCleanScaledInterMatrix(ds);
    }

    private float[][] makeCleanScaledInterMatrix(Dataset ds) {

        // height/weight chromosomes
        Map<Integer, Integer> indexToFilteredLength = calculateActualLengthForChromosomes(chromosomes, intraSubcompartments);
        Map<Integer, Integer> indexToCompressedLength = calculateCompressedLengthForChromosomes(chromosomes, intraSubcompartments);

        Pair<Integer, int[]> dimensions = calculateDimensionInterMatrix(chromosomes, indexToFilteredLength);
        Pair<Integer, int[]> compressedDimensions = calculateDimensionInterMatrix(chromosomes, indexToCompressedLength);

        float[][] interMatrix = new float[dimensions.getFirst()][compressedDimensions.getFirst()];
        int[] numCountsForCol = new int[compressedDimensions.getFirst()];
        for (int i = 0; i < chromosomes.length; i++) {
            Chromosome chr1 = chromosomes[i];

            for (int j = i; j < chromosomes.length; j++) {
                Chromosome chr2 = chromosomes[j];

                final MatrixZoomData zd = HiCFileTools.getMatrixZoomData(ds, chr1, chr2, resolution);

                fillInChromosomeRegion(interMatrix, numCountsForCol, ds, zd, chr1, dimensions.getSecond()[i], compressedDimensions.getSecond()[i], chr2, dimensions.getSecond()[j], compressedDimensions.getSecond()[j], i == j);
            }
        }

        FloatMatrixTools.scaleValuesByCount(interMatrix, numCountsForCol);
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
            if (isIntra) {
                allDataForRegion = HiCFileTools.getRealOEMatrixForChromosomeFloatMatrix(ds, zd, chr1, resolution, norm, threshold,
                        ExtractingOEDataUtils.ThresholdType.LOGEO, //LOG_OE_PLUS_AVG_BOUNDED_MADE_POS
                        true);
            } else {
                float[][] allDataForRegionMatrix = HiCFileTools.extractLocalBoundedRegioFloatMatrix(zd, 0,
                        lengthChr1, 0, lengthChr2, lengthChr1, lengthChr2, norm, isIntra);
                allDataForRegion = ExtractingOEDataUtils.logOEP1(allDataForRegionMatrix, zd.getAverageCount());
            }
        } catch (Exception e) {
            e.printStackTrace();
            System.exit(99);
        }

        if (allDataForRegion == null) {
            System.err.println("Missing Interchromosomal Data " + zd.getKey());
            System.exit(98);
        }

        for (int i = 0; i < allDataForRegion.length; i++) {
            for (int j = 0; j < allDataForRegion[0].length; j++) {
                if (Float.isNaN(allDataForRegion[i][j]) || Float.isInfinite(allDataForRegion[i][j]) || Math.abs(allDataForRegion[i][j]) < 1E-30) {
                    allDataForRegion[i][j] = 0;
                }
            }
        }

        Map<String, Float> allAreaBetweenClusters = new HashMap<>();
        Map<String, Float> allContactsBetweenClusters = new HashMap<>();

        for (SubcompartmentInterval interv1 : intervals1) {
            if (interv1.getWidthForResolution(resolution) >= minIntervalSizeAllowed) {
                Integer id1 = interv1.getClusterID();
                for (SubcompartmentInterval interv2 : intervals2) {
                    if (interv2.getWidthForResolution(resolution) >= minIntervalSizeAllowed) {
                        Integer id2 = interv2.getClusterID();
                        String regionKey = id1 + "-" + id2;

                        //updateCountsAndArea(allDataForRegion, interv1, interv2, regionKey, allContactsBetweenClusters, allAreaBetweenClusters);
                        updateNonZeroCountsAndArea(allDataForRegion, interv1, interv2, regionKey, allContactsBetweenClusters, allAreaBetweenClusters);

                    }
                }
            }
        }


        allDataForRegion = null;
        System.gc();

        Map<String, Float> densityBetweenClusters = getContactDensity(allContactsBetweenClusters, allAreaBetweenClusters);

        int internalOffset1 = offsetIndex1;
        int internalCompressedOffset1 = compressedOffsetIndex1;
        for (SubcompartmentInterval interv1 : intervals1) {
            Integer id1 = interv1.getClusterID();
            int numRows = interv1.getWidthForResolution(resolution);
            if (numRows >= minIntervalSizeAllowed) {
                int internalOffset2 = offsetIndex2;
                int internalCompressedOffset2 = compressedOffsetIndex2;
                for (SubcompartmentInterval interv2 : intervals2) {
                    Integer id2 = interv2.getClusterID();
                    int numCols = interv2.getWidthForResolution(resolution);
                    if (numCols >= minIntervalSizeAllowed) {
                        String regionKey = id1 + "-" + id2;
                        float density = densityBetweenClusters.get(regionKey);
                        //updateMasterMatrixWithRegionalDensities(matrix, density, interv1, internalOffset1, numRows, interv2, internalOffset2, numCols, isIntra);
                        updateMasterMatrixWithRegionalDensitiesCompressed(matrix, numCountsForCol, density, interv1, internalOffset1, internalCompressedOffset1, numRows, interv2, internalOffset2, internalCompressedOffset2, numCols, isIntra);
                        internalOffset2 += numCols;
                        internalCompressedOffset2 += 1;
                    }
                }
                internalOffset1 += numRows;
                internalCompressedOffset1 += 1;
            }
        }
    }

    private void updateCountsAndArea(float[][] allDataForRegion, SubcompartmentInterval interv1, SubcompartmentInterval interv2,
                                     String regionKey, Map<String, Float> allContactsBetweenClusters, Map<String, Float> allAreaBetweenClusters) {
        float countsBetweenClusters = getSumTotalCounts(allDataForRegion, interv1, interv2);
        float areaBetweenClusters = interv1.getWidthForResolution(resolution) * interv2.getWidthForResolution(resolution);

        if (allAreaBetweenClusters.containsKey(regionKey)) {
            allAreaBetweenClusters.put(regionKey, allAreaBetweenClusters.get(regionKey) + areaBetweenClusters);
            allContactsBetweenClusters.put(regionKey, allContactsBetweenClusters.get(regionKey) + countsBetweenClusters);
        } else {
            allAreaBetweenClusters.put(regionKey, areaBetweenClusters);
            allContactsBetweenClusters.put(regionKey, countsBetweenClusters);
        }
    }

    private void updateNonZeroCountsAndArea(float[][] allDataForRegion, SubcompartmentInterval interv1, SubcompartmentInterval interv2,
                                            String regionKey, Map<String, Float> allContactsBetweenClusters, Map<String, Float> allAreaBetweenClusters) {
        //float countsBetweenClusters = getSumTotalCounts(allDataForRegion, interv1, interv2);
        //int areaBetweenClusters = interv1.getWidthForResolution(resolution) * interv2.getWidthForResolution(resolution);
        Pair<Float, Float> areaAndCounts = getTotalNonZeroCounts(allDataForRegion, interv1, interv2);

        if (allAreaBetweenClusters.containsKey(regionKey)) {
            allAreaBetweenClusters.put(regionKey, allAreaBetweenClusters.get(regionKey) + areaAndCounts.getFirst());
            allContactsBetweenClusters.put(regionKey, allContactsBetweenClusters.get(regionKey) + areaAndCounts.getSecond());
        } else {
            allAreaBetweenClusters.put(regionKey, areaAndCounts.getFirst());
            allContactsBetweenClusters.put(regionKey, areaAndCounts.getSecond());
        }
    }

    private Map<String, Float> getContactDensity(Map<String, Float> contacts, Map<String, Float> area) {

        Map<String, Float> density = new HashMap<>();
        for (String key : area.keySet()) {
            density.put(key, contacts.get(key) / area.get(key));
        }

        return density;
    }

    private void updateMasterMatrixWithRegionalDensities(float[][] matrix, float density,
                                                         SubcompartmentInterval interv1, int offsetIndex1, int numRows,
                                                         SubcompartmentInterval interv2, int offsetIndex2, int numCols, boolean isIntra) {
        for (int i = 0; i < numRows; i++) {
            for (int j = 0; j < numCols; j++) {
                rowIndexToIntervalMap.put(offsetIndex1 + i, interv1);
                rowIndexToIntervalMap.put(offsetIndex2 + j, interv2);
                matrix[offsetIndex1 + i][offsetIndex2 + j] = density;
            }
        }

        if (!isIntra) {
            for (int i = 0; i < numRows; i++) {
                for (int j = 0; j < numCols; j++) {
                    matrix[offsetIndex2 + j][offsetIndex1 + i] = density;
                }
            }
        }
    }

    private void updateMasterMatrixWithRegionalDensitiesCompressed(float[][] matrix, int[] numCountsForCol, float density,
                                                                   SubcompartmentInterval interv1, int offsetIndex1, int offsetCompressedIndex1, int numRows,
                                                                   SubcompartmentInterval interv2, int offsetIndex2, int offsetCompressedIndex2, int numCols, boolean isIntra) {
        numCountsForCol[offsetCompressedIndex2] = numCols;
        //double scalar = Math.sqrt(numCols);
        for (int i = 0; i < numRows; i++) {
            rowIndexToIntervalMap.put(offsetIndex1 + i, interv1);
            matrix[offsetIndex1 + i][offsetCompressedIndex2] = density;
        }

        if (!isIntra) {
            numCountsForCol[offsetCompressedIndex1] = numRows;
            //scalar = Math.sqrt(numRows);
            for (int j = 0; j < numCols; j++) {
                rowIndexToIntervalMap.put(offsetIndex2 + j, interv2);
                matrix[offsetIndex2 + j][offsetCompressedIndex1] = density;
            }
        }
    }

    private Pair<Float, Float> getTotalNonZeroCounts(float[][] allDataForRegion, SubcompartmentInterval interv1, SubcompartmentInterval interv2) {
        int numNonZero = 0;
        float total = 0;
        int binXStart = interv1.getX1() / resolution;
        int binXEnd = interv1.getX2() / resolution;

        int binYStart = interv2.getX1() / resolution;
        int binYEnd = interv2.getX2() / resolution;

        //List<Float> values = new ArrayList<>();

        for (int i = binXStart; i < binXEnd; i++) {
            for (int j = binYStart; j < binYEnd; j++) {
                float val = allDataForRegion[i][j];
                if (val > 0) {
                    total += val;
                    numNonZero++;
                    //values.add(val);
                }
                //if(val > maxVal){
                //    maxVal = val;
                //}
            }
        }

        // try new scaling
        numNonZero = Math.max(1, numNonZero);
        float baseArea = (binXEnd - binXStart) * (binYEnd - binYStart);
        float area = (float) Math.sqrt(numNonZero * baseArea);
        return new Pair<>(area, total);

        //return new Pair<>((float)numNonZero, total);
        /*
        float medianVal = getMedian(values);
        //return new Pair<>(1f, medianVal);

        return new Pair<>(baseArea, numNonZero*medianVal);

         */
    }

    private float getMedian(List<Float> values) {
        if (values.size() > 0) {
            Collections.sort(values);
            int midpoint = values.size() / 2;
            if (values.size() % 2 == 0) {
                // even number of elements
                return (values.get(midpoint - 1) + values.get(midpoint)) / 2;
            } else {
                // odd number of elements
                return values.get(midpoint);
            }
        }

        return 0;
    }

    private float getSumTotalCounts(float[][] allDataForRegion, SubcompartmentInterval interv1, SubcompartmentInterval interv2) {
        float total = 0;
        int binXStart = interv1.getX1() / resolution;
        int binXEnd = interv1.getX2() / resolution;

        int binYStart = interv2.getX1() / resolution;
        int binYEnd = interv2.getX2() / resolution;

        for (int i = binXStart; i < binXEnd; i++) {
            for (int j = binYStart; j < binYEnd; j++) {
                total += allDataForRegion[i][j];
            }
        }
        return total;
    }

    public synchronized Pair<Double, int[]> processGWKmeansResult(Cluster[] clusters, GenomeWideList<SubcompartmentInterval> subcompartments) {

        Set<SubcompartmentInterval> subcompartmentIntervals = new HashSet<>();
        if (MixerGlobals.printVerboseComments) {
            System.out.println("GW Composite data vs clustered into " + clusters.length + " clusters");
        }

        double withinClusterSumOfSquares = 0;
        int genomewideCompartmentID = 0;

        int[] ids = new int[clusters.length];
        for (int z = 0; z < clusters.length; z++) {
            Cluster cluster = clusters[z];
            int currentClusterID = ++genomewideCompartmentID;
            ids[z] = currentClusterID;

            if (MixerGlobals.printVerboseComments) {
                System.out.println("Size of cluster " + currentClusterID + " - " + cluster.getMemberIndexes().length);
            }

            for (int i : cluster.getMemberIndexes()) {

                withinClusterSumOfSquares += ClusterTools.getVectorMSEDifference(cluster.getCenter(), gwCleanMatrix[i]);

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
                    } else {
                        System.err.println("is weird error?");
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

        return new Pair<>(withinClusterSumOfSquares, ids);
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
