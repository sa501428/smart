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
import mixer.commandline.utils.common.DoubleMatrixTools;
import mixer.commandline.utils.common.FloatMatrixTools;
import mixer.commandline.utils.drink.kmeansfloat.Cluster;
import mixer.commandline.utils.drink.kmeansfloat.ClusterTools;
import mixer.data.ChromosomeHandler;
import mixer.data.Dataset;
import mixer.data.HiCFileTools;
import mixer.data.MatrixZoomData;
import mixer.data.feature.GenomeWideList;
import mixer.windowui.NormalizationType;
import org.apache.commons.math.linear.RealMatrix;
import org.broad.igv.feature.Chromosome;
import org.broad.igv.util.Pair;

import java.io.File;
import java.util.*;

public class CompositeGenomeWideDensityMatrix {
    private final NormalizationType norm;
    private final int resolution;
    private final GenomeWideList<SubcompartmentInterval> intraSubcompartments;
    private final float[][] gwCleanMatrix;
    private final Map<Integer, SubcompartmentInterval> indexToIntervalMap = new HashMap<>();
    private final Chromosome[] chromosomes;
    private final float threshold;
    private final int minIntervalSizeAllowed;
    private final int OFFSET = 1;

    public CompositeGenomeWideDensityMatrix(ChromosomeHandler chromosomeHandler, Dataset ds, NormalizationType norm, int resolution,
                                            GenomeWideList<SubcompartmentInterval> intraSubcompartments, float oeThreshold, int minIntervalSizeAllowed) {
        this.minIntervalSizeAllowed = minIntervalSizeAllowed;
        this.norm = norm;
        this.resolution = resolution;
        this.intraSubcompartments = intraSubcompartments;
        threshold = oeThreshold;
        chromosomes = chromosomeHandler.getAutosomalChromosomesArray();
        float[][] tempCleanData = makeCleanScaledInterMatrix(ds);
        //gwCleanMatrix = MatrixTools.getMainAppendedDerivativeDownColumnV2(tempCleanData, threshold / 2, threshold);

        /*
        if (useNormalizationOfRows) {
            if (derivativeStatus == Drink.USE_ONLY_DERIVATIVE) {
                gwCleanMatrix = FloatMatrixTools.getNormalizedThresholdedByMedian(FloatMatrixTools.getRelevantDerivativeScaledPositive(tempCleanData, threshold / 2, threshold), threshold);
            } else if (derivativeStatus == Drink.IGNORE_DERIVATIVE) {
                gwCleanMatrix = FloatMatrixTools.getNormalizedThresholdedByMedian(tempCleanData, threshold);
            } else {
                gwCleanMatrix = FloatMatrixTools.getNormalizedThresholdedByMedian(FloatMatrixTools.getMainAppendedDerivativeScaledPosDownColumn(tempCleanData, threshold / 2, threshold), threshold);
            }
        } else {
            if (derivativeStatus == Drink.USE_ONLY_DERIVATIVE) {
                gwCleanMatrix = FloatMatrixTools.getRelevantDerivative(tempCleanData, threshold / 2, threshold);
            } else if (derivativeStatus == Drink.IGNORE_DERIVATIVE) {
                gwCleanMatrix = tempCleanData;
            } else {
                gwCleanMatrix = FloatMatrixTools.getMainAppendedDerivativeDownColumn(tempCleanData, threshold / 2, threshold);
            }
        }
        */
        gwCleanMatrix = FloatMatrixTools.getMainAppendedDerivativeDownColumn(tempCleanData, threshold, threshold);
        //gwCleanMatrix = FloatMatrixTools.inPlaceZscoreDownRows(FloatMatrixTools.getRelevantDerivative(tempCleanData, 1, 3), 3);

        //gwCleanMatrix = FloatMatrixTools.getSimpleAppendedDerivativeDownColumn(tempCleanData, 6);

    }

    private float[][] makeCleanScaledInterMatrix(Dataset ds) {

        // height/weight chromosomes
        Map<Integer, Integer> indexToFilteredLength = calculateActualLengthForChromosomes(chromosomes, intraSubcompartments);

        Pair<Integer, int[]> dimensions = calculateDimensionInterMatrix(chromosomes, indexToFilteredLength);

        float[][] interMatrix = new float[dimensions.getFirst()][dimensions.getFirst()];
        for (int i = 0; i < chromosomes.length; i++) {
            Chromosome chr1 = chromosomes[i];

            for (int j = i; j < chromosomes.length; j++) {
                Chromosome chr2 = chromosomes[j];

                final MatrixZoomData zd = HiCFileTools.getMatrixZoomData(ds, chr1, chr2, resolution);
                if (zd == null) continue;

                fillInChromosomeRegion(interMatrix, ds, zd, chr1, dimensions.getSecond()[i], chr2, dimensions.getSecond()[j], i == j);
            }
        }

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

    private void fillInChromosomeRegion(float[][] matrix, Dataset ds, MatrixZoomData zd, Chromosome chr1, int offsetIndex1,
                                        Chromosome chr2, int offsetIndex2, boolean isIntra) {

        int lengthChr1 = chr1.getLength() / resolution + 1;
        int lengthChr2 = chr2.getLength() / resolution + 1;
        List<SubcompartmentInterval> intervals1 = intraSubcompartments.getFeatures("" + chr1.getIndex());
        List<SubcompartmentInterval> intervals2 = intraSubcompartments.getFeatures("" + chr2.getIndex());

        if (intervals1.size() == 0 || intervals2.size() == 0) return;
        float[][] allDataForRegion = null;
        try {
            if (isIntra) {
                RealMatrix localizedRegionData = HiCFileTools.getRealOEMatrixForChromosome(ds, zd, chr1, resolution, norm, threshold, ExtractingOEDataUtils.ThresholdType.LOG_OE_PLUS_AVG_BOUNDED_MADE_POS, true);
                allDataForRegion = DoubleMatrixTools.convertToFloatMatrix(localizedRegionData.getData());
            } else {
                RealMatrix allDataForRegionMatrix = HiCFileTools.extractLocalBoundedRegion(zd, 0, lengthChr1, 0, lengthChr2, lengthChr1, lengthChr2, norm, isIntra);
                allDataForRegion = DoubleMatrixTools.convertToFloatMatrix(allDataForRegionMatrix.getData());
            }
        } catch (Exception e) {
            e.printStackTrace();
            System.exit(99);
        }

        for (int i = 0; i < allDataForRegion.length; i++) {
            for (int j = 0; j < allDataForRegion[0].length; j++) {
                if (Float.isNaN(allDataForRegion[i][j]) || Float.isInfinite(allDataForRegion[i][j]) || Math.abs(allDataForRegion[i][j]) < 1E-30) {
                    allDataForRegion[i][j] = 0;
                }
            }
        }

        if (allDataForRegion == null) {
            System.err.println("Missing Interchromosomal Data " + zd.getKey());
            return;
        }

        Map<String, Integer> allAreaBetweenClusters = new HashMap<>();
        Map<String, Float> allContactsBetweenClusters = new HashMap<>();

        for (SubcompartmentInterval interv1 : intervals1) {
            if (interv1.getWidthForResolution(resolution) >= minIntervalSizeAllowed) {
                Integer id1 = interv1.getClusterID();
                for (SubcompartmentInterval interv2 : intervals2) {
                    if (interv2.getWidthForResolution(resolution) >= minIntervalSizeAllowed) {
                        Integer id2 = interv2.getClusterID();
                        String regionKey = id1 + "-" + id2;

                        float countsBetweenClusters = getSumTotalCounts(allDataForRegion, interv1, interv2, OFFSET);
                        int areaBetweenClusters = (interv1.getWidthForResolution(resolution) - 2 * OFFSET) * (interv2.getWidthForResolution(resolution) - 2 * OFFSET);

                        if (allAreaBetweenClusters.containsKey(regionKey)) {
                            allAreaBetweenClusters.put(regionKey, allAreaBetweenClusters.get(regionKey) + areaBetweenClusters);
                            allContactsBetweenClusters.put(regionKey, allContactsBetweenClusters.get(regionKey) + countsBetweenClusters);
                        } else {
                            allAreaBetweenClusters.put(regionKey, areaBetweenClusters);
                            allContactsBetweenClusters.put(regionKey, countsBetweenClusters);
                        }
                    }
                }
            }
        }

        Map<String, Float> zscoreDensityBetweenClusters = getThresholdedZscoreOfDensity(allContactsBetweenClusters, allAreaBetweenClusters);

        int internalOffset1 = offsetIndex1;
        for (SubcompartmentInterval interv1 : intervals1) {
            Integer id1 = interv1.getClusterID();
            int numRows = interv1.getWidthForResolution(resolution);
            if (numRows >= minIntervalSizeAllowed) {
                int internalOffset2 = offsetIndex2;
                for (SubcompartmentInterval interv2 : intervals2) {
                    Integer id2 = interv2.getClusterID();
                    int numCols = interv2.getWidthForResolution(resolution);
                    if (numCols >= minIntervalSizeAllowed) {
                        String regionKey = id1 + "-" + id2;
                        float density = zscoreDensityBetweenClusters.get(regionKey);
                        updateMasterMatrixWithRegionalDensities(matrix, density, interv1, internalOffset1, numRows, interv2, internalOffset2, numCols, isIntra);
                        internalOffset2 += numCols;
                    }
                }
                internalOffset1 += numRows;
            }
        }
    }

    private Map<String, Float> getThresholdedZscoreOfDensity(Map<String, Float> contacts, Map<String, Integer> area) {

        double totalCounts = 0;
        long totalArea = 0;
        Map<String, Float> density = new HashMap<>();
        for (String key : area.keySet()) {
            density.put(key, contacts.get(key) / area.get(key));
            totalCounts += contacts.get(key);
            totalArea += area.get(key);
        }
        float avgDensity = (float) (totalCounts / totalArea);

        double stdDev = 0;
        for (String key : area.keySet()) {
            float x = density.get(key) - avgDensity;
            stdDev += (x * x) * area.get(key);
        }
        stdDev = (float) Math.sqrt(stdDev / totalArea);

        Map<String, Float> zScores = new HashMap<>();
        for (String key : area.keySet()) {
            double z = (density.get(key) - avgDensity) / stdDev;
            float zThresh = (float) Math.min(Math.max(z, -2f), 2f);
            zScores.put(key, zThresh);
        }

        return zScores;
    }

    private void updateMasterMatrixWithRegionalDensities(float[][] matrix, float density,
                                                         SubcompartmentInterval interv1, int offsetIndex1, int numRows,
                                                         SubcompartmentInterval interv2, int offsetIndex2, int numCols, boolean isIntra) {
        for (int i = 0; i < numRows; i++) {
            for (int j = 0; j < numCols; j++) {
                indexToIntervalMap.put(offsetIndex1 + i, interv1);
                indexToIntervalMap.put(offsetIndex2 + j, interv2);
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

    private float getSumTotalCounts(float[][] allDataForRegion, SubcompartmentInterval interv1, SubcompartmentInterval interv2, int offset) {
        float total = 0;
        int binXStart = interv1.getX1() / resolution;
        int binXEnd = interv1.getX2() / resolution;

        int binYStart = interv2.getX1() / resolution;
        int binYEnd = interv2.getX2() / resolution;

        for (int i = binXStart + offset; i < binXEnd - offset; i++) {
            for (int j = binYStart + offset; j < binYEnd - offset; j++) {
                try {
                    if (!Float.isNaN(allDataForRegion[i][j])) {
                        total += allDataForRegion[i][j];
                    }
                } catch (Exception e) {
                    System.err.println(binXStart + " - " + binXEnd);

                    System.err.println(binYStart + " - " + binYEnd);
                    System.err.println(i + " - " + j);
                    System.err.println(interv1.getChrIndex() + " - " + interv2.getChrIndex());
                    e.printStackTrace();
                    System.exit(-1);
                }
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

                    if (indexToIntervalMap.containsKey(i)) {
                        interv = indexToIntervalMap.get(i);
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

    public void exportData(File outputDirectory) {
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
