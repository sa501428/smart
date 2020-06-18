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

import mixer.commandline.utils.common.FloatMatrixTools;
import mixer.data.ChromosomeHandler;
import mixer.data.Dataset;
import mixer.data.HiCFileTools;
import mixer.data.MatrixZoomData;
import mixer.data.feature.GenomeWideList;
import mixer.windowui.NormalizationType;
import org.broad.igv.feature.Chromosome;
import org.broad.igv.util.Pair;

import java.io.File;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.Random;

public class DrinksMatrix extends CompositeGenomeWideDensityMatrix {

    private float threshold = 5f;

    public DrinksMatrix(ChromosomeHandler chromosomeHandler, Dataset ds, NormalizationType norm, int resolution, GenomeWideList<SubcompartmentInterval> intraSubcompartments, int minIntervalSizeAllowed, File outputDirectory, Random generator, String[] relativeTestFiles) {
        super(chromosomeHandler, ds, norm, resolution, intraSubcompartments, minIntervalSizeAllowed, outputDirectory, generator, relativeTestFiles);
    }

    float[][] makeCleanScaledInterMatrix(Dataset ds) {

        // height/weight chromosomes
        Map<Integer, Integer> indexToFilteredLength = calculateActualLengthForChromosomes(chromosomes, intraSubcompartments);
        Map<Integer, Integer> indexToCompressedLength = calculateCompressedLengthForChromosomes(chromosomes, intraSubcompartments);

        Pair<Integer, int[][]> dimensions = calculateDimensionInterMatrix(chromosomes, indexToFilteredLength);
        Pair<Integer, int[][]> compressedDimensions = calculateDimensionInterMatrix(chromosomes, indexToCompressedLength);

        System.out.println(".");

        float[][] interMatrix = new float[dimensions.getFirst()][compressedDimensions.getFirst()];
        int[] numCountsForCol = new int[compressedDimensions.getFirst()];
        for (int i = 0; i < chromosomes.length; i++) {
            Chromosome chr1 = chromosomes[i];

            for (int j = i; j < chromosomes.length; j++) {
                Chromosome chr2 = chromosomes[j];

                final MatrixZoomData zd = HiCFileTools.getMatrixZoomData(ds, chr1, chr2, resolution);

                fillInChromosomeRegion(interMatrix, numCountsForCol, ds, zd, chr1, dimensions.getSecond()[0][i], compressedDimensions.getSecond()[0][i], chr2, dimensions.getSecond()[0][j], compressedDimensions.getSecond()[0][j], i == j);
                System.out.print(".");
            }
        }
        System.out.println(".");

        //FloatMatrixTools.inPlaceZscoreDownCols(interMatrix);
        // todo check
        FloatMatrixTools.scaleValuesByCount(interMatrix, numCountsForCol);

        return interMatrix;
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

        FloatMatrixTools.cleanUpNansInfinitesNegatives(allDataForRegion);

        int internalOffset1 = offsetIndex1;
        int internalCompressedOffset1 = compressedOffsetIndex1;
        for (SubcompartmentInterval interv1 : intervals1) {
            int numRows = interv1.getWidthForResolution(resolution);
            if (numRows >= minIntervalSizeAllowed) {
                int internalOffset2 = offsetIndex2;
                int internalCompressedOffset2 = compressedOffsetIndex2;
                for (SubcompartmentInterval interv2 : intervals2) {
                    int numCols = interv2.getWidthForResolution(resolution);
                    if (numCols >= minIntervalSizeAllowed) {
                        //updateMasterMatrixWithRegionalDensities(matrix, density, interv1, internalOffset1, numRows, interv2, internalOffset2, numCols, isIntra);
                        copyValuesToArea(matrix, numCountsForCol, allDataForRegion, interv1, internalOffset1, internalCompressedOffset1, numRows, interv2, internalOffset2, internalCompressedOffset2, numCols, isIntra);
                        internalOffset2 += numCols;
                        internalCompressedOffset2 += 1;
                    }
                }
                internalOffset1 += numRows;
                internalCompressedOffset1 += 1;
            }
        }
    }

    private void copyValuesToArea(float[][] matrix, int[] numCountsForCol, float[][] dataForRegion,
                                  SubcompartmentInterval interv1, int offsetIndex1, int offsetCompressedIndex1, int numRows,
                                  SubcompartmentInterval interv2, int offsetIndex2, int offsetCompressedIndex2, int numCols, boolean isIntra) {

        float countsBetweenClusters = getSumTotalCounts(dataForRegion, interv1, interv2);
        float areaBetweenClusters = interv1.getWidthForResolution(resolution) * interv2.getWidthForResolution(resolution);
        float density = countsBetweenClusters / areaBetweenClusters;

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
}
