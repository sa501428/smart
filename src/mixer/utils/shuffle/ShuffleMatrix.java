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

package mixer.utils.shuffle;

import javastraw.featurelist.GenomeWideList;
import javastraw.reader.ChromosomeHandler;
import javastraw.reader.Dataset;
import javastraw.reader.HiCFileTools;
import javastraw.reader.MatrixZoomData;
import javastraw.reader.basics.Chromosome;
import javastraw.type.NormalizationType;
import mixer.utils.common.FloatMatrixTools;
import mixer.utils.common.Pair;
import mixer.utils.slice.structures.SubcompartmentInterval;

import java.util.HashMap;
import java.util.List;
import java.util.Map;

public class ShuffleMatrix {

    private final NormalizationType norm;
    private final int resolution;
    private final GenomeWideList<SubcompartmentInterval> intraSubcompartments;
    private final float[][] gwCleanMatrix;
    private final Map<Integer, SubcompartmentInterval> indexToInterval1Map = new HashMap<>();
    private final Map<Integer, SubcompartmentInterval> indexToInterval2Map = new HashMap<>();
    private final Chromosome[] rowsChromosomes;
    private final Chromosome[] colsChromosomes;
    private final int compressionFactor;

    public ShuffleMatrix(ChromosomeHandler chromosomeHandler, Dataset ds, NormalizationType norm, int resolution,
                         GenomeWideList<SubcompartmentInterval> intraSubcompartments,
                         InterMapType mapType, int compressionFactor) {
        this.norm = norm;
        this.resolution = resolution;
        this.intraSubcompartments = intraSubcompartments;
        this.compressionFactor = compressionFactor;

        switch (mapType) {
            case SKIP_BY_TWOS: // but start with CHR 1 separate
                rowsChromosomes = chromosomeHandler.splitAutosomesAndSkipByTwos().getFirst();
                colsChromosomes = chromosomeHandler.splitAutosomesAndSkipByTwos().getSecond();
                break;
            case FIRST_HALF_VS_SECOND_HALF:
                rowsChromosomes = chromosomeHandler.splitAutosomesIntoHalves().getFirst();
                colsChromosomes = chromosomeHandler.splitAutosomesIntoHalves().getSecond();
                break;
            case ODDS_VS_EVENS:
            default:
                rowsChromosomes = chromosomeHandler.extractOddOrEvenAutosomes(true);
                colsChromosomes = chromosomeHandler.extractOddOrEvenAutosomes(false);
                break;
        }

        gwCleanMatrix = makeCleanScaledInterMatrix(ds);
    }

    public enum InterMapType {ODDS_VS_EVENS, FIRST_HALF_VS_SECOND_HALF, SKIP_BY_TWOS}

    private float[][] makeCleanScaledInterMatrix(Dataset ds) {

        // assuming Odd vs Even
        // height chromosomes
        Pair<Integer, int[]> rowsDimension = calculateDimensionInterMatrix(rowsChromosomes);

        // width chromosomes
        Pair<Integer, int[]> colsDimension = calculateDimensionInterMatrix(colsChromosomes);

        float[][] interMatrix = new float[rowsDimension.getFirst()][colsDimension.getFirst()];
        for (int i = 0; i < rowsChromosomes.length; i++) {
            Chromosome chr1 = rowsChromosomes[i];

            for (int j = 0; j < colsChromosomes.length; j++) {
                Chromosome chr2 = colsChromosomes[j];

                if (chr1.getIndex() == chr2.getIndex()) continue;
                final MatrixZoomData zd = HiCFileTools.getMatrixZoomData(ds, chr1, chr2, resolution);
                if (zd == null) continue;

                // will need to flip across diagonal
                boolean needToFlip = chr2.getIndex() < chr1.getIndex();
                fillInInterChromosomeRegion(interMatrix, zd, chr1, rowsDimension.getSecond()[i], chr2, colsDimension.getSecond()[j], needToFlip);
            }
        }

        return interMatrix;
    }

    private Pair<Integer, int[]> calculateDimensionInterMatrix(Chromosome[] chromosomes) {
        int total = 0;
        int[] indices = new int[chromosomes.length];

        for (int i = 0; i < chromosomes.length; i++) {
            for (SubcompartmentInterval interval : intraSubcompartments.getFeatures("" + chromosomes[i].getIndex())) {
                total += interval.getWidthForResolution(resolution);
            }
            if (i < chromosomes.length - 1) {
                indices[i + 1] = total;
            }
        }

        return new Pair<>(total, indices);
    }

    private void fillInInterChromosomeRegion(float[][] matrix, MatrixZoomData zd, Chromosome chr1, int offsetIndex1,
                                             Chromosome chr2, int offsetIndex2, boolean needToFlip) {

        int chr1Index = chr1.getIndex();
        int chr2Index = chr2.getIndex();
        if (chr1Index == chr2Index) {
            System.err.println("Same chr " + chr1.getName());
            System.exit(989);
        }

        int lengthChr1 = (int) (chr1.getLength() / resolution);
        int lengthChr2 = (int) (chr2.getLength() / resolution);
        List<SubcompartmentInterval> intervals1 = intraSubcompartments.getFeatures("" + chr1.getIndex());
        List<SubcompartmentInterval> intervals2 = intraSubcompartments.getFeatures("" + chr2.getIndex());

        if (intervals1.size() == 0 || intervals2.size() == 0) return;
        float[][] allDataForRegion = null;
        try {
            if (needToFlip) {
                float[][] allDataForRegionMatrix = HiCFileTools.extractLocalBoundedRegionFloatMatrix(zd, 0, lengthChr2,
                        0, lengthChr1, lengthChr2, lengthChr1, norm, false);
                allDataForRegion = FloatMatrixTools.transpose(allDataForRegionMatrix);
            } else {
                allDataForRegion = HiCFileTools.extractLocalBoundedRegionFloatMatrix(zd, 0, lengthChr1,
                        0, lengthChr2, lengthChr1, lengthChr2, norm, false);
            }
        } catch (Exception e) {
            e.printStackTrace();
            System.exit(99);
        }

        if (allDataForRegion == null) {
            System.err.println("Missing Interchromosomal Data " + zd.getKey());
            return;
        }

        Map<String, Integer> allAreaBetweenClusters = new HashMap<>();
        Map<String, Float> allContactsBetweenClusters = new HashMap<>();

        for (SubcompartmentInterval interv1 : intervals1) {
            Integer id1 = interv1.getClusterID();
            for (SubcompartmentInterval interv2 : intervals2) {
                Integer id2 = interv2.getClusterID();
                String regionKey = id1 + "-" + id2;


            }
        }

        int internalOffset1 = offsetIndex1;
        for (SubcompartmentInterval interv1 : intervals1) {
            Integer id1 = interv1.getClusterID();
            int numRows = interv1.getWidthForResolution(resolution);

            int internalOffset2 = offsetIndex2;
            for (SubcompartmentInterval interv2 : intervals2) {
                Integer id2 = interv2.getClusterID();
                int numCols = interv2.getWidthForResolution(resolution);

            }
            internalOffset1 += numRows;
        }
    }

    private void updateMasterMatrixWithRegionalDensities(float[][] matrix, float density,
                                                         SubcompartmentInterval interv1, int offsetIndex1, int numRows,
                                                         SubcompartmentInterval interv2, int offsetIndex2, int numCols) {
        for (int i = 0; i < numRows; i++) {
            for (int j = 0; j < numCols; j++) {
                indexToInterval1Map.put(offsetIndex1 + i, interv1);
                indexToInterval2Map.put(offsetIndex2 + j, interv2);
                matrix[offsetIndex1 + i][offsetIndex2 + j] = density;
            }
        }
    }
}
