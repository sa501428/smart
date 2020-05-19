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
import java.util.*;

public class LinksMatrix extends CompositeGenomeWideDensityMatrix {
    public LinksMatrix(ChromosomeHandler chromosomeHandler, Dataset ds, NormalizationType norm, int resolution, GenomeWideList<SubcompartmentInterval> intraSubcompartments, int minIntervalSizeAllowed, File outputDirectory, Random generator, String[] referenceBedFiles) {
        super(chromosomeHandler, ds, norm, resolution, intraSubcompartments, minIntervalSizeAllowed, outputDirectory, generator, referenceBedFiles);
    }

    float[][] makeCleanScaledInterMatrix(Dataset ds) {

        Map<Integer, Integer> indexToFilteredLength = calculateActualLengthForChromosomes(chromosomes, intraSubcompartments);
        Pair<Integer, int[]> dimensions = calculateDimensionInterMatrix(chromosomes, indexToFilteredLength);
        float[][] interMatrix = new float[dimensions.getFirst()][dimensions.getFirst()];

        for (int i = 0; i < chromosomes.length; i++) {
            Chromosome chr1 = chromosomes[i];

            for (int j = i; j < chromosomes.length; j++) {
                Chromosome chr2 = chromosomes[j];

                final MatrixZoomData zd = HiCFileTools.getMatrixZoomData(ds, chr1, chr2, resolution);

                fillInChromosomeRegion(interMatrix, zd, chr1, dimensions.getSecond()[i], chr2, dimensions.getSecond()[j], i == j);
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

        return interMatrix;
    }


    private void fillInChromosomeRegion(float[][] matrix, MatrixZoomData zd, Chromosome chr1, int offsetIndex1,
                                        Chromosome chr2, int offsetIndex2, boolean isIntra) {

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

                FloatMatrixTools.cleanUpNansInfinitesNegatives(allDataForRegion);

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


}
