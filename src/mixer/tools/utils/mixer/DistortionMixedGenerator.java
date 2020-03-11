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

package mixer.tools.utils.mixer;

import mixer.data.Dataset;
import mixer.data.HiCFileTools;
import mixer.data.MatrixZoomData;
import mixer.tools.utils.common.MatrixTools;
import mixer.tools.utils.common.UNIXTools;
import mixer.tools.utils.mixer.grind.GrindUtils;
import mixer.windowui.NormalizationType;
import org.apache.commons.math.linear.RealMatrix;
import org.broad.igv.feature.Chromosome;
import org.broad.igv.util.Pair;

import java.io.File;
import java.util.Random;

public class DistortionMixedGenerator {

    private final int specificResolution;
    private final String prefixString = System.currentTimeMillis() + "_";
    private final Chromosome chromI, chromJ;
    private final Dataset ds, ds2;
    private final NormalizationType norm;
    private final File outFolder;
    private final Random generator = new Random(0);
    private Integer imgSliceWidth, imgHalfSliceWidth;
    private Integer numManipulations, numExamplesPerRegion;
    private String negPath, posPath;
    private int counter = 0, stride;

    // grind -k KR -r 5000,10000,25000,100000 --stride 3 -c 1,2,3 --dense-labels --distort <hic file> null <128,4,1000> <directory>
    public DistortionMixedGenerator(Dataset ds, Dataset ds2, Pair<Chromosome, Chromosome> chromosomePair, int x, int y, int z, int resolution, NormalizationType norm, int stride, File outputDirectory) {
        this.imgSliceWidth = x;
        imgHalfSliceWidth = x / 2;
        this.numManipulations = y;
        numExamplesPerRegion = z;
        this.specificResolution = resolution;
        chromI = chromosomePair.getFirst();
        chromJ = chromosomePair.getSecond();
        this.ds = ds;
        this.ds2 = ds2;
        this.norm = norm;
        outFolder = UNIXTools.makeDir(outputDirectory, chromI.getName() + "_" + chromJ.getName() + "_" + resolution);
        this.stride = stride;
    }

    public void makeExamples() {


        try {

            final MatrixZoomData zdIJA = HiCFileTools.getMatrixZoomData(ds, chromI, chromJ, specificResolution);
            final MatrixZoomData zdIJB = HiCFileTools.getMatrixZoomData(ds2, chromI, chromJ, specificResolution);
            if (zdIJA == null || zdIJB == null) return;

            updateLatestMainPaths();

            boolean isIntraChromosomal = chromI.getIndex() == chromJ.getIndex();
            System.out.println("Currently processing: " + chromI.getName() + " - " + chromJ.getName() +
                    " at specificResolution " + specificResolution);

            MatrixZoomData matrixZoomDataIA, matrixZoomDataJA, matrixZoomDataIB, matrixZoomDataJB;
            if (isIntraChromosomal) {
                iterateAcrossIntraChromosomalRegion(zdIJA, zdIJB, chromI, specificResolution);
            } else {
                matrixZoomDataIA = HiCFileTools.getMatrixZoomData(ds, chromI, chromI, specificResolution);
                matrixZoomDataJA = HiCFileTools.getMatrixZoomData(ds, chromJ, chromJ, specificResolution);
                matrixZoomDataIB = HiCFileTools.getMatrixZoomData(ds2, chromI, chromI, specificResolution);
                matrixZoomDataJB = HiCFileTools.getMatrixZoomData(ds2, chromJ, chromJ, specificResolution);
                if (matrixZoomDataIA == null || matrixZoomDataJA == null) return;
                if (matrixZoomDataIB == null || matrixZoomDataJB == null) return;

                iterateBetweenInterChromosomalRegions(matrixZoomDataIA, matrixZoomDataJA, zdIJA,
                        matrixZoomDataIB, matrixZoomDataJB, zdIJB, chromI, chromJ, specificResolution);
            }

        } catch (Exception ex) {
            ex.printStackTrace();
        }
    }

    private float[][] generateCompositeMatrixWithNansCleanedFromZDS(MatrixZoomData zd1, MatrixZoomData zd2, MatrixZoomData zd12,
                                                                    int box1RectUL, int box1RectLR, int box2RectUL, int box2RectLR,
                                                                    int imgHalfSliceWidth, NormalizationType norm) throws Exception {
        RealMatrix localizedRegionDataBox1 = HiCFileTools.extractLocalBoundedRegion(zd1,
                box1RectUL, box1RectLR, box1RectUL, box1RectLR, imgHalfSliceWidth, imgHalfSliceWidth, norm, true);
        if (GrindUtils.mapRegionIsProblematic(localizedRegionDataBox1, .3)) return null;
        RealMatrix localizedRegionDataBox2 = HiCFileTools.extractLocalBoundedRegion(zd2,
                box2RectUL, box2RectLR, box2RectUL, box2RectLR, imgHalfSliceWidth, imgHalfSliceWidth, norm, true);
        if (GrindUtils.mapRegionIsProblematic(localizedRegionDataBox2, .3)) return null;
        RealMatrix localizedRegionDataBox12 = HiCFileTools.extractLocalBoundedRegion(zd12,
                box1RectUL, box1RectLR, box2RectUL, box2RectLR, imgHalfSliceWidth, imgHalfSliceWidth, norm, false);

        return MatrixTools.generateCompositeMatrixWithNansCleaned(localizedRegionDataBox1, localizedRegionDataBox2, localizedRegionDataBox12);
    }

    private void iterateAcrossIntraChromosomalRegion(MatrixZoomData zdA, MatrixZoomData zdB, Chromosome chrom, int resolution) {

        // sliding along the diagonal
        int numberOfExamplesCounter = 0;
        int maxChrLength = (chrom.getLength() / resolution);

        double[] vectorA = ds.getNormalizationVector(zdA.getChr1Idx(), zdA.getZoom(), norm).getData();
        double[] vectorB = ds2.getNormalizationVector(zdB.getChr1Idx(), zdB.getZoom(), norm).getData();
        if (vectorA == null || vectorB == null) return;

        for (int posIndex1 = 0; posIndex1 < maxChrLength - imgSliceWidth; posIndex1 += stride) {
            if (regionOverlapsBadAreas(posIndex1, vectorA, imgHalfSliceWidth)) continue;
            if (regionOverlapsBadAreas(posIndex1, vectorB, imgHalfSliceWidth)) continue;
            int numTimesRegionUsed = 0;

            for (int posIndex2 = posIndex1 + imgHalfSliceWidth; posIndex2 < maxChrLength - imgHalfSliceWidth; posIndex2 += stride) {
                if (numTimesRegionUsed > numExamplesPerRegion) break;
                if (regionOverlapsBadAreas(posIndex2, vectorA, imgHalfSliceWidth)) continue;
                if (regionOverlapsBadAreas(posIndex2, vectorB, imgHalfSliceWidth)) continue;

                if (numberOfExamplesCounter > 1000) {
                    updateLatestMainPaths();
                    numberOfExamplesCounter = 0;
                }
                if (getTrainingDataAndSaveToFile(zdA, zdA, zdA, zdB, zdB, zdB, posIndex1, posIndex2, chrom.getName(), chrom.getName(),
                        posIndex2 == posIndex1 + imgHalfSliceWidth, imgHalfSliceWidth, norm, numManipulations,
                        posPath, negPath)) {
                    numberOfExamplesCounter++;
                    numTimesRegionUsed++;
                }
            }
        }
    }

    private boolean regionOverlapsBadAreas(int posIndex1, double[] vector, int distance) {
        for (int k = posIndex1; k < posIndex1 + distance; k++) {
            if (Double.isNaN(vector[k]) || Double.isInfinite(vector[k]) || Math.abs(vector[k]) < 1e-10) {
                return true;
            }
        }
        return false;
    }

    private void iterateBetweenInterChromosomalRegions(MatrixZoomData zdA1, MatrixZoomData zdA2, MatrixZoomData zdA12,
                                                       MatrixZoomData zdB1, MatrixZoomData zdB2, MatrixZoomData zdB12,
                                                       Chromosome chrom1, Chromosome chrom2, int resolution) {

        // iterating across both chromosomes
        int maxChrLength1 = (chrom1.getLength() / resolution);
        int maxChrLength2 = (chrom2.getLength() / resolution);
        int numberOfExamplesCounter = 0;

        double[] vectorA1 = ds.getNormalizationVector(zdA1.getChr1Idx(), zdA1.getZoom(), norm).getData();
        double[] vectorA2 = ds.getNormalizationVector(zdA2.getChr1Idx(), zdA2.getZoom(), norm).getData();
        double[] vectorB1 = ds2.getNormalizationVector(zdB1.getChr1Idx(), zdB1.getZoom(), norm).getData();
        double[] vectorB2 = ds2.getNormalizationVector(zdB2.getChr1Idx(), zdB2.getZoom(), norm).getData();
        if (vectorA1 == null || vectorA2 == null) return;
        if (vectorB1 == null || vectorB2 == null) return;

        for (int posIndex1 = 0; posIndex1 < maxChrLength1 - imgHalfSliceWidth; posIndex1 += stride) {
            if (regionOverlapsBadAreas(posIndex1, vectorA1, imgHalfSliceWidth)) continue;
            if (regionOverlapsBadAreas(posIndex1, vectorB1, imgHalfSliceWidth)) continue;
            int numTimesRegionUsed = 0;

            for (int posIndex2 = 0; posIndex2 < maxChrLength2 - imgHalfSliceWidth; posIndex2 += stride) {
                if (numTimesRegionUsed > numExamplesPerRegion) break;
                if (regionOverlapsBadAreas(posIndex2, vectorA2, imgHalfSliceWidth)) continue;
                if (regionOverlapsBadAreas(posIndex2, vectorB2, imgHalfSliceWidth)) continue;
                if (numberOfExamplesCounter > 1000) {
                    updateLatestMainPaths();
                    numberOfExamplesCounter = 0;
                }
                if (getTrainingDataAndSaveToFile(zdA1, zdA2, zdA12, zdB1, zdB2, zdB12, posIndex1, posIndex2, chrom1.getName(), chrom2.getName(),
                        false, imgHalfSliceWidth, norm, numManipulations,
                        posPath, negPath)) {
                    numberOfExamplesCounter++;
                    numTimesRegionUsed++;
                }
            }
        }
    }

    private void updateLatestMainPaths() {
        File newDir = UNIXTools.makeDir(outFolder, "negative_" + specificResolution + "_" + counter);
        negPath = newDir.getAbsolutePath();
        newDir = UNIXTools.makeDir(outFolder, "positive_" + specificResolution + "_" + counter);
        posPath = newDir.getAbsolutePath();
        counter++;
    }

    private float getRandScale() {
        return Math.max(generator.nextFloat(), .05f);
    }

    private boolean getTrainingDataAndSaveToFile(MatrixZoomData zdA1, MatrixZoomData zdA2, MatrixZoomData zdA12,
                                                 MatrixZoomData zdB1, MatrixZoomData zdB2, MatrixZoomData zdB12,
                                                 int box1XIndex, int box2XIndex, String chrom1Name, String chrom2Name,
                                                 boolean isContinuousRegion, int imgHalfSliceWidth, NormalizationType norm,
                                                 int numManipulations,
                                                 String posPath, String negPath) {

        int box1RectUL = box1XIndex;
        int box1RectLR = box1XIndex + imgHalfSliceWidth;

        int box2RectUL = box2XIndex;
        int box2RectLR = box2XIndex + imgHalfSliceWidth;

        try {
            float[][] compositeMatrixA = generateCompositeMatrixWithNansCleanedFromZDS(zdA1, zdA2, zdA12,
                    box1RectUL, box1RectLR, box2RectUL, box2RectLR, imgHalfSliceWidth, norm);
            float[][] compositeMatrixB = generateCompositeMatrixWithNansCleanedFromZDS(zdB1, zdB2, zdB12,
                    box1RectUL, box1RectLR, box2RectUL, box2RectLR, imgHalfSliceWidth, norm);
            if (compositeMatrixA == null || compositeMatrixB == null) return false;

            //if (!GrindUtils.isJustEmptyEnough(compositeMatrix)) return false;

            float[][] labelsMatrixA = GrindUtils.generateDefaultDistortionLabelsFile(compositeMatrixA.length, 4, isContinuousRegion);
            float[][] labelsMatrixB = GrindUtils.generateDefaultDistortionLabelsFile(compositeMatrixB.length, 4, isContinuousRegion);
            //GrindUtils.cleanUpLabelsMatrixBasedOnData(labelsMatrix, compositeMatrix);


            float[][] compositeMatrixAB = MatrixTools.add(compositeMatrixA, compositeMatrixB, 1f, getRandScale());
            float[][] labelsMatrixAB = MatrixTools.max(labelsMatrixA, labelsMatrixB);

            String filePrefix = prefixString + "orig_" + chrom1Name + "_" + box1XIndex + "_" + chrom2Name + "_" + box2XIndex + "_matrix";
            GrindUtils.saveGrindMatrixDataToFile(filePrefix, negPath, compositeMatrixAB, false);
            GrindUtils.saveGrindMatrixDataToFile(filePrefix + "_labels", negPath, labelsMatrixAB, false);

            for (int k = 0; k < numManipulations; k++) {
                Pair<float[][], float[][]> alteredMatrices = GrindUtils.randomlyManipulateMatrix(compositeMatrixB, labelsMatrixB, generator);
                compositeMatrixB = alteredMatrices.getFirst();
                labelsMatrixB = alteredMatrices.getSecond();

                if (k == 0 || k == (numManipulations - 1) || generator.nextBoolean()) {
                    compositeMatrixAB = MatrixTools.add(compositeMatrixA, compositeMatrixB, 1f, getRandScale());
                    labelsMatrixAB = MatrixTools.max(labelsMatrixA, labelsMatrixB);

                    filePrefix = prefixString + "dstrt_" + chrom1Name + "_" + box1XIndex + "_" + chrom2Name + "_" + box2XIndex + "_" + k + "_matrix";
                    GrindUtils.saveGrindMatrixDataToFile(filePrefix, posPath, compositeMatrixAB, false);
                    GrindUtils.saveGrindMatrixDataToFile(filePrefix + "_labels", posPath, labelsMatrixAB, false);
                }
            }
        } catch (Exception e) {
        }
        return true;
    }

}
