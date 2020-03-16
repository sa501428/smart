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

package mixer.commandline.utils.mixer;

import mixer.commandline.utils.common.MatrixTools;
import mixer.commandline.utils.common.UNIXTools;
import mixer.commandline.utils.grind.GrindUtils;
import mixer.data.Dataset;
import mixer.data.HiCFileTools;
import mixer.data.MatrixZoomData;
import mixer.windowui.NormalizationType;
import org.apache.commons.math.linear.RealMatrix;
import org.broad.igv.feature.Chromosome;
import org.broad.igv.util.Pair;

import java.io.File;
import java.util.Random;

public class DistortionGenerator {

    private final int specificResolution;
    private final String prefixString = System.currentTimeMillis() + "_";
    private final Chromosome chromI, chromJ;
    private final Dataset ds;
    private final NormalizationType norm;
    private final File outFolder;
    private final Random generator = new Random(0);
    private Integer imgSliceWidth, imgHalfSliceWidth;
    private Integer numManipulations, numExamplesPerRegion;
    private String negPath, posPath;
    private int counter = 0, stride;

    // grind -k KR -r 5000,10000,25000,100000 --stride 3 -c 1,2,3 --dense-labels --distort <hic file> null <128,4,1000> <directory>
    public DistortionGenerator(Dataset ds, Pair<Chromosome, Chromosome> chromosomePair, int x, int y, int z, int resolution, NormalizationType norm, int stride, File outputDirectory) {
        this.imgSliceWidth = x;
        imgHalfSliceWidth = x / 2;
        this.numManipulations = y;
        numExamplesPerRegion = z;
        this.specificResolution = resolution;
        chromI = chromosomePair.getFirst();
        chromJ = chromosomePair.getSecond();
        this.ds = ds;
        this.norm = norm;
        outFolder = UNIXTools.makeDir(outputDirectory, chromI.getName() + "_" + chromJ.getName() + "_" + resolution);
        this.stride = stride;
    }

    public void makeExamples() {


        try {

            final MatrixZoomData zd = HiCFileTools.getMatrixZoomData(ds, chromI, chromJ, specificResolution);
            if (zd == null) return;

            updateLatestMainPaths();

            boolean isIntraChromosomal = chromI.getIndex() == chromJ.getIndex();
            System.out.println("Currently processing: " + chromI.getName() + " - " + chromJ.getName() +
                    " at specificResolution " + specificResolution);

            MatrixZoomData matrixZoomDataI, matrixZoomDataJ;
            if (isIntraChromosomal) {
                iterateAcrossIntraChromosomalRegion(zd, chromI, specificResolution);
            } else {
                matrixZoomDataI = HiCFileTools.getMatrixZoomData(ds, chromI, chromI, specificResolution);
                matrixZoomDataJ = HiCFileTools.getMatrixZoomData(ds, chromJ, chromJ, specificResolution);
                if (matrixZoomDataI == null || matrixZoomDataJ == null) return;

                iterateBetweenInterChromosomalRegions(matrixZoomDataI, matrixZoomDataJ, zd, chromI, chromJ, specificResolution);
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

    private void iterateAcrossIntraChromosomalRegion(MatrixZoomData zd, Chromosome chrom, int resolution) {

        // sliding along the diagonal
        int numberOfExamplesCounter = 0;
        int maxChrLength = (chrom.getLength() / resolution);

        double[] vector = ds.getNormalizationVector(zd.getChr1Idx(), zd.getZoom(), norm).getData();
        if (vector == null) return;

        for (int posIndex1 = 0; posIndex1 < maxChrLength - imgSliceWidth; posIndex1 += stride) {
            if (regionOverlapsBadAreas(posIndex1, vector, imgHalfSliceWidth)) continue;
            int numTimesRegionUsed = 0;

            for (int posIndex2 = posIndex1 + imgHalfSliceWidth; posIndex2 < maxChrLength; posIndex2 += stride) {
                if (numTimesRegionUsed > numExamplesPerRegion) break;
                if (regionOverlapsBadAreas(posIndex2, vector, imgHalfSliceWidth)) continue;

                if (numberOfExamplesCounter > 1000) {
                    updateLatestMainPaths();
                    numberOfExamplesCounter = 0;
                }
                if (getTrainingDataAndSaveToFile(zd, zd, zd, posIndex1, posIndex2, chrom.getName(), chrom.getName(),
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

    private void iterateBetweenInterChromosomalRegions(MatrixZoomData zd1, MatrixZoomData zd2, MatrixZoomData zd12, Chromosome chrom1, Chromosome chrom2,
                                                       int resolution) {

        // iterating across both chromosomes
        int maxChrLength1 = (chrom1.getLength() / resolution);
        int maxChrLength2 = (chrom2.getLength() / resolution);
        int numberOfExamplesCounter = 0;

        double[] vector1 = ds.getNormalizationVector(zd1.getChr1Idx(), zd1.getZoom(), norm).getData();
        double[] vector2 = ds.getNormalizationVector(zd2.getChr1Idx(), zd2.getZoom(), norm).getData();
        if (vector1 == null || vector2 == null) return;

        for (int posIndex1 = 0; posIndex1 < maxChrLength1 - imgHalfSliceWidth; posIndex1 += stride) {
            if (regionOverlapsBadAreas(posIndex1, vector1, imgHalfSliceWidth)) continue;
            int numTimesRegionUsed = 0;

            for (int posIndex2 = 0; posIndex2 < maxChrLength2 - imgHalfSliceWidth; posIndex2 += stride) {
                if (numTimesRegionUsed > numExamplesPerRegion) break;
                if (regionOverlapsBadAreas(posIndex2, vector2, imgHalfSliceWidth)) continue;
                if (numberOfExamplesCounter > 1000) {
                    updateLatestMainPaths();
                    numberOfExamplesCounter = 0;
                }
                if (getTrainingDataAndSaveToFile(zd1, zd2, zd12, posIndex1, posIndex2, chrom1.getName(), chrom2.getName(),
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

    private boolean getTrainingDataAndSaveToFile(MatrixZoomData zd1, MatrixZoomData zd2, MatrixZoomData zd12,
                                                 int box1XIndex, int box2XIndex, String chrom1Name, String chrom2Name,
                                                 boolean isContinuousRegion, int imgHalfSliceWidth, NormalizationType norm,
                                                 int numManipulations,
                                                 String posPath, String negPath) {

        int box1RectUL = box1XIndex;
        int box1RectLR = box1XIndex + imgHalfSliceWidth;

        int box2RectUL = box2XIndex;
        int box2RectLR = box2XIndex + imgHalfSliceWidth;

        try {
            float[][] compositeMatrix = generateCompositeMatrixWithNansCleanedFromZDS(zd1, zd2, zd12,
                    box1RectUL, box1RectLR, box2RectUL, box2RectLR, imgHalfSliceWidth, norm);
            if (compositeMatrix == null) return false;


            //if (!GrindUtils.isJustEmptyEnough(compositeMatrix)) return false;

            float[][] labelsMatrix = GrindUtils.generateDefaultDistortionLabelsFile(compositeMatrix.length, 4, isContinuousRegion);
            //GrindUtils.cleanUpLabelsMatrixBasedOnData(labelsMatrix, compositeMatrix);

            String filePrefix = prefixString + "orig_" + chrom1Name + "_" + box1XIndex + "_" + chrom2Name + "_" + box2XIndex + "_matrix";
            GrindUtils.saveGrindMatrixDataToFile(filePrefix, negPath, compositeMatrix, false);
            GrindUtils.saveGrindMatrixDataToFile(filePrefix + "_labels", negPath, labelsMatrix, false);

            for (int k = 0; k < numManipulations; k++) {
                Pair<float[][], float[][]> alteredMatrices = GrindUtils.randomlyManipulateMatrix(compositeMatrix, labelsMatrix, generator);
                compositeMatrix = alteredMatrices.getFirst();
                labelsMatrix = alteredMatrices.getSecond();

                if (k == 0 || k == (numManipulations - 1) || generator.nextBoolean()) {
                    filePrefix = prefixString + "dstrt_" + chrom1Name + "_" + box1XIndex + "_" + chrom2Name + "_" + box2XIndex + "_" + k + "_matrix";
                    GrindUtils.saveGrindMatrixDataToFile(filePrefix, posPath, compositeMatrix, false);
                    GrindUtils.saveGrindMatrixDataToFile(filePrefix + "_labels", posPath, labelsMatrix, false);
                }
            }

        } catch (Exception e) {

        }
        return true;
    }

}
