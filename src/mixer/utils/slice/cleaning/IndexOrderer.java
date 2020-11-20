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

package mixer.utils.slice.cleaning;

import javastraw.reader.Dataset;
import javastraw.reader.ExpectedValueFunction;
import javastraw.reader.HiCFileTools;
import javastraw.reader.MatrixZoomData;
import javastraw.reader.basics.Block;
import javastraw.reader.basics.Chromosome;
import javastraw.reader.basics.ContactRecord;
import javastraw.type.NormalizationType;
import mixer.MixerGlobals;
import mixer.utils.similaritymeasures.RobustCorrelationSimilarity;
import mixer.utils.similaritymeasures.SimilarityMetric;

import java.io.IOException;
import java.util.*;

public class IndexOrderer {

    private final Map<Chromosome, int[]> chromToReorderedIndices = new HashMap<>();
    private final int DISTANCE = 5000000;
    private final int resolution;
    private final int minDistanceThreshold;
    private final int numColsToJoin;
    private final int IGNORE = -1;
    private final int DEFAULT = -5;
    private final int CHECK_VAL = -2;
    private final float CORR_MIN = 0.2f;
    private final float INCREMENT = .1f;
    private final Random generator = new Random(0);
    private final Map<Integer, Integer> indexToRearrangedLength = new HashMap<>();

    public IndexOrderer(Dataset ds, Chromosome[] chromosomes, int resolution, NormalizationType normalizationType,
                        int numColumnsToPutTogether, BadIndexFinder badIndexLocations, long seed) {
        this.resolution = resolution;
        minDistanceThreshold = DISTANCE / resolution;
        numColsToJoin = numColumnsToPutTogether;
        generator.setSeed(seed);
        for (Chromosome chrom : chromosomes) {
            final MatrixZoomData zd = HiCFileTools.getMatrixZoomData(ds, chrom, chrom, resolution);
            ExpectedValueFunction df = ds.getExpectedValuesOrExit(zd.getZoom(), normalizationType, chrom, true);
            try {
                float[][] matrix = extractObsOverExpBoundedRegion(zd, chrom, normalizationType, df);
                int[] newOrderIndexes = getNewOrderOfIndices(chrom, matrix, badIndexLocations.getBadIndices(chrom));
                chromToReorderedIndices.put(chrom, newOrderIndexes);
            } catch (Exception e) {
                e.printStackTrace();
            }
            System.out.print(".");
        }
    }

    public Map<Integer, Integer> getIndexToRearrangedLength() {
        return indexToRearrangedLength;
    }

    private static double getExpected(int dist, ExpectedValueFunction df, int chrIndex) {
        return df.getExpectedValue(chrIndex, dist);
    }

    public int[] get(Chromosome chrom) {
        return chromToReorderedIndices.get(chrom);
    }

    private int[] getNewOrderOfIndices(Chromosome chromosome, float[][] matrix, Set<Integer> badIndices) {
        int[] newIndexOrderAssignments = new int[matrix.length];
        Arrays.fill(newIndexOrderAssignments, DEFAULT);
        eraseTheRowsColumnsWeDontWant(badIndices, matrix, newIndexOrderAssignments);

        int gCounter = doFirstRoundOfAssignmentsByCentroids(matrix, newIndexOrderAssignments);
        gCounter = doSecondRoundOfAssignments(matrix, newIndexOrderAssignments, gCounter);
        indexToRearrangedLength.put(chromosome.getIndex(), gCounter);
        return newIndexOrderAssignments;
    }

    private void eraseTheRowsColumnsWeDontWant(Set<Integer> badIndices, float[][] matrix, int[] newIndexOrderAssignments) {
        Set<Integer> columnsToErase = new HashSet<>();
        for (int k : badIndices) {
            newIndexOrderAssignments[k] = IGNORE;
            columnsToErase.add(k / numColsToJoin);
        }

        for (int i = 0; i < matrix.length; i++) {
            for (int j : columnsToErase) {
                matrix[i][j] = Float.NaN;
            }
        }

        for (int i : badIndices) {
            Arrays.fill(matrix[i], Float.NaN);
        }
    }

    private int doFirstRoundOfAssignmentsByCentroids(float[][] matrix, int[] newIndexOrderAssignments) {

        int numCentroids = 10;

        float[][] centroids = new QuickCentroids(quickCleanMatrix(matrix, newIndexOrderAssignments), numCentroids, generator.nextLong()).generateCentroids();
        System.out.println("Planned centroids: " + numCentroids + " Actual centroids: " + centroids.length);
        SimilarityMetric corrMetric = RobustCorrelationSimilarity.SINGLETON;

        int vectorLength = newIndexOrderAssignments.length;
        int[] numDecentRelations = new int[numCentroids];
        float[][] correlationCentroidsWithData = new float[numCentroids][vectorLength];
        for (int k = 0; k < numCentroids; k++) {
            for (int z = 0; z < vectorLength; z++) {
                if (newIndexOrderAssignments[z] < CHECK_VAL) {
                    float corr = corrMetric.distance(centroids[k], matrix[z]);
                    correlationCentroidsWithData[k][z] = corr;
                    if (corr > CORR_MIN || corr < -CORR_MIN) {
                        numDecentRelations[k]++;
                    }
                } else {
                    correlationCentroidsWithData[k][z] = Float.NaN;
                }
            }
        }

        int maxIndex = 0;
        for (int k = 1; k < numDecentRelations.length; k++) {
            if (numDecentRelations[maxIndex] < numDecentRelations[k]) {
                maxIndex = k;
            }
        }

        int gCounter = doSequentialOrdering(correlationCentroidsWithData[maxIndex], newIndexOrderAssignments, 0);
        for (int c = 0; c < numCentroids; c++) {
            if (c == maxIndex) continue;
            gCounter = doSequentialOrdering(correlationCentroidsWithData[c],
                    newIndexOrderAssignments, gCounter);
        }

        return gCounter;
    }

    private float[][] quickCleanMatrix(float[][] matrix, int[] newIndexOrderAssignments) {
        List<Integer> actualIndices = new ArrayList<>();
        for (int z = 0; z < newIndexOrderAssignments.length; z++) {
            if (newIndexOrderAssignments[z] < CHECK_VAL) {
                actualIndices.add(z);
            }
        }

        float[][] tempCleanMatrix = new float[actualIndices.size()][matrix[0].length];
        for (int i = 0; i < actualIndices.size(); i++) {
            System.arraycopy(matrix[actualIndices.get(i)], 0, tempCleanMatrix[i], 0, tempCleanMatrix[i].length);
        }
        return tempCleanMatrix;
    }

    private int doSequentialOrdering(float[] correlationWithCentroid, int[] newIndexOrderAssignments, int startCounter) {
        int counter = startCounter;
        for (float cutoff = 1 - INCREMENT; cutoff >= CORR_MIN; cutoff -= INCREMENT) {
            for (int z = 0; z < correlationWithCentroid.length; z++) {
                if (newIndexOrderAssignments[z] < CHECK_VAL && correlationWithCentroid[z] > cutoff) {
                    newIndexOrderAssignments[z] = counter++;
                }
            }
        }

        counter = getUpdatedNoMixIndex(counter);

        for (float cutoff = CORR_MIN; cutoff < 1; cutoff += INCREMENT) {
            for (int z = 0; z < correlationWithCentroid.length; z++) {
                float corr = correlationWithCentroid[z];
                float cutoff1 = -cutoff;
                float cutoff2 = cutoff1 - INCREMENT;
                if (newIndexOrderAssignments[z] < CHECK_VAL && corr < cutoff1 && corr >= cutoff2) {
                    newIndexOrderAssignments[z] = counter++;
                }
            }
        }

        return getUpdatedNoMixIndex(counter);
    }

    private int getUpdatedNoMixIndex(int counter) {
        int temp = (counter / numColsToJoin);
        if (counter % numColsToJoin > 0) {
            temp++;
        }
        return temp * numColsToJoin;
    }

    private int doSecondRoundOfAssignments(float[][] matrix, int[] newIndexOrderAssignments, int startCounter) {
        int vectorLength = newIndexOrderAssignments.length;
        int numRoundsThatHappen = 0;
        int counter = startCounter;
        SimilarityMetric corrMetric = RobustCorrelationSimilarity.SINGLETON;
        for (int cI = 0; cI < vectorLength; cI++) {
            // handle stuff
            if (newIndexOrderAssignments[cI] < CHECK_VAL) {
                numRoundsThatHappen++;
                newIndexOrderAssignments[cI] = counter++;

                for (int z = cI + 1; z < vectorLength; z++) {
                    if (newIndexOrderAssignments[z] < CHECK_VAL) {
                        float val = corrMetric.distance(matrix[cI], matrix[z]);
                        if (val >= CORR_MIN) {
                            newIndexOrderAssignments[z] = counter++;
                        }
                    }
                }
                counter = getUpdatedNoMixIndex(counter);
            }
            // else it has already been handled
            // or is a bad index, so skip
        }
        if (MixerGlobals.printVerboseComments) {
            System.out.println("Num rounds " + numRoundsThatHappen);
        }
        return counter;
    }

    public float[][] extractObsOverExpBoundedRegion(MatrixZoomData zd, Chromosome chromosome,
                                                    NormalizationType normalizationType,
                                                    ExpectedValueFunction df) throws IOException {
        if (df == null) {
            System.err.println("DF is null");
            return null;
        }
        // numRows/numCols is just to ensure a set size in case bounds are approximate
        // left upper corner is reference for 0,0
        int maxBin = (int) (chromosome.getLength() / resolution) + 1;
        int maxCompressedBin = maxBin / numColsToJoin;
        if (maxBin % numColsToJoin > 0) maxCompressedBin += 1;

        List<Block> blocks = HiCFileTools.getAllRegionBlocks(zd, 0, maxBin, 0, maxBin,
                normalizationType, true);
        float[][] data = new float[maxBin][maxCompressedBin];

        if (blocks.size() > 0) {
            for (Block b : blocks) {
                if (b != null) {
                    for (ContactRecord rec : b.getContactRecords()) {

                        int dist = Math.abs(rec.getBinX() - rec.getBinY());
                        if (dist < minDistanceThreshold) {
                            addOEValInPosition(Float.NaN, rec, data);
                        } else {
                            double observed = rec.getCounts() + 1;
                            double expected = getExpected(dist, df, chromosome.getIndex()) + 1;
                            float answer = (float) (observed / expected);
                            if (Float.isNaN(answer) || Float.isInfinite(answer)) {
                                answer = Float.NaN;
                            }
                            addOEValInPosition(answer, rec, data);
                        }
                    }
                }
            }
        }
        // force cleanup
        blocks = null;
        return data;
    }

    private void addOEValInPosition(float oeVal, ContactRecord rec, float[][] data) {
        int rX = rec.getBinX();
        int rY = rec.getBinY() / numColsToJoin;
        data[rX][rY] += oeVal;

        rX = rec.getBinY();
        rY = rec.getBinX() / numColsToJoin;
        data[rX][rY] += oeVal;
    }
}
