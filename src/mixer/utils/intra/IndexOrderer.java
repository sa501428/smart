/*
 * The MIT License (MIT)
 *
 * Copyright (c) 2011-2022 Rice University, Baylor College of Medicine, Aiden Lab
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

package mixer.utils.intra;

import javastraw.expected.LogExpectedSpline;
import javastraw.reader.Dataset;
import javastraw.reader.basics.Chromosome;
import javastraw.reader.mzd.MatrixZoomData;
import javastraw.reader.type.NormalizationType;
import javastraw.tools.HiCFileTools;
import javastraw.tools.ParallelizationTools;
import mixer.SmartTools;
import mixer.utils.cleaning.SimilarityMatrixTools;
import mixer.utils.drive.BinMappings;
import mixer.utils.kmeans.QuickCentroids;
import mixer.utils.similaritymeasures.RobustCorrelationSimilarity;
import mixer.utils.similaritymeasures.RobustCosineSimilarity;
import mixer.utils.similaritymeasures.SimilarityMetric;
import mixer.utils.tracks.ColorMap;

import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.util.*;
import java.util.concurrent.atomic.AtomicInteger;

public class IndexOrderer {

    private static final int FIVE_MB = 5000000, FIFTY_MB = 50000000, HUNDRED_MB = 100000000;
    private static final int IGNORE = -1;
    private static final int DEFAULT = -5;
    private static final int CHECK_VAL = -2;
    private static final float CORR_MIN = 0.3f;
    public static final boolean SHOULD_FILTER_SINGLE_COLUMNS = false;

    public static BinMappings getInitialMappings(Dataset ds, Chromosome[] chromosomes,
                                                 int hires, Map<Integer, Set<Integer>> badIndices, NormalizationType norm,
                                                 long seed, File outputDirectory, Map<Integer, LogExpectedSpline> splines) {

        Random generator = new Random(seed);
        int[] offset = new int[]{0};
        BinMappings mappings = new BinMappings(hires, chromosomes);

        int lowRes = Math.max(hires, 100000);
        int resFactor = lowRes / hires;
        if (lowRes % hires != 0) {
            System.err.println(hires + "is not a factor of " + lowRes + ". Invalid resolutions.");
            System.exit(23);
        }
        for (Chromosome chrom : chromosomes) {
            final MatrixZoomData zd = HiCFileTools.getMatrixZoomData(ds, chrom, chrom, lowRes);
            if (zd != null) {
                try {

                    if (!splines.containsKey(chrom.getIndex())) {
                        LogExpectedSpline spline = new LogExpectedSpline(zd, norm, chrom, lowRes);
                        splines.put(chrom.getIndex(), spline);
                    }

                    float[][] matrix = OETools.getCleanOEMatrix(zd, chrom, lowRes, norm,
                            badIndices.get(chrom.getIndex()), resFactor, true,
                            true, splines.get(chrom.getIndex()));
                    int[] lowResNewOrderIndexes = getNewOrderOfIndices(chrom, matrix, badIndices.get(chrom.getIndex()),
                            offset, lowRes, generator.nextLong(), resFactor);
                    int[] newOrderIndexes = convertToHigherRes(lowResNewOrderIndexes, chrom, hires, resFactor);

                    mappings.putBinToProtoCluster(chrom, newOrderIndexes);
                } catch (Exception e) {
                    e.printStackTrace();
                }
            }
            System.out.print(".");
        }

        mappings.calculateGlobalIndices(chromosomes);

        writeOutInitialResults(outputDirectory, hires, chromosomes, mappings);

        return mappings;
    }

    public static float[][] quickCleanMatrix(float[][] matrix, int[] newIndexOrderAssignments) {
        List<Integer> actualIndices = new ArrayList<>();
        for (int z = 0; z < newIndexOrderAssignments.length; z++) {
            if (newIndexOrderAssignments[z] < CHECK_VAL && percentNaN(matrix[z]) < .7) {
                actualIndices.add(z);
            }
        }

        float[][] tempCleanMatrix = new float[actualIndices.size()][matrix[0].length];
        for (int i = 0; i < actualIndices.size(); i++) {
            System.arraycopy(matrix[actualIndices.get(i)], 0, tempCleanMatrix[i], 0, tempCleanMatrix[i].length);
        }
        if (SmartTools.printVerboseComments) {
            System.out.println("New clean matrix: " + tempCleanMatrix.length + " rows kept from " + matrix.length);
        }
        return tempCleanMatrix;
    }

    private static int[] convertToHigherRes(int[] lowResOrderIndexes, Chromosome chrom, int hires, int resFactor) {
        int hiResLength = (int) (chrom.getLength() / hires) + 1;
        int[] hiResOrderAssignments = new int[hiResLength];
        if ((hiResLength - 1) / resFactor >= lowResOrderIndexes.length) {
            System.err.println("chromosome lengths are off");
            System.exit(32);
        }

        for (int i = 0; i < hiResOrderAssignments.length; i++) {
            hiResOrderAssignments[i] = lowResOrderIndexes[i / resFactor];
        }
        return hiResOrderAssignments;
    }

    private static int[] getNewOrderOfIndices(Chromosome chromosome, float[][] oeMatrix1,
                                              Set<Integer> badIndices, int[] offset, int lowRes,
                                              long seed, int resFactor) {

        int[] newIndexOrderAssignments = generateNewAssignments(oeMatrix1.length, badIndices, resFactor);
        int numPotentialClusters = (int) (chromosome.getLength() / FIFTY_MB) + 5;

        float[][] matrixCorr1 = SimilarityMatrixTools.getSymmNonNanSimilarityMatrixWithMask(oeMatrix1,
                RobustCosineSimilarity.SINGLETON, newIndexOrderAssignments, CHECK_VAL);
        IntraMatrixCleaner.nanFillBadRowsColumns(badIndices, matrixCorr1, resFactor);
        IntraMatrixCleaner.nanFillNearDiagonal(matrixCorr1, FIVE_MB / lowRes);

        try {
            offset[0] = doAssignmentsByCorrWithCentroids(matrixCorr1, newIndexOrderAssignments, chromosome.getName(),
                    numPotentialClusters, seed, offset[0]);
            //indexToRearrangedLength.put(chromosome.getIndex(), gCounter);
        } catch (Exception e) {
            e.printStackTrace();
        }
        return newIndexOrderAssignments;
    }

    private static int doAssignmentsByCorrWithCentroids(float[][] matrix, int[] newIndexOrderAssignments, String chromName,
                                                        int numInitialClusters, long seed, int offset) {
        float[][] centroids = new QuickCentroids(quickCleanMatrix(matrix, newIndexOrderAssignments),
                numInitialClusters, seed, 20).generateCentroids(10, false);

        List<Integer> problemIndices = Collections.synchronizedList(new ArrayList<>());
        int[] clusterAssignment = new int[newIndexOrderAssignments.length];
        Arrays.fill(clusterAssignment, IGNORE);

        AtomicInteger currDataIndex = new AtomicInteger(0);
        ParallelizationTools.launchParallelizedCode(() -> {
            SimilarityMetric corrMetric = RobustCorrelationSimilarity.SINGLETON;
            int i = currDataIndex.getAndIncrement();
            while (i < (matrix).length) {
                if (newIndexOrderAssignments[i] < CHECK_VAL) {
                    int bestIndex = IGNORE;
                    float bestCorr = CORR_MIN;

                    for (int j = 0; j < centroids.length; j++) {
                        float corrVal = corrMetric.distance(centroids[j], matrix[i]);
                        if (corrVal > bestCorr) {
                            bestCorr = corrVal;
                            bestIndex = j;
                        }
                    }
                    if (bestIndex < 0) {
                        synchronized (problemIndices) {
                            problemIndices.add(bestIndex);
                        }
                    }
                    clusterAssignment[i] = bestIndex + offset;
                }
                i = currDataIndex.getAndIncrement();
            }
        });

        if (SmartTools.printVerboseComments) {
            synchronized (problemIndices) {
                double percentProblem = 100 * (problemIndices.size() + 0.0) / (matrix.length + 0.0);
                System.out.println("IndexOrderer problems: " + problemIndices.size() + " (" + percentProblem + " %)");
            }
        }

        if (SHOULD_FILTER_SINGLE_COLUMNS) {
            int filtered = 0;
            for (int z = 0; z < clusterAssignment.length; z++) {
                int zMinus1 = Math.max(0, z - 1);
                int zPlus1 = Math.min(z + 1, clusterAssignment.length - 1);
                if (clusterAssignment[z] != clusterAssignment[zMinus1]
                        && clusterAssignment[z] != clusterAssignment[zPlus1]) {
                    clusterAssignment[z] = IGNORE;
                    filtered++;
                }
            }
            if (SmartTools.printVerboseComments) {
                System.out.println("Post filtered: " + filtered);
            }
        }

        for (int i = 0; i < clusterAssignment.length; i++) {
            if (newIndexOrderAssignments[i] < CHECK_VAL) {
                newIndexOrderAssignments[i] = clusterAssignment[i];
            }
        }

        return offset + centroids.length;
    }

    private static int[] generateNewAssignments(int length, Set<Integer> badIndices, int resFactor) {
        int[] newIndexOrderAssignments = new int[length];
        Arrays.fill(newIndexOrderAssignments, DEFAULT);
        for (int k : badIndices) {
            newIndexOrderAssignments[k / resFactor] = IGNORE;
        }
        return newIndexOrderAssignments;
    }

    private static void writeOutInitialResults(File outputDirectory, int hires, Chromosome[] chromosomes,
                                               BinMappings mappings) {
        File problemFile = new File(outputDirectory, "problems.bed");
        File initFile = new File(outputDirectory, "initial_split.bed");

        try {
            final FileWriter fwProblem = new FileWriter(problemFile);
            final FileWriter fwInit = new FileWriter(initFile);

            for (Chromosome chrom : chromosomes) {
                int[] vals = mappings.getProtoclusterAssignments(chrom);
                for (int i = 0; i < vals.length; i++) {
                    if (vals[i] < 0) {
                        writeRegionToFile(fwProblem, chrom, i, vals[i], hires);
                    } else {
                        writeRegionToFile(fwInit, chrom, i, vals[i], hires);
                    }
                }
            }

            try {
                fwProblem.close();
                fwInit.close();
            } catch (IOException ww) {
                ww.printStackTrace();
            }
        } catch (IOException e) {
            e.printStackTrace();
        }
    }

    private static void writeRegionToFile(FileWriter fw, Chromosome chrom, int pos, int val, int resolution) throws IOException {
        int x1 = pos * resolution;
        int x2 = x1 + resolution;
        String bedLine = "chr" + chrom.getName() + "\t" + x1 + "\t" + x2 + "\t" + val + "\t" + val
                + "\t.\t" + x1 + "\t" + x2 + "\t" + ColorMap.getColorString(Math.abs(val));
        fw.write(bedLine + "\n");
    }

    public static float percentNaN(float[] array) {
        float counter = 0;
        for (float val : array) {
            if (Float.isNaN(val)) {
                counter++;
            }
        }
        return counter / array.length;
    }
}
