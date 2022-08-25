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

package mixer.utils.slice.cleaning;

import javastraw.reader.Dataset;
import javastraw.reader.basics.Chromosome;
import javastraw.reader.basics.ChromosomeHandler;
import javastraw.reader.mzd.MatrixZoomData;
import javastraw.reader.type.NormalizationType;
import javastraw.tools.HiCFileTools;
import mixer.MixerTools;
import mixer.clt.ParallelizedMixerTools;
import mixer.utils.common.ArrayTools;
import mixer.utils.similaritymeasures.RobustCorrelationSimilarity;
import mixer.utils.similaritymeasures.SimilarityMetric;
import mixer.utils.slice.drive.BinMappings;
import mixer.utils.slice.structures.ColorMap;

import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.util.*;
import java.util.concurrent.atomic.AtomicInteger;

public class IndexOrderer {

    private final Map<Chromosome, int[]> chromToReorderedIndices = new HashMap<>();
    private static final int DISTANCE = 5000000, FIVE_MB = 5000000, FIFTY_MB = 50000000;
    private static final int IGNORE = -1;
    private static final int DEFAULT = -5;
    private static final int CHECK_VAL = -2;
    private static final float CORR_MIN = 0.2f;

    public static BinMappings getInitialMappings(Dataset ds, ChromosomeHandler handler,
                                                 int resolution,
                                                 Map<Integer, Set<Integer>> badIndices, NormalizationType norm,
                                                 long seed, File outputDirectory) {
        Chromosome[] chromosomes = handler.getAutosomalChromosomesArray();
        Random generator = new Random(seed);
        int[] offset = new int[]{0};
        BinMappings mappings = new BinMappings(resolution);

        for (Chromosome chrom : chromosomes) {
            final MatrixZoomData zd = HiCFileTools.getMatrixZoomData(ds, chrom, chrom, resolution);
            if (zd != null) {
                try {
                    float[][] matrix = OETools.getOEMatrix(ds, zd, chrom, resolution, norm);
                    int[] newOrderIndexes = getNewOrderOfIndices(chrom, matrix, badIndices.get(chrom.getIndex()),
                            offset, resolution, generator.nextLong());

                    mappings.putBinToProtoCluster(chrom, newOrderIndexes);
                } catch (Exception e) {
                    e.printStackTrace();
                }
            }
            System.out.print(".");
        }

        mappings.calculateGlobalIndices(chromosomes);

        writeOutInitialResults(outputDirectory, resolution, chromosomes, mappings);

        return mappings;
    }

    public static float[][] quickCleanMatrix(float[][] matrix, int[] newIndexOrderAssignments) {
        List<Integer> actualIndices = new ArrayList<>();
        for (int z = 0; z < newIndexOrderAssignments.length; z++) {
            if (newIndexOrderAssignments[z] < CHECK_VAL && ArrayTools.percentNaN(matrix[z]) < .7) {
                actualIndices.add(z);
            }
        }

        float[][] tempCleanMatrix = new float[actualIndices.size()][matrix[0].length];
        for (int i = 0; i < actualIndices.size(); i++) {
            System.arraycopy(matrix[actualIndices.get(i)], 0, tempCleanMatrix[i], 0, tempCleanMatrix[i].length);
        }
        if (MixerTools.printVerboseComments) {
            System.out.println("New clean matrix: " + tempCleanMatrix.length + " rows kept from " + matrix.length);
        }
        return tempCleanMatrix;
    }

    public int[] get(Chromosome chrom) {
        return chromToReorderedIndices.get(chrom);
    }

    private static int[] getNewOrderOfIndices(Chromosome chromosome, float[][] oeMatrix1,
                                              Set<Integer> badIndices, int[] offset, int lowres,
                                              long seed) {

        IntraMatrixCleaner.oeClean(oeMatrix1, badIndices);
        int[] newIndexOrderAssignments = generateNewAssignments(oeMatrix1.length, badIndices);
        int numPotentialClusters = (int) (chromosome.getLength() / FIFTY_MB) + 7;

        float[][] matrixCorr1 = SimilarityMatrixTools.getSymmNonNanSimilarityMatrixWithMask(oeMatrix1,
                RobustCorrelationSimilarity.SINGLETON, newIndexOrderAssignments, CHECK_VAL);
        IntraMatrixCleaner.basicClean(matrixCorr1, badIndices, FIVE_MB / lowres);

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
        ParallelizedMixerTools.launchParallelizedCode(() -> {
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

        if (MixerTools.printVerboseComments) {
            synchronized (problemIndices) {
                double percentProblem = 100 * (problemIndices.size() + 0.0) / (matrix.length + 0.0);
                System.out.println("IndexOrderer problems: " + problemIndices.size() + " (" + percentProblem + " %)");
            }
        }

        for (int i = 0; i < clusterAssignment.length; i++) {
            if (newIndexOrderAssignments[i] < CHECK_VAL) {
                newIndexOrderAssignments[i] = clusterAssignment[i];
            }
        }

        return offset + centroids.length;
    }

    private static int[] generateNewAssignments(int length, Set<Integer> badIndices) {
        int[] newIndexOrderAssignments = new int[length];
        Arrays.fill(newIndexOrderAssignments, DEFAULT);
        for (int k : badIndices) {
            newIndexOrderAssignments[k] = IGNORE;
        }
        return newIndexOrderAssignments;
    }

    private static void writeOutInitialResults(File outputDirectory, int hires, Chromosome[] chromosomes, BinMappings mappings) {
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
}
