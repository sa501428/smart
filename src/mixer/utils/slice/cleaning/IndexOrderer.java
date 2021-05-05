/*
 * The MIT License (MIT)
 *
 * Copyright (c) 2011-2021 Rice University, Baylor College of Medicine, Aiden Lab
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
import javastraw.reader.mzd.MatrixZoomData;
import javastraw.reader.type.NormalizationType;
import javastraw.tools.ExtractingOEDataUtils;
import javastraw.tools.HiCFileTools;
import mixer.MixerGlobals;
import mixer.clt.ParallelizedMixerTools;
import mixer.utils.common.ArrayTools;
import mixer.utils.similaritymeasures.RobustCorrelationSimilarity;
import mixer.utils.similaritymeasures.SimilarityMetric;
import mixer.utils.slice.structures.SubcompartmentColors;

import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.util.*;
import java.util.concurrent.atomic.AtomicInteger;

public class IndexOrderer {

    private final Map<Chromosome, int[]> chromToReorderedIndices = new HashMap<>();
    private final int DISTANCE = 5000000, ONE_HUNDRED_KB = 100000, TEN_MB = 10000000;
    private final int IGNORE = -1;
    private final int DEFAULT = -5;
    private final int CHECK_VAL = -2;
    private final float CORR_MIN = 0.4f;
    private final Random generator = new Random(0);
    private final Map<Integer, Integer> indexToRearrangedLength = new HashMap<>();
    private final File problemFile, initFile;
    private final int resolution;

    public IndexOrderer(Dataset ds, Chromosome[] chromosomes, int resolution, NormalizationType normalizationType,
                        GWBadIndexFinder badIndexLocations, long seed, File outputDirectory,
                        int maxClusterSizeExpected) {
        this.resolution = resolution;
        problemFile = new File(outputDirectory, "problems.bed");
        initFile = new File(outputDirectory, "initial_split.bed");

        generator.setSeed(seed);
        int smoothingInterval = ONE_HUNDRED_KB / resolution;
        for (Chromosome chrom : chromosomes) {
            final MatrixZoomData zd = HiCFileTools.getMatrixZoomData(ds, chrom, chrom, resolution);
            try {
                float[][] matrix = HiCFileTools.getOEMatrixForChromosome(ds, zd, chrom, resolution,
                        normalizationType, 5f, ExtractingOEDataUtils.ThresholdType.TRUE_OE,
                        true, 1, 0);
                Set<Integer> badIndices = badIndexLocations.getBadIndices(chrom);

                matrix = IntraMatrixCleaner.cleanAndCompress(chrom, matrix, resolution, smoothingInterval, badIndices);
                int[] newOrderIndexes = getNewOrderOfIndices(chrom, matrix, badIndices);
                chromToReorderedIndices.put(chrom, newOrderIndexes);
            } catch (Exception e) {
                e.printStackTrace();
            }
            System.out.print(".");
        }

        writeOutInitialResults();
    }

    public Map<Integer, Integer> getIndexToRearrangedLength() {
        return indexToRearrangedLength;
    }

    public int[] get(Chromosome chrom) {
        return chromToReorderedIndices.get(chrom);
    }

    private int[] getNewOrderOfIndices(Chromosome chromosome, float[][] matrix, Set<Integer> badIndices) {
        int[] newIndexOrderAssignments = generateNewAssignments(matrix.length, badIndices);
        int numPotentialClusters = (int) (chromosome.getLength() / TEN_MB);
        //numPotentialClusters = Math.max(numPotentialClusters, maxClusterSizeExpected);

        int gCounter = doAssignmentsByCorrWithCentroids(matrix, newIndexOrderAssignments, chromosome.getName(),
                numPotentialClusters);
        indexToRearrangedLength.put(chromosome.getIndex(), gCounter);
        return newIndexOrderAssignments;
    }

    private int[] generateNewAssignments(int length, Set<Integer> badIndices) {
        int[] newIndexOrderAssignments = new int[length];
        Arrays.fill(newIndexOrderAssignments, DEFAULT);
        for (int k : badIndices) {
            newIndexOrderAssignments[k] = IGNORE;
        }
        return newIndexOrderAssignments;
    }

    private int doAssignmentsByCorrWithCentroids(float[][] matrix, int[] newIndexOrderAssignments, String chromName,
                                                 int numInitialClusters) {
        float[][] centroids = new QuickCentroids(quickCleanMatrix(matrix, newIndexOrderAssignments),
                numInitialClusters, generator.nextLong(), 100).generateCentroids(5);

        if (true || MixerGlobals.printVerboseComments) {
            System.out.println("IndexOrderer: num centroids (init " + numInitialClusters + ") for " + chromName + ": " + centroids.length);
        }

        List<Integer> problemIndices = Collections.synchronizedList(new ArrayList<>());

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
                    newIndexOrderAssignments[i] = bestIndex;
                }
                i = currDataIndex.getAndIncrement();
            }
        });


        synchronized (problemIndices) {
            double percentProblem = 100 * (problemIndices.size() + 0.0) / (matrix.length + 0.0);
            System.out.println("IndexOrderer problems: " + problemIndices.size() + " (" + percentProblem + " %)");
        }


        return centroids.length;
    }

    private float[][] quickCleanMatrix(float[][] matrix, int[] newIndexOrderAssignments) {
        List<Integer> actualIndices = new ArrayList<>();
        for (int z = 0; z < newIndexOrderAssignments.length; z++) {
            if (newIndexOrderAssignments[z] < CHECK_VAL && ArrayTools.percentNaN(matrix[z]) < .5) {
                actualIndices.add(z);
            }
        }

        float[][] tempCleanMatrix = new float[actualIndices.size()][matrix[0].length];
        for (int i = 0; i < actualIndices.size(); i++) {
            System.arraycopy(matrix[actualIndices.get(i)], 0, tempCleanMatrix[i], 0, tempCleanMatrix[i].length);
        }
        if (true || MixerGlobals.printVerboseComments) {
            System.out.println("New clean matrix: " + tempCleanMatrix.length + " rows kept from " + matrix.length);
        }
        return tempCleanMatrix;
    }

    private void writeOutInitialResults() {
        try {
            final FileWriter fwProblem = new FileWriter(problemFile);
            final FileWriter fwInit = new FileWriter(initFile);

            for (Chromosome chrom : chromToReorderedIndices.keySet()) {
                int[] vals = chromToReorderedIndices.get(chrom);
                for (int i = 0; i < vals.length; i++) {
                    if (vals[i] < 0) {
                        writeRegionToFile(fwProblem, chrom, i, vals[i]);
                    } else {
                        writeRegionToFile(fwInit, chrom, i, vals[i]);
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

    private void writeRegionToFile(FileWriter fw, Chromosome chrom, int pos, int val) throws IOException {

        int x1 = resolution * pos;
        int x2 = x1 + resolution;
        String bedLine = "chr" + chrom.getName() + "\t" + x1 + "\t" + x2 + "\t" + val + "\t" + val
                + "\t.\t" + x1 + "\t" + x2 + "\t" + SubcompartmentColors.getColorString(Math.abs(val));
        fw.write(bedLine + "\n");
    }
}
