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

import javastraw.tools.MatrixTools;
import mixer.MixerGlobals;
import mixer.utils.common.FloatMatrixTools;
import mixer.utils.common.ZScoreTools;
import mixer.utils.similaritymeasures.SimilarityMetric;
import mixer.utils.slice.structures.SubcompartmentInterval;
import tagbio.umap.Umap;

import java.io.File;
import java.util.*;

public class MatrixCleanerAndProjector {
    private final static float PERCENT_NAN_ALLOWED = .5f;
    public static int NUM_PER_CENTROID = 10;
    protected final File outputDirectory;
    protected float[][] data;
    protected final Random generator = new Random(0);
    private final SimilarityMetric metric;
    private static final float ZERO = 1e-10f;

    public MatrixCleanerAndProjector(float[][] data, long seed, File outputDirectory, SimilarityMetric metric) {
        this.metric = metric;
        this.outputDirectory = outputDirectory;
        generator.setSeed(seed);

        // after this, there should be no infinities, no negative numbers
        // just real numbers >= 0 or NaNs
        this.data = data;
        //simpleLogWithCleanup(this.data);
        System.out.println("matrix size " + data.length + " x " + data[0].length);
    }

    public static void simpleLogWithCleanup(float[][] matrix) {
        for (int i = 0; i < matrix.length; i++) {
            for (int j = 0; j < matrix[i].length; j++) {
                float val = matrix[i][j];
                if (!Float.isNaN(val)) {
                    val = (float) Math.log(val + 1);
                    if (Float.isInfinite(val)) {
                        matrix[i][j] = Float.NaN;
                    } else {
                        matrix[i][j] = val;
                    }
                }
            }
        }
    }

    public static Set<Integer> getBadIndices(float[][] matrix) {
        Set<Integer> badIndices = new HashSet<>();
        int[] numNans = getNumberOfNansOrZerosInRow(matrix);
        int maxBadEntriesAllowed = (int) Math.ceil(PERCENT_NAN_ALLOWED * matrix[0].length);
        for (int i = 0; i < numNans.length; i++) {
            if (numNans[i] > maxBadEntriesAllowed) {
                badIndices.add(i);
            }
        }

        return badIndices;
    }

    public static float[][] filterOutColumnsAndRowsNonSymmetricMatrix(float[][] interMatrix, Map<Integer, SubcompartmentInterval> original) {
        Set<Integer> badIndices = getBadIndices(interMatrix);
        if (badIndices.size() == 0) {
            return interMatrix;
        }

        if (MixerGlobals.printVerboseComments) {
            System.out.println("interMatrix.length " + interMatrix.length + " badIndices.size() " + badIndices.size());
        }

        int counter = 0;
        int[] newIndexToOrigIndex = new int[interMatrix.length - badIndices.size()];
        for (int i = 0; i < interMatrix.length; i++) {
            if (!badIndices.contains(i)) {
                newIndexToOrigIndex[counter++] = i;
            }
        }

        float[][] newMatrix = new float[newIndexToOrigIndex.length][interMatrix[0].length];
        Map<Integer, SubcompartmentInterval> newRowIndexToIntervalMap = new HashMap<>();
        for (int i = 0; i < newMatrix.length; i++) {
            int tempI = newIndexToOrigIndex[i];
            System.arraycopy(interMatrix[tempI], 0, newMatrix[i], 0, newMatrix[0].length);
            newRowIndexToIntervalMap.put(i, (SubcompartmentInterval) original.get(newIndexToOrigIndex[i]).deepClone());
        }

        original.clear();
        original.putAll(newRowIndexToIntervalMap);

        return newMatrix;
    }

    private static int[] getNumberOfNansOrZerosInRow(float[][] matrix) {
        int[] numInvalids = new int[matrix.length];
        for (int i = 0; i < matrix.length; i++) {
            for (int j = 0; j < matrix[i].length; j++) {
                float val = matrix[i][j];
                if (Float.isNaN(val) || val < ZERO) {
                    numInvalids[i]++;
                }
            }
        }
        return numInvalids;
    }

    public float[][] getCleanedSimilarityMatrix(Map<Integer, SubcompartmentInterval> rowIndexToIntervalMap,
                                                int[] weights) {
        data = filterOutColumnsAndRowsNonSymmetricMatrix(data, rowIndexToIntervalMap);
        if (MixerGlobals.printVerboseComments) {
            System.out.println("matrix size " + data.length + " x " + data[0].length);
        }

        //ZScoreTools.inPlaceRobustZscoreDownCol(data, weights);
        if (MixerGlobals.printVerboseComments) {
            File temp = new File(outputDirectory, "zscore_new.npy");
            FloatMatrixTools.saveMatrixTextNumpy(temp.getAbsolutePath(), data);
        }

        System.out.println("Generating similarity matrix");
        data = SimilarityMatrixTools.getNonNanSimilarityMatrix(data, metric, NUM_PER_CENTROID, generator.nextLong());

        if (MixerGlobals.printVerboseComments) {
            System.out.println("similarity matrix size " + data.length + " x " + data[0].length);
        }

        //ZScoreTools.inPlaceRobustZscoreDownCol(data);
        ZScoreTools.inPlaceZscoreDownCol(data);

        File temp = new File(outputDirectory, "data_new.npy");
        FloatMatrixTools.saveMatrixTextNumpy(temp.getAbsolutePath(), data);

        System.out.println("Running UMAP");
        runUmapAndSaveMatrices(data, outputDirectory, rowIndexToIntervalMap);
        System.out.println("Done running UMAP");

        return data;
    }

    private void runUmapAndSaveMatrices(float[][] data, File outputDirectory,
                                        Map<Integer, SubcompartmentInterval> rowIndexToIntervalMap) {
        int numThreads = Math.max(1, Runtime.getRuntime().availableProcessors());
        for (int dimensions = 2; dimensions < 4; dimensions++) {
            System.out.println("Running UMAP with " + dimensions + " dimensions");
            final Umap umap = new Umap();
            umap.setNumberComponents(dimensions);
            umap.setNumberNearestNeighbours(50); // 15 -> 50 for more global picture
            umap.setThreads(numThreads);
            umap.setMinDist(0.2f);  //0.1f -> 0.2f for more general features
            umap.setVerbose(false);
            umap.setSeed(generator.nextLong());
            final float[][] result = umap.fitTransform(data);
            File temp = new File(outputDirectory, "umap_" + dimensions + "d_embedding.npy");
            FloatMatrixTools.saveMatrixTextNumpy(temp.getAbsolutePath(), result);
            System.gc();
        }

        int[][] indices = new int[data.length][2];
        for (int i = 0; i < indices.length; i++) {
            SubcompartmentInterval interval = rowIndexToIntervalMap.get(i);
            indices[i][0] = interval.getChrIndex();
            indices[i][1] = interval.getX1();
        }
        File temp = new File(outputDirectory, "genome_indices.npy");
        MatrixTools.saveMatrixTextNumpy(temp.getAbsolutePath(), indices);
    }
}
