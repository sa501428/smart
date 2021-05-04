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

import mixer.MixerGlobals;
import mixer.utils.common.LogTools;
import mixer.utils.common.ZScoreTools;
import mixer.utils.similaritymeasures.SimilarityMetric;
import mixer.utils.slice.cleaning.utils.ColumnCleaner;
import mixer.utils.slice.cleaning.utils.RowCleaner;
import mixer.utils.slice.structures.SubcompartmentInterval;

import java.io.File;
import java.util.Map;
import java.util.Random;

public class SliceMatrixCleaner {

    public static int NUM_PER_CENTROID = 100;
    protected final File outputDirectory;
    protected float[][] data;
    protected final Random generator = new Random(0);

    public SliceMatrixCleaner(float[][] data, long seed, File outputDirectory, SimilarityMetric metric) {
        this.outputDirectory = outputDirectory;
        generator.setSeed(seed);

        // after this, there should be no infinities, no negative numbers
        // just real numbers >= 0 or NaNs
        this.data = data;
        System.out.println("Initial matrix size " + data.length + " x " + data[0].length);
    }

    public float[][] getCleanedSimilarityMatrix(Map<Integer, SubcompartmentInterval> rowIndexToIntervalMap,
                                                int[] weights) {
        LogTools.simpleLogWithCleanup(this.data, Float.NaN);
        System.out.println("Initial matrix size " + data.length + " x " + data[0].length);
        data = (new ColumnCleaner(data)).getCleanedData();
        System.out.println("Matrix size after column cleanup " + data.length + " x " + data[0].length);
        data = (new RowCleaner(data, rowIndexToIntervalMap)).getCleanedData();
        System.out.println("Matrix size after row cleanup " + data.length + " x " + data[0].length);

        if (MixerGlobals.printVerboseComments) {
            System.out.println("matrix size " + data.length + " x " + data[0].length);
        }

        //System.out.println("Generating similarity matrix");
        //data = SimilarityMatrixTools.getZscoredNonNanSimilarityMatrix(data, metric, NUM_PER_CENTROID, generator.nextLong());
        /*


        if (MixerGlobals.printVerboseComments) {
            System.out.println("similarity matrix size " + data.length + " x " + data[0].length);
        }

        File temp = new File(outputDirectory, "data_new.npy");
        FloatMatrixTools.saveMatrixTextNumpy(temp.getAbsolutePath(), data);

        /*
        System.out.println("Running UMAP");
        runUmapAndSaveMatrices(data, outputDirectory, rowIndexToIntervalMap);
        System.out.println("Done running UMAP");
        */

        //int[] newWeights = new int[data[0].length];
        //Arrays.fill(newWeights, 1);
        ZScoreTools.inPlaceZscoreDownCol(data);
        ZScoreTools.inPlaceScaleSqrtWeightCol(data, weights);
        //whitener
        //        zscore, scale

        return data;
    }

    /*
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

    public float[][] justRemoveBadRows(Set<Integer> badIndices, Map<Integer, SubcompartmentInterval> rowIndexToIntervalMap, int[] weights) {
        return filterOutColumnsAndRowsGivenBadIndices(badIndices, data, rowIndexToIntervalMap);
    }
    */
}
