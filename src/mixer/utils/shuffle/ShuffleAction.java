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

package mixer.utils.shuffle;

import javastraw.feature1D.GenomeWide1DList;
import javastraw.reader.Dataset;
import javastraw.reader.type.NormalizationType;
import javastraw.tools.ParallelizationTools;
import mixer.utils.common.FloatMatrixTools;
import mixer.utils.matrix.*;
import mixer.utils.tracks.SubcompartmentInterval;

import java.io.File;
import java.util.*;
import java.util.concurrent.atomic.AtomicInteger;

public class ShuffleAction {

    private final static int NUM_SCORES = 2;
    private final Dataset ds;
    private final NormalizationType norm;
    private final int compressionFactor;
    private final int resolution;
    private final int numRounds = 10;
    private final Partition.Type[] mapTypes;
    private final ScoreContainer scoreContainer;

    public ShuffleAction(Dataset ds, NormalizationType norm, int resolution, int compressionFactor,
                         Partition.Type[] maptypes) {
        this.resolution = resolution;
        this.compressionFactor = compressionFactor;
        this.ds = ds;
        this.norm = norm;
        this.mapTypes = maptypes;
        this.scoreContainer = new ScoreContainer(mapTypes.length, NUM_SCORES);
    }

    private static long[] getSeedsForRound(Random generator, int numRounds) {
        long[] seeds = new long[numRounds];
        for (int i = 0; i < numRounds; i++) {
            seeds[i] = generator.nextLong();
        }
        return seeds;
    }

    public void runInterAnalysis(GenomeWide1DList<SubcompartmentInterval> subcompartments, File outfolder,
                                 Random generator) {
        for (int y = 0; y < mapTypes.length; y++) {
            final HiCMatrix interMatrix = InterOnlyMatrix.getMatrix(ds, norm, resolution, mapTypes[y]);
            Map<Integer, List<Integer>> clusterToRowIndices = CHICTools.populateCluster(interMatrix.getRowChromosomes(),
                    interMatrix.getRowOffsets(), subcompartments, resolution);
            Map<Integer, List<Integer>> clusterToColIndices = CHICTools.populateCluster(interMatrix.getColChromosomes(),
                    interMatrix.getColOffsets(), subcompartments, resolution);
            shuffleMap(interMatrix, clusterToRowIndices, clusterToColIndices, outfolder, mapTypes[y].toString(), y, generator);
        }
        scoreContainer.calculateRatios();
    }

    public void savePlotsAndResults(File outfolder, String prefix) {
        String[] names;
        names = new String[mapTypes.length];
        for (int k = 0; k < mapTypes.length; k++) {
            names[k] = mapTypes[k].toString();
        }
        scoreContainer.savePlotsAndResults(outfolder, prefix, names);

    }

    public double getResult(int index) {
        return scoreContainer.getDirectScore(false, index);
    }

    private void shuffleMap(HiCMatrix interMatrix, Map<Integer, List<Integer>> clusterToRowIndices,
                            Map<Integer, List<Integer>> clusterToColIndices,
                            File outfolder, String name, int mapIndex, Random random) {

        final AggregateMatrix aggregate = new AggregateMatrix();
        AtomicInteger currRowIndex = new AtomicInteger(0);
        long[] seeds = getSeedsForRound(random, numRounds);

        final ShuffledIndices[] globalAllIndices = new ShuffledIndices[2];
        Random gen = new Random(random.nextLong());
        globalAllIndices[0] = getShuffledByClusterIndices(clusterToRowIndices, gen);
        globalAllIndices[1] = getShuffledByClusterIndices(clusterToColIndices, gen);

        ParallelizationTools.launchParallelizedCode(() -> {
            int k = currRowIndex.getAndIncrement();
            AggregateMatrix aggregateForThread = new AggregateMatrix();
            while (k < numRounds) {
                Random generator = new Random(seeds[k]);
                ShuffledIndices allRowIndices = getShuffledByClusterIndices(clusterToRowIndices, generator);
                ShuffledIndices allColIndices = getShuffledByClusterIndices(clusterToColIndices, generator);

                double[][] matrix = getShuffledMatrix(interMatrix, allRowIndices.allIndices, allColIndices.allIndices);
                FloatMatrixTools.log(matrix, 1);

                aggregateForThread.add(matrix);
                k = currRowIndex.getAndIncrement();
            }

            synchronized (aggregate) {
                if (aggregateForThread.hasData()) {
                    aggregate.add(aggregateForThread);
                }
            }
        });

        aggregate.scaleForNumberOfRounds(numRounds);
        if (outfolder != null) aggregate.export(outfolder, name);
        scoreContainer.updateAggregateScores(aggregate, globalAllIndices, mapIndex);
    }

    private ShuffledIndices getShuffledByClusterIndices(Map<Integer, List<Integer>> clusterToIndices,
                                                        Random generator) {
        List<Integer> allIndices = new ArrayList<>();
        List<Integer> order = new ArrayList<>(clusterToIndices.keySet());
        Collections.sort(order);

        int count = 0;
        List<Integer> boundaries = new ArrayList<>();
        boundaries.add(count);

        for (Integer clusterID : order) {
            List<Integer> indexList = new ArrayList<>(clusterToIndices.get(clusterID));
            Collections.shuffle(indexList, generator);
            int numToUse = (indexList.size() / compressionFactor) * compressionFactor;
            for (int z = 0; z < numToUse; z++) {
                allIndices.add(indexList.get(z));
            }
            count += (numToUse / compressionFactor);
            boundaries.add(count);
        }

        return new ShuffledIndices(allIndices, boundaries.toArray(new Integer[0]), order.toArray(new Integer[0]));
    }

    private double[][] getShuffledMatrix(HiCMatrix interMatrix, List<Integer> allRowIndices, List<Integer> allColIndices) {
        int numRows = allRowIndices.size() / compressionFactor;
        int numCols = allColIndices.size() / compressionFactor;
        int numRowsKept = numRows * compressionFactor;
        int numColsKept = numCols * compressionFactor;
        float[][] original = interMatrix.getMatrix();

        double[][] result = new double[numRows][numCols];
        for (int i = 0; i < numRowsKept; i++) {
            final int i0 = allRowIndices.get(i);
            for (int j = 0; j < numColsKept; j++) {
                final int j0 = allColIndices.get(j);
                result[i / compressionFactor][j / compressionFactor] += original[i0][j0];
            }
        }

        return result;
    }
}
