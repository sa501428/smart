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
import javastraw.reader.block.ContactRecord;
import javastraw.reader.type.NormalizationType;
import javastraw.tools.ParallelizationTools;
import mixer.utils.common.SimpleArray2DTools;
import mixer.utils.matrix.*;
import mixer.utils.tracks.SubcompartmentInterval;

import java.io.File;
import java.util.*;
import java.util.concurrent.atomic.AtomicInteger;

public class Shuffle {

    private final static int NUM_SCORES = 2;
    private final Dataset ds;
    private final NormalizationType norm;
    private final int compressionFactor;
    private final int resolution;
    private final int numRounds = 50;
    private final InterOnlyMatrix.InterMapType[] mapTypes = {InterOnlyMatrix.InterMapType.ODDS_VS_EVENS,
            InterOnlyMatrix.InterMapType.SKIP_BY_TWOS, InterOnlyMatrix.InterMapType.FIRST_HALF_VS_SECOND_HALF};
    private final ScoreContainer scoreContainer = new ScoreContainer(mapTypes.length, NUM_SCORES);

    public Shuffle(Dataset ds, NormalizationType norm, int resolution, int compressionFactor) {
        this.resolution = resolution;
        this.compressionFactor = compressionFactor;
        this.ds = ds;
        this.norm = norm;
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
                float[][] matrix = getShuffledMatrix(interMatrix, allRowIndices, allColIndices);
                SimpleArray2DTools.simpleLogWithCleanup(matrix, Float.NaN);
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
        aggregate.saveToPNG(outfolder, name);
        scoreContainer.updateAggregateScores(aggregate, globalAllIndices, mapIndex);
    }

    private ShuffledIndices getShuffledByClusterIndices(Map<Integer, List<Integer>> clusterToIndices,
                                                        Random generator) {
        List<Integer> allIndices = new ArrayList<>(clusterToIndices.size());
        List<Integer> order = new ArrayList<>(clusterToIndices.keySet());
        Collections.sort(order);

        int boundaryCount = 0;
        List<Integer> boundaries = new ArrayList<>();
        boundaries.add(boundaryCount);

        for (Integer index : order) {
            List<Integer> indexList = new ArrayList<>(clusterToIndices.get(index));
            Collections.shuffle(indexList, generator);
            int numToUse = (indexList.size() / compressionFactor) * compressionFactor;
            for (int z = 0; z < numToUse; z++) {
                allIndices.add(indexList.get(z));
            }
            boundaryCount += (numToUse / compressionFactor);
            boundaries.add(boundaryCount);
        }
        Integer[] output = new Integer[boundaries.size()];
        return new ShuffledIndices(allIndices, boundaries.toArray(output), compressionFactor);
    }

    private float[][] getShuffledMatrix(HiCMatrix interMatrix,
                                        ShuffledIndices rowIndices,
                                        ShuffledIndices colIndices) {
        Map<Integer, Integer> rowMap = rowIndices.allIndices;
        Map<Integer, Integer> colMap = colIndices.allIndices;
        int numRows = rowIndices.boundaries[rowIndices.boundaries.length - 1];
        int numCols = colIndices.boundaries[colIndices.boundaries.length - 1];

        float[][] result = new float[numRows][numCols];
        for (List<ContactRecord> records : interMatrix.getMatrix()) {
            for (ContactRecord record : records) {
                if (rowMap.containsKey(record.getBinX()) && colMap.containsKey(record.getBinY())) {
                    result[rowMap.get(record.getBinX())][colMap.get(record.getBinY())] += record.getCounts();
                }
            }
        }
        return result;
    }
}
