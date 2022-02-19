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
import javastraw.reader.basics.Chromosome;
import javastraw.reader.type.NormalizationType;
import mixer.clt.ParallelizedMixerTools;
import mixer.utils.common.FloatMatrixTools;
import mixer.utils.matrix.*;
import mixer.utils.shuffle.scoring.KLDivergenceScoring;
import mixer.utils.shuffle.scoring.VarianceScoring;
import mixer.utils.shuffle.stats.GenomeWideStatistics;
import mixer.utils.similaritymeasures.SimilarityMetric;
import mixer.utils.slice.cleaning.VectorOrganizer;
import mixer.utils.slice.structures.SliceUtils;
import mixer.utils.slice.structures.SubcompartmentInterval;

import java.io.File;
import java.util.*;
import java.util.concurrent.atomic.AtomicInteger;

public class ShuffleAction {

    private final static int NUM_SCORES = 2;
    //private final Random generator = new Random(0);
    private final Dataset ds;
    private final NormalizationType norm;
    private final int compressionFactor;
    private final int resolution;
    private final int numRounds = 50;
    private final InterOnlyMatrix.InterMapType[] mapTypes = {InterOnlyMatrix.InterMapType.ODDS_VS_EVENS,
            InterOnlyMatrix.InterMapType.SKIP_BY_TWOS, InterOnlyMatrix.InterMapType.FIRST_HALF_VS_SECOND_HALF};
    private final SimilarityMetric metric;
    private ScoreContainer scoreContainer = new ScoreContainer(mapTypes.length, NUM_SCORES);
    private Chromosome[] chromosomes;
    private HiCMatrix.INTRA_TYPE intraType = HiCMatrix.INTRA_TYPE.DEFAULT;
    private boolean isIntra = false;

    public ShuffleAction(Dataset ds, NormalizationType norm, int resolution, int compressionFactor,
                         SimilarityMetric metric) {
        this.resolution = resolution;
        this.compressionFactor = compressionFactor;
        this.ds = ds;
        this.norm = norm;
        this.metric = metric;
    }

    public ShuffleAction(Dataset ds, NormalizationType norm, int resolution, int compressionFactor,
                         SimilarityMetric metric, InterOnlyMatrix.INTRA_TYPE intraType) {
        this(ds, norm, resolution, compressionFactor, metric);
        isIntra = true;
        this.intraType = intraType;
        chromosomes = ds.getChromosomeHandler().getAutosomalChromosomesArray();
        scoreContainer = new ScoreContainer(chromosomes.length, 2);
    }

    public static void updateMatrixScores(double[][] scores, int k, float[][] matrix, Integer[] rowBounds, Integer[] colBounds,
                                          boolean isBaseline) {
        scores[0][k] = (new VarianceScoring(matrix, rowBounds, colBounds)).score(isBaseline);
        scores[1][k] = (new KLDivergenceScoring(matrix, rowBounds, colBounds, true)).score(isBaseline);
    }

    public void runGWStats(GenomeWide1DList<SubcompartmentInterval> subcompartments, File outfolder) {
        SliceUtils.collapseGWList(subcompartments);
        // todo, using only big size? todo sorting picture
        GenomeWideStatistics statistics = new GenomeWideStatistics(ds, resolution, norm, subcompartments);
        VectorOrganizer vo = new VectorOrganizer(statistics.getBasicResult());
        System.out.println("Interaction summary statistics saved");
    }

    public void runInterAnalysis(GenomeWide1DList<SubcompartmentInterval> subcompartments, File outfolder,
                                 Random generator) {
        for (int y = 0; y < mapTypes.length; y++) {
            final HiCMatrix interMatrix = InterOnlyMatrix.getMatrix(ds, norm, resolution, mapTypes[y], metric);
            Map<Integer, List<Integer>> clusterToRowIndices = populateCluster(interMatrix.getRowChromosomes(), interMatrix.getRowOffsets(), subcompartments);
            Map<Integer, List<Integer>> clusterToColIndices = populateCluster(interMatrix.getColChromosomes(), interMatrix.getColOffsets(), subcompartments);
            shuffleMap(interMatrix, clusterToRowIndices, clusterToColIndices, outfolder, mapTypes[y].toString(), y, generator);
        }
        scoreContainer.calculateRatios();
    }

    public void runIntraAnalysis(GenomeWide1DList<SubcompartmentInterval> subcompartments, File outfolder, Random generator) {
        for (int y = 0; y < chromosomes.length; y++) {
            final HiCMatrix matrix = new IntraOnlyMatrix(ds, norm, resolution, chromosomes[y], intraType, metric, compressionFactor);
            Map<Integer, List<Integer>> clusterToRowIndices = populateCluster(matrix.getRowChromosomes(), matrix.getRowOffsets(), subcompartments);
            shuffleMap(matrix, clusterToRowIndices, clusterToRowIndices, outfolder, chromosomes[y].getName(), y, generator);
        }
        scoreContainer.calculateRatios();
    }

    public void savePlotsAndResults(File outfolder, String prefix) {
        String[] names;
        if (isIntra) {
            names = new String[chromosomes.length];
            for (int k = 0; k < chromosomes.length; k++) {
                names[k] = chromosomes[k].getName();
            }
        } else {
            names = new String[mapTypes.length];
            for (int k = 0; k < mapTypes.length; k++) {
                names[k] = mapTypes[k].toString();
            }
        }
        scoreContainer.savePlotsAndResults(outfolder, prefix, names);

    }

    private void shuffleMap(HiCMatrix interMatrix, Map<Integer, List<Integer>> clusterToRowIndices,
                            Map<Integer, List<Integer>> clusterToColIndices,
                            File outfolder, String name, int mapIndex, Random generator) {
        boolean isBaseline = false;

        final AggregateMatrix aggregate = new AggregateMatrix(isBaseline);
        double[][] scoresBaselineForRound = new double[NUM_SCORES][numRounds];
        double[][] logScoresBaselineForRound = new double[NUM_SCORES][numRounds];
        double[][] scoresForRound = new double[NUM_SCORES][numRounds];
        double[][] logScoresForRound = new double[NUM_SCORES][numRounds];


        AtomicInteger currRowIndex = new AtomicInteger(0);


        final ShuffledIndices[] globalAllIndices = new ShuffledIndices[2];
        ParallelizedMixerTools.launchParallelizedCode(() -> {
            int k = currRowIndex.getAndIncrement();
            while (k < numRounds) {
                ShuffledIndices allRowIndices = getShuffledByClusterIndices(clusterToRowIndices, isBaseline, generator);
                ShuffledIndices allColIndices;
                if (isIntra) {
                    allColIndices = allRowIndices;
                } else {
                    allColIndices = getShuffledByClusterIndices(clusterToColIndices, isBaseline, generator);
                }
                if (k == 0) {
                    globalAllIndices[0] = getShuffledByClusterIndices(clusterToRowIndices, isBaseline, generator);
                    globalAllIndices[1] = getShuffledByClusterIndices(clusterToColIndices, isBaseline, generator);
                }

                float[][] matrix = getShuffledMatrix(interMatrix, allRowIndices.allIndices, allColIndices.allIndices);

                updateMatrixScores(scoresBaselineForRound, k, matrix, allRowIndices.boundaries, allColIndices.boundaries, true);
                updateMatrixScores(scoresForRound, k, matrix, allRowIndices.boundaries, allColIndices.boundaries, false);

                FloatMatrixTools.log(matrix, 1);
                updateMatrixScores(logScoresBaselineForRound, k, matrix, allRowIndices.boundaries, allColIndices.boundaries, true);
                updateMatrixScores(logScoresForRound, k, matrix, allRowIndices.boundaries, allColIndices.boundaries, false);

                aggregate.addBToA(matrix);
                k = currRowIndex.getAndIncrement();
            }
        });

        scoreContainer.calculateMeans(scoresForRound, logScoresForRound,
                scoresBaselineForRound, logScoresBaselineForRound, mapIndex);

        aggregate.scaleForNumberOfRounds(numRounds);
        aggregate.saveToPNG(outfolder, name);

        scoreContainer.updateAggregateScores(aggregate, globalAllIndices, mapIndex);
    }

    private ShuffledIndices getShuffledByClusterIndices(Map<Integer, List<Integer>> clusterToIndices,
                                                        boolean isBaseline, Random generator) {
        List<Integer> allIndices = new ArrayList<>();

        List<Integer> order = new ArrayList<>(clusterToIndices.keySet());
        Collections.sort(order);

        int count = 0;
        List<Integer> boundaries = new ArrayList<>();
        boundaries.add(count);

        for (Integer index : order) {
            List<Integer> indexList = clusterToIndices.get(index);
            Collections.shuffle(indexList, generator);
            int numToUse = (indexList.size() / compressionFactor) * compressionFactor;
            for (int z = 0; z < numToUse; z++) {
                allIndices.add(indexList.get(z));
            }
            count += (numToUse / compressionFactor);
            boundaries.add(count);
        }
        if (isBaseline) {
            Collections.shuffle(allIndices, generator);
            return new ShuffledIndices(allIndices, new Integer[]{0, count});
        }
        Integer[] output = new Integer[boundaries.size()];
        return new ShuffledIndices(allIndices, boundaries.toArray(output));
    }

    private float[][] getShuffledMatrix(HiCMatrix interMatrix, List<Integer> allRowIndices, List<Integer> allColIndices) {
        int numRows = allRowIndices.size() / compressionFactor;
        int numCols = allColIndices.size() / compressionFactor;
        int numRowsKept = numRows * compressionFactor;
        int numColsKept = numCols * compressionFactor;
        float[][] original = interMatrix.getMatrix();

        float[][] result = new float[numRows][numCols];
        for (int i = 0; i < numRowsKept; i++) {
            final int i0 = allRowIndices.get(i);
            for (int j = 0; j < numColsKept; j++) {
                final int j0 = allColIndices.get(j);
                result[i / compressionFactor][j / compressionFactor] += original[i0][j0];
            }
        }

        return result;
    }

    private Map<Integer, List<Integer>> populateCluster(Chromosome[] chromosomes, int[] offsets, GenomeWide1DList<SubcompartmentInterval> subcompartments) {
        Map<Integer, List<Integer>> clusterToIndices = new HashMap<>();

        for (int x = 0; x < chromosomes.length; x++) {
            Chromosome chrom = chromosomes[x];
            List<SubcompartmentInterval> intervalList = subcompartments.getFeatures("" + chrom.getIndex());
            for (SubcompartmentInterval interval : intervalList) {
                int xStart = interval.getX1() / resolution;
                int xEnd = interval.getX2() / resolution;
                int clusterID = interval.getClusterID();

                List<Integer> tempList = new ArrayList<>();
                for (int k = xStart; k < xEnd; k++) {
                    final int actualPosition = k + offsets[x];
                    tempList.add(actualPosition);
                }

                if (clusterToIndices.containsKey(clusterID)) {
                    clusterToIndices.get(clusterID).addAll(tempList);
                } else {
                    clusterToIndices.put(clusterID, tempList);
                }
            }
        }
        return clusterToIndices;
    }
}
