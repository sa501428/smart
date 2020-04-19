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

package mixer.commandline.utils.drink;

import mixer.MixerGlobals;
import mixer.commandline.utils.common.DoubleMatrixTools;
import mixer.commandline.utils.drink.kmeansfloat.ClusterTools;
import mixer.data.ChromosomeHandler;
import mixer.data.Dataset;
import mixer.data.HiCFileTools;
import mixer.data.feature.GenomeWideList;
import mixer.windowui.NormalizationType;
import org.apache.commons.math.linear.RealMatrix;
import org.broad.igv.feature.Chromosome;

import java.io.File;
import java.io.IOException;
import java.util.*;
import java.util.stream.IntStream;

import static java.util.Comparator.comparing;

public class DistanceSplitter {

    private final double maxPercentAllowBeZero = 0.75;
    private final List<Dataset> datasets;
    private final int numDatasets;
    private final ChromosomeHandler chromosomeHandler;
    private final int resolution;
    private final NormalizationType norm;
    private final List<GenomeWideList<SubcompartmentInterval>> comparativeSubcompartments = new ArrayList<>();
    private final boolean useStackingAlongRow;
    private final float threshold;
    int numNeighbors = 3;

    public DistanceSplitter(List<Dataset> datasetList, ChromosomeHandler chromosomeHandler,
                            int resolution, NormalizationType norm, float oeThreshold, boolean useStackingAlongRow) {
        this.datasets = datasetList;
        this.chromosomeHandler = chromosomeHandler;
        this.resolution = resolution;
        this.norm = norm;
        this.useStackingAlongRow = useStackingAlongRow;
        this.threshold = oeThreshold;

        if (useStackingAlongRow) {
            numDatasets = 1;
            comparativeSubcompartments.add(new GenomeWideList<>(chromosomeHandler));
        } else {
            numDatasets = datasets.size();
            for (int i = 0; i < numDatasets; i++) {
                comparativeSubcompartments.add(new GenomeWideList<>(chromosomeHandler));
            }
        }
    }

    public static int[] smallestNIndices(final double[] input, final int n) {
        return IntStream.range(0, input.length)
                .boxed()
                .sorted(comparing(i -> input[i]))
                .mapToInt(i -> i)
                .limit(n)
                .toArray();
    }

    public List<GenomeWideList<SubcompartmentInterval>> extractIntraSubcompartmentsTo(File outdir, double[] convolution) {

        for (final Chromosome chromosome : chromosomeHandler.getAutosomalChromosomesArray()) {
            try {
                List<double[][]> matrices = new ArrayList<>();

                for (Dataset ds : datasets) {
                    RealMatrix localizedRegionData = HiCFileTools.getRealOEMatrixForChromosome(ds, chromosome, resolution,
                            norm, threshold,
                            //ExtractingOEDataUtils.ThresholdType.LOG_OE_BOUNDED,
                            //ExtractingOEDataUtils.ThresholdType.LOGEO, //
                            ExtractingOEDataUtils.ThresholdType.TRUE_OE,
                            true);
                    if (localizedRegionData != null) {
                        matrices.add(localizedRegionData.getData());
                        if (MixerGlobals.printVerboseComments) {
                            DoubleMatrixTools.saveMatrixTextNumpy(new File(outdir, chromosome.getName() + "_matrix.npy").getAbsolutePath(), localizedRegionData.getData());
                        }
                    }
                }

                // can't assess vs non existent map
                if (matrices.size() != datasets.size() || matrices.size() < 1) continue;

                double[][] collapsedMatrix;
                if (useStackingAlongRow) {
                    collapsedMatrix = DoubleMatrixTools.stitchMultipleMatricesTogetherByColDim(matrices);
                } else {
                    collapsedMatrix = DoubleMatrixTools.stitchMultipleMatricesTogetherByRowDim(matrices);
                }

                DataCleanerV2 dataCleaner = new DataCleanerV2(matrices, collapsedMatrix, numDatasets,
                        maxPercentAllowBeZero, resolution, convolution);

                List<int[]> memberIndices = getGroupedRegions(dataCleaner);
                dataCleaner.postProcessSplitting(chromosome, memberIndices, comparativeSubcompartments);

            } catch (IOException e) {
                e.printStackTrace();
            }
        }
        return comparativeSubcompartments;
    }

    private List<int[]> getGroupedRegions(DataCleanerV2 cleanedData) {

        double[][] distanceMatrix = getDistancesFromVectors(cleanedData.getCleanedData());

        int[][] minIndices = getMinClosestIndices(distanceMatrix);

        List<Integer> newStartPoints = getBreakPoints(minIndices);

        return groupContiguousBlocks(newStartPoints, distanceMatrix.length);
    }

    private List<int[]> groupContiguousBlocks(List<Integer> newStartPoints, int length) {
        List<int[]> regions = new ArrayList<>();

        int numStarts = newStartPoints.size();

        for (int indx = 0; indx < numStarts - 1; indx++) {
            int currIndx = newStartPoints.get(indx);
            int nextIndx = newStartPoints.get(indx + 1);
            int[] region = generateRegionBlock(currIndx, nextIndx);
            regions.add(region);
        }

        // last region
        int[] region = generateRegionBlock(newStartPoints.get(numStarts - 1), length);
        regions.add(region);

        return regions;
    }

    private int[] generateRegionBlock(int currIndx, int nextIndx) {
        int[] block = new int[nextIndx - currIndx];
        int counter = 0;
        for (int i = currIndx; i < nextIndx; i++) {
            block[counter] = i;
            counter++;
        }

        return block;
    }

    private List<Integer> getBreakPoints(int[][] minIndices) {

        Set<Integer> newStartPoints = new HashSet<>();
        newStartPoints.add(0); // first node
        for (int i = 1; i < minIndices.length; i++) {
            boolean connectedToPrev = indexIsContainedInOrNear(i, i - 1, minIndices)
                    || indexIsContainedInOrNear(i - 1, i, minIndices);
            if (!connectedToPrev) {
                newStartPoints.add(i);
            }
        }

        List<Integer> sortedPoints = new ArrayList<>(newStartPoints);
        Collections.sort(sortedPoints);

        return sortedPoints;
    }

    // secondary contact
    private boolean indexIsContainedInOrNear(int node, int neighbor, int[][] minIndices) {
        for (int node2 : minIndices[neighbor]) {
            if (node == node2) {
                return true;
            }
            for (int node3 : minIndices[node2]) {
                if (node == node3) {
                    return true;
                }
            }
        }
        return false;
    }

    private double[][] getDistancesFromVectors(float[][] data) {

        int numVectors = data.length;
        double[][] distances = new double[numVectors][numVectors];

        for (int i = 0; i < numVectors; i++) {
            Arrays.fill(distances[i], Double.MAX_VALUE);
        }

        for (int i = 0; i < numVectors; i++) {
            for (int j = i + 1; j < numVectors; j++) {

                double val = ClusterTools.getDistance(data[i], data[j]);
                //double val = ClusterTools.getL1Distance(data[i], data[j]);
                distances[i][j] = val;
                distances[j][i] = val;
            }
        }

        return distances;
    }

    private int[][] getMinClosestIndices(double[][] distanceMatrix) {
        int[][] minIndices = new int[distanceMatrix.length][numNeighbors];
        for (int i = 0; i < distanceMatrix.length; i++) {
            int[] vals = smallestNIndices(distanceMatrix[i], numNeighbors);
            if (vals.length == numNeighbors) {
                minIndices[i] = vals;
            } else {
                System.err.println("Error encountered in splitter");
            }
        }
        return minIndices;
    }


}
