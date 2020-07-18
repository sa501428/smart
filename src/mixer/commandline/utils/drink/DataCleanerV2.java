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

import mixer.commandline.utils.drink.kmeansfloat.Cluster;
import mixer.commandline.utils.drink.obsolete.OldTools;
import mixer.data.feature.GenomeWideList;
import org.broad.igv.feature.Chromosome;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

public class DataCleanerV2 extends DataCleaner {

    private final List<Integer> dataSetSeparatingIndices = new ArrayList<>();
    private final int numDatasets;

    public DataCleanerV2(List<double[][]> matrices, double[][] data, int numDatasets, double maxPercentAllowedToBeZeroThreshold, int resolution, double[] convolution1d) {
        super(data, maxPercentAllowedToBeZeroThreshold, resolution, convolution1d);
        this.numDatasets = numDatasets;
        determineSeparatingIndices(matrices);
    }

    private int determineWhichDatasetThisBelongsTo(int originalRow) {
        for (int i = 0; i < numDatasets - 1; i++) {
            if (originalRow < dataSetSeparatingIndices.get(i + 1)) return i;
        }
        return numDatasets - 1;
    }

    private void determineSeparatingIndices(List<double[][]> data) {
        int rowOffSet = 0;
        for (int i = 0; i < numDatasets; i++) {
            dataSetSeparatingIndices.add(rowOffSet);
            rowOffSet += data.get(i).length;
        }
    }

    public synchronized List<Map<Integer, List<Integer>>> postProcessKmeansResultV2(Cluster[] clusters,
                                                                                    Map<Integer, float[]> idToCentroidMap) {

        List<Map<Integer, List<Integer>>> mapPosIndexToCluster = new ArrayList<>();
        for (int i = 0; i < numDatasets; i++) {
            mapPosIndexToCluster.add(new HashMap<>());
        }

        for (Cluster cluster : clusters) {
            int currentClusterID = initialClusterID.getAndIncrement();
            synchronized (idToCentroidMap) {
                idToCentroidMap.put(currentClusterID, cluster.getCenter());
            }

            for (int i : cluster.getMemberIndexes()) {
                int originalRow = getOriginalIndexRow(i);
                int datasetIndx = determineWhichDatasetThisBelongsTo(originalRow);
                int xIndex = originalRow - dataSetSeparatingIndices.get(datasetIndx);
                List<Integer> newList = new ArrayList<>();
                newList.add(currentClusterID);
                mapPosIndexToCluster.get(datasetIndx).put(xIndex, newList);
            }
        }

        return mapPosIndexToCluster;
    }

    public void postProcessSplitting(Chromosome chromosome, List<int[]> memberIndices,
                                     List<GenomeWideList<SubcompartmentInterval>> subcompartmentIntervals) {

        List<Map<Integer, Integer>> finalMemberIndices = verifyFinalMemberIndices(memberIndices, numDatasets);

        for (int mIdx = 0; mIdx < numDatasets; mIdx++) {
            Map<Integer, Integer> startAndEndMap = finalMemberIndices.get(mIdx);
            List<SubcompartmentInterval> intervals = new ArrayList<>();
            for (Integer start : startAndEndMap.keySet()) {
                Integer end = startAndEndMap.get(start);
                int metaID = initialClusterID.getAndIncrement();
                intervals.add(new SubcompartmentInterval(chromosome.getIndex(),
                        chromosome.getName(), start * getResolution(), end * getResolution(), metaID));
            }
            subcompartmentIntervals.get(mIdx).addAll(intervals);
        }
    }

    private List<Map<Integer, Integer>> verifyFinalMemberIndices(List<int[]> memberIndices, int numDatasets) {

        List<Map<Integer, Integer>> regionsForDataset = new ArrayList<>();
        for (int k = 0; k < numDatasets; k++) {
            regionsForDataset.add(new HashMap<Integer, Integer>());
        }

        for (int[] indices : memberIndices) {
            int xIndex1 = getOriginalIndexRow(indices[0]);
            int datasetIndex1 = determineWhichDatasetThisBelongsTo(xIndex1);
            xIndex1 = xIndex1 - dataSetSeparatingIndices.get(datasetIndex1);

            //int currIdx = xIndex;
            //int currDataset = datasetIndx;
            int datasetIndex2, xIndex2;
            int currXIdx = xIndex1;

            for (int k = 1; k < indices.length; k++) {
                xIndex2 = getOriginalIndexRow(indices[k]);
                datasetIndex2 = determineWhichDatasetThisBelongsTo(xIndex2);
                xIndex2 = xIndex2 - dataSetSeparatingIndices.get(datasetIndex2);

                if (datasetIndex1 == datasetIndex2 && currXIdx + 1 == xIndex2) {
                    // continuous
                    currXIdx = xIndex2;
                } else {
                    // break - shouldn't happen often
                    regionsForDataset.get(datasetIndex1).put(xIndex1, currXIdx + 1);

                    datasetIndex1 = datasetIndex2;
                    xIndex1 = xIndex2;
                    currXIdx = xIndex2;
                }
            }
            regionsForDataset.get(datasetIndex1).put(xIndex1, currXIdx + 1);
        }

        return regionsForDataset;
    }

    public float[][] getMatrixModIndicesOfColumns(int modIndx, int base) {
        return OldTools.getMatrixModIndicesOfColumns(getCleanedData(), modIndx, base);
    }

    public float[][] getHalfOfMatrix(boolean getFirstHalf) {
        return OldTools.getHalfOfMatrix(getCleanedData(), getFirstHalf);
    }
}
