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

package mixer.utils.slice;

import javastraw.feature1D.GenomeWideList;
import javastraw.reader.basics.ChromosomeHandler;
import mixer.utils.slice.kmeans.KmeansEvaluator;
import mixer.utils.slice.matrices.SliceMatrix;
import mixer.utils.slice.structures.ENCODESubcompartmentInterval;
import mixer.utils.slice.structures.SliceUtils;
import mixer.utils.slice.structures.SubcompartmentInterval;

import java.io.File;
import java.util.*;

public class EncodeExportUtils {
    public static void exportSubcompartments(SliceMatrix sliceMatrix, Map<Integer, List<List<Integer>>> kmeansIndicesMap,
                                             KmeansEvaluator evaluator, String prefix, File outputDirectory,
                                             ChromosomeHandler chromosomeHandler) {

        List<Integer> keys = new ArrayList<>(kmeansIndicesMap.keySet());
        Collections.sort(keys);

        int[] clusterSizes = convertToIntArray(keys);

        int bestIndex = getBestIndex(keys, evaluator);

        int[][] ids = new int[sliceMatrix.getData(false).length][keys.size()];
        for (int[] row : ids) {
            Arrays.fill(row, -1);
        }

        for (int index = 0; index < keys.size(); index++) {
            List<List<Integer>> colorList = kmeansIndicesMap.get(keys.get(index));
            for (int i = 0; i < colorList.size(); i++) {
                for (Integer val : colorList.get(i)) {
                    ids[val][index] = i;
                }
            }
        }

        Map<Integer, SubcompartmentInterval> map = sliceMatrix.getRowIndexToIntervalMap();
        GenomeWideList<SubcompartmentInterval> subcompartments = new GenomeWideList<>(chromosomeHandler);
        Set<SubcompartmentInterval> encodeSubcompartmentIntervals = new HashSet<>();
        for (int i = 0; i < ids.length; i++) {
            if (map.containsKey(i)) {
                SubcompartmentInterval interv = map.get(i);
                if (interv != null) {
                    encodeSubcompartmentIntervals.add(
                            generateENCODESubcompartment(interv, clusterSizes, ids[i], bestIndex));
                }
            }
        }
        subcompartments.addAll(new ArrayList<>(encodeSubcompartmentIntervals));
        SliceUtils.reSort(subcompartments);
        File outBedFile = new File(outputDirectory, prefix + "_subcompartment_clusters.bed");
        subcompartments.simpleExport(outBedFile);
    }

    private static ENCODESubcompartmentInterval generateENCODESubcompartment(SubcompartmentInterval interv,
                                                                             int[] clusterSizes, int[] id,
                                                                             int indx) {
        return new ENCODESubcompartmentInterval(interv.getChrIndex(), interv.getChrName(),
                interv.getX1(), interv.getX2(), id[indx], clusterSizes, id);
    }

    private static int[] convertToIntArray(List<Integer> keys) {
        int[] arr = new int[keys.size()];
        for (int i = 0; i < arr.length; i++) {
            arr[i] = keys.get(i);
        }
        return arr;
    }

    private static int getBestIndex(List<Integer> keys, KmeansEvaluator evaluator) {
        int bestIndex = 0;
        double bestScore = evaluator.getSilhouette(keys.get(bestIndex));
        for (int k = 0; k < keys.size(); k++) {
            Integer key = keys.get(k);
            if (bestScore < evaluator.getSilhouette(key)) {
                bestScore = evaluator.getSilhouette(key);
                bestIndex = k;
            }
        }
        return bestIndex;
    }
}
