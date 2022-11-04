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

package mixer.utils.tracks;

import javastraw.feature1D.GenomeWide1DList;

import java.util.HashMap;
import java.util.Map;

public class Concensus2DTools {

    public static void checkOverlapPerChrom(GenomeWide1DList<SubcompartmentInterval> file1, GenomeWide1DList<SubcompartmentInterval> file2) {
        Map<String, Map<Integer, Integer>> map1 = summarize(file1);
        Map<String, Map<Integer, Integer>> map2 = summarize(file2);
        for (String key : map1.keySet()) {
            int[][] summary = populateSummary(map1, map2, key);
            double totalDenom = sum(summary);
            double totalNum = 0;

            int minNum = Math.min(summary.length, summary[0].length);
            for (int q = 0; q < minNum; q++) {
                int[] coords = getMaxCoordinates(summary);
                if (coords != null) {
                    totalNum += summary[coords[0]][coords[1]];
                    clearSectionSlices(summary, coords);
                } else {
                    break;
                }
            }
            System.out.println("Accuracy for " + key + " " + (100 * totalNum / totalDenom));
        }
    }

    public static void checkOverlap(GenomeWide1DList<SubcompartmentInterval> file1,
                                    GenomeWide1DList<SubcompartmentInterval> file2) {

        int[][] summary = populateSummary(file1, file2);
        double totalDenom = sum(summary);
        double totalNum = 0;

        int minNum = Math.min(summary.length, summary[0].length);
        for (int q = 0; q < minNum; q++) {
            int[] coords = getMaxCoordinates(summary);
            if (coords != null) {
                totalNum += summary[coords[0]][coords[1]];
                clearSectionSlices(summary, coords);
            } else {
                break;
            }
        }

        System.out.println("Accuracy " + (100 * totalNum / totalDenom));
    }

    private static int sum(int[][] matrix) {
        int total = 0;
        for (int[] row : matrix) {
            for (int val : row) {
                total += val;
            }
        }
        return total;
    }

    private static int[][] populateSummary(GenomeWide1DList<SubcompartmentInterval> file1, GenomeWide1DList<SubcompartmentInterval> file2) {

        Map<String, Map<Integer, Integer>> map1 = summarize(file1);
        Map<String, Map<Integer, Integer>> map2 = summarize(file2);

        int n = Math.max(getMaxId(map1), getMaxId(map2)) + 1;
        int[][] counts = new int[n][n];
        for (String key : map1.keySet()) {
            populateContingencyTable(map1, map2, key, counts);
        }

        return counts;
    }

    private static int[][] populateSummary(Map<String, Map<Integer, Integer>> map1,
                                           Map<String, Map<Integer, Integer>> map2,
                                           String key) {
        int n = Math.max(getMaxId(map1), getMaxId(map2)) + 1;
        int[][] counts = new int[n][n];
        populateContingencyTable(map1, map2, key, counts);
        return counts;
    }

    private static void populateContingencyTable(Map<String, Map<Integer, Integer>> map1, Map<String, Map<Integer, Integer>> map2, String key, int[][] counts) {
        Map<Integer, Integer> mapping1 = map1.get(key);
        Map<Integer, Integer> mapping2 = map2.get(key);
        for (Integer pos : mapping1.keySet()) {
            if (mapping2.containsKey(pos)) {
                if (mapping1.get(pos) >= 0 && mapping2.get(pos) >= 0)
                    counts[mapping1.get(pos)][mapping2.get(pos)]++;
            }
        }
    }

    private static int getMaxId(Map<String, Map<Integer, Integer>> map) {
        int maxVal = 0;
        for (Map<Integer, Integer> mapping : map.values()) {
            for (Integer val : mapping.values()) {
                maxVal = Math.max(maxVal, val);
            }
        }
        return maxVal;
    }

    private static Map<String, Map<Integer, Integer>> summarize(GenomeWide1DList<SubcompartmentInterval> file1) {
        Map<String, Map<Integer, Integer>> summary = new HashMap<>();
        file1.processLists((s, list) -> {
            Map<Integer, Integer> indexToId = new HashMap<>();
            for (SubcompartmentInterval interval : list) {
                indexToId.put(interval.getX1(), interval.getClusterID());
            }
            summary.put(s, indexToId);
        });
        return summary;
    }

    private static void clearSectionSlices(int[][] tensor, int[] coords) {
        for (int j = 0; j < tensor[0].length; j++) {
            tensor[coords[0]][j] = -1;
        }
        for (int i = 0; i < tensor.length; i++) {
            tensor[i][coords[1]] = -1;
        }
    }

    public static int[] getMaxCoordinates(int[][] tensor) {
        int maxVal = 0;
        int[] coords = null;
        for (int i = 0; i < tensor.length; i++) {
            for (int j = 0; j < tensor.length; j++) {
                if (tensor[i][j] > maxVal) {
                    maxVal = tensor[i][j];
                    coords = new int[]{i, j};
                }
            }
        }
        return coords;
    }

    private static String makeKey(int[] coords) {
        return coords[0] + "_" + coords[1] + "_" + coords[2];
    }

    private static String makeKey(int i, int j, int k) {
        return i + "_" + j + "_" + k;
    }

}
