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
import javastraw.reader.basics.Chromosome;
import mixer.utils.slice.structures.SubcompartmentInterval;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

public class CHICTools {
    public static Map<Integer, List<Integer>> populateCluster(Chromosome[] chromosomes, int[] offsets,
                                                              GenomeWide1DList<SubcompartmentInterval> subcompartments,
                                                              int resolution) {
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
