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

package mixer.utils.drive;

import javastraw.feature1D.GenomeWide1DList;
import javastraw.reader.basics.Chromosome;
import mixer.utils.tracks.Concensus2DTools;
import mixer.utils.tracks.SubcompartmentInterval;

import java.util.Arrays;
import java.util.Map;

public class BedFileMappings extends BinMappings {


    public BedFileMappings(Chromosome[] chromosomes, int resolution,
                           GenomeWide1DList<SubcompartmentInterval> clusters) {
        super(resolution, chromosomes);

        Map<String, Map<Integer, Integer>> map1 = Concensus2DTools.summarize(clusters);
        int numClusters = Concensus2DTools.getMaxId(map1) + 1;

        int counter = 0;
        for (Chromosome chrom : chromosomes) {
            int length = (int) (chrom.getLength() / resolution) + 1;
            int[] binToProtocluster = new int[length];
            Arrays.fill(binToProtocluster, IGNORE);
            fillInClusterAssigments(binToProtocluster, map1.get("" + chrom.getIndex()), counter, resolution);
            putBinToProtoCluster(chrom, binToProtocluster);
            counter += numClusters;
        }

        calculateGlobalIndices(chromosomes);
    }

    private void fillInClusterAssigments(int[] binToProtocluster, Map<Integer, Integer> posToID,
                                         int counter, int resolution) {
        for (Integer pos : posToID.keySet()) {
            int id = posToID.get(pos);
            binToProtocluster[pos / resolution] = id + counter;
        }
    }
}
