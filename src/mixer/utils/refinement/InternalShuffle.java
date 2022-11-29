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

package mixer.utils.refinement;

import javastraw.feature1D.GenomeWide1DList;
import javastraw.reader.Dataset;
import javastraw.reader.basics.ChromosomeHandler;
import javastraw.reader.type.NormalizationType;
import mixer.utils.BedTools;
import mixer.utils.shuffle.Partition;
import mixer.utils.shuffle.ShuffleAction;
import mixer.utils.tracks.SliceUtils;
import mixer.utils.tracks.SubcompartmentInterval;

import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.Random;

public class InternalShuffle {
    public static Map<Integer, GenomeWide1DList<SubcompartmentInterval>> determineBest(Map<Integer, List<String>> allBedFiles,
                                                                                       int resolution,
                                                                                       ChromosomeHandler handler,
                                                                                       Dataset ds,
                                                                                       NormalizationType norm) {

        int compressionFactor = 1600000 / resolution;
        Random generator = new Random(0);

        Partition.Type[] mapTypes = {Partition.Type.ODDS_VS_EVENS};

        Map<Integer, GenomeWide1DList<SubcompartmentInterval>> bestClusterings = new HashMap<>();

        for (int k : allBedFiles.keySet()) {
            List<String> bedFiles = allBedFiles.get(k);
            double bestScore = 0;
            GenomeWide1DList<SubcompartmentInterval> bestClustering = null;
            for (String bedFile : bedFiles) {
                GenomeWide1DList<SubcompartmentInterval> subcompartments =
                        BedTools.loadBedFileAtResolution(handler, bedFile, resolution);
                ShuffleAction matrix = new ShuffleAction(ds, norm, resolution, compressionFactor, mapTypes);
                matrix.runInterAnalysis(subcompartments, null, generator);
                double score = matrix.getResult(0);
                //System.out.println(score);
                if (score > bestScore) {
                    bestScore = score;
                    bestClustering = subcompartments;
                }
            }
            if (bestClustering != null) {
                SliceUtils.collapseGWList(bestClustering);
                bestClusterings.put(k, bestClustering);
            }
        }

        return bestClusterings;
    }

    public static Map<Integer, GenomeWide1DList<SubcompartmentInterval>> getDefault(Map<Integer, List<String>> allBedFiles,
                                                                                    int resolution,
                                                                                    ChromosomeHandler handler) {

        Map<Integer, GenomeWide1DList<SubcompartmentInterval>> bestClusterings = new HashMap<>();
        for (int k : allBedFiles.keySet()) {
            List<String> bedFiles = allBedFiles.get(k);
            if (bedFiles.size() > 1) {
                System.err.println("too many bed files created - internal warning");
            }
            GenomeWide1DList<SubcompartmentInterval> bestClustering = BedTools.loadBedFileAtResolution(handler,
                    bedFiles.get(0), resolution);
            SliceUtils.collapseGWList(bestClustering);
            bestClusterings.put(k, bestClustering);
        }
        return bestClusterings;
    }
}
