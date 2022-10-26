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
import javastraw.reader.block.ContactRecord;
import javastraw.reader.mzd.MatrixZoomData;
import javastraw.reader.type.NormalizationType;
import javastraw.tools.HiCFileTools;
import javastraw.tools.MatrixTools;
import javastraw.tools.ParallelizationTools;
import mixer.utils.common.FloatMatrixTools;
import mixer.utils.tracks.SubcompartmentInterval;

import java.io.File;
import java.io.FileWriter;
import java.util.*;
import java.util.concurrent.atomic.AtomicInteger;

public class GenomeWideStatistics {
    private final Dataset ds;
    private final int resolution;
    private final NormalizationType norm;
    private final Chromosome[] chromosomes;
    private final Map<Integer, Integer> clusterToFIdxMap;
    private final Map<Integer, Map<Integer, Integer>> chromToIndexToID;
    private final GenomeWide1DList<SubcompartmentInterval> subcompartments;
    private final double[][][] density;
    private final double klScoreBaseline, klScoreShuffle, klScoreShuffleSymm;
    private final double varScoreBaseline, varScoreShuffle, varScoreShuffleSymm;
    private final int n;
    private final BinSplitter splitter;

    public GenomeWideStatistics(Dataset ds, int resolution, NormalizationType norm,
                                GenomeWide1DList<SubcompartmentInterval> subcompartments, int numSplits) {
        this.ds = ds;
        this.resolution = resolution;
        this.norm = norm;
        this.subcompartments = subcompartments;
        this.splitter = new BinSplitter(numSplits);
        chromosomes = ds.getChromosomeHandler().getAutosomalChromosomesArray();
        clusterToFIdxMap = makeClusterToFIdxMap(subcompartments);
        chromToIndexToID = makeChromToIndexToIDMap();
        n = clusterToFIdxMap.keySet().size();
        double[][][] total = new double[splitter.getNumGroups()][n][n];
        long[][][] areasGW = new long[splitter.getNumGroups()][n][n];
        populateStatistics(total, areasGW);
        density = TensorTools.divide(total, areasGW);
        //inPlaceLog(density);
        varScoreBaseline = Scores.getVarScore(areasGW, density, true, false);
        klScoreBaseline = Scores.getKLScore(areasGW, density, true, true, false);

        varScoreShuffle = Scores.getVarScore(areasGW, density, false, false);
        klScoreShuffle = Scores.getKLScore(areasGW, density, true, false, false);

        varScoreShuffleSymm = Scores.getVarScore(areasGW, density, false, true);
        klScoreShuffleSymm = Scores.getKLScore(areasGW, density, true, false, true);
    }

    private void inPlaceLog(double[][][] a) {
        for (int i = 0; i < a.length; i++) {
            for (int j = 0; j < a[i].length; j++) {
                for (int k = 0; k < a[i][j].length; k++) {
                    a[i][j][k] = Math.log(1 + a[i][j][k]);
                }
            }
        }
    }

    private void populateStatistics(double[][][] allTotals, long[][][] allAreas) {

        AtomicInteger index = new AtomicInteger(0);
        ParallelizationTools.launchParallelizedCode(() -> {
            double[][][] totals = new double[splitter.getNumGroups()][n][n];
            long[][][] areas = new long[splitter.getNumGroups()][n][n];

            int i = index.getAndIncrement();
            while (i < chromosomes.length) {
                Chromosome chrom1 = chromosomes[i];
                Map<Integer, Integer> binToID1 = chromToIndexToID.get(chrom1.getIndex());

                for (int j = i + 1; j < chromosomes.length; j++) {
                    Chromosome chrom2 = chromosomes[j];
                    Map<Integer, Integer> binToID2 = chromToIndexToID.get(chrom2.getIndex());

                    final MatrixZoomData zd = HiCFileTools.getMatrixZoomData(ds, chrom1, chrom2, resolution);
                    if (zd != null) {
                        populateCounts(zd.getNormalizedIterator(norm), binToID1, binToID2, totals);
                        populateAreas(binToID1, binToID2, areas);
                        System.out.print(".");
                    }
                }
                System.out.println(".");
                i = index.getAndIncrement();
            }

            synchronized (allAreas) {
                TensorTools.addBtoA(allAreas, areas);
                TensorTools.addBtoA(allTotals, totals);
            }
        });
    }

    private Map<Integer, Map<Integer, Integer>> makeChromToIndexToIDMap() {
        Map<Integer, Map<Integer, Integer>> chromToIndexToID = new HashMap<>(chromosomes.length);
        for (Chromosome chromosome : chromosomes) {
            int numEntries = (int) ((chromosome.getLength() / resolution) + 1);
            Map<Integer, Integer> indexToID = new HashMap<>(numEntries);
            List<SubcompartmentInterval> intervals = subcompartments.getFeatures("" + chromosome.getIndex());
            for (SubcompartmentInterval interval : intervals) {
                int realID = clusterToFIdxMap.get(interval.getClusterID());
                int xStart = interval.getX1() / resolution;
                int xEnd = interval.getX2() / resolution;
                for (int r = xStart; r < xEnd; r++) {
                    indexToID.put(r, realID);
                }
            }
            chromToIndexToID.put(chromosome.getIndex(), indexToID);
        }
        return chromToIndexToID;
    }

    private Map<Integer, Integer> makeClusterToFIdxMap(GenomeWide1DList<SubcompartmentInterval> subcompartments) {
        Set<Integer> clusterIDs = new HashSet<>();
        for (Chromosome chromosome : chromosomes) {
            for (SubcompartmentInterval interval : subcompartments.getFeatures("" + chromosome.getIndex())) {
                clusterIDs.add(interval.getClusterID());
            }
        }

        List<Integer> ids = new ArrayList<>(clusterIDs);
        Collections.sort(ids);
        Map<Integer, Integer> clusterIDToFIdx = new HashMap<>(ids.size());
        for (int i = 0; i < ids.size(); i++) {
            clusterIDToFIdx.put(ids.get(i), i);
        }
        return clusterIDToFIdx;
    }

    private void populateCounts(Iterator<ContactRecord> iterator, Map<Integer, Integer> binToID1,
                                Map<Integer, Integer> binToID2, double[][][] counts) {
        while (iterator.hasNext()) {
            ContactRecord record = iterator.next();
            if (record.getCounts() > 0
                    && binToID1.containsKey(record.getBinX())
                    && binToID2.containsKey(record.getBinY())) {
                int id1 = binToID1.get(record.getBinX());
                int id2 = binToID2.get(record.getBinY());
                int sectionID = splitter.getSectionID(record.getBinX(), record.getBinY());
                counts[sectionID][id1][id2] += record.getCounts(); //todo explore Math.log(1+record.getCounts());
            }
        }
    }

    private void populateAreas(Map<Integer, Integer> binToID1,
                               Map<Integer, Integer> binToID2, long[][][] areas) {
        for (int bin1 : binToID1.keySet()) {
            int id1 = binToID1.get(bin1);
            for (int bin2 : binToID2.keySet()) {
                int id2 = binToID2.get(bin2);
                int sectionID = splitter.getSectionID(bin1, bin2);
                areas[sectionID][id1][id2]++;
            }
        }
    }

    public void writeToFile(File outfolder, String filename) {
        try {
            FileWriter myWriter = new FileWriter(new File(outfolder, filename + "_stats.txt"));
            myWriter.write("Variance Score: " + printDivision(varScoreBaseline, varScoreShuffle) + "\n");
            myWriter.write("Variance (Symm) Score: " + printDivision(varScoreBaseline, varScoreShuffleSymm) + "\n");
            myWriter.write("KL Divergence Score: " + printDivision(klScoreBaseline, klScoreShuffle) + "\n");
            myWriter.write("KL Divergence (Symm) Score: " + printDivision(klScoreBaseline, klScoreShuffleSymm) + "\n");
            myWriter.close();
        } catch (Exception e) {
            e.printStackTrace();
        }

        float[][] flattenedDensity = splitter.flatten(density);
        MatrixTools.saveMatrixTextNumpy(new File(outfolder, filename + "_density.npy").getAbsolutePath(),
                flattenedDensity);
        FloatMatrixTools.saveMatrixToPNG(new File(outfolder, filename + "_density.png"), flattenedDensity, false);
    }

    private String printDivision(double a, double b) {
        return roundI(a / b) + " (" + a + ")/(" + b + ")";
    }

    private int roundI(double v) {
        return (int) v;
    }
}
