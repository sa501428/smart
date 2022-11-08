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

package mixer.utils.cv;

import javastraw.expected.Welford;
import javastraw.feature1D.GenomeWide1DList;
import javastraw.reader.Dataset;
import javastraw.reader.basics.Chromosome;
import javastraw.reader.block.ContactRecord;
import javastraw.reader.mzd.MatrixZoomData;
import javastraw.reader.type.NormalizationType;
import javastraw.tools.HiCFileTools;
import javastraw.tools.MatrixTools;
import mixer.utils.tracks.SubcompartmentInterval;

import java.util.*;

public class GenomeWideStatistics {
    private final Dataset ds;
    private final int resolution;
    private final NormalizationType norm;
    private final Chromosome[] chromosomes;
    private final Map<Integer, Integer> clusterToFIdxMap;
    private final Map<Integer, Map<Integer, Integer>> chromToIndexToID;
    private final GenomeWide1DList<SubcompartmentInterval> subcompartments;
    private final int n;
    private final double[] cvs;

    public GenomeWideStatistics(Dataset ds, int resolution, NormalizationType norm,
                                GenomeWide1DList<SubcompartmentInterval> subcompartments, int numSplits) {
        this.ds = ds;
        this.resolution = resolution;
        this.norm = norm;
        this.subcompartments = subcompartments;
        chromosomes = ds.getChromosomeHandler().getAutosomalChromosomesArray();
        clusterToFIdxMap = makeClusterToFIdxMap(subcompartments);
        chromToIndexToID = makeChromToIndexToIDMap();
        n = clusterToFIdxMap.keySet().size();


        Welford[][] welfords = populateStatistics();

        Welford finalWelford = new Welford();
        List<Double> values = new ArrayList<>();
        for (int i = 0; i < n; i++) {
            for (int j = 0; j < n; j++) {
                if (welfords[i][j].getCounts() > 3) {
                    double cv = welfords[i][j].getStdDev() / welfords[i][j].getMean();
                    values.add(cv);
                    finalWelford.addValue(cv);
                }
            }
        }
        cvs = toArray(values);
        System.out.println("STD of CVs = " + finalWelford.getStdDev());
    }

    public static double[][] divide(double[][] a, long[][] b) {
        double[][] result = new double[a.length][a[0].length];
        for (int i = 0; i < a.length; i++) {
            for (int j = 0; j < a[i].length; j++) {
                if (b[i][j] > 0) {
                    result[i][j] = a[i][j] / b[i][j];
                }
            }
        }
        return result;
    }

    private double[] toArray(List<Double> values) {
        double[] answer = new double[values.size()];
        int counter = 0;
        for (Double val : values) {
            answer[counter++] = val;
        }
        return answer;
    }

    private Welford[][] initializeWelfords(int n) {
        Welford[][] result = new Welford[n][n];
        for (int i = 0; i < n; i++) {
            for (int j = 0; j < n; j++) {
                result[i][j] = new Welford();
            }
        }
        return result;
    }

    private Welford[][] populateStatistics() {

        Welford[][] welfords = initializeWelfords(n);

        for (int i = 0; i < chromosomes.length; i++) {
            Chromosome chrom1 = chromosomes[i];
            Map<Integer, Integer> binToID1 = chromToIndexToID.get(chrom1.getIndex());

            for (int j = i + 1; j < chromosomes.length; j++) {
                Chromosome chrom2 = chromosomes[j];
                Map<Integer, Integer> binToID2 = chromToIndexToID.get(chrom2.getIndex());

                final MatrixZoomData zd = HiCFileTools.getMatrixZoomData(ds, chrom1, chrom2, resolution);
                if (zd != null) {
                    double[][] totals = new double[n][n];
                    long[][] areas = new long[n][n];

                    populateCounts(zd.getNormalizedIterator(norm), binToID1, binToID2, totals);
                    populateAreas(binToID1, binToID2, areas);

                    double[][] density = divide(totals, areas);
                    update(welfords, density);

                    System.out.print(".");
                }
            }
        }
        return welfords;
    }

    private void update(Welford[][] welfords, double[][] density) {
        for (int i = 0; i < n; i++) {
            for (int j = 0; j < n; j++) {
                welfords[i][j].addValue(density[i][j]);
            }
        }
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
                                Map<Integer, Integer> binToID2, double[][] counts) {
        while (iterator.hasNext()) {
            ContactRecord record = iterator.next();
            if (record.getCounts() > 0
                    && binToID1.containsKey(record.getBinX())
                    && binToID2.containsKey(record.getBinY())) {
                int id1 = binToID1.get(record.getBinX());
                int id2 = binToID2.get(record.getBinY());
                counts[id1][id2] += (record.getCounts()); //todo explore Math.log(1+record.getCounts());
            }
        }
    }

    private void populateAreas(Map<Integer, Integer> binToID1,
                               Map<Integer, Integer> binToID2, long[][] areas) {
        for (int bin1 : binToID1.keySet()) {
            int id1 = binToID1.get(bin1);
            for (int bin2 : binToID2.keySet()) {
                int id2 = binToID2.get(bin2);
                areas[id1][id2]++;
            }
        }
    }

    public void export(String prefix) {
        MatrixTools.saveMatrixTextNumpy(prefix + "_coeff_variances.npy", cvs);
    }
}
