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

package mixer.utils.slice;

import javastraw.featurelist.GenomeWideList;
import javastraw.reader.ChromosomeHandler;
import javastraw.reader.Dataset;
import javastraw.reader.basics.Chromosome;
import javastraw.type.NormalizationType;
import mixer.MixerGlobals;
import mixer.utils.common.FloatMatrixTools;
import mixer.utils.common.Pair;
import mixer.utils.slice.kmeansfloat.Cluster;
import mixer.utils.slice.kmeansfloat.ClusterTools;

import java.io.File;
import java.util.*;

public abstract class CompositeGenomeWideDensityMatrix {
    protected final NormalizationType norm;
    protected final int resolution;
    protected final Map<Integer, SubcompartmentInterval> rowIndexToIntervalMap = new HashMap<>();
    protected final Chromosome[] chromosomes;
    protected final Random generator;
    protected final File outputDirectory;
    private final float[][] gwCleanMatrix;
    private final List<Map<Integer, Map<Integer, Integer>>> chrIndxTorowIndexToGoldIDMapList = new ArrayList<>();
    protected GenomewideBadIndexFinder badIndexLocations;

    public CompositeGenomeWideDensityMatrix(ChromosomeHandler chromosomeHandler, Dataset ds, NormalizationType norm, int resolution,
                                            File outputDirectory, Random generator, String[] relativeTestFiles) {
        this.norm = norm;
        this.resolution = resolution;
        this.outputDirectory = outputDirectory;
        this.generator = generator;

        chrIndxTorowIndexToGoldIDMapList.clear();
        if (relativeTestFiles != null) {
            for (String filename : relativeTestFiles) {
                chrIndxTorowIndexToGoldIDMapList.add(SliceUtils.createGoldStandardLookup(filename, resolution, chromosomeHandler));
            }
        }

        chromosomes = chromosomeHandler.getAutosomalChromosomesArray();
        badIndexLocations = new GenomewideBadIndexFinder(ds, chromosomes, resolution, norm);
        gwCleanMatrix = makeCleanScaledInterMatrix(ds);
    }

    abstract float[][] makeCleanScaledInterMatrix(Dataset ds);

    public synchronized Pair<Double, List<int[][]>> processGWKmeansResult(Cluster[] clusters, GenomeWideList<SubcompartmentInterval> subcompartments) {

        Set<SubcompartmentInterval> subcompartmentIntervals = new HashSet<>();
        if (MixerGlobals.printVerboseComments) {
            System.out.println("GW Composite data vs clustered into " + clusters.length + " clusters");
        }

        double withinClusterSumOfSquares = 0;
        int genomewideCompartmentID = 0;

        int[][] ids = new int[1][clusters.length];
        int[][] idsForIndex = new int[chrIndxTorowIndexToGoldIDMapList.size() + 1][gwCleanMatrix.length];
        for (int q = 0; q < idsForIndex.length; q++) {
            Arrays.fill(idsForIndex[q], -1);
        }

        for (int z = 0; z < clusters.length; z++) {
            Cluster cluster = clusters[z];
            int currentClusterID = ++genomewideCompartmentID;
            ids[0][z] = currentClusterID;

            if (MixerGlobals.printVerboseComments) {
                System.out.println("Size of cluster " + currentClusterID + " - " + cluster.getMemberIndexes().length);
            }

            for (int i : cluster.getMemberIndexes()) {
                withinClusterSumOfSquares += ClusterTools.getNonNanVectorSumOfSquares(cluster.getCenter(), gwCleanMatrix[i]);

                try {
                    SubcompartmentInterval interv;

                    if (rowIndexToIntervalMap.containsKey(i)) {
                        interv = rowIndexToIntervalMap.get(i);
                        if (interv == null) continue; // probably a zero row

                        int chrIndex = interv.getChrIndex();
                        String chrName = interv.getChrName();
                        int x1 = interv.getX1();
                        int x2 = interv.getX2();

                        subcompartmentIntervals.add(
                                new SubcompartmentInterval(chrIndex, chrName, x1, x2, currentClusterID));

                        idsForIndex[chrIndxTorowIndexToGoldIDMapList.size()][i] = currentClusterID;

                        for (int q = 0; q < chrIndxTorowIndexToGoldIDMapList.size(); q++) {
                            Map<Integer, Map<Integer, Integer>> map = chrIndxTorowIndexToGoldIDMapList.get(q);
                            if (map.containsKey(chrIndex)) {
                                if (map.get(chrIndex).containsKey(x1)) {
                                    idsForIndex[q][i] = map.get(chrIndex).get(x1);
                                }
                            }
                        }

                    } else {
                        System.err.println("********* is weird error?");
                    }

                } catch (Exception e) {
                    e.printStackTrace();
                    System.exit(87);
                }
            }
        }

        if (MixerGlobals.printVerboseComments) {
            System.out.println("Final WCSS " + withinClusterSumOfSquares);
        }

        subcompartments.addAll(new ArrayList<>(subcompartmentIntervals));
        SliceUtils.reSort(subcompartments);

        List<int[][]> outputs = new ArrayList<>();
        outputs.add(ids);
        outputs.add(idsForIndex);

        return new Pair<>(withinClusterSumOfSquares, outputs);
    }

    protected Pair<Integer, int[][]> calculateDimensionInterMatrix(Chromosome[] chromosomes, Map<Integer, Integer> indexToFilteredLength) {
        int total = 0;
        int[][] indices = new int[2][chromosomes.length];

        for (int i = 0; i < chromosomes.length; i++) {
            int val = indexToFilteredLength.get(chromosomes[i].getIndex());
            total += val;
            if (i < chromosomes.length - 1) {
                indices[0][i + 1] = total;
            }
            indices[1][i] = val;
        }

        return new Pair<>(total, indices);
    }

    public float[][] getCleanedData() {
        return gwCleanMatrix;
    }

    public int getLength() {
        return gwCleanMatrix.length;
    }

    public int getWidth() {
        return gwCleanMatrix[0].length;
    }

    public void exportData() {
        System.out.println(getLength() + " -v- " + getWidth());
        FloatMatrixTools.saveMatrixTextNumpy(new File(outputDirectory, "data_matrix.npy").getAbsolutePath(), getCleanedData());
    }

    public void appendDataAlongExistingRows(CompositeGenomeWideDensityMatrix additionalData) {
        if (getLength() != additionalData.getLength()) {
            System.err.println("***************************************\n" +
                    "Dimension mismatch: " + getLength() + " != " + additionalData.getLength());
        } else {

        }
    }

    public GenomewideBadIndexFinder getBadIndices() {
        return badIndexLocations;
    }
}
