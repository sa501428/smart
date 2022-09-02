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
import javastraw.tools.MatrixTools;
import mixer.utils.cleaning.EmptyRowCleaner;
import mixer.utils.common.ZScoreTools;
import mixer.utils.tracks.SliceUtils;
import mixer.utils.tracks.SubcompartmentInterval;
import robust.concurrent.kmeans.clustering.Cluster;

import java.io.File;
import java.util.*;

public class MatrixAndWeight {
    public float[][] matrix;
    public int[] weights;
    private final Map<Integer, SubcompartmentInterval> map = new HashMap<>();
    private final Mappings mappings;


    public MatrixAndWeight(float[][] interMatrix, int[] weights, Mappings mappings) {
        this.matrix = interMatrix;
        this.weights = weights;
        this.mappings = mappings;
        if (mappings != null) populateRowIndexToIntervalMap(mappings);
    }

    public void inPlaceScaleSqrtWeightCol() {
        ZScoreTools.inPlaceScaleSqrtWeightCol(matrix, weights);
    }

    public void export(File outputDirectory, String stem) {
        String path1 = new File(outputDirectory, stem + ".matrix.npy").getAbsolutePath();
        MatrixTools.saveMatrixTextNumpy(path1, matrix);

        String path2 = new File(outputDirectory, stem + ".weights.npy").getAbsolutePath();
        MatrixTools.saveMatrixTextNumpy(path2, weights);
    }

    public void processKMeansClusteringResult(Cluster[] clusters, GenomeWide1DList<SubcompartmentInterval> subcompartments) {

        Set<SubcompartmentInterval> subcompartmentIntervals = new HashSet<>();

        int genomewideCompartmentID = 0;
        for (Cluster cluster : clusters) {
            int currentClusterID = ++genomewideCompartmentID;
            for (int i : cluster.getMemberIndexes()) {
                if (map.containsKey(i)) {
                    SubcompartmentInterval interv = map.get(i);
                    if (interv != null) {
                        subcompartmentIntervals.add(generateNewSubcompartment(interv, currentClusterID));
                    }
                }
            }
        }

        subcompartments.addAll(new ArrayList<>(subcompartmentIntervals));
        SliceUtils.reSort(subcompartments);
    }

    protected SubcompartmentInterval generateNewSubcompartment(SubcompartmentInterval interv, int currentClusterID) {
        SubcompartmentInterval newInterv = (SubcompartmentInterval) interv.deepClone();
        newInterv.setClusterID(currentClusterID);
        return newInterv;
    }

    private void populateRowIndexToIntervalMap(Mappings mappings) {
        int resolution = mappings.getResolution();
        Chromosome[] chromosomes = mappings.getChromosomes();
        for (Chromosome chromosome : chromosomes) {
            int maxGenomeLen = (int) chromosome.getLength();
            int[] globalIndices = mappings.getGlobalIndex(chromosome);
            for (int i = 0; i < globalIndices.length; i++) {
                if (globalIndices[i] > -1) {
                    int coord = globalIndices[i];
                    int x1 = i * resolution;
                    int x2 = Math.min(x1 + resolution, maxGenomeLen);
                    SubcompartmentInterval newRInterval = new SubcompartmentInterval(chromosome, x1, x2, -1);
                    map.put(coord, newRInterval);
                }
            }
        }
    }

    public void removeAllZeroRows() {
        matrix = EmptyRowCleaner.cleanUpMatrix(matrix, mappings);
    }

    public int[] getSumOfAllLoci(Chromosome[] chromosomes) {
        int[] totalLoci = new int[mappings.getNumCols()];
        for (Chromosome chromosome : chromosomes) {
            int[] row = mappings.getDistributionForChrom(chromosome);
            for (int z = 0; z < row.length; z++) {
                totalLoci[z] += row[z];
            }
        }
        return totalLoci;
    }

    public void updateWeights(Chromosome[] chromosomes) {
        int[] totalDistribution = getSumOfAllLoci(chromosomes);
        System.arraycopy(totalDistribution, 0, weights, 0, mappings.getNumCols());
    }
}

