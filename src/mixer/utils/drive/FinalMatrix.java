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
import javastraw.reader.basics.ChromosomeHandler;
import javastraw.tools.MatrixTools;
import mixer.utils.cleaning.NaNRowCleaner;
import mixer.utils.common.ZScoreTools;
import mixer.utils.tracks.SliceUtils;
import mixer.utils.tracks.SubcompartmentInterval;
import robust.concurrent.kmeans.clustering.Cluster;

import java.io.File;
import java.util.ArrayList;
import java.util.HashSet;
import java.util.Map;
import java.util.Set;

public class FinalMatrix {
    private final Map<Integer, SubcompartmentInterval> map;
    private final Mappings mappings;
    public float[][] matrix;
    public int[] weights;

    public FinalMatrix(float[][] matrix, int[] weights, Mappings mappings, Map<Integer, SubcompartmentInterval> map) {
        this.matrix = matrix;
        this.weights = weights;
        this.mappings = mappings;
        this.map = map;
    }

    public void removeAllNanRows() {
        matrix = NaNRowCleaner.cleanUpMatrix(matrix, mappings);
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

    public GenomeWide1DList<SubcompartmentInterval> getClusteringResult(int[] assignments, ChromosomeHandler handler) {
        GenomeWide1DList<SubcompartmentInterval> subcompartments = new GenomeWide1DList<>(handler);

        Set<SubcompartmentInterval> subcompartmentIntervals = new HashSet<>();
        for (int i = 0; i < assignments.length; i++) {
            if (assignments[i] > -1) {
                int currentClusterID = assignments[i];
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
        return subcompartments;
    }

    protected SubcompartmentInterval generateNewSubcompartment(SubcompartmentInterval interv, int currentClusterID) {
        SubcompartmentInterval newInterv = (SubcompartmentInterval) interv.deepClone();
        newInterv.setClusterID(currentClusterID);
        return newInterv;
    }

    public int getNumRows() {
        return matrix.length;
    }

    public int getNumCols() {
        return matrix[0].length;
    }

    public boolean notEmpty() {
        return matrix.length > 10 && matrix[0].length > 2;
    }


}
