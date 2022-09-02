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

package mixer.utils.magic;

import javastraw.reader.Dataset;
import javastraw.reader.basics.Chromosome;
import javastraw.reader.basics.ChromosomeHandler;
import javastraw.reader.type.NormalizationType;
import javastraw.tools.MatrixTools;
import mixer.utils.bed.BedFileMappings;
import mixer.utils.common.ZScoreTools;
import mixer.utils.drive.MatrixAndWeight;
import mixer.utils.tracks.SubcompartmentInterval;
import mixer.utils.translocations.InterChromosomeRegion;

import java.io.File;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

public class MagicMatrix extends DriveMatrix {

    private int numRows, numCols;
    private final float[][] matrix;
    private int[] weights;
    private final Map<Integer, SubcompartmentInterval> rowIndexToIntervalMap;
    private final List<InterChromosomeRegion> regionsToIgnore;


    public MagicMatrix(Dataset ds, ChromosomeHandler chromosomeHandler, int resolution,
                       NormalizationType norm, File outputDirectory, long seed, BedFileMappings mappings,
                       List<InterChromosomeRegion> regionsToIgnore, boolean doScale, boolean useZscore) {


        this.regionsToIgnore = regionsToIgnore;
        this.rowIndexToIntervalMap = createIndexToIntervalMap(chromosomeHandler, resolution);
        matrix = null; // todo populateMatrix(ds, chromosomeHandler, resolution, norm, mappings, doScale);
        System.out.println("MAGIC matrix loaded");

        if (useZscore) {
            ZScoreTools.inPlaceZscorePositivesDownColAndSetZeroToNan(matrix);
            System.out.println("MAGIC matrix zScored");
        }

        inPlaceScaleSqrtWeightCol();
        System.out.println("final magic matrix num rows: " + matrix.length);
    }

    private Map<Integer, SubcompartmentInterval> createIndexToIntervalMap(ChromosomeHandler handler,
                                                                          int resolution) {
        Map<Integer, SubcompartmentInterval> map = new HashMap<>();
        int counter = 0;
        Chromosome[] chromosomes = handler.getAutosomalChromosomesArray();
        for (Chromosome chromosome : chromosomes) {
            int maxGenomeLen = (int) chromosome.getLength();
            int chrBinLength = (int) (chromosome.getLength() / resolution + 1);
            for (int i = 0; i < chrBinLength; i++) {
                int x1 = i * resolution;
                int x2 = Math.min(x1 + resolution, maxGenomeLen);
                SubcompartmentInterval newRInterval = new SubcompartmentInterval(chromosome, x1, x2, counter);
                map.put(counter, newRInterval);
                counter++;
            }
        }
        return map;
    }




    public void export(String path) {
        MatrixTools.saveMatrixTextNumpy(path, matrix);
    }

    @Override
    public MatrixAndWeight getData() {
        return new MatrixAndWeight(matrix, weights, null);
    }

    @Override
    public Map<Integer, SubcompartmentInterval> getRowIndexToIntervalMap() {
        return rowIndexToIntervalMap;
    }

    @Override
    public void inPlaceScaleSqrtWeightCol() {
        ZScoreTools.inPlaceScaleSqrtWeightCol(matrix, weights);
    }


}
