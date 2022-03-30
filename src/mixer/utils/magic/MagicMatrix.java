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
import javastraw.reader.Matrix;
import javastraw.reader.basics.Chromosome;
import javastraw.reader.basics.ChromosomeHandler;
import javastraw.reader.block.ContactRecord;
import javastraw.reader.mzd.MatrixZoomData;
import javastraw.reader.type.HiCZoom;
import javastraw.reader.type.NormalizationType;
import javastraw.tools.MatrixTools;
import mixer.utils.drive.DriveMatrix;
import mixer.utils.slice.structures.SubcompartmentInterval;

import java.io.File;
import java.util.HashMap;
import java.util.Iterator;
import java.util.Map;

public class MagicMatrix extends DriveMatrix {

    private final float[][] matrix;
    private final Map<Integer, SubcompartmentInterval> rowIndexToIntervalMap;


    public MagicMatrix(Dataset ds, ChromosomeHandler chromosomeHandler, int resolution,
                       NormalizationType norm, File outputDirectory, long seed, int[] offset,
                       Map<Integer, Integer> binToClusterID, int numRows, int numCols) {
        matrix = new float[numRows][numCols];
        rowIndexToIntervalMap = createIndexToIntervalMap(chromosomeHandler, resolution);
        populateMatrix(ds, chromosomeHandler, resolution, norm, offset, binToClusterID);
    }

    private Map<Integer, SubcompartmentInterval> createIndexToIntervalMap(ChromosomeHandler handler,
                                                                          int resolution) {
        Map<Integer, SubcompartmentInterval> map = new HashMap<>();
        int counter = 0;
        Chromosome[] chromosomes = handler.getAutosomalChromosomesArray();
        for (int q = 0; q < chromosomes.length; q++) {
            Chromosome chromosome = chromosomes[q];
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


    private void populateMatrix(Dataset ds, ChromosomeHandler handler, int resolution,
                                NormalizationType norm, int[] offset, Map<Integer, Integer> binToClusterID) {
        // incorporate norm
        Chromosome[] chromosomes = handler.getAutosomalChromosomesArray();
        for (int i = 0; i < chromosomes.length; i++) {
            for (int j = i + 1; j < chromosomes.length; j++) {
                Matrix m1 = ds.getMatrix(chromosomes[i], chromosomes[j]);
                if (m1 == null) continue;
                MatrixZoomData zd = m1.getZoomData(new HiCZoom(resolution));
                if (zd == null) continue;
                populateMatrixFromIterator(zd.getDirectIterator(), offset[i], offset[j], binToClusterID);
            }
        }
    }

    private void populateMatrixFromIterator(Iterator<ContactRecord> iterator, int rowOffset, int colOffset,
                                            Map<Integer, Integer> binToClusterID) {
        while (iterator.hasNext()) {
            ContactRecord cr = iterator.next();

            if (cr.getCounts() > 0) {
                int r = cr.getBinX() + rowOffset;
                int c = cr.getBinY() + colOffset;
                if (binToClusterID.containsKey(c)) {
                    matrix[r][binToClusterID.get(c)] += cr.getCounts();
                }

                if (binToClusterID.containsKey(r)) {
                    matrix[c][binToClusterID.get(r)] += cr.getCounts();
                }
            }
        }
    }

    public void export(String path) {
        MatrixTools.saveMatrixTextNumpy(path, matrix);
    }

    @Override
    public float[][] getData(boolean getCorrelationMatrix) {
        return matrix;
    }

    @Override
    public Map<Integer, SubcompartmentInterval> getRowIndexToIntervalMap() {
        return rowIndexToIntervalMap;
    }
}
