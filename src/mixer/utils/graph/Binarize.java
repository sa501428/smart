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

package mixer.utils.graph;

import javastraw.expected.LogExpectedSpline;
import javastraw.reader.Dataset;
import javastraw.reader.basics.Chromosome;
import javastraw.reader.mzd.Matrix;
import javastraw.reader.mzd.MatrixZoomData;
import javastraw.reader.norm.NormalizationPicker;
import javastraw.reader.type.HiCZoom;
import javastraw.reader.type.NormalizationType;
import javastraw.tools.HiCFileTools;
import javastraw.tools.MatrixTools;
import javastraw.tools.ParallelizationTools;
import mixer.utils.intra.OETools;

import java.util.HashSet;
import java.util.concurrent.atomic.AtomicInteger;

public class Binarize {
    private static void dumpMatrices() {

        String filename = "/Users/muhammad/Desktop/hicfiles/gm12878_rh14_30.hic";
        Dataset ds = HiCFileTools.extractDatasetForCLT(filename, false, false, false);
        NormalizationType norm = NormalizationPicker.getFirstValidNormInThisOrder(ds, new String[]{"KR", "SCALE", "VC", "VC_SQRT", "NONE"});
        System.out.println("Norm being used: " + norm.getLabel());
        int resolution = 50000;
        Chromosome[] chromosomes = ds.getChromosomeHandler().getChromosomeArrayWithoutAllByAll();

        for (Chromosome chromosome : chromosomes) {
            Matrix matrix = ds.getMatrix(chromosome, chromosome);
            if (matrix == null) continue;
            MatrixZoomData zd = matrix.getZoomData(new HiCZoom(resolution));
            if (zd == null) continue;
            LogExpectedSpline spline = new LogExpectedSpline(zd, norm, chromosome, resolution);

            float[][] data = OETools.getCleanOEMatrix(zd, chromosome, resolution, norm,
                    new HashSet<>(), 1, false, false, spline);

            MatrixTools.saveMatrixTextNumpy(chromosome.getName() + "_raw.npy", data);

            binarizeForThresholdLessThan(data, 2);

            MatrixTools.saveMatrixTextNumpy(chromosome.getName() + "_binary.npy", data);

            /*
            float[][] data2 = multiply(data, data);
            data = null;

            MatrixTools.saveMatrixTextNumpy(chromosome.getName()+"_bsquared.npy", data);

             */


        }
    }

    private static void binarizeForThresholdLessThan(float[][] data, int threshold) {
        for (float[] row : data) {
            for (int k = 0; k < row.length; k++) {
                if (row[k] > threshold) {
                    row[k] = 1;
                } else {
                    row[k] = 0;
                }
            }
        }
    }

    public static float[][] multiply(float[][] a, float[][] b) {
        int numARows = a.length;
        int numACols = a[0].length;
        int numBCols = b[0].length;

        float[][] result = new float[numARows][numBCols];
        AtomicInteger currRowIndex = new AtomicInteger(0);
        ParallelizationTools.launchParallelizedCode(() -> {
            int i = currRowIndex.getAndIncrement();
            while (i < numARows) {
                for (int j = 0; j < numBCols; j++) {
                    for (int k = 0; k < numACols; k++) {
                        result[i][j] += a[i][k] * b[k][j];
                    }
                }
            }
        });
        return result;
    }
}
