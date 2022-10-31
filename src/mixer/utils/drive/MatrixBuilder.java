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

import javastraw.reader.Dataset;
import javastraw.reader.basics.Chromosome;
import javastraw.reader.block.ContactRecord;
import javastraw.reader.mzd.Matrix;
import javastraw.reader.mzd.MatrixZoomData;
import javastraw.reader.type.HiCZoom;
import javastraw.reader.type.NormalizationType;
import mixer.SmartTools;
import mixer.utils.intra.OETools;
import mixer.utils.translocations.SimpleTranslocationFinder;

import java.io.File;
import java.util.Iterator;
import java.util.List;

public class MatrixBuilder {

    private static final int FIVE_MB = 5000000;

    public static MatrixAndWeight populateMatrix(Dataset ds, Chromosome[] chromosomes, int resolution,
                                                 NormalizationType interNorm,
                                                 NormalizationType intraNorm,
                                                 Mappings mappings,
                                                 SimpleTranslocationFinder translocations,
                                                 File outputDirectory) {
        int numRows = mappings.getNumRows();
        int numCols = mappings.getNumCols();
        int[] weights = new int[numCols];
        float[][] matrix = new float[numRows][numCols];
        float[][] intra = new float[numRows][numCols];
        float[][] matrix2 = new float[numRows][numCols];
        float[][] counts = new float[numRows][numCols];

        System.out.println(".");
        if (SmartTools.printVerboseComments) {
            mappings.printStatus();
        }

        for (int i = 0; i < chromosomes.length; i++) {
            for (int j = i; j < chromosomes.length; j++) {

                // INTRA region
                if (i == j) {
                    fillInNans(matrix, mappings, chromosomes[i], chromosomes[i]);
                    fillInNans(matrix2, mappings, chromosomes[i], chromosomes[i]);

                    Matrix m1 = ds.getMatrix(chromosomes[i], chromosomes[j]);
                    if (m1 != null) {
                        MatrixZoomData zd = m1.getZoomData(new HiCZoom(resolution));
                        if (zd != null) {
                            populateIntraMatrix(intra, counts, zd, intraNorm, mappings, chromosomes[i],
                                    resolution);
                        }
                    }
                } else {
                    fillInNans(intra, mappings, chromosomes[i], chromosomes[j]);
                    if (translocations.contains(chromosomes[i], chromosomes[j])) {
                        fillInNans(matrix, mappings, chromosomes[i], chromosomes[j]);
                        fillInNans(matrix2, mappings, chromosomes[i], chromosomes[j]);
                        continue;
                    }

                    Matrix m1 = ds.getMatrix(chromosomes[i], chromosomes[j]);
                    if (m1 != null) {
                        MatrixZoomData zd = m1.getZoomData(new HiCZoom(resolution));
                        if (zd != null) {
                            double[] norm1R = getNorm(ds, chromosomes[i], resolution, interNorm);
                            double[] norm1C = getNorm(ds, chromosomes[j], resolution, interNorm);
                            double[] norm2R = getNorm(ds, chromosomes[i], resolution, intraNorm);
                            double[] norm2C = getNorm(ds, chromosomes[j], resolution, intraNorm);
                            populateTwoMatriceFromIterator(matrix, matrix2,
                                    norm1R, norm1C, norm2R, norm2C,
                                    zd.getDirectIterator(), mappings, chromosomes[i], chromosomes[j]);
                        }
                    }
                }

                System.out.print(".");
            }
            System.out.println(".");
        }

        return new MatrixAndWeight(matrix, normalize(intra, counts), matrix2, weights, mappings);
    }

    private static double[] getNorm(Dataset ds, Chromosome chromosome, int resolution, NormalizationType norm) {
        return ds.getNormalizationVector(chromosome.getIndex(), new HiCZoom(resolution), norm).getData().getValues().get(0);
    }

    private static void fillInNans(float[][] matrix, Mappings mappings, Chromosome c1, Chromosome c2) {
        if (mappings.contains(c1) && mappings.contains(c2)) {
            int[] binToClusterID1 = mappings.getProtocluster(c1);
            int[] binToClusterID2 = mappings.getProtocluster(c2);
            int[] binToGlobalIndex1 = mappings.getGlobalIndex(c1);
            int[] binToGlobalIndex2 = mappings.getGlobalIndex(c2);

            for (int r = 0; r < binToClusterID1.length; r++) {
                for (int c = 0; c < binToClusterID2.length; c++) {
                    if (binToClusterID1[r] > -1 && binToClusterID2[c] > -1) {
                        matrix[binToGlobalIndex1[r]][binToClusterID2[c]] = Float.NaN;
                        matrix[binToGlobalIndex2[c]][binToClusterID1[r]] = Float.NaN;
                    }
                }
            }
        } else {
            System.err.println("Error with reading from " + c1.getName() + " " + c2.getName());
        }
    }


    // todo move to mapping?
    private static void populateMatrixFromIterator(float[][] matrix, Iterator<ContactRecord> iterator,
                                                   Mappings mappings, Chromosome c1, Chromosome c2) {

        if (mappings.contains(c1) && mappings.contains(c2)) {
            int[] binToClusterID1 = mappings.getProtocluster(c1);
            int[] binToClusterID2 = mappings.getProtocluster(c2);
            int[] binToGlobalIndex1 = mappings.getGlobalIndex(c1);
            int[] binToGlobalIndex2 = mappings.getGlobalIndex(c2);

            while (iterator.hasNext()) {
                ContactRecord cr = iterator.next();
                if (cr.getCounts() > 0) {
                    int r = cr.getBinX();
                    int c = cr.getBinY();
                    if (binToClusterID1[r] > -1 && binToClusterID2[c] > -1) {
                        matrix[binToGlobalIndex1[r]][binToClusterID2[c]] += cr.getCounts();
                        matrix[binToGlobalIndex2[c]][binToClusterID1[r]] += cr.getCounts();
                    }
                }
            }
        } else {
            System.err.println("Error with reading from " + c1.getName() + " " + c2.getName());
        }
    }

    private static void populateTwoMatriceFromIterator(float[][] matrix1, float[][] matrix2,
                                                       double[] norm1R, double[] norm1C,
                                                       double[] norm2R, double[] norm2C,
                                                       Iterator<ContactRecord> iterator,
                                                       Mappings mappings, Chromosome c1, Chromosome c2) {

        if (mappings.contains(c1) && mappings.contains(c2)) {
            int[] binToClusterID1 = mappings.getProtocluster(c1);
            int[] binToClusterID2 = mappings.getProtocluster(c2);
            int[] binToGlobalIndex1 = mappings.getGlobalIndex(c1);
            int[] binToGlobalIndex2 = mappings.getGlobalIndex(c2);

            while (iterator.hasNext()) {
                ContactRecord cr = iterator.next();
                if (cr.getCounts() > 0) {
                    int r = cr.getBinX();
                    int c = cr.getBinY();
                    addValueToMatrix(matrix1, norm1R, norm1C, binToClusterID1, binToClusterID2,
                            binToGlobalIndex1, binToGlobalIndex2, r, c, cr.getCounts());
                    addValueToMatrix(matrix2, norm2R, norm2C, binToClusterID1, binToClusterID2,
                            binToGlobalIndex1, binToGlobalIndex2, r, c, cr.getCounts());
                }
            }
        } else {
            System.err.println("Error with reading from " + c1.getName() + " " + c2.getName());
        }
    }

    private static void addValueToMatrix(float[][] matrix2, double[] normR, double[] normC,
                                         int[] binToClusterID1, int[] binToClusterID2,
                                         int[] binToGlobalIndex1, int[] binToGlobalIndex2,
                                         int r, int c, float counts) {
        if (normR[r] > 0 && normC[c] > 0) {
            if (binToClusterID1[r] > -1 && binToClusterID2[c] > -1) {
                float val = (float) (counts / (normR[r] * normC[c]));
                matrix2[binToGlobalIndex1[r]][binToClusterID2[c]] += val;
                matrix2[binToGlobalIndex2[c]][binToClusterID1[r]] += val;
            }
        }
    }

    private static void populateIntraMatrix(float[][] matrix, float[][] counts, MatrixZoomData zd, NormalizationType intraNorm,
                                            Mappings mappings, Chromosome chromosome, int resolution) {
        List<ContactRecord> filteredContacts = OETools.filter(resolution, zd.getNormalizedIterator(intraNorm));

        LogExpectedSubset expected = new LogExpectedSubset(filteredContacts, chromosome, resolution);

        if (mappings.contains(chromosome)) {
            int[] binToClusterID = mappings.getProtocluster(chromosome);
            int[] binToGlobalIndex = mappings.getGlobalIndex(chromosome);
            for (ContactRecord cr : filteredContacts) {
                if (expected.isInInterval(cr)) {
                    int r = cr.getBinX();
                    int c = cr.getBinY();
                    if (binToClusterID[r] > -1 && binToClusterID[c] > -1) {
                        float z = expected.getZscoreForObservedUncompressedBin(cr);
                        if (Math.abs(z) < 5) {
                            matrix[binToGlobalIndex[r]][binToClusterID[c]] += z;
                            matrix[binToGlobalIndex[c]][binToClusterID[r]] += z;
                            counts[binToGlobalIndex[r]][binToClusterID[c]]++;
                            counts[binToGlobalIndex[c]][binToClusterID[r]]++;
                        }
                    }
                }
            }
        } else {
            System.err.println("Error with intra reading from " + chromosome.getName());
        }

        filteredContacts.clear();
        filteredContacts = null;
    }

    private static float[][] normalize(float[][] intra, float[][] counts) {
        for (int i = 0; i < intra.length; i++) {
            for (int j = 0; j < intra[i].length; j++) {
                if (counts[i][j] > 10) {
                    intra[i][j] = intra[i][j] / counts[i][j];
                } else {
                    intra[i][j] = Float.NaN;
                }
            }
        }
        return intra;
    }
}
