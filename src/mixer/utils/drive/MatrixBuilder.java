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
import mixer.utils.translocations.TranslocationSet;

import java.util.Iterator;
import java.util.List;

public class MatrixBuilder {

    public static MatrixAndWeight populateMatrix(Dataset ds, Chromosome[] chromosomes, int resolution,
                                                 NormalizationType interNorm,
                                                 NormalizationType intraNorm,
                                                 Mappings mappings,
                                                 TranslocationSet translocations,
                                                 boolean fillInIntraMatrix) {
        int numRows = mappings.getNumRows();
        int numCols = mappings.getNumCols();
        int[] weights = new int[numCols];
        float[][] inter = new float[numRows][numCols];
        float[][] intra = new float[numRows][numCols];
        float[][] counts = new float[numRows][numCols];

        System.out.println(".");
        if (SmartTools.printVerboseComments) {
            mappings.printStatus();
        }

        for (int i = 0; i < chromosomes.length; i++) {
            for (int j = i; j < chromosomes.length; j++) {
                if (i == j) { // INTRA region
                    fillInNans(inter, mappings, chromosomes[i], chromosomes[i]);
                    if (fillInIntraMatrix) {
                        Matrix m1 = ds.getMatrix(chromosomes[i], chromosomes[j]);
                        if (m1 != null) {
                            MatrixZoomData zd = m1.getZoomData(new HiCZoom(resolution));
                            if (zd != null) {
                                populateIntraMatrix(intra, counts, zd, intraNorm, mappings, chromosomes[i],
                                        resolution);
                            }
                        }
                    }
                } else {
                    fillInNans(intra, mappings, chromosomes[i], chromosomes[j]);
                    if (translocations.contains(chromosomes[i], chromosomes[j])) {
                        fillInNans(inter, mappings, chromosomes[i], chromosomes[j]);
                        continue;
                    }

                    Matrix m1 = ds.getMatrix(chromosomes[i], chromosomes[j]);
                    if (m1 != null) {
                        MatrixZoomData zd = m1.getZoomData(new HiCZoom(resolution));
                        if (zd != null) {
                            populateFromIterator(inter,
                                    getIterator(zd, interNorm), mappings, chromosomes[i], chromosomes[j]);
                        }
                    }
                }
                System.out.print(".");
            }
            System.out.println(".");
        }

        if (fillInIntraMatrix) {
            return new MatrixAndWeight(inter, normalize(intra, counts), weights, mappings);
        }
        return new MatrixAndWeight(inter, intra, weights, mappings);
    }

    private static Iterator<ContactRecord> getIterator(MatrixZoomData zd, NormalizationType norm) {
        if (norm.getLabel().equalsIgnoreCase("none")) {
            return zd.getDirectIterator();
        }
        return zd.getNormalizedIterator(norm);
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

    private static void populateFromIterator(float[][] matrix1,
                                             Iterator<ContactRecord> iterator,
                                             Mappings mappings, Chromosome c1, Chromosome c2) {

        if (mappings.contains(c1) && mappings.contains(c2) && iterator != null) {
            int[] binToClusterID1 = mappings.getProtocluster(c1);
            int[] binToClusterID2 = mappings.getProtocluster(c2);
            int[] binToGlobalIndex1 = mappings.getGlobalIndex(c1);
            int[] binToGlobalIndex2 = mappings.getGlobalIndex(c2);

            while (iterator.hasNext()) {
                ContactRecord cr = iterator.next();
                if (cr.getCounts() > 0) {
                    int r = cr.getBinX();
                    int c = cr.getBinY();
                    addValueToMatrix(matrix1, binToClusterID1, binToClusterID2,
                            binToGlobalIndex1, binToGlobalIndex2, r, c, cr.getCounts());
                }
            }
        } else {
            System.err.println("Error with reading from " + c1.getName() + " " + c2.getName());
        }
    }

    private static void addValueToMatrix(float[][] matrix2,
                                         int[] binToClusterID1, int[] binToClusterID2,
                                         int[] binToGlobalIndex1, int[] binToGlobalIndex2,
                                         int r, int c, float counts) {
        if (counts > 0) {
            if (binToClusterID1[r] > -1 && binToClusterID2[c] > -1) {
                matrix2[binToGlobalIndex1[r]][binToClusterID2[c]] += counts;
                matrix2[binToGlobalIndex2[c]][binToClusterID1[r]] += counts;
            }
        }
    }

    private static void populateIntraMatrix(float[][] matrix, float[][] counts, MatrixZoomData zd, NormalizationType intraNorm,
                                            Mappings mappings, Chromosome chromosome, int resolution) {
        List<ContactRecord> filteredContacts = OETools.filter(resolution, getIterator(zd, intraNorm));

        if (filteredContacts.size() > 1) {
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
        } else {
            System.err.println("No intra data from " + chromosome.getName());
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
