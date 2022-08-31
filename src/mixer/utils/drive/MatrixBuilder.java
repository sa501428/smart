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
import mixer.MixerTools;
import mixer.utils.translocations.SimpleTranslocationFinder;

import java.io.File;
import java.util.HashMap;
import java.util.Iterator;
import java.util.Map;

public class MatrixBuilder {
    public static MatrixAndWeight populateMatrix(Dataset ds, Chromosome[] chromosomes, int resolution,
                                                 NormalizationType norm, Mappings mappings,
                                                 SimpleTranslocationFinder translocations,
                                                 File outputDirectory, boolean useNone) {
        int numRows = mappings.getNumRows();
        int numCols = mappings.getNumCols();
        int[] weights = new int[numCols];
        float[][] matrix = new float[numRows][numCols];

        Map<Integer, int[]> genomewideDistributionForChrom = new HashMap<>();
        for (Chromosome chromosome : chromosomes) {
            genomewideDistributionForChrom.put(chromosome.getIndex(), new int[numCols]);
        }

        System.out.println(".");
        if (MixerTools.printVerboseComments) {
            mappings.printStatus();
        }

        for (int i = 0; i < chromosomes.length; i++) {
            fillInNans(matrix, mappings, chromosomes[i], chromosomes[i]);

            for (int j = i + 1; j < chromosomes.length; j++) {
                if (translocations.contains(chromosomes[i], chromosomes[j])) {
                    fillInNans(matrix, mappings, chromosomes[i], chromosomes[j]);
                    continue;
                }

                Matrix m1 = ds.getMatrix(chromosomes[i], chromosomes[j]);
                if (m1 != null) {
                    MatrixZoomData zd = m1.getZoomData(new HiCZoom(resolution));
                    if (zd != null) {
                        if (useNone || norm.getLabel().equalsIgnoreCase("none")) {
                            populateMatrixFromIterator(matrix, zd.getDirectIterator(), mappings, chromosomes[i], chromosomes[j]);
                        } else {
                            populateMatrixFromIterator(matrix, zd.getNormalizedIterator(norm), mappings, chromosomes[i], chromosomes[j]);
                        }

                        updateNumberOfLoci(genomewideDistributionForChrom.get(chromosomes[i].getIndex()),
                                mappings.getDistributionForChrom(chromosomes[j]));
                        updateNumberOfLoci(genomewideDistributionForChrom.get(chromosomes[j].getIndex()),
                                mappings.getDistributionForChrom(chromosomes[i]));
                    }
                }

                System.out.print(".");
            }
            System.out.println(".");
        }

        return new MatrixAndWeight(matrix, weights, mappings);
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

    private static void updateNumberOfLoci(int[] totalLoci, int[] lociForRegion) {
        for (int i = 0; i < totalLoci.length; i++) {
            totalLoci[i] += lociForRegion[i];
        }
    }

}
