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

import javastraw.reader.Dataset;
import javastraw.reader.ExpectedValueFunction;
import javastraw.reader.HiCFileTools;
import javastraw.reader.MatrixZoomData;
import javastraw.reader.basics.Block;
import javastraw.reader.basics.Chromosome;
import javastraw.reader.basics.ContactRecord;
import javastraw.type.NormalizationType;

import java.io.IOException;
import java.util.*;

public class IndexOrderer {

    private final Map<Chromosome, int[]> chromToReorderedIndices = new HashMap<>();
    private final int DISTANCE = 5000000;
    private final int resolution;
    private final int minDistanceThreshold;
    private final int numColsToJoin;

    public IndexOrderer(Dataset ds, Chromosome[] chromosomes, int resolution, NormalizationType normalizationType,
                        int numColumnsToPutTogether, GenomewideBadIndexFinder badIndexLocations) {
        this.resolution = resolution;
        minDistanceThreshold = DISTANCE / resolution;
        numColsToJoin = numColumnsToPutTogether;
        for (Chromosome chrom : chromosomes) {
            final MatrixZoomData zd = HiCFileTools.getMatrixZoomData(ds, chrom, chrom, resolution);
            ExpectedValueFunction df = ds.getExpectedValuesOrExit(zd.getZoom(), normalizationType, chrom, true);
            try {
                float[][] matrix = extractObsOverExpBoundedRegion(zd, chrom, normalizationType, df);
                int[] newOrderIndexes = getNewOrderOfIndices(matrix, badIndexLocations.getBadIndices(chrom));
                chromToReorderedIndices.put(chrom, newOrderIndexes);
            } catch (Exception e) {
                e.printStackTrace();
            }
            System.out.print(".");
        }
    }

    private static double getExpected(int dist, ExpectedValueFunction df, int chrIndex) {
        return df.getExpectedValue(chrIndex, dist);
    }

    public int[] get(Chromosome chrom) {
        return chromToReorderedIndices.get(chrom);
    }

    private int[] getNewOrderOfIndices(float[][] matrix, Set<Integer> badIndices) {
        int[] indices = new int[matrix.length];
        Arrays.fill(indices, -5);

        Set<Integer> columnsToErase = new HashSet<>();
        for (int k : badIndices) {
            indices[k] = -1;
            columnsToErase.add(k / numColsToJoin);
        }

        // erase columns that include these regions
        for (int i = 0; i < matrix.length; i++) {
            for (int j : columnsToErase) {
                matrix[i][j] = Float.NaN;
            }
        }

        int numRoundsThatHappen = 0;
        int counter = 0;
        for (int cI = 0; cI < indices.length; cI++) {
            // handle stuff
            if (indices[cI] < -2) {
                numRoundsThatHappen++;
                indices[cI] = counter++;
                for (int z = cI + 1; z < indices.length; z++) {
                    if (indices[z] < -2) {
                        float val = CorrelationTools.getNonNanPearsonCorrelation(matrix[cI], matrix[z]);
                        if (val >= .5) {
                            indices[z] = counter++;
                        }
                    }
                }
            }
            // else it has already been handled
            // or is a bad index, so skip
        }
        System.out.println("Num rounds " + numRoundsThatHappen);
        return indices;
    }

    public float[][] extractObsOverExpBoundedRegion(MatrixZoomData zd, Chromosome chromosome,
                                                    NormalizationType normalizationType,
                                                    ExpectedValueFunction df) throws IOException {
        if (df == null) {
            System.err.println("DF is null");
            return null;
        }
        // numRows/numCols is just to ensure a set size in case bounds are approximate
        // left upper corner is reference for 0,0
        int maxBin = (int) (chromosome.getLength() / resolution) + 1;
        int maxCompressedBin = maxBin / numColsToJoin;
        if (maxBin % numColsToJoin > 0) maxCompressedBin += 1;

        List<Block> blocks = HiCFileTools.getAllRegionBlocks(zd, 0, maxBin, 0, maxBin,
                normalizationType, true);
        float[][] data = new float[maxBin][maxCompressedBin];

        if (blocks.size() > 0) {
            for (Block b : blocks) {
                if (b != null) {
                    for (ContactRecord rec : b.getContactRecords()) {

                        int dist = Math.abs(rec.getBinX() - rec.getBinY());
                        if (dist < minDistanceThreshold) {
                            placeOEValInPosition(Float.NaN, rec, data);
                        } else {
                            double observed = rec.getCounts() + 1;
                            double expected = getExpected(dist, df, chromosome.getIndex()) + 1;
                            float answer = (float) (observed / expected);
                            // answer = Math.log(observed / expected);

                            if (Float.isNaN(answer) || Float.isInfinite(answer)) {
                                answer = Float.NaN;
                            } else if (answer < 1) {
                                answer = 0;
                            }


                            placeOEValInPosition(answer, rec, data);
                        }
                    }
                }
            }
        }
        // force cleanup
        blocks = null;
        return data;
    }

    private void placeOEValInPosition(float oeVal, ContactRecord rec, float[][] data) {
        int rX = rec.getBinX();
        int rY = rec.getBinY() / numColsToJoin;
        data[rX][rY] += oeVal;

        rX = rec.getBinY();
        rY = rec.getBinX() / numColsToJoin;
        data[rX][rY] += oeVal;
    }
}
