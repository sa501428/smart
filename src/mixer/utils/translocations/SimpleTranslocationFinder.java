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

package mixer.utils.translocations;

import javastraw.reader.Dataset;
import javastraw.reader.basics.Chromosome;
import javastraw.reader.block.ContactRecord;
import javastraw.reader.expected.ExpectedValueFunction;
import javastraw.reader.mzd.Matrix;
import javastraw.reader.mzd.MatrixZoomData;
import javastraw.reader.type.HiCZoom;
import javastraw.reader.type.NormalizationType;
import mixer.algos.Slice;

import java.io.File;
import java.util.Iterator;
import java.util.List;
import java.util.Map;
import java.util.Set;

public class SimpleTranslocationFinder {

    private final static int minToBeTranslocation = 10;
    private static final int distance = 1000000;

    public static TranslocationSet find(Dataset ds, NormalizationType[] norms, File outputDirectory,
                                        Map<Integer, Set<Integer>> badIndices, int hiRes) {
        TranslocationSet tSet = new TranslocationSet();

        Chromosome[] chroms = ds.getChromosomeHandler().getAutosomalChromosomesArray();
        int lowRes = Math.max(getLowestResolution(ds, 1000000), hiRes);
        NormalizationType norm = norms[Slice.INTRA_SCALE_INDEX];
        int factor = lowRes / hiRes;

        for (int i = 0; i < chroms.length; i++) {
            for (int j = i + 1; j < chroms.length; j++) {
                if (hasTranslocation(chroms[i], chroms[j], ds, norm, badIndices, factor, lowRes)) {
                    tSet.add(chroms[i], chroms[j]);
                    System.out.println("Potential translocation at " + chroms[i].getName() + " - " + chroms[j].getName());
                }
                System.out.print(".");
            }
            System.out.println(".");
        }
        return tSet;
    }

    /**
     * @return preferred resolution if available, otherwise the lowest resolution present
     */
    public static int getLowestResolution(Dataset ds, int preferredResolution) {
        List<HiCZoom> zooms = ds.getBpZooms();
        int maxResolution = zooms.get(0).getBinSize();
        for (HiCZoom zoom : zooms) {
            if (zoom.getBinSize() > maxResolution) {
                maxResolution = zoom.getBinSize();
            }
            if (zoom.getBinSize() == preferredResolution) return preferredResolution;
        }
        return maxResolution;
    }

    private static boolean hasTranslocation(Chromosome chrom1, Chromosome chrom2, Dataset ds, NormalizationType norm,
                                            Map<Integer, Set<Integer>> badIndices, int factor, int lowRes) {
        Matrix matrix = ds.getMatrix(chrom1, chrom2);
        if (matrix != null) {
            MatrixZoomData zd = matrix.getZoomData(new HiCZoom(lowRes));
            if (zd != null) {
                double translocationCutoff = Double.MAX_VALUE;
                try {
                    translocationCutoff = getCutoff(ds.getExpectedValues(new HiCZoom(lowRes), norm, false),
                            chrom1, chrom2, distance / lowRes);
                } catch (Exception e) {
                    System.err.println("Expected missing; skipping translocation check");
                    return false;
                }

                Set<Integer> badSet1 = badIndices.get(chrom1.getIndex());
                Set<Integer> badSet2 = badIndices.get(chrom2.getIndex());

                int count = 0;
                Iterator<ContactRecord> iterator = zd.getNormalizedIterator(norm);
                while (iterator.hasNext()) {
                    ContactRecord record = iterator.next();
                    if (record.getCounts() > translocationCutoff) {
                        if (isGoodRow(record.getBinX(), badSet1, factor)
                                && isGoodRow(record.getBinY(), badSet2, factor)) {
                            if (count++ > minToBeTranslocation) {
                                return true;
                            }
                        }
                    }
                }
            }
        }
        return false;
    }

    private static boolean isGoodRow(int bin, Set<Integer> badSet, int factor) {
        int binH = bin * factor;
        for (int i = 0; i < factor; i++) {
            if (badSet.contains(binH + i)) return false;
        }
        return true;
    }

    private static double getCutoff(ExpectedValueFunction expectedValues, Chromosome chrom1, Chromosome chrom2, int dIndex) {
        double val1 = expectedValues.getExpectedValuesWithNormalization(chrom1.getIndex()).getValues().get(0)[dIndex];
        double val2 = expectedValues.getExpectedValuesWithNormalization(chrom2.getIndex()).getValues().get(0)[dIndex];
        return Math.min(val1, val2);
    }
}
