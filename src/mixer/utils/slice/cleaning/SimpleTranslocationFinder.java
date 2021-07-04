/*
 * The MIT License (MIT)
 *
 * Copyright (c) 2011-2021 Rice University, Baylor College of Medicine, Aiden Lab
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

package mixer.utils.slice.cleaning;

import javastraw.reader.Dataset;
import javastraw.reader.basics.Chromosome;
import javastraw.reader.block.Block;
import javastraw.reader.block.ContactRecord;
import javastraw.reader.mzd.MatrixZoomData;
import javastraw.reader.type.NormalizationType;
import javastraw.tools.HiCFileTools;
import mixer.algos.Slice;
import mixer.utils.slice.structures.HiCInterTools;
import org.apache.commons.math3.stat.descriptive.DescriptiveStatistics;

import java.util.ArrayList;
import java.util.List;
import java.util.Set;

public class SimpleTranslocationFinder {

    private final List<TranslocationSet> translocations = new ArrayList<>();

    public SimpleTranslocationFinder(List<Dataset> datasets,
                                     Chromosome[] chroms,
                                     List<NormalizationType[]> normalizationTypes, GWBadIndexFinder badIndexFinder) {

        for (int z = 0; z < datasets.size(); z++) {

            Dataset ds = datasets.get(z);
            int lowestResZoom = HiCInterTools.getLowestResolution(ds);
            NormalizationType normGW = normalizationTypes.get(z)[Slice.GW_SCALE_INDEX];
            NormalizationType normNone = ds.getNormalizationHandler().getNormTypeFromString("NONE");
            NormalizationType normIntra = normalizationTypes.get(z)[Slice.INTRA_SCALE_INDEX];


            TranslocationSet translocationSet = new TranslocationSet();

            for (int i = 0; i < chroms.length; i++) {
                Chromosome chr1 = chroms[i];
                for (int j = i + 1; j < chroms.length; j++) {
                    Chromosome chr2 = chroms[j];
                    final MatrixZoomData zd = HiCFileTools.getMatrixZoomData(ds, chr1, chr2, lowestResZoom);
                    int lengthChr1 = (int) (chr1.getLength() / lowestResZoom + 1);
                    int lengthChr2 = (int) (chr2.getLength() / lowestResZoom + 1);

                    Set<Integer> badIndices1 = badIndexFinder.getBadGenomePositionsAtResolution(chr1, lowestResZoom);
                    Set<Integer> badIndices2 = badIndexFinder.getBadGenomePositionsAtResolution(chr2, lowestResZoom);


                    try {
                        List<Block> blocks = HiCFileTools.getAllRegionBlocks(zd, 0, lengthChr1, 0,
                                lengthChr2, normGW, false);
                        translocationSet.put(chr1, chr2, getMean(blocks, badIndices1, badIndices2));
                    } catch (Exception e) {
                        System.err.println(chr1.getName() + " - " + chr2.getName());
                        e.printStackTrace();
                    }
                }
            }
            translocationSet.determineTranslocations();
        }
    }

    private static float getMean(List<Block> blocks, Set<Integer> badIndices1, Set<Integer> badIndices2) {
        //float max = 0;
        DescriptiveStatistics stats0 = new DescriptiveStatistics();
        for (Block b : blocks) {
            if (b != null) {
                for (ContactRecord cr : b.getContactRecords()) {
                    if (isValidEntry(cr, badIndices1, badIndices2)) {
                        //max = Math.max(max, cr.getCounts());
                        stats0.addValue(cr.getCounts());
                    }
                }
            }
        }
        //return (float) (stats.getPercentile(0.99) - 0*stats.getPercentile(0.05));
        /*
        double minCutoff = stats0.getPercentile(0.5);
        double maxCutoff = stats0.getPercentile(0.9);

        DescriptiveStatistics stats2 = new DescriptiveStatistics();
        for (Block b : blocks) {
            if (b != null) {
                for (ContactRecord cr : b.getContactRecords()) {
                    float val = cr.getCounts();
                    if(!Float.isNaN(val)) {
                        //max = Math.max(max, cr.getCounts());
                        if(val > minCutoff && val < maxCutoff) {
                            stats2.addValue(val);
                        }
                    }
                }
            }
        }
        */

        return (float) stats0.getMean();
    }

    private static boolean isValidEntry(ContactRecord cr, Set<Integer> badIndices1, Set<Integer> badIndices2) {
        if (Float.isNaN(cr.getCounts())) return false;
        if (badIndices1.contains(cr.getBinX())) return false;
        return !badIndices2.contains(cr.getBinY());
    }

    public TranslocationSet getSet(int index) {
        return translocations.get(index);
    }
}
