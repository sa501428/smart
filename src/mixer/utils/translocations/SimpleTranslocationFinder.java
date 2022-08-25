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

import javastraw.feature2D.Feature2DList;
import javastraw.reader.Dataset;
import javastraw.reader.basics.Chromosome;
import javastraw.reader.block.Block;
import javastraw.reader.block.ContactRecord;
import javastraw.reader.mzd.MatrixZoomData;
import javastraw.reader.type.NormalizationType;
import javastraw.tools.HiCFileTools;
import mixer.MixerTools;
import mixer.algos.Slice;
import mixer.utils.slice.cleaning.BadIndexFinder;
import mixer.utils.slice.cleaning.TranslocationSet;
import org.apache.commons.math3.stat.descriptive.DescriptiveStatistics;

import java.awt.*;
import java.io.File;
import java.util.ArrayList;
import java.util.HashSet;
import java.util.List;
import java.util.Set;

public class SimpleTranslocationFinder {

    private final List<TranslocationSet> translocations = new ArrayList<>();
    private final File outputDirectory;
    private final int MAX_CUTOFF = 6;

    public SimpleTranslocationFinder(List<Dataset> datasets,
                                     Chromosome[] chroms,
                                     List<NormalizationType[]> normalizationTypes,
                                     BadIndexFinder badIndexFinder, File outputDirectory) {
        this.outputDirectory = outputDirectory;

        for (int z = 0; z < datasets.size(); z++) {

            Dataset ds = datasets.get(z);
            int lowestResZoom = 1000000;// HiCInterTools.getLowestResolution(ds);
            NormalizationType normGW = normalizationTypes.get(z)[Slice.GW_SCALE_INDEX];
            NormalizationType normNone = ds.getNormalizationHandler().getNormTypeFromString("NONE");
            NormalizationType normIntra = normalizationTypes.get(z)[Slice.INTRA_SCALE_INDEX];


            //TranslocationSet translocationSet = new TranslocationSet();
            DescriptiveStatistics gwStats = initialPass(chroms, ds, lowestResZoom, normGW, badIndexFinder);
            double mean = gwStats.getMean();
            double std = gwStats.getStandardDeviation();
            System.out.println("mu " + mean + " sigma " + std);
            TranslocationSet ts = secondPass(chroms, ds, lowestResZoom, normGW, badIndexFinder,
                    mean, std);

            translocations.add(ts);
            //translocationSet.determineTranslocations();
        }
    }

    private static void updateStatsFirstPass(DescriptiveStatistics stats, List<Block> blocks, Set<Integer> badIndices1, Set<Integer> badIndices2) {
        for (Block b : blocks) {
            if (b != null) {
                for (ContactRecord cr : b.getContactRecords()) {
                    if (isValidEntry(cr, badIndices1, badIndices2)) {
                        //max = Math.max(max, cr.getCounts());
                        stats.addValue(Math.log(1 + cr.getCounts()));
                    }
                }
            }
        }
    }

    private static boolean isValidEntry(ContactRecord cr, Set<Integer> badIndices1, Set<Integer> badIndices2) {
        if (Float.isNaN(cr.getCounts())) return false;
        if (badIndices1.contains(cr.getBinX())) return false;
        return !badIndices2.contains(cr.getBinY());
        //return true;
    }

    private DescriptiveStatistics initialPass(Chromosome[] chroms, Dataset ds, int lowestResZoom,
                                              NormalizationType normGW, BadIndexFinder badIndexFinder) {
        DescriptiveStatistics gwStats = new DescriptiveStatistics();
        for (int i = 0; i < chroms.length; i++) {
            Chromosome chr1 = chroms[i];
            for (int j = i + 1; j < chroms.length; j++) {
                Chromosome chr2 = chroms[j];
                final MatrixZoomData zd = HiCFileTools.getMatrixZoomData(ds, chr1, chr2, lowestResZoom);
                int lengthChr1 = (int) (chr1.getLength() / lowestResZoom + 1);
                int lengthChr2 = (int) (chr2.getLength() / lowestResZoom + 1);

                Set<Integer> badIndices1 = new HashSet<>(); // todo  badIndexFinder.getBadGenomePositionsAtResolution(chr1, lowestResZoom);
                Set<Integer> badIndices2 = new HashSet<>(); // todo badIndexFinder.getBadGenomePositionsAtResolution(chr2, lowestResZoom);

                try {
                    List<Block> blocks = HiCFileTools.getAllRegionBlocks(zd, 0, lengthChr1, 0,
                            lengthChr2, normGW, false);
                    updateStatsFirstPass(gwStats, blocks, badIndices1, badIndices2);
                    //translocationSet.put(chr1, chr2, getMean(blocks, badIndices1, badIndices2));
                } catch (Exception e) {
                    System.err.println(chr1.getName() + " - " + chr2.getName());
                    e.printStackTrace();
                }
            }
        }
        return gwStats;
    }

    private TranslocationSet secondPass(Chromosome[] chroms, Dataset ds, int lowestResZoom,
                                        NormalizationType normGW, BadIndexFinder badIndexFinder,
                                        double mean, double stdDev) {
        TranslocationSet translocationSet = new TranslocationSet();
        Feature2DList feature2DList = new Feature2DList();
        Feature2DList allFeature2DList = new Feature2DList();
        for (int i = 0; i < chroms.length; i++) {
            Chromosome chr1 = chroms[i];
            for (int j = i + 1; j < chroms.length; j++) {
                Chromosome chr2 = chroms[j];
                final MatrixZoomData zd = HiCFileTools.getMatrixZoomData(ds, chr1, chr2, lowestResZoom);
                int lengthChr1 = (int) (chr1.getLength() / lowestResZoom + 1);
                int lengthChr2 = (int) (chr2.getLength() / lowestResZoom + 1);

                Set<Integer> badIndices1 = new HashSet<>(); // todo  badIndexFinder.getBadGenomePositionsAtResolution(chr1, lowestResZoom);
                Set<Integer> badIndices2 = new HashSet<>(); // todo  badIndexFinder.getBadGenomePositionsAtResolution(chr2, lowestResZoom);

                try {
                    List<Block> blocks = HiCFileTools.getAllRegionBlocks(zd, 0, lengthChr1, 0,
                            lengthChr2, normGW, false);
                    List<Rectangle> bounds = new ArrayList<>();

                    /*
                    Feature2DList feature2DListTemp = findTranslocations(mean, stdDev, blocks, badIndices1, badIndices2, chr1, chr2,
                            lowestResZoom, allFeature2DList, bounds);
                    feature2DList.add(feature2DListTemp);
                    if (bounds.size() > 0) {
                        translocationSet.put(chr1, chr2, bounds);
                    }
                    */
                } catch (Exception e) {
                    System.err.println(chr1.getName() + " - " + chr2.getName());
                    e.printStackTrace();
                }
            }
        }

        if (true || MixerTools.printVerboseComments) {
            File outfile = new File(outputDirectory, "suspected_translocations_" + lowestResZoom + ".bedpe");
            feature2DList.exportFeatureList(outfile, false, Feature2DList.ListFormat.NA);
            outfile = new File(outputDirectory, "all_suspected_translocations_" + lowestResZoom + ".bedpe");
            allFeature2DList.exportFeatureList(outfile, false, Feature2DList.ListFormat.NA);
        }
        return translocationSet;
    }

    public TranslocationSet getSet(int index) {
        return translocations.get(index);
    }
}
