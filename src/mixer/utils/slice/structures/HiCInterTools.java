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

package mixer.utils.slice.structures;

import javastraw.reader.ChromosomeHandler;
import javastraw.reader.Dataset;
import javastraw.reader.HiCFileTools;
import javastraw.reader.MatrixZoomData;
import javastraw.reader.basics.Block;
import javastraw.reader.basics.Chromosome;
import javastraw.reader.basics.ContactRecord;
import javastraw.type.HiCZoom;
import javastraw.type.NormalizationType;

import java.util.List;

public class HiCInterTools {

    private static int getIdealCompression(double genomelength, double resolution, double counts) {
        double x = genomelength / resolution;
        return (int) Math.ceil(1.5 * x * x / counts);
    }

    public static int calculateIdealWidth(Dataset ds, int resolution) {

        ChromosomeHandler handler = ds.getChromosomeHandler();
        Chromosome[] chroms = handler.getChromosomeArrayWithoutAllByAll();
        long genomelength = 0L;
        for (Chromosome chrom : chroms) {
            genomelength += chrom.getLength();
        }

        long totalCounts = 0L;
        int lowestResZoom = getLowestResolution(ds);
        NormalizationType normNone = ds.getNormalizationHandler().getNormTypeFromString("NONE");

        for (int i = 0; i < chroms.length; i++) {
            Chromosome chr1 = chroms[i];
            for (int j = i + 1; j < chroms.length; j++) {
                Chromosome chr2 = chroms[j];
                final MatrixZoomData zd = HiCFileTools.getMatrixZoomData(ds, chr1, chr2, lowestResZoom);
                int lengthChr1 = (int) (chr1.getLength() / resolution + 1);
                int lengthChr2 = (int) (chr2.getLength() / resolution + 1);

                try {
                    List<Block> blocks = HiCFileTools.getAllRegionBlocks(zd, 0, lengthChr1, 0,
                            lengthChr2, normNone, false);
                    totalCounts += getTotalCounts(blocks);
                } catch (Exception e) {
                    System.err.println(chr1.getName() + " - " + chr2.getName());
                    e.printStackTrace();
                }
            }
        }

        System.out.println("Total counts: " + totalCounts);

        return getIdealCompression(genomelength, resolution, totalCounts);
    }

    private static int getLowestResolution(Dataset ds) {
        List<HiCZoom> zooms = ds.getBpZooms();
        int maxResolution = zooms.get(0).getBinSize();
        for (HiCZoom zoom : zooms) {
            if (zoom.getBinSize() > maxResolution) {
                maxResolution = zoom.getBinSize();
            }
        }
        return maxResolution;
    }

    private static long getTotalCounts(List<Block> blocks) {
        double total = 0;

        for (Block b : blocks) {
            if (b != null) {
                for (ContactRecord cr : b.getContactRecords()) {
                    total += cr.getCounts();
                }
            }
        }

        return (long) total;
    }
}
