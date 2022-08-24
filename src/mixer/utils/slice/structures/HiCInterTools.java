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

package mixer.utils.slice.structures;

import javastraw.reader.Dataset;
import javastraw.reader.basics.Chromosome;
import javastraw.reader.basics.ChromosomeHandler;
import javastraw.reader.block.ContactRecord;
import javastraw.reader.mzd.MatrixZoomData;
import javastraw.reader.type.HiCZoom;
import javastraw.reader.type.NormalizationType;
import javastraw.tools.HiCFileTools;

import java.util.Iterator;
import java.util.List;

public class HiCInterTools {

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
            for (int j = i + 1; j < chroms.length; j++) {
                final MatrixZoomData zd = HiCFileTools.getMatrixZoomData(ds, chroms[i], chroms[j], lowestResZoom);
                Iterator<ContactRecord> iterator = zd.getDirectIterator();
                while (iterator.hasNext()) {
                    ContactRecord cr = iterator.next();
                    if (cr.getCounts() > 0) {
                        totalCounts += cr.getCounts();
                    }
                }
            }
        }
        System.out.println("Total counts: " + totalCounts);

        return getIdealCompression(genomelength, resolution, totalCounts);
    }

    private static int getIdealCompression(double genomelength, double resolution, double counts) {
        double x = genomelength / resolution;
        return (int) Math.ceil(20 * x * x / counts);
    }

    public static int getLowestResolution(Dataset ds) {
        List<HiCZoom> zooms = ds.getBpZooms();
        int maxResolution = zooms.get(0).getBinSize();
        for (HiCZoom zoom : zooms) {
            if (zoom.getBinSize() > maxResolution) {
                maxResolution = zoom.getBinSize();
            }
        }
        return maxResolution;
    }
}
