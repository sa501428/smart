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

package mixer.utils.matrix;

import javastraw.reader.Dataset;
import javastraw.reader.basics.Chromosome;
import javastraw.reader.block.ContactRecord;
import javastraw.reader.mzd.MatrixZoomData;
import javastraw.reader.type.NormalizationType;
import javastraw.tools.HiCFileTools;
import javastraw.tools.ParallelizationTools;

import java.util.ArrayList;
import java.util.List;
import java.util.concurrent.atomic.AtomicInteger;

abstract public class HiCMatrix {

    protected final NormalizationType norm;
    protected final int resolution;
    List<List<ContactRecord>> allRecords;
    protected final Chromosome[] rowsChromosomes;
    protected final Chromosome[] colsChromosomes;
    protected final Dimension rowsDimension, colsDimension;

    public HiCMatrix(Dataset ds, NormalizationType norm, int resolution,
                     Chromosome[] rowsChromosomes, Chromosome[] colsChromosomes) {
        this.norm = norm;
        this.resolution = resolution;
        this.rowsChromosomes = rowsChromosomes;
        this.colsChromosomes = colsChromosomes;
        rowsDimension = new Dimension(rowsChromosomes, resolution);
        colsDimension = new Dimension(colsChromosomes, resolution);
        allRecords = makeCleanScaledInterMatrix(ds);
    }

    private List<List<ContactRecord>> makeCleanScaledInterMatrix(Dataset ds) {

        final List<List<ContactRecord>> allRecords = new ArrayList<>(rowsDimension.length * colsDimension.length);

        AtomicInteger index = new AtomicInteger(0);
        ParallelizationTools.launchParallelizedCode(() -> {
            int i = index.getAndIncrement();
            while (i < rowsChromosomes.length) {
                Chromosome chr1 = rowsChromosomes[i];
                for (int j = 0; j < colsChromosomes.length; j++) {
                    Chromosome chr2 = colsChromosomes[j];

                    final MatrixZoomData zd = HiCFileTools.getMatrixZoomData(ds, chr1, chr2, resolution);
                    if (zd == null) continue;

                    boolean needToFlip = chr2.getIndex() < chr1.getIndex();
                    List<ContactRecord> records = getRecords(ds, zd, chr1, rowsDimension.offset[i],
                            chr2, colsDimension.offset[j], needToFlip);
                    synchronized (allRecords) {
                        allRecords.add(records);
                    }
                }
                System.out.print(".");
                i = index.getAndIncrement();
            }
        });
        System.out.println(".");

        return allRecords;
    }

    abstract protected List<ContactRecord> getRecords(Dataset ds, MatrixZoomData zd,
                                                      Chromosome chr1, int offsetIndex1,
                                                      Chromosome chr2, int offsetIndex2, boolean needToFlip);

    public Chromosome[] getRowChromosomes() {
        return rowsChromosomes;
    }

    public Chromosome[] getColChromosomes() {
        return colsChromosomes;
    }

    public List<List<ContactRecord>> getMatrix() {
        return allRecords;
    }

    public int[] getRowOffsets() {
        return rowsDimension.offset;
    }

    public int[] getColOffsets() {
        return colsDimension.offset;
    }
}
