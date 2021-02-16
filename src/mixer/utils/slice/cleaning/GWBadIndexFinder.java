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
import javastraw.reader.HiCFileTools;
import javastraw.reader.MatrixZoomData;
import javastraw.reader.basics.Block;
import javastraw.reader.basics.Chromosome;
import javastraw.type.NormalizationType;

import java.io.IOException;
import java.util.List;

public class GWBadIndexFinder extends BadIndexFinder {
    private final GWRegionStatistics stats = new GWRegionStatistics();

    public GWBadIndexFinder(List<Dataset> datasets, Chromosome[] chromosomes, int resolution, NormalizationType[] norms) {
        super(datasets, chromosomes, resolution, norms);
    }

    @Override
    protected void determineBadIndicesForRegion(Chromosome chr1, Chromosome chr2, MatrixZoomData zd,
                                                int dIndex) throws IOException {
        int lengthChr1 = (int) (chr1.getLength() / resolution + 1);
        int lengthChr2 = (int) (chr2.getLength() / resolution + 1);

        List<Block> blocks = HiCFileTools.getAllRegionBlocks(zd, 0, lengthChr1, 0, lengthChr2,
                norms[dIndex], false);
        stats.update(chr1.getIndex(), chr2.getIndex(), lengthChr1, lengthChr2, blocks);
    }

    @Override
    protected void postprocess(Chromosome[] chromosomes) {
        stats.postprocess();

        for (Chromosome chromosome : chromosomes) {
            getBadCoverageIndicesByZscore(chromosome,
                    stats.getSums(chromosome), stats.getSumMean(), stats.getSumStd(), false);
            getBadCoverageIndicesByZscore(chromosome,
                    stats.getNonZeros(chromosome), stats.getNonZeroMean(), stats.getNonZeroStd(), true);
        }
    }

}
