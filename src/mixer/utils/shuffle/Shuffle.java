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

package mixer.utils.shuffle;

import javastraw.feature1D.GenomeWide1DList;
import javastraw.reader.Dataset;
import javastraw.reader.type.NormalizationType;
import mixer.utils.tracks.SliceUtils;
import mixer.utils.tracks.SubcompartmentInterval;

import java.io.File;

public class Shuffle {

    private final Dataset ds;
    private final NormalizationType norm;
    private final int resolution;

    public Shuffle(Dataset ds, NormalizationType norm, int resolution) {
        this.ds = ds;
        this.norm = norm;
        this.resolution = resolution;
    }

    public void runGWStats(GenomeWide1DList<SubcompartmentInterval> subcompartments, File outfolder, String prefix) {
        SliceUtils.collapseGWList(subcompartments);
        GenomeWideStatistics statistics = new GenomeWideStatistics(ds, resolution, norm, subcompartments);
        statistics.writeToFile(outfolder, prefix);
        System.out.println("Interaction summary statistics saved");
    }
}
