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

package mixer.utils.cleaning;

import mixer.utils.drive.MatrixAndWeight;
import mixer.utils.tracks.SubcompartmentInterval;

import java.io.File;
import java.util.Map;
import java.util.Random;

public class SliceMatrixCleaner {
    protected final File outputDirectory;
    protected float[][] data;
    protected final Random generator = new Random(2352);
    protected int resolution;

    public SliceMatrixCleaner(float[][] data, long seed, File outputDirectory, int resolution) {
        this.outputDirectory = outputDirectory;
        generator.setSeed(seed);
        this.resolution = resolution;
        this.data = data;
    }

    public MatrixAndWeight getCleanFilteredZscoredMatrix(Map<Integer, SubcompartmentInterval> rowIndexToIntervalMap,
                                                         int[] weights) {

        //MatrixAndWeight mw = (new ColumnCleaner(data, weights)).getCleanedData();
        //data = mw.matrix;
        //weights = mw.weights;
        //System.out.println("Matrix size after column cleanup " + mw.matrix.length + " x " + mw.matrix[0].length);

        data = (new RowCleaner(data, rowIndexToIntervalMap, weights)).getCleanedData(resolution, outputDirectory).matrix;

        //ZScoreTools.inPlaceZscoreDownCol(data);

        return new MatrixAndWeight(data, weights, null);
    }

}
