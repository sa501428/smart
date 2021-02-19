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

package mixer.utils.slice.cleaning;

import javastraw.reader.basics.Chromosome;

public class NearDiagonalTrim {

    private static final int DISTANCE_CUTOFF = 3000000;

    public static void trim(Chromosome chrom, float[][] matrix, int resolution) {
        if (chrom.getLength() > 5 * DISTANCE_CUTOFF) {
            trimDiagonalWithin(matrix, resolution);
        }
    }

    private static void trimDiagonalWithin(float[][] data, int resolution) {
        int pixelDistance = DISTANCE_CUTOFF / resolution;
        for (int i = 0; i < data.length; i++) {
            int limit = Math.min(data[i].length, i + pixelDistance);
            for (int j = i; j < limit; j++) {
                data[i][j] = Float.NaN;
                data[j][i] = Float.NaN;
            }
        }
    }
}
