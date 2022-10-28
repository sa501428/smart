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

package mixer.utils.transform;

import javastraw.expected.Welford;
import javastraw.expected.Zscore;
import mixer.utils.drive.Mappings;

public class LogisticTransform {
    public static void transform(float[][] matrix, Mappings mappings) {
        for (int i = 0; i < matrix.length; i++) {
            normalizeRegion0(matrix[i]);
        }
    }

    private static void normalizeRegion0(float[] row) {
        Welford welford = new Welford();
        for (float val : row) {
            if (val > 0) {
                welford.addValue(val);
            }
        }
        if (welford.getCounts() > 2) {
            Zscore zscore = welford.getZscore();
            for (int c = 0; c < row.length; c++) {
                row[c] = (float) zscore.getZscore(row[c]); // transform
            }
        }
    }

    private static int transform(double v) {
        return (int) Math.round(logisticFunction(v));
    }

    private static double logisticFunction(double x) {
        return 4 * ((1 / (1 + Math.exp(-x))) - 0.5);
    }
}
