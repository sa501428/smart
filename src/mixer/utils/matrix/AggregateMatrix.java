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

import javastraw.tools.UNIXTools;
import mixer.utils.common.FloatMatrixTools;

import java.io.File;

public class AggregateMatrix {

    private float[][] aggregate = null;

    private static void addBToA(float[][] a, float[][] b) {
        if (a.length == b.length && a[0].length == b[0].length) {
            for (int i = 0; i < a.length; i++) {
                for (int j = 0; j < a[i].length; j++) {
                    a[i][j] += b[i][j];
                }
            }
        } else {
            System.err.println("dimensions incorrect " + a.length + "==" + b.length
                    + "; " + a[0].length + "==" + b[0].length);
        }
    }

    private static void divideBy(float[][] matrix, float scalar) {
        for (int i = 0; i < matrix.length; i++) {
            for (int j = 0; j < matrix[i].length; j++) {
                matrix[i][j] /= scalar;
            }
        }
    }

    public synchronized void addBToA(float[][] matrix) {
        if (aggregate == null) {
            aggregate = new float[matrix.length][matrix[0].length];
        }
        addBToA(aggregate, matrix);
    }

    public void scaleForNumberOfRounds(int numRounds) {
        divideBy(aggregate, numRounds);
    }

    private static void saveToPNG(float[][] matrix, File outfolder, String name) {
        File mapLogFile = new File(outfolder, name + ".png");
        FloatMatrixTools.saveMatrixToPNG(mapLogFile, matrix, false);
    }

    public void saveToPNG(File outfolder, String name) {
        File outfolder2 = new File(outfolder, name);
        UNIXTools.makeDir(outfolder2);
        saveToPNG(aggregate, outfolder2, name);
    }

    public float[][] getMatrixCopy() {
        return FloatMatrixTools.deepClone(aggregate);
    }

    public boolean hasData() {
        return aggregate != null;
    }
}
