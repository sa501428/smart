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

package mixer.utils.slice.matrices;

import javastraw.tools.MatrixTools;
import mixer.utils.common.ZScoreTools;

import java.io.File;

public class MatrixAndWeight {
    public float[][] matrix;
    public int[] weights;

    public MatrixAndWeight(float[][] interMatrix, int[] weights) {
        this.matrix = interMatrix;
        this.weights = weights;
    }

    public void inPlaceScaleSqrtWeightCol() {
        ZScoreTools.inPlaceScaleSqrtWeightCol(matrix, weights);
    }

    public void export(File outputDirectory, String stem) {
        String path1 = new File(outputDirectory, stem + ".matrix.npy").getAbsolutePath();
        MatrixTools.saveMatrixTextNumpy(path1, matrix);

        String path2 = new File(outputDirectory, stem + ".weights.npy").getAbsolutePath();
        MatrixTools.saveMatrixTextNumpy(path2, weights);
    }
}
