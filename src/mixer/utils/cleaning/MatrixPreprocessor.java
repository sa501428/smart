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

import javastraw.reader.basics.Chromosome;
import mixer.utils.drive.FinalMatrix;
import mixer.utils.drive.MatrixAndWeight;

public class MatrixPreprocessor {

    private static final int ZSCORE_LIMIT = 3, MAX_ZSCORE_UPPER_LIMIT = 5;

    public static FinalMatrix preprocess(MatrixAndWeight matrix, Chromosome[] chromosomes,
                                         boolean includeIntra, boolean useLog, boolean useBothNorms,
                                         boolean appendIntra, boolean useRowZ, boolean shouldRegularize) {
        matrix.updateWeights(chromosomes);
        matrix.divideColumnsByWeights(useBothNorms);

        if (shouldRegularize) {
            matrix.applyLog(useBothNorms);
            matrix.removeHighGlobalThresh(useBothNorms, MAX_ZSCORE_UPPER_LIMIT);
            matrix.setZerosToNan(useBothNorms);
            matrix.regularize(useBothNorms, ZSCORE_LIMIT);
            matrix.applyExpm1(useBothNorms);
        }

        if (useRowZ) {
            matrix.zscoreByRows(ZSCORE_LIMIT, useBothNorms);
        } else {
            matrix.zscoreByCols(ZSCORE_LIMIT, useBothNorms);
        }

        if (includeIntra && !appendIntra) {
            matrix.putIntraIntoMainMatrix();
        }
        FinalMatrix result = matrix.getFinalMatrix(useBothNorms, includeIntra && appendIntra);
        result.removeAllNanRows();
        /*
        if (false) {
            SimpleMatrixAndWeight mw = SimilarityMatrixTools.getCompressedCosineSimilarityMatrix(matrix.matrix,
                    50, 0);
            matrix.matrix = mw.matrix;
            matrix.weights = mw.weights;
        } */
        return result;
    }
}
