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

package mixer.utils.slice.gmm;

import org.apache.commons.math3.linear.CholeskyDecomposition;
import org.apache.commons.math3.linear.RealMatrix;

public class CovarianceMatrixInverseAndDeterminant {
    public double determinant;
    public RealMatrix inverse;

    private CovarianceMatrixInverseAndDeterminant(RealMatrix covMatrix) {
        CholeskyDecomposition chol = new CholeskyDecomposition(covMatrix);
        determinant = chol.getDeterminant();
        inverse = chol.getSolver().getInverse();
    }

    public static CovarianceMatrixInverseAndDeterminant[] convertToArray(RealMatrix[] covMatrices) {
        CovarianceMatrixInverseAndDeterminant[] covs = new CovarianceMatrixInverseAndDeterminant[covMatrices.length];
        for (int k = 0; k < covs.length; k++) {
            covs[k] = new CovarianceMatrixInverseAndDeterminant(covMatrices[k]);
        }
        return covs;
    }
}
