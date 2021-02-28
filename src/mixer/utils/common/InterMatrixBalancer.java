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

package mixer.utils.common;

import org.apache.commons.math.stat.StatUtils;

import java.util.Arrays;

public class InterMatrixBalancer {
    private final float[][] matrix;
    private final int totSize;
    private final int colOffset;
    private final boolean isSymmetric;
    private float[] norm = null;

    public InterMatrixBalancer(float[][] matrix, boolean isSymmetric) {
        this.matrix = matrix;
        this.isSymmetric = isSymmetric;
        if (isSymmetric) {
            totSize = matrix.length;
            colOffset = 0;
        } else {
            totSize = matrix.length + matrix[0].length;
            colOffset = matrix.length;
        }
    }

    public float[][] getNormalizedMatrix() {
        if (norm == null) {
            norm = computeKR();
        }

        if (norm != null) {
            //norm = normalizeVectorByScaleFactor(norm);
            float[][] normedMatrix = normalizeMatrix();
            return normedMatrix;
        }
        return null;
    }

    private float[][] normalizeMatrix() {
        float[][] result = new float[matrix.length][matrix[0].length];

        for (int i = 0; i < matrix.length; i++) {
            for (int j = 0; j < matrix[i].length; j++) {
                result[i][j] = matrix[i][j] / (norm[i] * norm[j + colOffset]);
            }
        }

        return result;
    }


    private double[] computeKRNormVector(int[] offset, double[] x0) {

        double tol = 0.000001;
        double delta = 0.1;
        int n = x0.length;
        double[] e = new double[n];
        Arrays.fill(e, 1);

        double g = 0.9;
        double etamax = 0.1;
        double eta = etamax;

        double rt = Math.pow(tol, 2);

        double[] v = sparseMultiplyFromContactRecords(offset, x0);
        double[] rk = new double[v.length];
        for (int i = 0; i < v.length; i++) {
            v[i] *= x0[i];
            rk[i] = 1 - v[i];
        }
        double rho_km1 = 0;
        for (double aRk : rk) {
            rho_km1 += aRk * aRk;
        }
        double rout = rho_km1;
        double rold = rout;
        int MVP = 0;  // We'll count matrix vector products.

        int not_changing = 0;
        while (rout > rt && not_changing < 100) {    // Outer iteration
            int k = 0;
            double[] y = deepClone(e);
            double[] ynew = new double[e.length];
            double[] Z = new double[e.length];
            double[] p = new double[e.length];
            double[] w = new double[e.length];
            double alpha;
            double beta;
            double gamma;
            double rho_km2 = rho_km1;


            double innertol = Math.max(Math.pow(eta, 2) * rout, rt);
            while (rho_km1 > innertol) {   // Inner iteration by CG
                k++;

                if (k == 1) {
                    rho_km1 = 0;
                    for (int i = 0; i < Z.length; i++) {
                        double rkVal = rk[i];
                        double zVal = rkVal / v[i];
                        Z[i] = zVal;
                        rho_km1 += rkVal * zVal;
                    }
                    p = deepClone(Z);

                } else {
                    beta = rho_km1 / rho_km2;
                    for (int q = 0; q < p.length; q++) {
                        p[q] *= beta;
                    }
                    for (int i = 0; i < p.length; i++) {
                        p[i] += Z[i];
                    }
                }
                double[] tmp = new double[e.length];
                for (int i = 0; i < tmp.length; i++) {
                    tmp[i] = x0[i] * p[i];
                }
                tmp = sparseMultiplyFromContactRecords(offset, tmp);
                alpha = 0;
                // Update search direction efficiently.
                for (int i = 0; i < tmp.length; i++) {
                    double pVal = p[i];
                    double wVal = (x0[i] * tmp[i] + v[i] * pVal);
                    w[i] = wVal;
                    alpha += pVal * wVal;
                }
                alpha = rho_km1 / alpha;
                double minynew = Double.MAX_VALUE;
                // Test distance to boundary of cone.
                for (int i = 0; i < p.length; i++) {
                    double yVal = y[i] + alpha * p[i];
                    ynew[i] = yVal;
                    if (yVal < minynew) {
                        minynew = yVal;
                    }
                }
                if (minynew <= delta) {
                    if (delta == 0) break;     // break out of inner loop?
                    gamma = Double.MAX_VALUE;
                    for (int i = 0; i < ynew.length; i++) {
                        double pVal = p[i];
                        if (alpha * pVal < 0) {
                            double yVal = y[i];
                            if ((delta - yVal) / (alpha * pVal) < gamma) {
                                gamma = ((delta - yVal) / (alpha * pVal));
                            }
                        }
                    }
                    for (int i = 0; i < y.length; i++) {
                        y[i] += gamma * alpha * p[i];
                    }
                    break;   // break out of inner loop?
                }
                rho_km2 = rho_km1;
                rho_km1 = 0;
                y = deepClone(ynew);
                for (int i = 0; i < y.length; i++) {
                    rk[i] += -alpha * w[i];
                    double rkVal = rk[i];
                    Z[i] = rkVal / v[i];
                    rho_km1 += rkVal * Z[i];
                }

            } // end inner loop
            for (int i = 0; i < x0.length; i++) {
                x0[i] *= y[i];
            }
            v = sparseMultiplyFromContactRecords(offset, x0);
            rho_km1 = 0;
            for (int i = 0; i < v.length; i++) {
                v[i] *= x0[i];
                double rkVal = 1 - v[i];
                rk[i] = rkVal;

                rho_km1 += rkVal * rkVal;
            }
            if (Math.abs(rho_km1 - rout) < 0.000001 || Double.isInfinite(rho_km1)) {
                not_changing++;
            }
            rout = rho_km1;
            MVP = MVP + k + 1;
            //  Update inner iteration stopping criterion.
            double rat = rout / rold;
            rold = rout;
            double r_norm = Math.sqrt(rout);
            double eta_o = eta;
            eta = g * rat;
            if (g * Math.pow(eta_o, 2) > 0.1) {
                eta = Math.max(eta, g * Math.pow(eta_o, 2));
            }
            eta = Math.max(Math.min(eta, etamax), 0.5 * tol / r_norm);
        }
        if (not_changing >= 100) {
            return null;
        }
        return x0;
    }

    private double[] deepClone(double[] e) {
        double[] clone = new double[e.length];
        System.arraycopy(e, 0, clone, 0, e.length);
        return clone;
    }

    private double[] sparseMultiplyFromContactRecords(int[] offset, double[] vector) {
        double[] result = new double[vector.length];

        for (int i = 0; i < matrix.length; i++) {
            int row0 = i;
            for (int j = 0; j < matrix[i].length; j++) {
                int col0 = j + colOffset;
                float value = matrix[i][j];

                int row = offset[row0];
                int col = offset[col0];

                if (row != -1 && col != -1) {
                    result[row] += vector[col] * value;
                    if (row != col) {
                        result[col] += vector[row] * value;
                    }
                }
            }
        }

        return result;
    }


    private float[] computeKR() {

        boolean recalculate = true;
        int[] offset = getOffset(0);
        float[] kr = null;
        int iteration = 1;

        while (recalculate && iteration <= 20) {
            // create new matrix indices upon every iteration, because we've thrown out rows
            // newSize is size of new sparse matrix (non-sparse rows)
            int newSize = 0;
            for (int offset1 : offset) {
                if (offset1 != -1) newSize++;
            }

            // initialize x0 for call the compute KR norm
            double[] x0 = new double[newSize];
            Arrays.fill(x0, 1);

            x0 = computeKRNormVector(offset, x0);

            // assume all went well and we don't need to recalculate
            recalculate = false;
            int rowsTossed = 0;

            if (x0 == null || iteration == 5) {
                // if x0 is no good, throw out some percentage of rows and reset the offset array that gives those rows
                recalculate = true;
                if (iteration < 5) {
                    offset = getOffset(iteration);
                } else {
                    offset = getOffset(10);
                }
                //   System.out.print(" " + iteration + "%");
            } else {
                // otherwise, check to be sure there are no tiny KR values
                // create true KR vector
                kr = new float[totSize];
                int krIndex = 0;
                for (int offset1 : offset) {
                    if (offset1 == -1) {
                        kr[krIndex++] = Float.NaN;
                    } else {
                        kr[krIndex++] = (float) (1.0f / x0[offset1]);
                    }
                }
                // find scaling factor
                double mySum = getSumFactor(kr);

                // if any values are too small, recalculate.  set those rows to be thrown out and reset the offset
                // note that if no rows are thrown out, the offset should not change
                int index = 0;
                for (int i = 0; i < kr.length; i++) {
                    if (kr[i] * mySum < 0.01) {
                        offset[i] = -1;
                        rowsTossed++;
                        recalculate = true;
                    } else {
                        if (offset[i] != -1) {
                            offset[i] = index++;
                        }
                    }
                }
                // if (recalculate) System.out.print(" " + rowsTossed);
            }
            iteration++;
            System.gc();
        }
        if (iteration > 6 && recalculate) {
            return null; // error did not converge
        }

        return kr;
    }

    private int[] getOffset(double percent) {
        double[] rowSums = new double[totSize];

        for (int i = 0; i < matrix.length; i++) {
            int x = i;
            for (int j = 0; j < matrix[i].length; j++) {
                int y = j + colOffset;
                float value = matrix[i][j];

                rowSums[x] += value;
                if (x != y) {
                    rowSums[y] += value;
                }
            }
        }

        double thresh = 0;
        if (percent > 0) {
            // Get percent threshold from positive row sums (nonzero)
            int j = 0;
            for (double sum : rowSums) {
                if (sum != 0) {
                    j++;
                }
            }
            double[] posRowSums = new double[j];
            j = 0;
            for (double sum : rowSums) {
                if (sum != 0) {
                    posRowSums[j++] = sum;
                }
            }
            thresh = StatUtils.percentile(posRowSums, percent);
        }

        int[] offset = new int[rowSums.length];
        int index = 0;
        for (int i = 0; i < rowSums.length; i++) {
            if (rowSums[i] <= thresh) {
                offset[i] = -1;
            } else {
                offset[i] = index++;
            }
        }

        return offset;

    }

    public double getSumFactor(float[] norm) {
        double[] normMatrixSums = getNormMatrixSumFactor(norm);
        return Math.sqrt(normMatrixSums[0] / normMatrixSums[1]);
    }

    public double[] getNormMatrixSumFactor(float[] norm) {
        double matrix_sum = 0;
        double norm_sum = 0;

        for (int i = 0; i < matrix.length; i++) {
            int x = i;
            for (int j = 0; j < matrix[i].length; j++) {
                int y = j + colOffset;
                float value = matrix[i][j];

                double valX = norm[x];
                double valY = norm[y];
                if (!Double.isNaN(valX) && !Double.isNaN(valY) && valX > 0 && valY > 0) {
                    // want total sum of matrix, not just upper triangle
                    if (x == y) {
                        norm_sum += value / (valX * valY);
                        matrix_sum += value;
                    } else {
                        norm_sum += 2 * value / (valX * valY);
                        matrix_sum += 2 * value;
                    }
                }
            }
        }
        return new double[]{norm_sum, matrix_sum};
    }

    private float[] normalizeVectorByScaleFactor(float[] newNormVector) {

        /*
        if(matrix.length == matrix[0].length && newNormVector.length < 30){
            return scalePerDensestCluster(newNormVector);
        }
        */

        for (int k = 0; k < newNormVector.length; k++) {
            float kVal = newNormVector[k];
            if (kVal <= 0 || Double.isNaN(kVal)) {
                newNormVector[k] = Float.NaN;
            } else {
                newNormVector[k] = 1.f / kVal;
            }
        }

        double normalizedSumTotal = 0, sumTotal = 0;

        for (int i = 0; i < matrix.length; i++) {
            int x = i;
            for (int j = 0; j < matrix[i].length; j++) {
                int y = j + colOffset;
                float value = matrix[i][j];

                double valX = newNormVector[x];
                double valY = newNormVector[y];

                if (!Double.isNaN(valX) && !Double.isNaN(valY)) {
                    double normalizedValue = value / (valX * valY);
                    normalizedSumTotal += normalizedValue;
                    sumTotal += value;
                    if (x != y) {
                        normalizedSumTotal += normalizedValue;
                        sumTotal += value;
                    }
                }
            }
        }

        return scaleVectorBy(newNormVector, Math.sqrt(normalizedSumTotal / sumTotal));
    }

    private float[] scalePerDensestCluster(float[] newNormVector) {
        int maxValIndex = 0;
        for (int i = 1; i < matrix.length; i++) {
            if (matrix[i][i] > matrix[maxValIndex][maxValIndex]) {
                maxValIndex = i;
            }
        }

        float denom = newNormVector[maxValIndex] * newNormVector[colOffset + maxValIndex];
        return scaleVectorBy(newNormVector, 1.0 / Math.sqrt(denom));
    }

    private float[] scaleVectorBy(float[] newNormVector, double scaleFactor) {
        for (int q = 0; q < newNormVector.length; q++) {
            newNormVector[q] *= scaleFactor;
        }
        return newNormVector;
    }

}
