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

package mixer.utils.shuffle.scoring;

import java.util.HashMap;
import java.util.Map;

public class KLDivergenceScoring extends ShuffleScore {
    private final boolean matrixIsP;

    public KLDivergenceScoring(float[][] matrix, Integer[] rowBounds, Integer[] colBounds, boolean matrixIsP) {
        super(matrix, rowBounds, colBounds);
        this.matrixIsP = matrixIsP;
    }

    @Override
    protected float score(Integer[] rBounds, Integer[] cBounds) {
        double sumTotal = 0;
        Map<String, Double> regionalAverage = new HashMap<>();
        for (int rI = 0; rI < rBounds.length - 1; rI++) {
            for (int cI = 0; cI < cBounds.length - 1; cI++) {
                double sum = 0;
                int numVals = 0;
                for (int i = rBounds[rI]; i < rBounds[rI + 1]; i++) {
                    for (int j = cBounds[cI]; j < cBounds[cI + 1]; j++) {
                        sum += matrix[i][j];
                        numVals++;
                    }
                }
                sumTotal += sum;
                regionalAverage.put(getKey(rI, cI), sum / numVals);
            }
        }

        double klDivergence = 0;
        for (int rI = 0; rI < rBounds.length - 1; rI++) {
            for (int cI = 0; cI < cBounds.length - 1; cI++) {
                double q = regionalAverage.get(getKey(rI, cI)) / sumTotal;
                for (int i = rBounds[rI]; i < rBounds[rI + 1]; i++) {
                    for (int j = cBounds[cI]; j < cBounds[cI + 1]; j++) {
                        double p = matrix[i][j] / sumTotal;
                        if (matrixIsP) {
                            klDivergence += p * Math.log(p / q);
                        } else {
                            klDivergence += q * Math.log(q / p);
                        }
                    }
                }
            }
        }

        return (float) klDivergence;
    }

    private String getKey(int rI, int cI) {
        return rI + "_" + cI;
    }
}
