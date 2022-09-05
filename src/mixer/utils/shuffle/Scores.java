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

public class Scores {
    public static double getVarScore(long[][][] areas, double[][][] counts, boolean isBaseline) {

        long totalArea = getTotalArea(areas);
        double totalWeightedSum = getWeightedSum(counts, areas);
        double[][] mu = getMean(isBaseline, totalWeightedSum, totalArea, areas, counts);

        double sumOfSquareErr = 0;

        for (int i = 0; i < areas.length; i++) {
            for (int j = 0; j < areas[i].length; j++) {
                for (int k = 0; k < areas[i][j].length; k++) {
                    double v = counts[i][j][k] - mu[j][k];
                    sumOfSquareErr += (v * v) * areas[i][j][k];
                }
            }
        }

        return (float) (sumOfSquareErr / totalArea);
    }

    public static double getKLScore(long[][][] areas, double[][][] counts, boolean matrixIsP, boolean isBaseline) {
        long totalArea = getTotalArea(areas);
        double totalWeightedSum = getWeightedSum(counts, areas);
        double[][] mu = getMean(isBaseline, totalWeightedSum, totalArea, areas, counts);

        double klDivergence = 0;
        for (int i = 0; i < areas.length; i++) {
            for (int j = 0; j < areas[i].length; j++) {
                for (int k = 0; k < areas[i][j].length; k++) {
                    if (areas[i][j][k] > 0) {
                        double q = mu[j][k] / totalWeightedSum;
                        double p = counts[i][j][k] / totalWeightedSum;
                        if (matrixIsP) {
                            klDivergence += (p * Math.log(p / q)) * areas[i][j][k];
                        } else {
                            klDivergence += (q * Math.log(q / p)) * areas[i][j][k];
                        }
                    }
                }
            }
        }

        return klDivergence;
    }

    private static double[][] getMean(boolean isBaseline, double totalWeightedSum, long totalArea,
                                      long[][][] areas, double[][][] counts) {
        int n = areas[0].length;
        double[][] result = new double[n][n];
        if (isBaseline) {
            for (int i = 0; i < n; i++) {
                for (int j = 0; j < n; j++) {
                    result[i][j] = totalWeightedSum / totalArea;
                }
            }
        } else {
            double[][] denom = new double[n][n];

            for (int m = 0; m < areas.length; m++) {
                for (int i = 0; i < n; i++) {
                    for (int j = 0; j < n; j++) {
                        result[i][j] += counts[m][i][j] * areas[m][i][j];
                        denom[i][j] += areas[m][i][j];
                    }
                }
            }

            for (int i = 0; i < n; i++) {
                for (int j = 0; j < n; j++) {
                    result[i][j] /= denom[i][j];
                }
            }
        }
        return result;
    }

    private static double weightedAverage(double d1, long a1, double d2, long a2) {
        if ((a1 + a2) > 0) return (d1 * a1 + d2 * a2) / (a1 + a2);
        return 0;
    }


    private static double getWeightedSum(double[][][] densityMatrix, long[][][] areas) {
        double weightedSum = 0;
        for (int i = 0; i < areas.length; i++) {
            for (int j = 0; j < areas[i].length; j++) {
                for (int k = 0; k < areas[i][j].length; k++) {
                    weightedSum += densityMatrix[i][j][k] * areas[i][j][k];
                }
            }
        }
        return weightedSum;
    }

    private static long getTotalArea(long[][][] areas) {
        long totalArea = 0;
        for (long[][] area : areas) {
            for (long[] longs : area) {
                for (long aLong : longs) {
                    totalArea += aLong;
                }
            }
        }
        return totalArea;
    }
}
