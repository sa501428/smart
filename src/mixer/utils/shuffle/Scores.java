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
    public static double getVarScore(long[][] areas, double[][] density) {

        long totalArea = getTotalArea(areas);
        double totalWeightedSum = getWeightedSum(density, areas);
        double mu = totalWeightedSum / totalArea;

        double sumOfSquareErr = 0;

        for (int i = 0; i < areas.length; i++) {
            for (int j = 0; j < areas[i].length; j++) {
                double v = density[i][j] - mu;
                sumOfSquareErr += (v * v) * areas[i][j];
            }
        }

        return (float) (sumOfSquareErr / totalArea);
    }

    private static double getWeightedSum(double[][] densityMatrix, long[][] areas) {
        double weightedSum = 0;
        for (int i = 0; i < areas.length; i++) {
            for (int j = 0; j < areas[i].length; j++) {
                weightedSum += densityMatrix[i][j] * areas[i][j];
            }
        }
        return weightedSum;
    }

    private static long getTotalArea(long[][] areas) {
        long totalArea = 0;
        for (int i = 0; i < areas.length; i++) {
            for (int j = 0; j < areas[i].length; j++) {
                totalArea += areas[i][j];
            }
        }
        return totalArea;
    }

    public static double getKLScore(long[][] areas, double[][] density, boolean matrixIsP) {
        long totalArea = getTotalArea(areas);
        double totalWeightedSum = getWeightedSum(density, areas);
        double mu = totalWeightedSum / totalArea;
        double q = mu / totalWeightedSum;

        double klDivergence = 0;
        for (int i = 0; i < areas.length; i++) {
            for (int j = 0; j < areas[i].length; j++) {
                if (areas[i][j] > 0) {
                    double p = density[i][j] / totalWeightedSum;
                    if (matrixIsP) {
                        klDivergence += (p * Math.log(p / q)) * areas[i][j];
                    } else {
                        klDivergence += (q * Math.log(q / p)) * areas[i][j];
                    }
                }
            }
        }

        return (float) klDivergence;
    }


}
