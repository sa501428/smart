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

package mixer.utils.slice.kmeans;

import javastraw.tools.MatrixTools;

import java.io.File;
import java.util.Arrays;

public class KmeansEvaluator {

    private final double[][] iterToWcssAicBic;
    private static final int NUM_VALUES = 5;
    private static final int K_INDEX = 0;
    private static final int SUM_SQUARES_INDEX = 1;
    private static final int AIC_INDEX = 2;
    private static final int BIC_INDEX = 3;
    private static final int S_INDEX = 4;


    public KmeansEvaluator(int numClusterSizes) {
        iterToWcssAicBic = new double[NUM_VALUES][numClusterSizes];
        for (double[] row : iterToWcssAicBic) {
            Arrays.fill(row, Double.MAX_VALUE);
        }

    }

    public double getWCSS(int index) {
        return iterToWcssAicBic[SUM_SQUARES_INDEX][index];
    }

    public double getSilhouette(int index) {
        return iterToWcssAicBic[S_INDEX][index];
    }

    public void setMseAicBicValues(int z, int numRows, int numColumns, KmeansResult result) {
        int numClusters = result.getNumActualClusters();
        double sumOfSquares = result.getWithinClusterSumOfSquares();
        double silhouette = result.getSilhouette();

        iterToWcssAicBic[K_INDEX][z] = numClusters;
        iterToWcssAicBic[SUM_SQUARES_INDEX][z] = sumOfSquares;
        iterToWcssAicBic[AIC_INDEX][z] = sumOfSquares + 2 * numColumns * numClusters;
        iterToWcssAicBic[BIC_INDEX][z] = sumOfSquares + 0.5 * numColumns * numClusters * Math.log(numRows);
        iterToWcssAicBic[S_INDEX][z] = silhouette;
        System.out.println("Fin Silhouette for c=" + numClusters + " " + silhouette);
    }

    public void export(File outputDirectory, String kstem) {
        String outIterPath = new File(outputDirectory, kstem + "_cluster_size_WCSS_AIC_BIC.npy").getAbsolutePath();
        MatrixTools.saveMatrixTextNumpy(outIterPath, iterToWcssAicBic);
    }
}
