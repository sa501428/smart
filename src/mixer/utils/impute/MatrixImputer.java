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

package mixer.utils.impute;

import javastraw.tools.MatrixTools;
import mixer.utils.common.ArrayTools;

import java.util.Random;

public class MatrixImputer {

    private static final Random generator = new Random();

    public static float[][] imputeUntilNoNansOnlyNN(float[][] data, long seed) {
        generator.setSeed(seed);
        //int[] allnumClusters = ;
        //float[][] imputed = MatrixColImputer.impute(MatrixTools.deepClone(data));
        float[][] imputed = MatrixTools.deepClone(data);
        int numClusters = 1280;
        boolean useKmedians = false;

        /*
        for(int numClusters : new int[]{200, 200, 100, 100, 50, 50, 40, 40, 30, 30, 20, 20, 10, 10}){
            fillInImputedMatrix(imputed, numClusters);
            if(!checkIfHasNan(imputed)) break;
            imputed = FloatMatrixTools.transpose(imputed);
            fillInImputedMatrix(imputed, imputed.length/4);
            imputed = FloatMatrixTools.transpose(imputed);
            if(!checkIfHasNan(imputed)) break;
        }
        */

        while (checkIfHasNan(imputed) && numClusters >= 10) {
            boolean nothingChanged = fillInImputedMatrix(imputed, numClusters, false);
            if (nothingChanged) numClusters /= 2;
        }

        numClusters = 160;
        while (checkIfHasNan(imputed) && numClusters >= 10) {
            boolean nothingChanged = fillInImputedMatrix(imputed, numClusters, true);
            if (nothingChanged) numClusters /= 2;
        }

        return imputed;
    }

    private static boolean checkIfHasNan(float[][] data) {
        for (int i = 0; i < data.length; i++) {
            for (int j = 0; j < data[i].length; j++) {
                if (Float.isNaN(data[i][j])) {
                    return true;
                }
            }
        }
        return false;
    }

    private static boolean fillInImputedMatrix(float[][] imputed, int numClusters, boolean useKmedians) {
        System.out.println("Imputing... (k=" + numClusters + ")");
        int cont1 = numNans("Starting Nans: ", imputed);
        CentroidImputer.updateBasedOnCentroids(imputed, numClusters, generator, useKmedians);
        System.out.println(".");
        int cont2 = numNans("Done imputing; Nans: ", imputed);
        return cont1 == cont2;
    }

    private static int numNans(String message, float[][] data) {
        int counter = 0;
        int[] numNanRows = new int[data.length];
        int[] numNanCols = new int[data[0].length];

        for (int i = 0; i < data.length; i++) {
            for (int j = 0; j < data[i].length; j++) {
                if (Float.isNaN(data[i][j])) {
                    counter++;
                    numNanRows[i]++;
                    numNanCols[j]++;
                }
            }
        }
        int maxNanInRows = ArrayTools.max(numNanRows);
        int maxNanInCols = ArrayTools.max(numNanCols);
        int avgNanInRows = ArrayTools.mean(numNanRows);
        int avgNanInCols = ArrayTools.mean(numNanCols);

        System.out.println(message + counter + "/(" + maxNanInRows + " : " + avgNanInRows + ")/(" +
                maxNanInCols + " : " + avgNanInCols + ")");
        return counter;
    }
}
