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

package mixer.utils.aba;

import javastraw.tools.HiCFileTools;
import javastraw.tools.MatrixTools;
import mixer.MixerGlobals;
import org.apache.commons.math.linear.RealMatrix;

import java.io.File;
import java.text.SimpleDateFormat;
import java.util.Date;

public class ABADataStack {

    private static final Object key = new Object();
    // genome wide variables
    private static boolean genomeWideVariablesNotSet = true;
    private static RealMatrix gwAggregateMatrix;
    private static RealMatrix gwAggSupDiagNormedMatrix;
    private static RealMatrix gwAggDiagNormedMatrix;
    // saving data variables
    private static File dataDirectory;

    // chr variables
    private RealMatrix aggregateMatrix;
    private RealMatrix aggSupDiagNormedMatrix;
    private RealMatrix aggDiagNormedMatrix;

    /**
     * class for saving data from chromosme wide run of ABA, keeps static class to store genome wide data
     *
     * @param n            width of matrix
     * @param outputFolder location for saving data
     * @param customPrefix optional file/folder prefix
     */
    public ABADataStack(int n, File outputFolder, String customPrefix) {
        aggregateMatrix = MatrixTools.cleanArray2DMatrix(n, n);
        aggSupDiagNormedMatrix = MatrixTools.cleanArray2DMatrix(n, n);
        aggDiagNormedMatrix = MatrixTools.cleanArray2DMatrix(n, n);

        initializeGenomeWideVariables(n);
        initializeDataSaveFolder(outputFolder, customPrefix);
    }

    /**
     * Ensure that directory for saving exists
     *
     * @param outputFolderDirectory to directory
     * @param prefix                of files to be saved
     */
    public static void initializeDataSaveFolder(File outputFolderDirectory, String prefix) {
        if (prefix.length() < 1) {// no preference specified
            dataDirectory = new File(outputFolderDirectory,
                    new SimpleDateFormat("yyyy.MM.dd.HH.mm").format(new Date()));
        } else {
            dataDirectory = new File(outputFolderDirectory, prefix);
        }
        dataDirectory = HiCFileTools.createValidDirectory(dataDirectory.getAbsolutePath());
    }

    private static void initializeGenomeWideVariables(int n) {
        if (genomeWideVariablesNotSet) {
            gwAggregateMatrix = MatrixTools.cleanArray2DMatrix(n, n);
            gwAggSupDiagNormedMatrix = MatrixTools.cleanArray2DMatrix(n, n);
            gwAggDiagNormedMatrix = MatrixTools.cleanArray2DMatrix(n, n);
            genomeWideVariablesNotSet = false;
        }
    }

    public static void exportGenomeWideData(Integer numRegions, boolean saveAllData) {
        double gwNPeaksUsedInv = 1. / numRegions;
        gwAggSupDiagNormedMatrix = gwAggSupDiagNormedMatrix.scalarMultiply(gwNPeaksUsedInv);
        gwAggDiagNormedMatrix = gwAggDiagNormedMatrix.scalarMultiply(gwNPeaksUsedInv);

        RealMatrix[] matrices = {gwAggregateMatrix, gwAggSupDiagNormedMatrix, gwAggDiagNormedMatrix};
        String[] titles = {"ABA", "supDiagNormedABA", "diagNormedABA"};

        saveDataSet("gw", matrices, titles);
    }

    private static void saveDataSet(String prefix,
                                    RealMatrix[] matrices,
                                    String[] titles) {

        File subFolder = HiCFileTools.createValidDirectory(new File(dataDirectory, prefix).getAbsolutePath());
        if (MixerGlobals.printVerboseComments) {
            System.out.println("Saving chr " + prefix + " data to " + subFolder);
        }

        for (int i = 0; i < matrices.length; i++) {
            MatrixTools.saveMatrixTextNumpy(
                    (new File(subFolder, titles[i] + ".npy")).getAbsolutePath(),
                    matrices[i].getData());
        }
    }

    public static void clearAllData() {
        dataDirectory = null;
        genomeWideVariablesNotSet = true;
        gwAggregateMatrix = null;
        gwAggSupDiagNormedMatrix = null;
        gwAggDiagNormedMatrix = null;
    }

    public void addData(RealMatrix newData) {
        MatrixTools.cleanUpNaNs(newData);
        aggregateMatrix = aggregateMatrix.add(newData);
        aggSupDiagNormedMatrix = aggSupDiagNormedMatrix.add(supDiagNormalization(newData));
        aggDiagNormedMatrix = aggDiagNormedMatrix.add(diagNormalization(newData));
    }

    public synchronized void updateGenomeWideData() {
        synchronized (key) {
            gwAggregateMatrix = gwAggregateMatrix.add(aggregateMatrix);
            gwAggSupDiagNormedMatrix = gwAggSupDiagNormedMatrix.add(aggSupDiagNormedMatrix);
            gwAggDiagNormedMatrix = gwAggDiagNormedMatrix.add(aggDiagNormedMatrix);
        }
    }

    private RealMatrix diagNormalization(RealMatrix matrix) {
        return matrix.copy().scalarMultiply(1. / Math.max(1., diagonalMean(matrix.getData())));
    }

    private RealMatrix supDiagNormalization(RealMatrix matrix) {
        return matrix.copy().scalarMultiply(1. / Math.max(1., supDiagonalMean(matrix.getData())));
    }

    private double supDiagonalMean(double[][] data) {
        int counter = 0;
        double accum = 0;
        for (int i = 0; i < data.length - 1; i++) {
            int j = i + 1;
            if (data[i][j] > 0) {
                accum += data[i][j];
                counter++;
            }
        }
        return accum / Math.max(counter, 1);
    }

    private double diagonalMean(double[][] data) {
        int counter = 0;
        double accum = 0;
        for (int i = 0; i < data.length; i++) {
            if (data[i][i] > 0) {
                accum += data[i][i];
                counter++;
            }
        }
        return accum / Math.max(counter, 1);
    }
}
