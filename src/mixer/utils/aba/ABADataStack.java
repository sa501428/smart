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
    private static RealMatrix gwAPAMatrix;
    private static RealMatrix gwMeanNormedAPAMatrix;
    private static RealMatrix gwMedianNormedAPAMatrix;
    // saving data variables
    private static File dataDirectory;

    // chr variables
    private RealMatrix APAMatrix;
    private RealMatrix meanNormedAPAMatrix;
    private RealMatrix medianNormedAPAMatrix;

    /**
     * class for saving data from chromosme wide run of APA, keeps static class to store genomide data
     *
     * @param n            width of matrix
     * @param outputFolder location for saving data
     * @param customPrefix optional file/folder prefix
     */
    public ABADataStack(int n, File outputFolder, String customPrefix) {
        APAMatrix = MatrixTools.cleanArray2DMatrix(n, n);
        meanNormedAPAMatrix = MatrixTools.cleanArray2DMatrix(n, n);
        medianNormedAPAMatrix = MatrixTools.cleanArray2DMatrix(n, n);

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
            gwAPAMatrix = MatrixTools.cleanArray2DMatrix(n, n);
            gwMeanNormedAPAMatrix = MatrixTools.cleanArray2DMatrix(n, n);
            gwMedianNormedAPAMatrix = MatrixTools.cleanArray2DMatrix(n, n);
            genomeWideVariablesNotSet = false;
        }
    }

    public static void exportGenomeWideData(Integer numRegions, boolean saveAllData) {
        double gwNPeaksUsedInv = 1. / numRegions;
        gwMeanNormedAPAMatrix = gwMeanNormedAPAMatrix.scalarMultiply(gwNPeaksUsedInv);
        gwMedianNormedAPAMatrix = gwMedianNormedAPAMatrix.scalarMultiply(gwNPeaksUsedInv);

        RealMatrix[] matrices = {gwAPAMatrix, gwMeanNormedAPAMatrix, gwMedianNormedAPAMatrix};
        String[] titles = {"APA", "normedAPA", "centerNormedAPA"};

        saveDataSet("gw", matrices, titles);
    }

    private static void saveDataSet(String prefix,
                                    RealMatrix[] apaMatrices,
                                    String[] apaDataTitles) {

        File subFolder = HiCFileTools.createValidDirectory(new File(dataDirectory, prefix).getAbsolutePath());
        if (MixerGlobals.printVerboseComments) {
            System.out.println("Saving chr " + prefix + " data to " + subFolder);
        }

        for (int i = 0; i < apaMatrices.length; i++) {
            MatrixTools.saveMatrixTextNumpy(
                    (new File(subFolder, apaDataTitles[i] + ".npy")).getAbsolutePath(),
                    apaMatrices[i].getData());
        }
    }

    public static void clearAllData() {
        dataDirectory = null;
        genomeWideVariablesNotSet = true;
        gwAPAMatrix = null;
        gwMeanNormedAPAMatrix = null;
        gwMedianNormedAPAMatrix = null;
    }

    public void addData(RealMatrix newData) {
        MatrixTools.cleanUpNaNs(newData);
        APAMatrix = APAMatrix.add(newData);
        meanNormedAPAMatrix = meanNormedAPAMatrix.add(standardNormalization(newData));
        medianNormedAPAMatrix = medianNormedAPAMatrix.add(centerNormalization(newData));
    }

    private RealMatrix centerNormalization(RealMatrix newData) {
        return newData;
    }

    private RealMatrix standardNormalization(RealMatrix newData) {
        return newData;
    }

    public synchronized void updateGenomeWideData() {
        synchronized (key) {
            gwAPAMatrix = gwAPAMatrix.add(APAMatrix);
            gwMeanNormedAPAMatrix = gwMeanNormedAPAMatrix.add(meanNormedAPAMatrix);
            gwMedianNormedAPAMatrix = gwMedianNormedAPAMatrix.add(medianNormedAPAMatrix);
        }
    }

    public void exportDataSet(String subFolderName, Integer[] peakNumbers) {
        double nPeaksUsedInv = 1. / peakNumbers[0];
        meanNormedAPAMatrix = meanNormedAPAMatrix.scalarMultiply(nPeaksUsedInv);
        medianNormedAPAMatrix = medianNormedAPAMatrix.scalarMultiply(nPeaksUsedInv);

        RealMatrix[] matrices = {APAMatrix, meanNormedAPAMatrix, medianNormedAPAMatrix};
        String[] titles = {"APA", "meanNormedAPA", "medianNormedAPA"};

        saveDataSet(subFolderName, matrices, titles);
    }

    public void thresholdPlots(int val) {
        MatrixTools.thresholdValues(APAMatrix, val);
    }
}
