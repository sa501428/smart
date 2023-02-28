/*
 * The MIT License (MIT)
 *
 * Copyright (c) 2011-2023 Rice University, Baylor College of Medicine, Aiden Lab
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

package mixer.utils.eig;

import javastraw.feature1D.GenomeWide1DList;
import javastraw.reader.basics.ChromosomeHandler;
import mixer.utils.cleaning.SimilarityMatrixTools;
import mixer.utils.drive.FinalMatrix;
import mixer.utils.similaritymeasures.RobustCorrelationSimilarity;
import mixer.utils.tracks.EigenvectorInterval;
import org.apache.commons.math3.linear.EigenDecomposition;
import org.apache.commons.math3.linear.MatrixUtils;
import org.apache.commons.math3.linear.RealMatrix;

import java.io.File;
import java.util.List;
import java.util.Map;

public class EigenvectorTools {
    public static void runEigenvectorAnalysis(String prefix, Map<Integer, List<String>> bedFiles,
                                              ChromosomeHandler handler, FinalMatrix matrix) {
        float[] eigenvector = run(matrix);
        GenomeWide1DList<EigenvectorInterval> result = matrix.processEigenvectorResult(eigenvector, handler);
        File outBedFile = new File("EIG_" + prefix + ".bedgraph");
        result.simpleExport(outBedFile);
    }

    public static float[] run(FinalMatrix matrix) {
        float[][] array = SimilarityMatrixTools.getSymmetricDistanceMatrix(matrix.matrix,
                RobustCorrelationSimilarity.SINGLETON);
        RealMatrix realMatrix = MatrixUtils.createRealMatrix(convert(array));
        //RealMatrix covarianceMatrix = new Covariance(realMatrix).getCovarianceMatrix();
        EigenDecomposition ed = new EigenDecomposition(realMatrix);
        return convert(ed.getEigenvector(0).toArray());
    }

    public static float[] convert(double[] row) {
        float[] result = new float[row.length];
        for (int i = 0; i < result.length; i++) {
            result[i] = (float) row[i];
        }
        return result;
    }

    public static double[][] convert(float[][] matrix) {
        double[][] result = new double[matrix.length][matrix[0].length];
        for (int i = 0; i < result.length; i++) {
            for (int j = 0; j < result[i].length; j++) {
                result[i][j] = matrix[i][j];
            }
        }
        return result;
    }
}
