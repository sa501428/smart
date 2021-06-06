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

import javastraw.reader.basics.Chromosome;
import mixer.utils.common.FloatMatrixTools;
import mixer.utils.common.LogTools;
import mixer.utils.slice.matrices.Dimension;

import java.util.ArrayList;
import java.util.List;

public class GMMClusterCaster {
    public static float[][] cast(float[][] interMatrix, Chromosome[] chromosomes,
                                 Dimension dimensions, Dimension compressedDimensions) {
        List<float[][]> output = new ArrayList<>();

        LogTools.simpleLogWithCleanup(interMatrix, Float.NaN);

        for (int i = 0; i < chromosomes.length; i++) {
            System.out.println("chrom " + chromosomes[i].getName());
            float[][] extraction = getDataJustExcludingChromosome(interMatrix, i, compressedDimensions);
            output.add(runGMM(extraction));
        }

        return FloatMatrixTools.concatenateAll(output);
    }

    private static float[][] getDataJustExcludingChromosome(float[][] interMatrix, int k, Dimension compressedDimensions) {
        float[][] subMatrix = new float[interMatrix.length][compressedDimensions.interval[k]];
        for (int i = 0; i < subMatrix.length; i++) {
            System.arraycopy(interMatrix[i], compressedDimensions.offset[k],
                    subMatrix[i], 0, compressedDimensions.interval[k]);
        }
        return subMatrix;
    }

    private static float[][] runGMM(float[][] extraction) {
        for (int z = 10; z > 4; z--) {
            try {
                GaussianMixtureModels gmm = new GaussianMixtureModels(extraction,
                        z, 50, false);
                gmm.fit();
                gmm.trimEmptyClusters();
                return FloatMatrixTools.convert(gmm.nanPredict(extraction, true));
            } catch (Exception e) {
                System.err.println("Z " + z + " - failed");
                System.err.println(e.getLocalizedMessage());
                e.printStackTrace();
            }
        }
        System.exit(73);
        return null;
    }
}
