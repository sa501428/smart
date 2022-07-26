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

import javastraw.feature1D.GenomeWide1DList;
import javastraw.reader.basics.ChromosomeHandler;
import mixer.utils.slice.matrices.CompositeGenomeWideMatrix;
import mixer.utils.slice.structures.SubcompartmentInterval;

import java.util.List;
import java.util.Map;

public class GenomeWideGMMRunner {
    private final CompositeGenomeWideMatrix interMatrix;
    private final float[][] data;
    private final ChromosomeHandler chromosomeHandler;

    public GenomeWideGMMRunner(ChromosomeHandler chromosomeHandler, CompositeGenomeWideMatrix interMatrix) {
        this.interMatrix = interMatrix;
        data = interMatrix.getData(true);
        this.chromosomeHandler = chromosomeHandler;
    }

    public void launch(int numClusters, long seed, Map<Integer, GenomeWide1DList<SubcompartmentInterval>> results,
                       List<List<Integer>> startingIndices, boolean useRobustGMM) {
        if (data.length > 0 && data[0].length > 0) {
            GaussianMixtureModels gmm = new GaussianMixtureModels(data,
                    numClusters, 100, startingIndices, useRobustGMM);
            try {
                gmm.fit();
                int[] result = gmm.predict();
                populateMap(result, results, numClusters);
            } catch (Exception g) {
                System.err.println("Unable to run GMM with " + numClusters + " clusters on data.\n" +
                        "Matrix may be singular or some other issue encountered");
            }
        }
    }

    private void populateMap(int[] result, Map<Integer, GenomeWide1DList<SubcompartmentInterval>> results, int numClusters) {
        GenomeWide1DList<SubcompartmentInterval> finalCompartments = new GenomeWide1DList<>(chromosomeHandler);
        getCounts(result, numClusters);
        interMatrix.processGMMClusteringResult(result, finalCompartments);
        results.put(numClusters, finalCompartments);
    }

    private int[] getCounts(int[] result, int numClusters) {
        int[] counts = new int[numClusters];
        for (int val : result) {
            counts[val]++;
        }
        for (int count : counts) {
            if (count < 5) {
                System.err.println("GMM clustering may not have worked for " + numClusters + " clusters; count: " + count);
            }
        }
        return counts;
    }


}
