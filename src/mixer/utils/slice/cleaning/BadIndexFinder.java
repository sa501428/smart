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

package mixer.utils.slice.cleaning;

import javastraw.reader.Dataset;
import javastraw.reader.basics.Chromosome;
import javastraw.reader.basics.ChromosomeHandler;
import javastraw.reader.norm.NormalizationVector;
import javastraw.reader.type.HiCZoom;
import javastraw.reader.type.NormalizationHandler;
import javastraw.reader.type.NormalizationType;
import org.apache.commons.math3.stat.descriptive.DescriptiveStatistics;

import java.util.HashMap;
import java.util.HashSet;
import java.util.Map;
import java.util.Set;

public class BadIndexFinder {
    protected static final float ZSCORE_MIN_NONZERO_COVERAGE = -2;
    protected static final NormalizationType norm = NormalizationHandler.VC;
    private static final float MIN_NORM_VAL = 0.1f;

    public static Map<Integer, Set<Integer>> getBadIndices(Dataset dataset,
                                                           ChromosomeHandler handler, int resolution) {
        Map<Integer, Set<Integer>> badIndices = new HashMap<>();
        Chromosome[] chromosomes = handler.getAutosomalChromosomesArray();
        for (Chromosome chromosome : chromosomes) {
            badIndices.put(chromosome.getIndex(),
                    updateCoverageStats(dataset.getNormalizationVector(chromosome.getIndex(),
                            new HiCZoom(resolution), NormalizationHandler.VC)));
        }
        return badIndices;
    }

    private static Set<Integer> updateCoverageStats(NormalizationVector normalizationVector) {
        double[] vector = normalizationVector.getData().getValues().get(0);
        double[] muAndStd = getMuAndStd(vector);

        Set<Integer> values = new HashSet<>();
        for (int i = 0; i < vector.length; i++) {
            if (vector[i] > MIN_NORM_VAL) {
                double val = Math.log(vector[i]);
                double z = (val - muAndStd[0]) / muAndStd[1];
                if (z < ZSCORE_MIN_NONZERO_COVERAGE) {
                    values.add(i);
                }
            } else {
                values.add(i);
            }
        }
        return values;
    }

    private static double[] getMuAndStd(double[] vector) {
        DescriptiveStatistics stats = new DescriptiveStatistics();
        for (double val : vector) {
            if (val > MIN_NORM_VAL) {
                stats.addValue(Math.log(val));
            }
        }
        double mu = stats.getMean();
        double std = stats.getStandardDeviation();
        stats.clear();
        return new double[]{mu, std};
    }
}
