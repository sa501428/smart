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

package mixer.utils.cleaning;

import javastraw.expected.Welford;
import javastraw.expected.Zscore;
import javastraw.reader.Dataset;
import javastraw.reader.basics.Chromosome;
import javastraw.reader.norm.NormalizationVector;
import javastraw.reader.type.HiCZoom;
import javastraw.reader.type.NormalizationHandler;
import javastraw.reader.type.NormalizationType;

import java.util.HashMap;
import java.util.HashSet;
import java.util.Map;
import java.util.Set;

public class BadIndexFinder {
    protected static final NormalizationType VC = NormalizationHandler.VC;
    private static final int ZSCORE_MIN_NONZERO_COVERAGE = -3;
    private static final int ZSCORE_MAX_NONZERO_COVERAGE = 5;

    public static Map<Integer, Set<Integer>> getBadIndices(Dataset dataset, Chromosome[] chromosomes,
                                                           int resolution, NormalizationType norm) {
        Map<Integer, Set<Integer>> badIndices = new HashMap<>();

        for (Chromosome chromosome : chromosomes) {
            Set<Integer> badValues = getBadCoverageRowsFromNorm(dataset, VC, chromosome, resolution);
            badValues.addAll(getBadCoverageRowsFromNorm(dataset, norm, chromosome, resolution));
            badIndices.put(chromosome.getIndex(), badValues);
        }
        return badIndices;
    }

    private static Set<Integer> getBadCoverageRowsFromNorm(Dataset ds, NormalizationType norm, Chromosome chromosome,
                                                           int resolution) {
        double[] vector;
        try {
            NormalizationVector nv = ds.getNormalizationVector(chromosome.getIndex(), new HiCZoom(resolution), norm);
            vector = nv.getData().getValues().get(0);
        } catch (Exception e) {
            e.printStackTrace();
            System.err.println("\n" + norm.getLabel() + " vector missing " + chromosome.getName() + "\n");
            return new HashSet<>();
        }
        if (vector.length < 2) return new HashSet<>();

        Zscore zLog = getZscoreOfLogs(vector);

        int realLength = (int) (1 + (chromosome.getLength() / resolution));
        Set<Integer> badIndices = new HashSet<>();
        for (int i = 0; i < realLength; i++) {
            if (vector[i] > 0) {
                double z = zLog.getZscore(Math.log(vector[i]));
                if (z < ZSCORE_MIN_NONZERO_COVERAGE || z > ZSCORE_MAX_NONZERO_COVERAGE) {
                    badIndices.add(i);
                }
            } else {
                badIndices.add(i);
            }
        }
        return badIndices;
    }

    private static Zscore getZscoreOfLogs(double[] vector) {
        Welford welford = new Welford();
        for (double val : vector) {
            if (val > 0) {
                welford.addValue(Math.log(val));
            }
        }
        return welford.getZscore();
    }
}
