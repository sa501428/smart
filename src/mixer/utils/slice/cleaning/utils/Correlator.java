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

package mixer.utils.slice.cleaning.utils;

import java.util.ArrayList;
import java.util.HashSet;
import java.util.List;
import java.util.Set;

public class Correlator {
    private final float threshold = 0.5f;
    private final int cutoff;
    private final float[][] matrix;
    private final List<VectorGroup> groups = new ArrayList<>();

    public Correlator(float[][] matrix, boolean onlyUseCorr) {
        this.matrix = matrix;
        if (onlyUseCorr) {
            cutoff = 2;
        } else {
            cutoff = 20;
        }
    }

    public Set<Integer> getOutlierIndices() {
        for (int k = 0; k < matrix.length; k++) {
            boolean unassigned = true;
            for (VectorGroup group : groups) {
                if (group.corr(matrix[k]) > threshold) {
                    group.append(k);
                    unassigned = false;
                    break;
                }
            }
            if (unassigned) {
                groups.add(new VectorGroup(k, matrix[k]));
            }
        }

        Set<Integer> outlierIndices = new HashSet<>();
        for (VectorGroup group : groups) {
            if (group.size() < cutoff) {
                outlierIndices.addAll(group.getMembers());
            }
        }

        return outlierIndices;
    }
}
