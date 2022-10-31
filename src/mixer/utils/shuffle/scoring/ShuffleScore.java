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

package mixer.utils.shuffle.scoring;

import mixer.utils.matrix.ShuffledIndices;

import java.util.Map;

public abstract class ShuffleScore {
    protected final float[][] matrix;
    protected final boolean useSymmetry;
    protected final Integer[] rBounds, cBounds;
    protected final Integer[] rIDs, cIDs;

    public ShuffleScore(float[][] matrix, ShuffledIndices rBounds, ShuffledIndices cBounds, boolean useSymmetry) {
        this.matrix = matrix;
        this.useSymmetry = useSymmetry;
        this.rBounds = rBounds.boundaries;
        this.cBounds = cBounds.boundaries;
        this.rIDs = rBounds.ids;
        this.cIDs = cBounds.ids;
    }

    public abstract double score();

    protected String getKey(int rI, int cI) {
        if (useSymmetry) {
            int id1 = rIDs[rI];
            int id2 = cIDs[cI];
            if (id1 <= id2) return id1 + "_" + id2;
            return id2 + "_" + id1;
        }
        return rIDs[rI] + "_" + cIDs[cI];
    }

    protected long populateSumMap(Map<String, Double> sumMap, Map<String, Long> numRegionMap) {
        long totalNumElements = 0;
        sumMap.clear();
        numRegionMap.clear();
        for (int rI = 0; rI < rBounds.length - 1; rI++) {
            for (int cI = 0; cI < cBounds.length - 1; cI++) {
                double localSum = 0;
                long localNumElements = 0;
                for (int i = rBounds[rI]; i < rBounds[rI + 1]; i++) {
                    for (int j = cBounds[cI]; j < cBounds[cI + 1]; j++) {
                        localSum += matrix[i][j];
                        localNumElements++;
                    }
                }

                String key = getKey(rI, cI);
                if (sumMap.containsKey(key)) {
                    sumMap.put(key, sumMap.get(key) + localSum);
                    numRegionMap.put(key, numRegionMap.get(key) + localNumElements);
                } else {
                    sumMap.put(key, localSum);
                    numRegionMap.put(key, localNumElements);
                }
                totalNumElements += localNumElements;
            }
        }
        return totalNumElements;
    }
}
