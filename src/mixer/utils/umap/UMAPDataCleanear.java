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

package mixer.utils.umap;

import java.util.Arrays;

public class UMAPDataCleanear {

    final int width = 2000;
    final int height;
    final float[][] points;
    final int[][] ids;
    final int numPlots;
    private final boolean useOnlyFullAgreement = false;

    public UMAPDataCleanear(float[][] initialPoints, int[][] initialIDs) {
        numPlots = initialIDs[0].length;
        boolean[] shouldUse = getUseStatus(initialIDs);
        int realLength = getRealLength(shouldUse);

        // initial bounds
        float minX = getMin(initialPoints, shouldUse, 0);
        float maxX = getMax(initialPoints, shouldUse, 0, minX);
        float minY = getMin(initialPoints, shouldUse, 1);
        float maxY = getMax(initialPoints, shouldUse, 1, minY);

        // expand window a little bit
        float widthX0 = (maxX - minX);
        float widthY0 = (maxY - minY);
        minX -= 0.025 * widthX0;
        maxX += 0.025 * widthX0;
        minY -= 0.025 * widthY0;
        maxY += 0.025 * widthY0;
        float widthX = (maxX - minX);
        float widthY = (maxY - minY);

        // clean up data and transform points to fit image
        // maintain aspect ratio; rotate if need be
        points = new float[realLength][2];
        ids = getRealIDs(realLength, initialIDs, shouldUse);

        int counter = 0;
        if (widthX >= widthY) {
            height = (int) Math.ceil((widthY / widthX) * width);
            for (int k = 0; k < initialPoints.length; k++) {
                if (shouldUse[k]) {
                    points[counter][0] = transform(initialPoints[k][0], minX, maxX, width);
                    points[counter][1] = transform(initialPoints[k][1], minY, maxY, height);
                    counter++;
                }
            }
        } else {
            height = (int) Math.ceil((widthX / widthY) * width);
            for (int k = 0; k < initialPoints.length; k++) {
                if (shouldUse[k]) {
                    points[counter][0] = transform(initialPoints[k][1], minY, maxY, width);
                    points[counter][1] = transform(initialPoints[k][0], minX, maxX, height);
                    counter++;
                }
            }
        }
    }

    private int[][] getRealIDs(int realLength, int[][] initialIDs, boolean[] shouldUse) {
        int[][] realIDs = new int[realLength][numPlots];
        int counter = 0;
        for (int k = 0; k < shouldUse.length; k++) {
            if (shouldUse[k]) {
                for (int z = 0; z < numPlots; z++) {
                    realIDs[counter][z] = initialIDs[k][z];
                }
                counter++;
            }
        }
        return realIDs;
    }

    private float transform(float val, float minVal, float maxVal, float size) {
        return size * (val - minVal) / (maxVal - minVal);
    }

    private int getRealLength(boolean[] shouldUse) {
        int counter = 0;
        for (boolean val : shouldUse) {
            if (val) {
                counter++;
            }
        }
        return counter;
    }

    private float getMin(float[][] initialPoints, boolean[] shouldUse, int i) {
        float minVal = Float.MAX_VALUE;
        for (int k = 0; k < initialPoints.length; k++) {
            if (shouldUse[k]) {
                minVal = Math.min(minVal, initialPoints[k][i]);
            }
        }
        return minVal;
    }

    private float getMax(float[][] initialPoints, boolean[] shouldUse, int i, float startVal) {
        float maxVal = startVal;
        for (int k = 0; k < initialPoints.length; k++) {
            if (shouldUse[k]) {
                maxVal = Math.max(maxVal, initialPoints[k][i]);
            }
        }
        return maxVal;
    }

    private boolean[] getUseStatus(int[][] ids) {
        int[] statusVal = new int[ids.length];
        Arrays.fill(statusVal, 1);
        for (int i = 0; i < ids.length; i++) {
            for (int j = 0; j < ids[i].length; j++) {
                statusVal[i] *= ids[i][j];
            }
        }

        boolean[] status = new boolean[ids.length];
        for (int k = 0; k < statusVal.length; k++) {
            status[k] = statusVal[k] > 0;
        }

        if (useOnlyFullAgreement) {
            for (int k = 0; k < statusVal.length; k++) {
                int baseline = ids[k][0];
                boolean allIdentical = true;
                for (int z = 1; z < numPlots; z++) {
                    allIdentical &= baseline == ids[k][z];
                }
                status[k] &= allIdentical;
            }
        }

        return status;
    }
}
