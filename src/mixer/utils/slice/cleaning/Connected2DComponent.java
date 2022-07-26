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

package mixer.utils.slice.cleaning;

import javastraw.feature2D.Feature2D;
import javastraw.feature2D.Feature2DList;
import javastraw.reader.basics.Chromosome;

import java.awt.*;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;

public class Connected2DComponent {
    private static final int BRIDGE_SIZE = 10; // 10MB
    private static final int MIN_AREA = 20; // 20MB^2

    public static Feature2DList generateConnectedComponentFeature(Chromosome chrom1, Chromosome chrom2,
                                                                  int resolution, List<int[]> feature2DS,
                                                                  List<Rectangle> bounds) {
        int maxRow = getMaxInDim(feature2DS, 0);
        int maxCol = getMaxInDim(feature2DS, 1);
        boolean[][] bitmap = generateBitmap(maxRow, maxCol, feature2DS);

        List<int[]> boundedBoxes = getBoundedBoxes(bitmap);
        Feature2DList feature2DList = new Feature2DList();
        for (int[] box : boundedBoxes) {
            if (size(box) > MIN_AREA) {
                feature2DList.add(chrom1.getIndex(), chrom2.getIndex(),
                        new Feature2D(Feature2D.FeatureType.NONE,
                                chrom1.getName(), (long) box[0] * resolution, (long) box[1] * resolution,
                                chrom2.getName(), (long) box[2] * resolution, (long) box[3] * resolution,
                                Color.blue, new HashMap<>()));
                bounds.add(new Rectangle(box[0] * resolution, box[2] * resolution,
                        (box[1] - box[0]) * resolution, (box[3] - box[2]) * resolution));
            }
        }
        return feature2DList;
    }

    private static int size(int[] box) {
        return (box[1] - box[0]) * (box[3] - box[2]);
    }


    private static List<int[]> getBoundedBoxes(boolean[][] bitmap) {
        List<int[]> boxes = new ArrayList<>();
        while (hasProblems(bitmap)) {
            boxes.add(grabBoxFrom(bitmap));
        }
        return boxes;
    }

    private static int[] grabBoxFrom(boolean[][] bitmap) {
        int[] position = getFirstPosition(bitmap);
        if (position == null) {
            System.err.println("Invalid bitmap situation");
            System.exit(987);
        }
        int[] currentBound = new int[]{position[0], position[0] + 1, position[1], position[1] + 1};

        boolean isChanging = true;
        while (isChanging) {
            isChanging = expandBounds(bitmap, currentBound);
        }

        eraseBoundsContent(bitmap, currentBound);
        return currentBound;
    }

    private static boolean expandBounds(boolean[][] bitmap, int[] currentBound) {
        boolean isChanging = false;
        if (canExpandEast(bitmap, currentBound)) {
            currentBound[3] += 1;
            isChanging = true;
        }
        if (canExpandWest(bitmap, currentBound)) {
            currentBound[2] -= 1;
            isChanging = true;
        }
        if (canExpandSouth(bitmap, currentBound)) {
            currentBound[1] += 1;
            isChanging = true;
        }
        return isChanging;
    }

    private static boolean canExpandEast(boolean[][] bitmap, int[] currentBound) {
        for (int i = currentBound[0]; i < currentBound[1] + 1; i++) {
            for (int j = currentBound[3]; j < currentBound[3] + 1; j++) {
                try {
                    if (bitmap[i][j]) return true;
                } catch (Exception ignored) {
                }
            }
        }
        return false;
    }

    private static boolean canExpandSouth(boolean[][] bitmap, int[] currentBound) {
        for (int i = currentBound[1]; i < currentBound[1] + 1; i++) {
            for (int j = currentBound[2] - 1; j < currentBound[3] + 1; j++) {
                try {
                    if (bitmap[i][j]) return true;
                } catch (Exception ignored) {
                }
            }
        }
        return false;
    }

    private static boolean canExpandWest(boolean[][] bitmap, int[] currentBound) {
        for (int i = currentBound[0]; i < currentBound[1] + 1; i++) {
            for (int j = currentBound[2] - 1; j < currentBound[2]; j++) {
                try {
                    if (bitmap[i][j]) return true;
                } catch (Exception ignored) {
                }
            }
        }
        return false;
    }

    private static void eraseBoundsContent(boolean[][] bitmap, int[] currentBound) {
        for (int i = currentBound[0]; i < currentBound[1]; i++) {
            for (int j = currentBound[2]; j < currentBound[3]; j++) {
                bitmap[i][j] = false;
            }
        }
    }

    private static int[] getFirstPosition(boolean[][] bitmap) {
        for (int i = 0; i < bitmap.length; i++) {
            for (int j = 0; j < bitmap[i].length; j++) {
                if (bitmap[i][j]) {
                    return new int[]{i, j};
                }
            }
        }
        return null;
    }

    private static boolean hasProblems(boolean[][] bitmap) {
        for (boolean[] arr : bitmap) {
            for (boolean val : arr) {
                if (val) return true;
            }
        }
        return false;
    }

    private static boolean[][] generateBitmap(int maxRow, int maxCol, List<int[]> feature2DS) {
        boolean[][] bitmap = new boolean[maxRow][maxCol];
        for (int[] point : feature2DS) {
            int minR = Math.max(point[0] - 1, 0);
            int maxR = Math.min(point[0] + 1, maxRow - 1);
            int minC = Math.max(point[1] - 1, 0);
            int maxC = Math.min(point[1] + 1, maxCol - 1);

            for (int i = minR; i <= maxR; i++) {
                for (int j = minC; j <= maxC; j++) {
                    bitmap[i][j] = true;
                }
            }
        }

        spikeInBridges(bitmap, BRIDGE_SIZE);

        return bitmap;
    }

    private static void spikeInBridges(boolean[][] bitmap, int bridge) {
        bridgeRows(bitmap, bridge);
        bridgeCols(bitmap, bridge);
    }

    private static void bridgeCols(boolean[][] bitmap, int bridge) {
        for (int i = 0; i < bitmap.length - 1; i++) {
            for (int j = 0; j < bitmap[i].length; j++) {
                if (bitmap[i][j] && !bitmap[i + 1][j]) {
                    attemptBridgeCol(bitmap, i, j, bridge);
                }
            }
        }
    }

    private static void bridgeRows(boolean[][] bitmap, int bridge) {
        for (int i = 0; i < bitmap.length; i++) {
            for (int j = 0; j < bitmap[i].length - 1; j++) {
                if (bitmap[i][j] && !bitmap[i][j + 1]) {
                    attemptBridgeRow(bitmap, i, j, bridge);
                }
            }
        }
    }

    private static void attemptBridgeRow(boolean[][] bitmap, int i0, int j0, int bridge) {
        int maxCol = Math.min(j0 + bridge, bitmap[i0].length - 1);
        for (int j1 = maxCol; j1 > j0; j1--) {
            if (bitmap[i0][j1]) {
                for (int j2 = j0 + 1; j2 < j1; j2++) {
                    bitmap[i0][j2] = true;
                }
                break;
            }
        }
    }

    private static void attemptBridgeCol(boolean[][] bitmap, int i0, int j0, int bridge) {
        int maxRow = Math.min(i0 + bridge, bitmap.length - 1);
        for (int i1 = maxRow; i1 > i0; i1--) {
            if (bitmap[i1][j0]) {
                for (int i2 = i0 + 1; i2 < i1; i2++) {
                    bitmap[i2][j0] = true;
                }
                break;
            }
        }
    }

    private static int getMaxInDim(List<int[]> feature2DS, int index) {
        int maxVal = feature2DS.get(0)[index];
        for (int[] point : feature2DS) {
            maxVal = Math.max(maxVal, point[index]);
        }
        return maxVal + 1;
    }

}
