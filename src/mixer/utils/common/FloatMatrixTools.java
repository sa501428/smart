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

package mixer.utils.common;

import javax.imageio.ImageIO;
import java.awt.*;
import java.awt.image.BufferedImage;
import java.io.File;


/**
 * Helper methods to handle matrix operations
 */
@SuppressWarnings("ForLoopReplaceableByForEach")
public class FloatMatrixTools {

    public static float[][] convert(double[][] matrix) {
        float[][] result = new float[matrix.length][matrix[0].length];
        for (int i = 0; i < result.length; i++) {
            for (int j = 0; j < result[i].length; j++) {
                result[i][j] = (float) matrix[i][j];
            }
        }
        //System.out.println("Converted");
        return result;
    }

    public static float[][] transpose(float[][] matrix) {
        float[][] result = new float[matrix[0].length][matrix.length];
        for (int i = 0; i < matrix.length; i++) {
            for (int j = 0; j < matrix[i].length; j++) {
                result[j][i] = matrix[i][j];
            }
        }
        return result;
    }

    public static float[][] cleanUpMatrix(float[][] matrix, boolean shouldZeroNans) {
        for (int r = 0; r < matrix.length; r++) {
            for (int c = 0; c < matrix[r].length; c++) {
                float val = matrix[r][c];
                if (Float.isInfinite(val) || Math.abs(val) < 1E-10) {
                    matrix[r][c] = 0;
                }
                if (shouldZeroNans && Float.isNaN(val)) {
                    matrix[r][c] = 0;
                }
            }
        }
        return matrix;
    }

    private static final int MIN_PIXEL_WIDTH = 50;

    public static void saveMatrixToPNG(File file, float[][] matrix0, boolean useLog) {
        float[][] matrix = matrix0;
        if (matrix.length < MIN_PIXEL_WIDTH || matrix[0].length < MIN_PIXEL_WIDTH) {
            matrix = expand(matrix0);
        }

        double range = getMaxVal(matrix);
        double minVal = 0;
        if (useLog) {
            minVal = Math.log(1 + getMinVal(matrix));
            range = Math.log(1 + range) - minVal;
        }

        BufferedImage image = new BufferedImage(matrix[0].length, matrix.length, BufferedImage.TYPE_INT_RGB);
        for (int i = 0; i < matrix.length; i++) {
            for (int j = 0; j < matrix[0].length; j++) {
                if (useLog) {
                    image.setRGB(j, i, mixColors((Math.log(1 + matrix[i][j]) - minVal) / range));
                } else {
                    image.setRGB(j, i, mixColors(matrix[i][j] / range));
                }
            }
        }

        try {
            ImageIO.write(image, "png", file);
        } catch (Exception e) {
            e.printStackTrace();
        }
    }

    private static float[][] expand(float[][] input) {
        float[][] result = new float[MIN_PIXEL_WIDTH * input.length][MIN_PIXEL_WIDTH * input[0].length];
        for (int i = 0; i < result.length; i++) {
            for (int j = 0; j < result[0].length; j++) {
                result[i][j] = input[i / MIN_PIXEL_WIDTH][j / MIN_PIXEL_WIDTH];
            }
        }
        return result;
    }

    public static float[][] deepClone(float[][] data) {
        float[][] copy = new float[data.length][data[0].length];
        for (int i = 0; i < data.length; i++) {
            System.arraycopy(data[i], 0, copy[i], 0, data[i].length);
        }
        return copy;
    }

    private static float getMaxVal(float[][] matrix) {
        float maxVal = matrix[0][0];
        for (int i = 0; i < matrix.length; i++) {
            for (int j = 0; j < matrix[i].length; j++) {
                if (matrix[i][j] > maxVal) {
                    maxVal = matrix[i][j];
                }
            }
        }
        return maxVal;
    }

    private static float getMinVal(float[][] matrix) {
        float minVal = matrix[0][0];
        for (int i = 0; i < matrix.length; i++) {
            for (int j = 0; j < matrix[i].length; j++) {
                if (matrix[i][j] < minVal) {
                    minVal = matrix[i][j];
                }
            }
        }
        return minVal;
    }

    public static int mixColors(double ratio) {
        int color1 = Color.WHITE.getRGB();
        int color2 = Color.RED.getRGB();
        return maskColors(ratio, color1, color2);
    }

    private static int maskColors(double ratio, int color1, int color2) {
        int mask1 = 0x00ff00ff;
        int mask2 = 0xff00ff00;

        int f2 = (int) (256 * ratio);
        int f1 = 256 - f2;

        return (((((color1 & mask1) * f1) + ((color2 & mask1) * f2)) >> 8) & mask1)
                | (((((color1 & mask2) * f1) + ((color2 & mask2) * f2)) >> 8) & mask2);
    }
}