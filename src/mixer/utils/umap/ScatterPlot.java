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

import mixer.utils.slice.structures.SubcompartmentColors;

import javax.imageio.ImageIO;
import java.awt.*;
import java.awt.image.BufferedImage;
import java.io.File;

public class ScatterPlot {
    private static final Color BACKGROUND_COLOR = Color.BLACK; //Color.WHITE;
    private final int circleOffset = 3;
    private final int circleWidth = circleOffset * 2;
    private final DataCleaner data;

    public ScatterPlot(float[][] points, int[][] ids) {
        data = new DataCleaner(points, ids);
    }

    public void plot(String absPath) {
        for (int i = 0; i < data.numPlots; i++) {
            plotIndividualMap(absPath + "_" + i + ".png", i);
        }
        plotOrderedMap(absPath + "_order.png");
    }

    private void plotOrderedMap(String absPath) {
        BufferedImage image = new BufferedImage(data.width, data.height, BufferedImage.TYPE_INT_ARGB);
        Graphics2D g = (Graphics2D) image.getGraphics();
        colorBackground(g);

        double maxLength = data.points.length;
        for (int i = 0; i < data.points.length; i++) {
            double p = ((double) i) / maxLength;
            drawGradientCircle(g, data.points[i], p);
        }

        save(image, absPath);
    }

    private void plotIndividualMap(String absPath, int mapIndex) {
        BufferedImage image = new BufferedImage(data.width, data.height, BufferedImage.TYPE_INT_ARGB);
        Graphics2D g = (Graphics2D) image.getGraphics();
        colorBackground(g);

        for (int i = 0; i < data.points.length; i++) {
            drawCircle(g, data.points[i], data.ids[i][mapIndex]);
        }

        save(image, absPath);
    }

    private void colorBackground(Graphics2D g) {
        g.setColor(BACKGROUND_COLOR);
        g.fillRect(0, 0, data.width, data.height);
    }

    private void save(BufferedImage image, String absPath) {
        File file = new File(absPath);
        try {
            ImageIO.write(image, "png", file);
        } catch (Exception e) {
            e.printStackTrace();
        }
    }

    private void drawCircle(Graphics g, float[] point, int id) {
        Color color = SubcompartmentColors.getColorWithAlpha(id, 128);
        g.setColor(color);
        g.fillOval((int) (point[0] - circleOffset), (int) (point[1] - circleOffset), circleWidth, circleWidth);
    }

    private void drawGradientCircle(Graphics2D g, float[] point, double p) {
        Color color = getGradientColor(p);
        g.setColor(color);
        g.fillOval((int) (point[0] - circleOffset), (int) (point[1] - circleOffset), circleWidth, circleWidth);
    }

    private Color getGradientColor(double p) {
        int[] color1, color2;
        double p2;
        if (p < .5) {
            // R -> G
            color1 = new int[]{255, 0, 0};
            color2 = new int[]{0, 255, 0};
            p2 = p * 2;
        } else {
            // G -> B
            color1 = new int[]{0, 255, 0};
            color2 = new int[]{0, 0, 255};
            p2 = 2 * p - 1;
        }
        return mix(color1, color2, p2);
    }

    private Color mix(int[] color1, int[] color2, double p) {
        int R = (int) (color1[0] * p + color2[0] * (1 - p));
        int G = (int) (color1[1] * p + color2[1] * (1 - p));
        int B = (int) (color1[2] * p + color2[2] * (1 - p));
        return new Color(R, G, B, 128);
    }
}
