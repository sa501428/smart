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

package mixer.utils.impute;

import tagbio.umap.Umap;
import tagbio.umap.metric.ManhattanMetric;

public class UmapTransform {

    public static float[][] run(float[][] points, int dim) {
        if (points[0].length > dim) {
            final Umap umap = new Umap();
            umap.setNumberComponents(dim);
            umap.setNumberNearestNeighbours(50); // 50 // 15 -> 50 for more global picture
            umap.setThreads(10);
            umap.setMinDist(0.5f); // 0.2 ->0.8 -> 0.5  //0.1f -> 0.2f for more general features
            umap.setVerbose(false);
            umap.setSeed(1L);
            umap.setMetric(ManhattanMetric.SINGLETON);
            return umap.fitTransform(points);
        }
        return points;
    }
}