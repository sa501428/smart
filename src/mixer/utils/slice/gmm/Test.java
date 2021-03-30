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

import javastraw.tools.MatrixTools;

import java.util.ArrayList;
import java.util.List;
import java.util.Random;

public class Test {

    private static final float[][] data = new float[1000][2];
    private static final int[] id = new int[1000];
    private static final Random generator = new Random(5);

    public static void main(String[] args) {

        List<List<Integer>> startingIndices = new ArrayList<>();
        startingIndices.add(new ArrayList<>());
        startingIndices.add(new ArrayList<>());
        startingIndices.add(new ArrayList<>());

        for (int k = 0; k < data.length; k++) {
            if (p50()) {
                data[k][0] = (float) (generator.nextGaussian() * 15 + 1); // 15 3
                data[k][1] = (float) (generator.nextGaussian() * 15 - 2); // 15 3
                id[k] = 0;
                updateGuesses(startingIndices, 0, 1, 2, k);
            } else if (p50()) {
                data[k][0] = (float) (generator.nextGaussian() * 3 + 12);
                data[k][1] = (float) (generator.nextGaussian() * 3 - 8);
                id[k] = 1;
                updateGuesses(startingIndices, 1, 0, 2, k);
            } else {
                data[k][0] = (float) (generator.nextGaussian() * 2 - 7);
                data[k][1] = (float) (generator.nextGaussian() * 2 + 10);
                id[k] = 2;
                updateGuesses(startingIndices, 2, 1, 0, k);
            }
        }

        MatrixTools.saveMatrixTextNumpy("/Users/mshamim/Desktop/testgmm/actual_data.npy", data);
        MatrixTools.saveMatrixTextNumpy("/Users/mshamim/Desktop/testgmm/actual_ids.npy", id);
        for (int i = 0; i < startingIndices.size(); i++) {
            int m = startingIndices.get(i).size();
            int[] sIDs = new int[m];
            for (int k = 0; k < m; k++) {
                sIDs[k] = startingIndices.get(i).get(k);
            }
            MatrixTools.saveMatrixTextNumpy("/Users/mshamim/Desktop/testgmm/init_ids_" + i + ".npy", sIDs);
        }
        SimpleScatterPlot plotter = new SimpleScatterPlot(data, id);
        plotter.plot("/Users/mshamim/Desktop/testgmm/actual");

        plotter = new SimpleScatterPlot(data, startingIndices);
        plotter.plot("/Users/mshamim/Desktop/testgmm/initial");


        GaussianMixtureModels gmm = new GaussianMixtureModels(data, 3, 20, startingIndices);
        gmm.fit();
        int[] result = gmm.predict();

        plotter = new SimpleScatterPlot(data, result);
        plotter.plot("/Users/mshamim/Desktop/testgmm/gmm_result");

    }

    private static void updateGuesses(List<List<Integer>> startingIndices, int good, int bad1, int bad2, int index) {
        if (p75()) {
            startingIndices.get(good).add(index);
        } else if (p50()) {
            startingIndices.get(bad1).add(index);
        } else {
            startingIndices.get(bad2).add(index);
        }
    }

    private static boolean p50() {
        return generator.nextBoolean();
    }

    private static boolean p75() {
        return generator.nextBoolean() || generator.nextBoolean();
    }

    private static boolean p25() {
        return generator.nextBoolean() && generator.nextBoolean();
    }
}
