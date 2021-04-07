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

import java.util.*;

public class Test {

    private static float[][] data;
    private static int[] id;
    private static final Random generator = new Random(5);
    private static int NUM_CLUSTERS = 3;

    public static void main(String[] args) {

        List<List<Integer>> startingIndices = generateData();

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
        SimpleScatterPlot plotter = new SimpleScatterPlot(data);
        plotter.plot(id, "/Users/mshamim/Desktop/testgmm/actual");
        plotter.plot(startingIndices, "/Users/mshamim/Desktop/testgmm/initial");

        putInNans(data, .2); // .2 .1

        //approximateCorrection(data);

        System.out.println("Starting GMM");
        long start = System.nanoTime();
        GaussianMixtureModels gmm = new GaussianMixtureModels(data, NUM_CLUSTERS, 20, startingIndices);
        try {
            gmm.fit();
            int[] result = gmm.predict();
            long end = System.nanoTime();
            System.out.println("Elapsed time: " + (end - start) * 1e-9);
            plotter.plot(result, "/Users/mshamim/Desktop/testgmm/gmm_nan_result");
        } catch (GMMException g) {
            System.err.println("Unable to run GMM on data; matrices are likely singular or some other issue encountered");
        }
    }

    private static void putInNans(float[][] data, double percent) {
        int numCols = data[0].length;
        int numToNanifyInRow = (int) (percent * numCols);
        System.out.println("Num nans per row " + numToNanifyInRow);
        for (int i = 0; i < data.length; i++) {
            Set<Integer> indices = getRandomSetIndices(numToNanifyInRow, numCols);
            for (int j : indices) {
                data[i][j] = Float.NaN;
            }
            if (i < 10) {
                System.out.println(Arrays.toString(data[i]));
            }
        }
    }

    private static Set<Integer> getRandomSetIndices(int size, int maxVal) {
        Set<Integer> vals = new HashSet<>();
        while (vals.size() < size) {
            vals.add(generator.nextInt(maxVal));
        }
        return vals;
    }

    private static List<List<Integer>> generateData() {

        NUM_CLUSTERS = 6;
        int columnNums = 40;//100; // 2000
        int numPoints = 300000; //50k

        data = new float[numPoints][columnNums];
        id = new int[numPoints];
        float[][] weights = generateRandomMatrix(NUM_CLUSTERS, 6); //6
        float[][] offsets = generateRandomMatrix(NUM_CLUSTERS, 6);
        List<List<Integer>> startingIndices = new ArrayList<>();
        for (int s = 0; s < NUM_CLUSTERS; s++) {
            startingIndices.add(new ArrayList<>());
        }

        for (int k = 0; k < data.length; k++) {
            if (k % 10 == 0) {
                id[k] = 0;
            } else if (k % 4 == 0) {
                id[k] = 1;
            } else if (k % 5 == 0) {
                id[k] = 2;
            } else if (k % 3 == 0) {
                id[k] = 3;
            } else if (k % 2 == 0) {
                id[k] = 4;
            } else {
                id[k] = 5;
            }
            fillInRow(data, weights[id[k]], offsets[id[k]], startingIndices, k);
        }
        return startingIndices;
    }

    private static float[][] generateRandomMatrix(int numRows, int numCols) {
        float[][] result = new float[numRows][numCols];
        for (int i = 0; i < numRows; i++) {
            for (int j = 0; j < numCols; j++) {
                if (p50()) {
                    result[i][j] = 10 * (generator.nextFloat() - .5f);
                } else {
                    result[i][j] = (float) (6 * generator.nextGaussian());
                }
            }
        }
        return result;
    }

    private static void fillInRow(float[][] data, float[] weight, float[] offset, List<List<Integer>> startingIndices, int k) {
        for (int z = 0; z < weight.length; z++) {
            data[k][z] = (float) (generator.nextGaussian() * weight[z] + offset[z]);
        }
        for (int z = weight.length; z < data[k].length; z++) {
            data[k][z] = (float) generator.nextGaussian() + data[k][z % weight.length];
        }
        for (int z = 0; z < data[k].length; z++) {
            data[k][z] += 5 * (generator.nextFloat() - .5f);
        }

        if (p75()) {
            startingIndices.get(id[k]).add(k);
        } else {
            startingIndices.get(generator.nextInt(NUM_CLUSTERS)).add(k);
        }
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

    private static List<List<Integer>> generateData0() {
        List<List<Integer>> startingIndices = new ArrayList<>();
        startingIndices.add(new ArrayList<>());
        startingIndices.add(new ArrayList<>());
        startingIndices.add(new ArrayList<>());
        NUM_CLUSTERS = 3;
        data = new float[1000][2];
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
        return startingIndices;
    }
}
