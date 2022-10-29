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

package mixer.utils.shuffle;

public class BinSplitter {
    private final int N;

    public BinSplitter(int n) {
        this.N = n;
    }

    public int getSectionID(int x, int y) {
        return (N * (x % N)) + (y % N);
    }

    public int getNumGroups() {
        return N * N;
    }

    public float[][] flatten(double[][][] input) {
        int n = input[0].length;
        float[][] matrix = new float[N * n][N * n];
        for (int i = 0; i < n; i++) {
            for (int j = 0; j < n; j++) {
                for (int ni = 0; ni < N; ni++) {
                    for (int nj = 0; nj < N; nj++) {
                        int id = getSectionID(ni, nj);
                        matrix[(N * i) + ni][(N * j) + nj] = (float) input[id][i][j];
                    }
                }
            }
        }
        return matrix;
    }
}
