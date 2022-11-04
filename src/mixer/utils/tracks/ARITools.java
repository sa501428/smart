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

package mixer.utils.tracks;

import javastraw.feature1D.GenomeWide1DList;

import java.util.Arrays;

import static mixer.utils.tracks.Concensus2DTools.populateSummary;

public class ARITools {

    public static double getARI(GenomeWide1DList<SubcompartmentInterval> file1, GenomeWide1DList<SubcompartmentInterval> file2) {
        int[][] summary = populateSummary(file1, file2);

        for (int[] row : summary) {
            System.out.println(Arrays.toString(row));
        }

        int[] rowsSums = sumRows(summary);
        int[] colsSums = sumCols(summary);
        int n = Concensus2DTools.sum(summary);

        double nijTerm = getMatrixTerm(summary);
        double aTerm = getArrayTerm(rowsSums);
        double bTerm = getArrayTerm(colsSums);
        double nTerm = comb2(n);

        double abnTerm = ((aTerm * bTerm) / nTerm);
        double num = nijTerm - abnTerm;
        double denom = 0.5 * (aTerm + bTerm) - abnTerm;

        return num / denom;
    }

    private static double getArrayTerm(int[] array) {
        double total = 0;
        for (int val : array) {
            total += comb2(val);
        }
        return total;
    }

    private static double getMatrixTerm(int[][] summary) {
        double total = 0;
        for (int[] row : summary) {
            for (int val : row) {
                total += comb2(val);
            }
        }
        return total;
    }

    private static int[] sumRows(int[][] summary) {
        int[] sums = new int[summary.length];
        for (int i = 0; i < summary.length; i++) {
            for (int j = 0; j < summary[i].length; j++) {
                sums[i] += summary[i][j];
            }
        }
        return sums;
    }

    private static int[] sumCols(int[][] summary) {
        int[] sums = new int[summary[0].length];
        for (int i = 0; i < summary.length; i++) {
            for (int j = 0; j < summary[i].length; j++) {
                sums[j] += summary[i][j];
            }
        }
        return sums;
    }

    private static double comb2(int n) {
        return (double) (n * (n - 1)) / 2.0;
    }
}
