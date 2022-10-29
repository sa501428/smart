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

package mixer.utils.drive;

import javastraw.expected.ExpectedUtils;
import javastraw.expected.WelfordArray;
import javastraw.expected.ZScoreArray;
import javastraw.reader.basics.Chromosome;
import javastraw.reader.block.ContactRecord;

import java.util.List;

public class LogExpectedSubset {
    private final float minX, maxX;
    private final ZScoreArray zscores;

    public LogExpectedSubset(List<ContactRecord> records, Chromosome chrom, int res) {
        int maxBin = (int) (chrom.getLength() / res) + 1;
        int[] countsPerDist = getCountsPerDist(records, maxBin);
        int[] minMaxDist = setMinMaxFromCountsNeeded(countsPerDist, 50);
        minX = minMaxDist[0] + 1;
        maxX = minMaxDist[1] - 1;

        WelfordArray distribution = getExpectedDistribution(records, maxBin);
        zscores = distribution.getZscores();
    }

    public static double logp1(double x) {
        return Math.log(1.0 + x);
    }

    public static int getDist(ContactRecord record) {
        return Math.abs(record.getBinX() - record.getBinY());
    }

    private WelfordArray getExpectedDistribution(List<ContactRecord> records, int maxBin) {
        WelfordArray distribution = new WelfordArray(logp1i(maxBin) + 1);
        for (ContactRecord record : records) {
            int dist = getDist(record);
            distribution.addValue(logp1i(dist), logp1(record.getCounts()));
        }
        return distribution;
    }

    public boolean isInInterval(ContactRecord cr) {
        return isInInterval(getDist(cr));
    }

    public boolean isInInterval(int dist) {
        return dist >= minX && dist <= maxX;
    }

    public float getZscoreForObservedUncompressedBin(ContactRecord cr) {
        return (float) zscores.getZscore(logp1i(getDist(cr)), logp1(cr.getCounts()));
    }

    private int[] setMinMaxFromCountsNeeded(int[] counts, int cutoff) {
        int minDist = counts.length;
        int maxDist = 0;
        for (int i = 0; i < counts.length; i++) {
            if (counts[i] < cutoff) {
                counts[i] = 0;
            } else {
                // valid dist
                minDist = Math.min(minDist, i);
                maxDist = Math.max(maxDist, i);
            }
        }

        for (int i = minDist; i < maxDist; i++) {
            if (counts[i] < cutoff) {
                // should not happen; means there's a gap/dip and we should lower the maxDist
                maxDist = i - 1;
                break;
            }
        }

        return new int[]{minDist, maxDist};
    }

    private int[] getCountsPerDist(List<ContactRecord> records, int maxBin) {
        int[] counts = new int[maxBin];
        for (ContactRecord record : records) {
            counts[ExpectedUtils.getDist(record)]++;
        }
        return counts;
    }

    public int logp1i(int x) {
        return (int) Math.log(1 + x);
    }
}
