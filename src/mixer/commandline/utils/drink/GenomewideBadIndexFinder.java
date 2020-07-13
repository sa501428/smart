/*
 * The MIT License (MIT)
 *
 * Copyright (c) 2011-2020 Rice University, Baylor College of Medicine, Aiden Lab
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

package mixer.commandline.utils.drink;

import mixer.commandline.utils.common.FloatMatrixTools;
import mixer.data.Dataset;
import mixer.data.HiCFileTools;
import mixer.data.MatrixZoomData;
import mixer.windowui.NormalizationType;
import org.broad.igv.feature.Chromosome;
import org.broad.igv.util.Pair;

import java.io.IOException;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Map;
import java.util.Set;

public class GenomewideBadIndexFinder {
	
	private static final float ZSCORE_THRESHOLD_HIGHER = -1.28f; // bottom 10% dropped
	private static final float ZSCORE_THRESHOLD_LOWER = -3; // bottom 0.15% dropped
	private final int resolution;
	private final Map<Integer, Set<Integer>> badIndices = new HashMap<>();
	private final Map<Integer, Set<Integer>> worstIndices = new HashMap<>();
	private NormalizationType norm;
	
	public GenomewideBadIndexFinder(Dataset ds, Chromosome[] chromosomes, int resolution, NormalizationType norm) {
		this.resolution = resolution;
		this.norm = norm;
		createInternalBadList(ds, chromosomes);
	}
	
	/**
	 * @param matrix
	 * @return
	 */
	private static Pair<int[], int[]> getNumberOfNonZeros(float[][] matrix) {
		int[] numNonZerosRows = new int[matrix.length];
		int[] numNonZerosCols = new int[matrix[0].length];
		
		for (int i = 0; i < matrix.length; i++) {
			for (int j = 0; j < matrix[i].length; j++) {
				if (matrix[i][j] > 1e-1) { // not zero
					numNonZerosRows[i]++;
					numNonZerosCols[j]++;
				}
			}
		}
		return new Pair<>(numNonZerosRows, numNonZerosCols);
	}
	
	private void createInternalBadList(Dataset ds, Chromosome[] chromosomes) {
		
		for (Chromosome chr1 : chromosomes) {
			badIndices.put(chr1.getIndex(), new HashSet<Integer>());
			worstIndices.put(chr1.getIndex(), new HashSet<Integer>());
		}
		
		for (int i = 0; i < chromosomes.length; i++) {
			Chromosome chr1 = chromosomes[i];
			for (int j = i; j < chromosomes.length; j++) {
				Chromosome chr2 = chromosomes[j];
				final MatrixZoomData zd = HiCFileTools.getMatrixZoomData(ds, chr1, chr2, resolution);
				determineBadIndicesForRegion(chr1, chr2, zd, i == j);
			}
		}
	}
	
	private void determineBadIndicesForRegion(Chromosome chr1, Chromosome chr2, MatrixZoomData zd, boolean isIntra) {
		
		int lengthChr1 = chr1.getLength() / resolution + 1;
		int lengthChr2 = chr2.getLength() / resolution + 1;
		
		float[][] allDataForRegion = null;
		try {
			allDataForRegion = HiCFileTools.extractLocalBoundedRegionFloatMatrix(zd, 0,
					lengthChr1, 0, lengthChr2, lengthChr1, lengthChr2, norm, isIntra);
		} catch (IOException e) {
			e.printStackTrace();
		}
		
		if (allDataForRegion == null) {
			System.err.println("Missing data for " + zd.getKey());
			System.exit(98);
		}
		
		FloatMatrixTools.cleanUpNansInfinitesNegatives(allDataForRegion);
		
		determineBadIndicesForRegion(chr1, chr2, allDataForRegion, isIntra);
	}
	
	private void determineBadIndicesForRegion(Chromosome chr1, Chromosome chr2, float[][] data, boolean isIntra) {
		Pair<int[], int[]> results = getNumberOfNonZeros(data);
		removeSparserIndicesZscoredCount(chr1, results.getFirst());
		if (!isIntra) {
			removeSparserIndicesZscoredCount(chr2, results.getSecond());
		}
	}
	
	private void removeSparserIndicesZscoredCount(Chromosome chr1, int[] numNonZeros) {
		float mean = getNonZeroMean(numNonZeros);
		float stdDev = getNonZeroStd(numNonZeros, mean);
		getBadIndicesByZscore(chr1, numNonZeros, mean, stdDev);
	}
	
	private void getBadIndicesByZscore(Chromosome chr1, int[] numNonZeros, float mean, float stdDev) {
		for (int k = 0; k < numNonZeros.length; k++) {
			if (numNonZeros[k] > 0) {
				float zval = (numNonZeros[k] - mean) / stdDev;
				if (zval < ZSCORE_THRESHOLD_LOWER) {
					badIndices.get(chr1.getIndex()).add(k);
					worstIndices.get(chr1.getIndex()).add(k);
				} else if (zval < ZSCORE_THRESHOLD_HIGHER) {
					badIndices.get(chr1.getIndex()).add(k);
				}
			} else {
				badIndices.get(chr1.getIndex()).add(k);
				worstIndices.get(chr1.getIndex()).add(k);
			}
		}
	}
	
	private float getNonZeroStd(int[] numNonZeros, float mean) {
		int count = 0;
		double total = 0;
		for (int val : numNonZeros) {
			if (val > 0) {
				float diff = val - mean;
				total += (diff * diff);
				count++;
			}
		}
		return (float) Math.sqrt(total / count);
	}
	
	/**
	 * @param numNonZeros
	 * @return
	 */
	private float getNonZeroMean(int[] numNonZeros) {
		int count = 0;
		float total = 0;
		for (int val : numNonZeros) {
			if (val > 0) {
				total += val;
				count++;
			}
		}
		return total / count;
	}
	
	public Set<Integer> getBadIndices(Chromosome chrom) {
		return badIndices.get(chrom.getIndex());
	}
	
	public Set<Integer> getWorstIndices(Chromosome chrom) {
		return worstIndices.get(chrom.getIndex());
	}
}
