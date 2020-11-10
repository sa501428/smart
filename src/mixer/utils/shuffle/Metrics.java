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

package mixer.utils.shuffle;

import mixer.utils.slice.kmeansfloat.ClusterTools;

public class Metrics {

    public static void fillEntries(int i, int colStart, float[][] matrix, float[][] centroids, float[][] result,
                                   Type type, boolean isSymmetric) {
        switch (type) {
            case CORRELATION:
            case PRE_Z_CORRELATION:
                for (int j = colStart; j < centroids.length; j++) {
                    float val = getCorrFromCosineStyleSimilarity(matrix[i], centroids[j]);
                    result[i][j] = arctanh(val);
                }
                break;
            case MSE:
            case PRE_Z_MSE:
                for (int j = colStart; j < centroids.length; j++) {
                    result[i][j] = (float) ClusterTools.getNonNanMeanSquaredError(matrix[i], centroids[j]);
                }
                break;
            case MAE:
            case PRE_Z_MAE:
                for (int j = colStart; j < centroids.length; j++) {
                    result[i][j] = (float) ClusterTools.getNonNanMeanAbsoluteError(matrix[i], centroids[j]);
                }
                break;
            case EARTHMOVER:
                for (int j = colStart; j < centroids.length; j++) {
                    result[i][j] = (float) nonNanEarthMoversDistance(matrix[i], centroids[j]);
                }
                break;
            case KL:
                for (int j = colStart; j < centroids.length; j++) {
                    result[i][j] = (float) nonNanKLDivergence(matrix[i], centroids[j]);
                }
                break;
            case JS:
                for (int j = colStart; j < centroids.length; j++) {
                    result[i][j] = (float) nonNanJSDistance(matrix[i], centroids[j]);
                }
                break;
            case COSINE:
            case PRE_Z_COSINE:
            default:
                for (int j = colStart; j < centroids.length; j++) {
                    float val = cosineSimilarity(matrix[i], centroids[j]);
                    result[i][j] = arctanh(val);
                }
                break;
        }

        if (isSymmetric) {
            for (int j = colStart; j < centroids.length; j++) {
                result[j][i] = result[i][j];
            }
        }
    }

    private static float arctanh(float x) {
        float val = Math.max(x, -.99f);
        val = Math.min(val, .99f);
        val = (float) (Math.log(1 + val) - Math.log(1 - val)) / 2;
        if (Float.isInfinite(val)) {
            val = Float.NaN;
        }
        return val;
    }

    private static float cosineSimilarity(float[] vectorA, float[] vectorB) {
        double dotProduct = 0.0;
        double normA = 0.0;
        double normB = 0.0;
        for (int i = 0; i < vectorA.length; i++) {
            boolean entryIsBad = Float.isNaN(vectorA[i]) || Float.isNaN(vectorB[i]);
            if (!entryIsBad) {
                dotProduct += vectorA[i] * vectorB[i];
                normA += vectorA[i] * vectorA[i];
                normB += vectorB[i] * vectorB[i];
            }
        }
        return (float) (dotProduct / (Math.sqrt(normA) * Math.sqrt(normB)));
    }

    public static float getCorrFromCosineStyleSimilarity(float[] vectorA, float[] vectorB) {
        int counter = 0;
        double sumA = 0;
        double sumB = 0;
        for (int i = 0; i < vectorA.length; i++) {
            boolean entryIsBad = Float.isNaN(vectorA[i]) || Float.isNaN(vectorB[i]);
            if (!entryIsBad) {
                sumA += vectorA[i];
                sumB += vectorB[i];
                counter++;
            }
        }
        double avgA = sumA / counter;
        double avgB = sumB / counter;

        double dotProduct = 0.0;
        double normA = 0.0;
        double normB = 0.0;
        for (int i = 0; i < vectorA.length; i++) {
            boolean entryIsBad = Float.isNaN(vectorA[i]) || Float.isNaN(vectorB[i]);
            if (!entryIsBad) {
                double tempA = vectorA[i] - avgA;
                double tempB = vectorB[i] - avgB;

                dotProduct += tempA * tempB;
                normA += tempA * tempA;
                normB += tempB * tempB;
            }
        }
        return (float) (dotProduct / (Math.sqrt(normA) * Math.sqrt(normB)));
    }

    public static double nonNanEarthMoversDistance(float[] a, float[] b) {
        double sumA = 0, sumB = 0;
        boolean[] isValid = new boolean[a.length];
        for (int k = 0; k < a.length; k++) {
            isValid[k] = !(Float.isNaN(a[k]) || Float.isNaN(b[k]));
            if (isValid[k]) {
                sumA += a[k];
                sumB += b[k];
            }
        }
        sumA = sumA / 100;
        sumB = sumB / 100;

        double lastDistance = 0;
        double totalDistance = 0;
        for (int i = 0; i < a.length; i++) {
            if (isValid[i]) {
                double tempA = a[i] / sumA;
                double tempB = b[i] / sumB;
                final double currentDistance = (tempA + lastDistance) - tempB;
                totalDistance += Math.abs(currentDistance);
                lastDistance = currentDistance;
            }
        }
        return totalDistance;
    }

    public static double nonNanKLDivergence(float[] a, float[] centroid) {
        double sumA = 0, sumB = 0;
        boolean[] isValid = new boolean[a.length];
        for (int k = 0; k < a.length; k++) {
            isValid[k] = !(Float.isNaN(a[k]) || Float.isNaN(centroid[k]));
            if (isValid[k]) {
                sumA += a[k];
                sumB += centroid[k];
            }
        }

        double totalDistance = 0;
        for (int i = 0; i < a.length; i++) {
            if (isValid[i]) {
                double p = a[i] / sumA;
                double q = centroid[i] / sumB;
                totalDistance += (p * Math.log(p / q));
            }
        }
        return totalDistance;
    }

    public static double nonNanJSDistance(float[] a, float[] b) {
        double sumA = 0, sumB = 0;
        boolean[] isValid = new boolean[a.length];
        float[] average = new float[a.length];
        for (int k = 0; k < a.length; k++) {
            isValid[k] = !(Float.isNaN(a[k]) || Float.isNaN(b[k]));
            if (isValid[k]) {
                sumA += a[k];
                sumB += b[k];
                average[k] = a[k] + b[k];
            }
        }
        double sumC = sumA + sumB;

        double distP = 0, distQ = 0;
        for (int i = 0; i < a.length; i++) {
            if (isValid[i]) {
                double p = a[i] / sumA;
                double q = b[i] / sumB;
                double m = average[i] / sumC;
                distP += (p * Math.log(p / m));
                distQ += (q * Math.log(q / m));
            }
        }
        return Math.sqrt(distP / 2 + distQ / 2);
    }

    public enum Type {
        MSE, MAE, COSINE, CORRELATION,
        EARTHMOVER, KL, JS, PRE_Z_MSE, PRE_Z_MAE, PRE_Z_COSINE, PRE_Z_CORRELATION
    }
}
