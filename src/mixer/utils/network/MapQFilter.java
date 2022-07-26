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

package mixer.utils.network;

import javastraw.feature2D.Feature2D;
import javastraw.feature2D.Feature2DList;
import javastraw.reader.Dataset;
import javastraw.reader.basics.Chromosome;
import javastraw.reader.basics.ChromosomeHandler;
import javastraw.reader.norm.NormalizationVector;
import javastraw.reader.type.NormalizationType;
import javastraw.tools.HiCFileTools;
import mixer.MixerGlobals;

import java.io.IOException;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

public class MapQFilter {

    private final static int NEIGHBOR_DISTANCE = 5;

    public static void removeLowMapQFeatures(Feature2DList list, final int resolution,
                                             final Dataset ds, final ChromosomeHandler chromosomeHandler,
                                             final NormalizationType norm) {

        final Map<String, Integer> chrNameToIndex = new HashMap<>();
        for (Chromosome chr : chromosomeHandler.getChromosomeArray()) {
            chrNameToIndex.put(Feature2DList.getKey(chr, chr), chr.getIndex());
        }
        if (MixerGlobals.printVerboseComments) {
            System.out.println("Initial: " + list.getNumTotalFeatures());
        }
        list.filterLists((chr, feature2DList) -> {
            try {
                return removeLowMapQ(resolution, chrNameToIndex.get(chr), ds, feature2DList, norm);
            } catch (Exception e) {
                System.err.println("Unable to remove low mapQ entries for " + chr);
                //e.printStackTrace();
            }
            return new ArrayList<>();
        });


    }

    private static List<Feature2D> removeLowMapQ(int res, int chrIndex, Dataset ds, List<Feature2D> list, NormalizationType norm) throws IOException {

        List<Feature2D> features = new ArrayList<>();
        NormalizationVector normVectorContainer = ds.getNormalizationVector(chrIndex, ds.getZoomForBPResolution(res),
                norm);
        if (normVectorContainer == null) {
            HiCFileTools.triggerNormError(norm);
        } else {
            double[] normalizationVector = normVectorContainer.getData().getValues().get(0);
            for (Feature2D feature : list) {
                int index1 = (int) (feature.getStart1() / res);
                int index2 = (int) (feature.getStart2() / res);
                if (nearbyValuesClear(normalizationVector, index1) && nearbyValuesClear(normalizationVector, index2)) {
                    features.add(feature);
                }
            }
        }

        return features;
    }

    private static boolean nearbyValuesClear(double[] normalizationVector, int index) {
        for (int i = index - NEIGHBOR_DISTANCE; i <= index + NEIGHBOR_DISTANCE; i++) {
            if (Double.isNaN(normalizationVector[i]))
                return false;
        }
        return true;
    }
}
