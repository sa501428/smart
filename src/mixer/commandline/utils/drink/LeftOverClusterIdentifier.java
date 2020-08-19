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

import mixer.commandline.handling.AggregateProcessing;
import mixer.commandline.utils.common.DoubleMatrixTools;
import mixer.commandline.utils.drink.kmeansfloat.ClusterTools;
import mixer.data.ChromosomeHandler;
import mixer.data.Dataset;
import mixer.data.HiCFileTools;
import mixer.data.MatrixZoomData;
import mixer.data.feature.GenomeWideList;
import mixer.windowui.NormalizationType;
import org.apache.commons.math.linear.RealMatrix;
import org.broad.igv.feature.Chromosome;

import java.util.*;

public class LeftOverClusterIdentifier {
    public static void identify(ChromosomeHandler chromosomeHandler, Dataset ds, NormalizationType norm, int resolution,
                                Map<Integer, GenomeWideList<SubcompartmentInterval>> results,
                                GenomewideBadIndexFinder badIndexFinder,
                                float threshold) {
    
        for (Chromosome chr1 : chromosomeHandler.getAutosomalChromosomesArray()) {
            final MatrixZoomData zd = HiCFileTools.getMatrixZoomData(ds, chr1, chr1, resolution);
            if (zd == null) continue;
        
            float[][] allDataForRegion = null;
            try {
                RealMatrix localizedRegionData = HiCFileTools.getRealOEMatrixForChromosome(ds, zd, chr1, resolution,
                        norm, threshold,
                        AggregateProcessing.afterThresholdType,
                        true);
            
                allDataForRegion = DoubleMatrixTools.convertToFloatMatrix(
                        DoubleMatrixTools.cleanUpMatrix(localizedRegionData.getData()));
            
            } catch (Exception e) {
                e.printStackTrace();
                System.exit(99);
            }
        
            if (allDataForRegion == null) {
                System.err.println("Missing Data " + zd.getKey());
                return;
            }
        
            Set<Integer> worstIndices = badIndexFinder.getWorstIndices(chr1);
            Set<Integer> indicesMissing = new HashSet<>();
        
        
            for (int k = 0; k < chr1.getLength() / resolution + 1; k++) {
                if (worstIndices.contains(k)) {
                    continue;
                }
            
                indicesMissing.add(k);
            }
        
            // do it this way because of additional internal filter
            for (SubcompartmentInterval handledInterval : results.get(2).getFeatures("" + chr1.getIndex())) {
                int binXStart = handledInterval.getX1() / resolution;
                int binXEnd = handledInterval.getX2() / resolution;

                for (int j = binXStart; j < binXEnd; j++) {
                    indicesMissing.remove(j);
                }
            }

            if (indicesMissing.size() > 0) {
                for (Integer key : results.keySet()) {
                    GenomeWideList<SubcompartmentInterval> listForKey = results.get(key);
                    Map<Integer, float[]> cIDToCenter = getClusterCenters(allDataForRegion, listForKey.getFeatures("" + chr1.getIndex()), resolution);

                    List<SubcompartmentInterval> newlyAssignedSubcompartments = getNewlyAssignedCompartments(chr1, cIDToCenter, indicesMissing, allDataForRegion, resolution);

                    listForKey.addAll(newlyAssignedSubcompartments);

                    results.put(key, listForKey);
                }
            }
            System.out.print(".");
        }
        System.out.println(".");
    }

    private static Map<Integer, float[]> getClusterCenters(float[][] allDataForRegion, List<SubcompartmentInterval> intervals, int resolution) {

        Map<Integer, float[]> cIDToCenter = new HashMap<>();
        Map<Integer, Integer> cIDToSize = new HashMap<>();


        for (SubcompartmentInterval interval : intervals) {
            int binXStart = interval.getX1() / resolution;
            int binXEnd = interval.getX2() / resolution;
            int cID = interval.getClusterID();
            float[] total = new float[allDataForRegion[0].length];

            for (int i = binXStart; i < binXEnd; i++) {
                for (int j = 0; j < allDataForRegion[i].length; j++) {
                    total[j] += allDataForRegion[i][j];
                }
            }

            if (cIDToSize.containsKey(cID)) {
                cIDToSize.put(cID, cIDToSize.get(cID) + binXEnd - binXStart);
                float[] vec = cIDToCenter.get(cID);
                for (int j = 0; j < vec.length; j++) {
                    total[j] += vec[j];
                }
            } else {
                cIDToSize.put(cID, binXEnd - binXStart);
            }
            cIDToCenter.put(cID, total);
        }

        for (Integer key : cIDToCenter.keySet()) {
            cIDToCenter.put(key, ClusterTools.normalize(cIDToCenter.get(key), cIDToSize.get(key)));
        }

        return cIDToCenter;
    }


    private static List<SubcompartmentInterval> getNewlyAssignedCompartments(Chromosome chromosome, Map<Integer, float[]> cIDToCenter,
                                                                             Set<Integer> indicesMissing, float[][] allDataForRegion, int resolution) {

        List<SubcompartmentInterval> intervals = new ArrayList<>();

        for (int indx : indicesMissing) {

            int metaID = getClosestClusterID(allDataForRegion[indx], cIDToCenter);

            int x1 = indx * resolution;
            int x2 = x1 + resolution;

            intervals.add(new SubcompartmentInterval(chromosome.getIndex(), chromosome.getName(), x1, x2, metaID));
        }
        return intervals;
    }

    private static int getClosestClusterID(float[] vector, Map<Integer, float[]> cIDToCenter) {
        int currID = Integer.MAX_VALUE;
        double overallDistance = Double.MAX_VALUE;
        boolean nothingChanged = true;

        for (Integer key : cIDToCenter.keySet()) {
            double newDistance;
            if (AggregateProcessing.useL1Norm) {
                newDistance = ClusterTools.getL1Distance(cIDToCenter.get(key), vector);
            } else {
                newDistance = ClusterTools.getDistance(cIDToCenter.get(key), vector);
            }

            if (newDistance < overallDistance) {
                overallDistance = newDistance;
                currID = key;
                nothingChanged = false;
            }
        }
        if (nothingChanged) {
            System.err.println("Error 787 " + overallDistance + " - " + cIDToCenter.keySet());
        }
        return currID;
    }
}
