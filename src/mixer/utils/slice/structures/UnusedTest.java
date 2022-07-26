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

package mixer.utils.slice.structures;

import javastraw.feature1D.GenomeWide1DList;
import javastraw.reader.basics.ChromosomeHandler;

import java.io.File;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

public class UnusedTest {
    public static void createCommonList(String locationHuntley, String locationSNIPER, String locationSCI, ChromosomeHandler handler) {

        final int resolution = 100000;

        GenomeWide1DList<SubcompartmentInterval> subcHuntley = SliceUtils.loadFromSubcompartmentBEDFile(handler, locationHuntley);
        SliceUtils.splitGWList(subcHuntley, resolution);

        GenomeWide1DList<SubcompartmentInterval> subcSNIPER = SliceUtils.loadFromSubcompartmentBEDFile(handler, locationSNIPER);
        SliceUtils.splitGWList(subcSNIPER, resolution);

        GenomeWide1DList<SubcompartmentInterval> subcSCI = SliceUtils.loadFromSubcompartmentBEDFile(handler, locationSCI);
        SliceUtils.splitGWList(subcSCI, resolution);

        GenomeWide1DList<SubcompartmentInterval> resultList = SliceUtils.loadFromSubcompartmentBEDFile(handler, locationHuntley);

        resultList.filterLists((key, features) -> {
            if (features.size() > 0) {
                List<SubcompartmentInterval> featureHuntley = subcHuntley.getFeatures(key);
                List<SubcompartmentInterval> featureSNIPER = subcSNIPER.getFeatures(key);
                List<SubcompartmentInterval> featureSCI = subcSCI.getFeatures(key);

                Map<Integer, List<Integer>> positionToID = new HashMap<>();

                updateMapWithIDsFromList(positionToID, featureHuntley);
                updateMapWithIDsFromList(positionToID, featureSNIPER);
                updateMapWithIDsFromList(positionToID, featureSCI);

                List<SubcompartmentInterval> result = new ArrayList<>();

                int chrIndex = features.get(0).getChrIndex();
                String chrName = features.get(0).getChrName();

                for (Integer x1 : positionToID.keySet()) {
                    int winnerID = getModeID(positionToID.get(x1));
                    if (winnerID > 0) {
                        result.add(new SubcompartmentInterval(chrIndex, chrName, x1, x1 + resolution, winnerID));
                    }
                }

                return result;
            }
            return new ArrayList<>();
        });

        SliceUtils.collapseGWList(resultList);

        resultList.simpleExport(new File(locationHuntley + "_New_gold_standard.bed"));
    }

    public static void createUnanimousList(String locationHuntley, String locationSNIPER, String locationSCI, ChromosomeHandler handler) {

        final int resolution = 100000;

        GenomeWide1DList<SubcompartmentInterval> subcHuntley = SliceUtils.loadFromSubcompartmentBEDFile(handler, locationHuntley);
        SliceUtils.splitGWList(subcHuntley, resolution);

        GenomeWide1DList<SubcompartmentInterval> subcSNIPER = SliceUtils.loadFromSubcompartmentBEDFile(handler, locationSNIPER);
        SliceUtils.splitGWList(subcSNIPER, resolution);

        GenomeWide1DList<SubcompartmentInterval> subcSCI = SliceUtils.loadFromSubcompartmentBEDFile(handler, locationSCI);
        SliceUtils.splitGWList(subcSCI, resolution);

        GenomeWide1DList<SubcompartmentInterval> resultList = SliceUtils.loadFromSubcompartmentBEDFile(handler, locationHuntley);

        resultList.filterLists((key, features) -> {
            if (features.size() > 0) {
                List<SubcompartmentInterval> featureHuntley = subcHuntley.getFeatures(key);
                List<SubcompartmentInterval> featureSNIPER = subcSNIPER.getFeatures(key);
                List<SubcompartmentInterval> featureSCI = subcSCI.getFeatures(key);

                Map<Integer, List<Integer>> positionToID = new HashMap<>();

                updateMapWithIDsFromList(positionToID, featureHuntley);
                updateMapWithIDsFromList(positionToID, featureSNIPER);
                updateMapWithIDsFromList(positionToID, featureSCI);

                List<SubcompartmentInterval> result = new ArrayList<>();

                int chrIndex = features.get(0).getChrIndex();
                String chrName = features.get(0).getChrName();


                for (Integer x1 : positionToID.keySet()) {
                    int winnerID = getUnanimousID(positionToID.get(x1));
                    if (winnerID > 0) {
                        result.add(new SubcompartmentInterval(chrIndex, chrName, x1, x1 + resolution, winnerID));
                    }
                }

                return result;
            }
            return new ArrayList<>();
        });

        SliceUtils.collapseGWList(resultList);

        resultList.simpleExport(new File(locationHuntley + "_New_STRICT_gold_standard.bed"));
    }


    public static void createCommonListType2(String locationHuntley, String locationSNIPER, String locationSCI, ChromosomeHandler handler) {

        final int resolution = 100000;

        GenomeWide1DList<SubcompartmentInterval> subcHuntley = SliceUtils.loadFromSubcompartmentBEDFile(handler, locationHuntley);
        SliceUtils.splitGWList(subcHuntley, resolution);

        GenomeWide1DList<SubcompartmentInterval> subcSNIPER = SliceUtils.loadFromSubcompartmentBEDFile(handler, locationSNIPER);
        SliceUtils.splitGWList(subcSNIPER, resolution);

        GenomeWide1DList<SubcompartmentInterval> subcSCI = SliceUtils.loadFromSubcompartmentBEDFile(handler, locationSCI);
        SliceUtils.splitGWList(subcSCI, resolution);

        GenomeWide1DList<SubcompartmentInterval> resultList = SliceUtils.loadFromSubcompartmentBEDFile(handler, locationHuntley);

        resultList.filterLists((key, features) -> {
            if (features.size() > 0) {
                List<SubcompartmentInterval> featureHuntley = subcHuntley.getFeatures(key);
                List<SubcompartmentInterval> featureSNIPER = subcSNIPER.getFeatures(key);
                List<SubcompartmentInterval> featureSCI = subcSCI.getFeatures(key);

                Map<Integer, List<Integer>> positionToIDHuntleySCI = new HashMap<>();
                Map<Integer, List<Integer>> positionToIDSNIPERSCI = new HashMap<>();
                Map<Integer, List<Integer>> allPositionToID = new HashMap<>();

                updateMapWithIDsFromList(positionToIDHuntleySCI, featureHuntley);
                updateMapWithIDsFromList(positionToIDHuntleySCI, featureSCI);

                updateMapWithIDsFromList(positionToIDSNIPERSCI, featureSNIPER);
                updateMapWithIDsFromList(positionToIDSNIPERSCI, featureSCI);

                updateMapWithIDsFromList(allPositionToID, featureHuntley);
                updateMapWithIDsFromList(allPositionToID, featureSNIPER);
                updateMapWithIDsFromList(allPositionToID, featureSCI);

                List<SubcompartmentInterval> result = new ArrayList<>();

                int chrIndex = features.get(0).getChrIndex();
                String chrName = features.get(0).getChrName();

                for (Integer x1 : allPositionToID.keySet()) {
                    int winnerID = getModeID(positionToIDHuntleySCI.get(x1));
                    if (winnerID > 0) {
                        result.add(new SubcompartmentInterval(chrIndex, chrName, x1, x1 + resolution, winnerID));
                    } else {
                        winnerID = getModeID(positionToIDSNIPERSCI.get(x1));
                        if (winnerID > 0) {
                            result.add(new SubcompartmentInterval(chrIndex, chrName, x1, x1 + resolution, winnerID));
                        }
                    }
                }

                return result;
            }
            return new ArrayList<>();
        });

        SliceUtils.collapseGWList(resultList);

        resultList.simpleExport(new File(locationHuntley + "_New_gold_standard_type2.bed"));
    }


    private static int getModeID(List<Integer> ids) {
        if (ids == null) return -2;
        if (ids.size() == 2) {
            if (ids.get(0).equals(ids.get(1))) {
                return ids.get(0);
            }
        } else if (ids.size() == 3) {
            if (ids.get(0).equals(ids.get(1)) || ids.get(0).equals(ids.get(2))) {
                return ids.get(0);
            } else if (ids.get(1).equals(ids.get(2))) {
                return ids.get(1);
            }
        }
        return -1;
    }

    private static int getUnanimousID(List<Integer> ids) {
        if (ids == null) return -2;
        if (ids.size() == 3) {
            if (ids.get(0).equals(ids.get(1)) && ids.get(0).equals(ids.get(2))) {
                return ids.get(0);
            }
        }
        return -1;
    }

    private static void updateMapWithIDsFromList(Map<Integer, List<Integer>> positionToID, List<SubcompartmentInterval> features) {
        for (SubcompartmentInterval interval : features) {
            if (positionToID.containsKey(interval.getX1())) {
                positionToID.get(interval.getX1()).add(interval.getClusterID());
            } else {
                List<Integer> ids = new ArrayList<>();
                ids.add(interval.getClusterID());
                positionToID.put(interval.getX1(), ids);
            }
        }
    }
}
