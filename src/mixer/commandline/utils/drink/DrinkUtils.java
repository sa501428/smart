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

import mixer.MixerGlobals;
import mixer.data.ChromosomeHandler;
import mixer.data.HiCFileTools;
import mixer.data.feature.FeatureFilter;
import mixer.data.feature.FeatureFunction;
import mixer.data.feature.GenomeWideList;
import org.broad.igv.Globals;
import org.broad.igv.feature.Chromosome;
import org.broad.igv.util.ParsingUtils;

import java.io.BufferedReader;
import java.io.File;
import java.io.IOException;
import java.io.InputStreamReader;
import java.util.*;
import java.util.concurrent.atomic.AtomicInteger;

public class DrinkUtils {

    public static void reSort(GenomeWideList<SubcompartmentInterval> subcompartments) {
        subcompartments.filterLists(new FeatureFilter<SubcompartmentInterval>() {
            @Override
            public List<SubcompartmentInterval> filter(String chr, List<SubcompartmentInterval> featureList) {
                Collections.sort(featureList);
                return featureList;
            }
        });
    }

    public static void collapseGWList(GenomeWideList<SubcompartmentInterval> intraSubcompartments) {
        intraSubcompartments.filterLists(new FeatureFilter<SubcompartmentInterval>() {
            @Override
            public List<SubcompartmentInterval> filter(String chr, List<SubcompartmentInterval> featureList) {
                return collapseSubcompartmentIntervals(featureList);
            }
        });
    }

    private static List<SubcompartmentInterval> collapseSubcompartmentIntervals(List<SubcompartmentInterval> intervals) {
        if (intervals.size() > 0) {

            Collections.sort(intervals);
            SubcompartmentInterval collapsedInterval = (SubcompartmentInterval) intervals.get(0).deepClone();

            Set<SubcompartmentInterval> newIntervals = new HashSet<>();
            for (SubcompartmentInterval nextInterval : intervals) {
                if (collapsedInterval.overlapsWith(nextInterval)) {
                    collapsedInterval = collapsedInterval.absorbAndReturnNewInterval(nextInterval);
                } else {
                    newIntervals.add(collapsedInterval);
                    collapsedInterval = (SubcompartmentInterval) nextInterval.deepClone();
                }
            }
            newIntervals.add(collapsedInterval);

            List<SubcompartmentInterval> newIntervalsSorted = new ArrayList<>(newIntervals);
            Collections.sort(newIntervalsSorted);

            return newIntervalsSorted;
        }
        return intervals;
    }

    public static void splitGWList(GenomeWideList<SubcompartmentInterval> intraSubcompartments, int width) {
        intraSubcompartments.filterLists(new FeatureFilter<SubcompartmentInterval>() {
            @Override
            public List<SubcompartmentInterval> filter(String chr, List<SubcompartmentInterval> featureList) {
                return splitSubcompartmentIntervals(featureList, width);
            }
        });
    }

    private static List<SubcompartmentInterval> splitSubcompartmentIntervals(List<SubcompartmentInterval> intervals, int width) {
        if (intervals.size() > 0) {

            Collections.sort(intervals);

            Set<SubcompartmentInterval> newIntervals = new HashSet<>();
            for (SubcompartmentInterval currInterval : intervals) {
                newIntervals.addAll(currInterval.splitByWidth(width));
            }

            List<SubcompartmentInterval> newIntervalsSorted = new ArrayList<>(newIntervals);
            Collections.sort(newIntervalsSorted);

            return newIntervalsSorted;
        }
        return intervals;
    }


    public static String cleanUpPath(String filePath) {
        String[] breakUpFileName = filePath.split("/");
        return breakUpFileName[breakUpFileName.length - 1].replaceAll(".hic", "");
    }

    /*
    public static void saveFileBeforeAndAfterCollapsing(GenomeWideList<SubcompartmentInterval> subcompartments,
                                                        File outputDirectory, String preCollapsingFileName,
                                                        String postCollapsingFileName) {
        File outputFile = new File(outputDirectory, preCollapsingFileName);
        subcompartments.simpleExport(outputFile);

        collapseGWList(subcompartments);

        File outputFile2 = new File(outputDirectory, postCollapsingFileName);
        subcompartments.simpleExport(outputFile2);
    }

    /*
        // find the most common
        List<SubcompartmentInterval> frequentFliers = new ArrayList();
        for(Integer x : allSubcompartmentIntervalsMap.keySet()) {

            Map<Integer, Integer> counts = new HashMap<>();
            for (Pair<Integer, Integer> pair :allSubcompartmentIntervalsMap.get(x)){
                int value = pair.getValue();
                if(counts.containsKey(value)){
                    counts.put(value,counts.get(value)+1);
                }
                else {
                    counts.put(value,1);
                }
            }
            int maxFrequency = Ints.max(Ints.toArray(counts.values()));
            if(maxFrequency > 1){
                int commonClusterID = -1;
                for(Integer clusterID : counts.keySet()){
                    if(counts.get(clusterID) >= maxFrequency){
                        commonClusterID = clusterID;
                        break;
                    }
                }

                int x2 = x + getResolution();
                frequentFliers.add(new SubcompartmentInterval(chromosome.getIndex(), chromosome.getName(), x, x2, commonClusterID));
            }
        }
        mostFrequentSubcompartment.addAll(frequentFliers);
        */

    public static void readInFileAndCollapse(String location) {

        ChromosomeHandler handler = HiCFileTools.loadChromosomes("hg19");
        GenomeWideList<SubcompartmentInterval> subc = loadFromSubcompartmentBEDFile(handler, location);
        DrinkUtils.collapseGWList(subc);
        subc.simpleExport(new File(location + "_collapsed.bed"));
    }


    public static void readInFileAndSplitByResolutionLevel(String location) {

        ChromosomeHandler handler = HiCFileTools.loadChromosomes("hg19");
        GenomeWideList<SubcompartmentInterval> subc = loadFromSubcompartmentBEDFile(handler, location);
        DrinkUtils.splitGWList(subc, 100000);
        subc.simpleExport(new File(location + "_split.bed"));

    }

    public static void createCommonList(String locationHuntley, String locationSNIPER, String locationSCI) {

        final int resolution = 100000;

        ChromosomeHandler handler = HiCFileTools.loadChromosomes("hg19");
        GenomeWideList<SubcompartmentInterval> subcHuntley = loadFromSubcompartmentBEDFile(handler, locationHuntley);
        splitGWList(subcHuntley, resolution);

        GenomeWideList<SubcompartmentInterval> subcSNIPER = loadFromSubcompartmentBEDFile(handler, locationSNIPER);
        splitGWList(subcSNIPER, resolution);

        GenomeWideList<SubcompartmentInterval> subcSCI = loadFromSubcompartmentBEDFile(handler, locationSCI);
        splitGWList(subcSCI, resolution);

        GenomeWideList<SubcompartmentInterval> resultList = loadFromSubcompartmentBEDFile(handler, locationHuntley);

        resultList.filterLists(new FeatureFilter<SubcompartmentInterval>() {
            @Override
            public List<SubcompartmentInterval> filter(String key, List<SubcompartmentInterval> features) {
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
            }
        });

        collapseGWList(resultList);

        resultList.simpleExport(new File(locationHuntley + "_New_gold_standard.bed"));
    }

    public static void createUnanimousList(String locationHuntley, String locationSNIPER, String locationSCI) {

        final int resolution = 100000;

        ChromosomeHandler handler = HiCFileTools.loadChromosomes("hg19");
        GenomeWideList<SubcompartmentInterval> subcHuntley = loadFromSubcompartmentBEDFile(handler, locationHuntley);
        splitGWList(subcHuntley, resolution);

        GenomeWideList<SubcompartmentInterval> subcSNIPER = loadFromSubcompartmentBEDFile(handler, locationSNIPER);
        splitGWList(subcSNIPER, resolution);

        GenomeWideList<SubcompartmentInterval> subcSCI = loadFromSubcompartmentBEDFile(handler, locationSCI);
        splitGWList(subcSCI, resolution);

        GenomeWideList<SubcompartmentInterval> resultList = loadFromSubcompartmentBEDFile(handler, locationHuntley);

        resultList.filterLists(new FeatureFilter<SubcompartmentInterval>() {
            @Override
            public List<SubcompartmentInterval> filter(String key, List<SubcompartmentInterval> features) {
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
            }
        });

        collapseGWList(resultList);

        resultList.simpleExport(new File(locationHuntley + "_New_STRICT_gold_standard.bed"));
    }


    public static void createCommonListType2(String locationHuntley, String locationSNIPER, String locationSCI) {

        final int resolution = 100000;

        ChromosomeHandler handler = HiCFileTools.loadChromosomes("hg19");
        GenomeWideList<SubcompartmentInterval> subcHuntley = loadFromSubcompartmentBEDFile(handler, locationHuntley);
        splitGWList(subcHuntley, resolution);

        GenomeWideList<SubcompartmentInterval> subcSNIPER = loadFromSubcompartmentBEDFile(handler, locationSNIPER);
        splitGWList(subcSNIPER, resolution);

        GenomeWideList<SubcompartmentInterval> subcSCI = loadFromSubcompartmentBEDFile(handler, locationSCI);
        splitGWList(subcSCI, resolution);

        GenomeWideList<SubcompartmentInterval> resultList = loadFromSubcompartmentBEDFile(handler, locationHuntley);

        resultList.filterLists(new FeatureFilter<SubcompartmentInterval>() {
            @Override
            public List<SubcompartmentInterval> filter(String key, List<SubcompartmentInterval> features) {
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
            }
        });

        collapseGWList(resultList);

        resultList.simpleExport(new File(locationHuntley + "_New_gold_standard_type2.bed"));
    }

    private static int getModeID(List<Integer> ids) {
        if (ids == null) return -2;
        if (ids.size() == 2) {
            if (ids.get(0) == ids.get(1)) {
                return ids.get(0);
            }
        } else if (ids.size() == 3) {
            if (ids.get(0) == ids.get(1) || ids.get(0) == ids.get(2)) {
                return ids.get(0);
            } else if (ids.get(1) == ids.get(2)) {
                return ids.get(1);
            }
        }
        return -1;
    }

    private static int getUnanimousID(List<Integer> ids) {
        if (ids == null) return -2;
        if (ids.size() == 3) {
            if (ids.get(0) == ids.get(1) && ids.get(0) == ids.get(2)) {
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


    /**
     * @param handler
     * @param bedFilePath
     * @return List of motif anchors from the provided bed file
     */
    public static GenomeWideList<SubcompartmentInterval> loadFromSubcompartmentBEDFile(ChromosomeHandler handler, String bedFilePath) {
        List<SubcompartmentInterval> anchors = new ArrayList<>();

        try {
            //BufferedReader br = ParsingUtils.openBufferedReader(bedFilePath);
            BufferedReader br = new BufferedReader(new InputStreamReader(ParsingUtils.openInputStream(bedFilePath)), MixerGlobals.bufferSize);
            anchors.addAll(parseSubcompartmentBEDFile(br, handler));
        } catch (IOException ec) {
            ec.printStackTrace();
        }

        return new GenomeWideList<>(handler, anchors);
    }

    /**
     * Methods for handling BED Files
     */

    /**
     * Helper function for actually parsing BED file
     * Ignores any attributes beyond the third column (i.e. just chr and positions are read)
     *
     * @param bufferedReader
     * @param handler
     * @return list of motifs
     * @throws IOException
     */
    private static List<SubcompartmentInterval> parseSubcompartmentBEDFile(BufferedReader bufferedReader, ChromosomeHandler handler) throws IOException {
        Set<SubcompartmentInterval> anchors = new HashSet<>();
        String nextLine;

        Map<String, Integer> allIdsToIntId = new HashMap<>();
        // 1 - A1, 2 - A2, 3 - B1, 4 - B2, 5 - B3, 6 - B4
        int counter = 7;

        int errorCount = 0;
        while ((nextLine = bufferedReader.readLine()) != null) {
            String[] tokens = Globals.tabPattern.split(nextLine);

            if (tokens.length > 3 && tokens[3].equalsIgnoreCase("NA")) {
                continue;
            }

            if (tokens[0].startsWith("chr") && tokens.length > 2) {
                // valid line
                String chr1Name = tokens[0];
                int start1 = Integer.parseInt(tokens[1]);
                int end1 = Integer.parseInt(tokens[2]);
                String id = tokens[3].toUpperCase();

                if (!allIdsToIntId.containsKey(id)) {
                    int newID = getKnownSubcompartmentType(id);
                    if (newID > 0) {
                        allIdsToIntId.put(id, newID);
                    } else {
                        allIdsToIntId.put(id, counter);
                        System.out.println(id + "  " + counter);
                        counter++;
                    }
                }
                int val = allIdsToIntId.get(id);


                Chromosome chr = handler.getChromosomeFromName(chr1Name);
                if (chr == null) {
                    if (errorCount < 10) {
                        System.out.println("Skipping line: " + nextLine);
                    } else if (errorCount == 10) {
                        System.err.println("Maximum error count exceeded.  Further errors will not be logged");
                    }

                    errorCount++;
                    continue;
                }

                anchors.add(new SubcompartmentInterval(chr.getIndex(), chr.getName(), start1, end1, val));
            }
        }
        if (anchors.size() < 1) System.err.println("BED File empty - file may have problems or error was encountered");
        bufferedReader.close();
        return new ArrayList<>(anchors);
    }

    private static Integer getKnownSubcompartmentType(String upperCaseID) {
        switch (upperCaseID) {
            case "A1":
            case "C1":
                return 1;
            case "A2":
            case "C2":
                return 2;
            case "B1":
            case "C3":
                return 3;
            case "B2":
            case "C4":
                return 4;
            case "B3":
            case "C5":
                return 5;
            case "B4":
                return 6;
            default:
                return -1;
        }
    }

    public static GenomeWideList<SubcompartmentInterval> redoAllIds(GenomeWideList<SubcompartmentInterval> intraSubcompartments) {
        AtomicInteger newIds = new AtomicInteger(1);
        intraSubcompartments.filterLists(new FeatureFilter<SubcompartmentInterval>() {
            @Override
            public List<SubcompartmentInterval> filter(String chr, List<SubcompartmentInterval> featureList) {
                List<SubcompartmentInterval> newIdIntervals = new ArrayList<>();

                for (SubcompartmentInterval interval : featureList) {
                    SubcompartmentInterval interval2 = (SubcompartmentInterval) interval.deepClone();
                    interval2.setClusterID(newIds.getAndIncrement());
                    newIdIntervals.add(interval2);
                }

                return newIdIntervals;
            }
        });
        return intraSubcompartments;
    }

    public static Map<Integer, Map<Integer, Integer>> createGoldStandardLookup(String locationHuntley) {
        Map<Integer, Map<Integer, Integer>> goldenMap = new HashMap<>();
        int resolution = 100000;

        ChromosomeHandler handler = HiCFileTools.loadChromosomes("hg19");
        GenomeWideList<SubcompartmentInterval> subcHuntley = loadFromSubcompartmentBEDFile(handler, locationHuntley);
        splitGWList(subcHuntley, resolution);

        subcHuntley.processLists(new FeatureFunction<SubcompartmentInterval>() {
            @Override
            public void process(String chr, List<SubcompartmentInterval> featureList) {
                if (featureList.size() > 0) {
                    int chrIndex = featureList.get(0).getChrIndex();
                    Map<Integer, Integer> indxToId = new HashMap<>();
                    for (SubcompartmentInterval interval : featureList) {
                        indxToId.put(interval.getX1(), interval.getClusterID());
                    }
                    goldenMap.put(chrIndex, indxToId);
                }
            }
        });
        return goldenMap;
    }
}
