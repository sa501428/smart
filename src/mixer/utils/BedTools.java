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

package mixer.utils;

import javastraw.feature1D.GenomeWide1DList;
import javastraw.reader.basics.Chromosome;
import javastraw.reader.basics.ChromosomeHandler;
import mixer.SmartTools;
import mixer.utils.tracks.SubcompartmentInterval;

import java.io.BufferedReader;
import java.io.FileReader;
import java.io.IOException;
import java.util.*;
import java.util.concurrent.atomic.AtomicInteger;
import java.util.regex.Pattern;

public class BedTools {

    public static GenomeWide1DList<SubcompartmentInterval> loadBedFile(ChromosomeHandler handler, String bedFilePath) {
        List<SubcompartmentInterval> anchors = new ArrayList<>();

        try {
            //BufferedReader br = ParsingUtils.openBufferedReader(bedFilePath);
            anchors.addAll(parseSubcompartmentBEDFile(bedFilePath, handler));
        } catch (IOException ec) {
            ec.printStackTrace();
        }

        return new GenomeWide1DList<>(handler, anchors);
    }

    /**
     * Helper function for actually parsing BED file
     * Ignores any attributes beyond the 5th column (i.e. just chr, positions, and subcompartment are read)
     *
     * @return list of motifs
     */
    private static List<SubcompartmentInterval> parseSubcompartmentBEDFile(String bedFilePath, ChromosomeHandler handler) throws IOException {
        BufferedReader bufferedReader = new BufferedReader(new FileReader(bedFilePath), SmartTools.bufferSize);

        final AtomicInteger idCounter = new AtomicInteger(0);

        Set<SubcompartmentInterval> anchors = new HashSet<>();
        String nextLine;
        Map<String, Integer> idToVal = new HashMap<>();

        int errorCount = 0;
        while ((nextLine = bufferedReader.readLine()) != null) {
            String[] tokens = Pattern.compile("\t").split(nextLine);

            if (tokens.length > 3 && tokens[3].equalsIgnoreCase("NA")) {
                continue;
            }

            if (tokens[0].startsWith("chr") && tokens.length > 4) {
                // valid line
                String chr1Name = tokens[0];
                int start1 = Integer.parseInt(tokens[1]);
                int end1 = Integer.parseInt(tokens[2]);
                String id = tokens[3].toUpperCase();

                int val = updateMapAndConfirmVal(id, idToVal, idCounter);

                Chromosome chrom = handler.getChromosomeFromName(chr1Name);
                if (chrom == null) {
                    if (errorCount < 3) {
                        System.out.println("Skipping line: " + nextLine);
                    } else if (errorCount == 3) {
                        System.err.println("Maximum error count exceeded.  Further errors will not be logged");
                        System.err.println("If chromosomes were intentionally filtered, ignore this error.");
                    }
                    errorCount++;
                    continue;
                }

                anchors.add(new SubcompartmentInterval(chrom, start1, end1, val, id));
            } else {
                System.out.println("Skipping line: " + nextLine);
            }
        }
        if (anchors.size() < 1) System.err.println("BED File empty - file may have problems or error was encountered");
        bufferedReader.close();
        return new ArrayList<>(anchors);
    }

    private static int updateMapAndConfirmVal(String id, Map<String, Integer> idToVal, AtomicInteger idCounter) {
        if (!idToVal.containsKey(id)) {
            idToVal.put(id, idCounter.incrementAndGet());
        }
        return idToVal.get(id);
    }


    public static void splitGWList(GenomeWide1DList<SubcompartmentInterval> intraSubcompartments, int width) {
        intraSubcompartments.filterLists((chr, featureList) -> splitSubcompartmentIntervals(featureList, width));
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

    public static GenomeWide1DList<SubcompartmentInterval> loadBedFileAtResolution(ChromosomeHandler handler, String path, int resolution) {
        GenomeWide1DList<SubcompartmentInterval> output = loadBedFile(handler, path);
        splitGWList(output, resolution);
        return output;
    }
}
