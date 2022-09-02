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

package mixer.utils.bed;

import javastraw.feature1D.GenomeWide1DList;
import javastraw.reader.basics.Chromosome;
import javastraw.reader.basics.ChromosomeHandler;
import mixer.SmartTools;
import mixer.utils.tracks.SubcompartmentInterval;

import java.io.BufferedReader;
import java.io.FileReader;
import java.io.IOException;
import java.util.ArrayList;
import java.util.HashSet;
import java.util.List;
import java.util.Set;
import java.util.regex.Pattern;

public class SubcompartmentBedFile {

    private final GenomeWide1DList<SubcompartmentInterval> featureGenomeWide1DList;

    public SubcompartmentBedFile(String bedpath, ChromosomeHandler handler) {
        featureGenomeWide1DList = populateBedFile(handler, bedpath);
    }

    public static GenomeWide1DList<SubcompartmentInterval> populateBedFile(ChromosomeHandler handler, String bedFilePath) {
        List<SubcompartmentInterval> anchors = new ArrayList<>();

        try {
            anchors.addAll(parseBEDFile(bedFilePath, handler));
        } catch (IOException ec) {
            ec.printStackTrace();
        }

        return new GenomeWide1DList<>(handler, anchors);
    }

    private static List<SubcompartmentInterval> parseBEDFile(String bedFilePath, ChromosomeHandler handler) throws IOException {
        BufferedReader bufferedReader = new BufferedReader(new FileReader(bedFilePath), SmartTools.bufferSize);

        Set<SubcompartmentInterval> positions = new HashSet<>();
        String nextLine;

        int errorCount = 0;
        while ((nextLine = bufferedReader.readLine()) != null) {
            String[] tokens = Pattern.compile("\t").split(nextLine);

            if (nextLine.startsWith("#")) {
                continue;
            }

            if (tokens.length > 4) {
                // valid line
                String chr1Name = tokens[0];
                int start1 = Integer.parseInt(tokens[1]);
                int end1 = Integer.parseInt(tokens[2]);
                int id = Integer.parseInt(tokens[4]);

                Chromosome chr = handler.getChromosomeFromName(chr1Name);
                if (chr == null) {
                    if (errorCount < 3) {
                        System.err.println("Skipping line: " + nextLine);
                    } else if (errorCount == 3) {
                        System.err.println("Maximum error count exceeded.  Further errors will not be logged");
                        System.err.println("If chromosomes were intentionally filtered, ignore this error.");
                    }
                    errorCount++;
                    continue;
                }
                positions.add(new SubcompartmentInterval(chr, start1, end1, id));
            }
        }
        if (positions.size() < 1)
            System.err.println("BED File empty - file may have problems or error was encountered");
        bufferedReader.close();
        return new ArrayList<>(positions);
    }

    public List<SubcompartmentInterval> get(int index) {
        return featureGenomeWide1DList.getFeatures("" + index);
    }

    public int getNumTotalFeatures() {
        return featureGenomeWide1DList.size();
    }
}
