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

package mixer.algos;

import mixer.MixerGlobals;
import mixer.clt.CommandLineParserForMixer;
import mixer.clt.MixerCLT;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.FileReader;
import java.io.FileWriter;
import java.util.HashMap;
import java.util.Map;
import java.util.regex.Pattern;

public class Rename extends MixerCLT {
    Map<String, String> oldNameToNewName = new HashMap<>();
    Map<String, String> oldNameToNewColor = new HashMap<>();
    private String inputBedFile, outputBedFile;
    private String[] changes;

    public Rename() {
        super("rename <old_name,new_name,R,G,B;A,B,0,0,0;B,A,1,1,1> <subcompartment.bed> <new_subcompartment.bed>");
    }

    @Override
    protected void readMixerArguments(String[] args, CommandLineParserForMixer mixerParser) {
        if (args.length != 4) {
            printUsageAndExit(51);
        }

        changes = args[1].split(";");
        inputBedFile = args[2];
        outputBedFile = args[3];
    }

    @Override
    public void run() {
        populateChanges();
        try {
            parseAndWrite();
        } catch (Exception e) {
            e.printStackTrace();
        }
    }

    private void populateChanges() {
        for (String change : changes) {
            String[] values = change.split(",");
            String oldID = values[0].toUpperCase();
            String newID = values[1].toUpperCase();
            String rgb = values[2] + "," + values[3] + "," + values[4];
            oldNameToNewName.put(oldID, newID);
            oldNameToNewColor.put(oldID, rgb);
        }
    }

    private void parseAndWrite() throws Exception {

        BufferedReader bufferedReader = new BufferedReader(new FileReader(inputBedFile), MixerGlobals.bufferSize);
        BufferedWriter bufferedWriter = new BufferedWriter(new FileWriter(outputBedFile));

        String nextLine;
        while ((nextLine = bufferedReader.readLine()) != null) {
            String[] tokens = Pattern.compile("\t").split(nextLine);

            if (nextLine.startsWith("#")) {
                bufferedWriter.write(nextLine);
                continue;
            }

            if (tokens.length > 3 && tokens[3].equalsIgnoreCase("NA")) {
                continue;
            }

            if (tokens[0].startsWith("chr") && tokens.length > 2) {
                String oldID = tokens[3].toUpperCase();
                String newID = oldNameToNewName.get(oldID);
                String newColor = oldNameToNewColor.get(oldID);
                String outLine = generateOutputLine(tokens, newID, newColor);
                bufferedWriter.write(outLine);
            }
        }
    }

    private String generateOutputLine(String[] tokens, String newID, String newColor) {
        return tokens[0] + "\t" + tokens[1] + "\t" + tokens[2] +
                "\t" + newID + "\t" + newID + "\t.\t" +
                tokens[1] + "\t" + tokens[2] + newColor;
    }
}
