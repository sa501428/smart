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

import java.io.File;
import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Paths;

@SuppressWarnings("unused")
public class EntropyCalculations {

    double entropyRatios;
    double entropyLogRatios;
    double entropyPerPixel;
    double entropyLogPerPixel;

    public EntropyCalculations(File shuffleFile, File baselineFile, File shuffleLogFile, File baselineLogFile,
                               float[][] shuffleM) {

        try {
            entropyRatios = getFileSize(shuffleFile) / getFileSize(baselineFile);
            entropyLogRatios = getFileSize(shuffleLogFile) / getFileSize(baselineLogFile);
            entropyPerPixel = getFileSize(shuffleFile) / (shuffleM.length * shuffleM[0].length);
            entropyLogPerPixel = getFileSize(shuffleLogFile) / (shuffleM.length * shuffleM[0].length);
        } catch (Exception e) {
            e.printStackTrace();
        }

        //myWriter.write(mapTypes[y].toString() + " pngEntropyRatio      : " + entropyRatios[y] + "\n");
        //myWriter.write(mapTypes[y].toString() + " pngLogEntropyRatio   : " + entropyLogRatios[y] + "\n");
        //myWriter.write(mapTypes[y].toString() + " pngEntropyPerPixel   : " + entropyPerPixel[y] + "\n");
        //myWriter.write(mapTypes[y].toString() + " pngLogEntropyPerPixel: " + entropyLogPerPixel[y] + "\n\n");


        //System.out.println("Scores: " + sumN + "  baseline: " + sumD + "  ratio:" + (sumN / sumD));
        //TTest test = new TTest();
        //System.out.println("Ttest " + test.pairedTTest(baseline, scores) / 2);
        //System.out.println("Ttest " + test.pairedTTest(baseline, scores, 0.1));
    }

    private double getFileSize(File fileName) throws IOException {
        return (double) Files.size(Paths.get(fileName.getAbsolutePath()));
    }
}
