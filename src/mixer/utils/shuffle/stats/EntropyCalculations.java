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

package mixer.utils.shuffle.stats;

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
    }

    private double getFileSize(File fileName) throws IOException {
        return (double) Files.size(Paths.get(fileName.getAbsolutePath()));
    }
}
