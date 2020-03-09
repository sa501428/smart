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

package mixer.tools.clt;

import mixer.data.ChromosomeHandler;
import mixer.data.Dataset;
import mixer.data.HiCFileTools;
import mixer.data.Matrix;
import mixer.windowui.NormalizationHandler;
import mixer.windowui.NormalizationType;
import org.broad.igv.feature.Chromosome;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

/**
 * Created by muhammadsaadshamim on 9/21/15.
 */
public abstract class MixerCLT {

    protected NormalizationType norm = NormalizationHandler.KR;
    protected List<String> givenChromosomes = null; //TODO set to protected
    protected static int numCPUThreads = 1;
    private static String usage;
    protected Dataset dataset = null;

    protected MixerCLT(String usage) {
        setUsage(usage);
    }

    protected int determineHowManyChromosomesWillActuallyRun(Dataset ds, ChromosomeHandler chromosomeHandler) {
        int maxProgressStatus = 0;
        for (Chromosome chr : chromosomeHandler.getChromosomeArrayWithoutAllByAll()) {
            Matrix matrix = ds.getMatrix(chr, chr);
            if (matrix == null) continue;
            maxProgressStatus++;
        }
        return maxProgressStatus;
    }

    public void readArguments(String[] args, CommandLineParserForMixer parser) {
        CommandLineParserForMixer mixerParser = parser;
        assessIfChromosomesHaveBeenSpecified(mixerParser);
        readMixerArguments(args, mixerParser);
    }

    protected void updateNumberOfCPUThreads(CommandLineParserForMixer mixerParser) {
        int numThreads = mixerParser.getNumThreads();
        if (numThreads > 0) {
            numCPUThreads = numThreads;
        } else if (numThreads < 0) {
            numCPUThreads = Runtime.getRuntime().availableProcessors();
        } else {
            numCPUThreads = 1;
        }
        System.out.println("Using " + numCPUThreads + " CPU thread(s)");
    }

    protected abstract void readMixerArguments(String[] args, CommandLineParserForMixer mixerParser);

    public static String[] splitToList(String nextLine) {
        return nextLine.trim().split("\\s+");
    }

    private void assessIfChromosomesHaveBeenSpecified(CommandLineParserForMixer mixerParser) {
        List<String> possibleChromosomes = mixerParser.getChromosomeListOption();
        if (possibleChromosomes != null && possibleChromosomes.size() > 0) {
            givenChromosomes = new ArrayList<>(possibleChromosomes);
        }
    }

    public abstract void run();

    private void setUsage(String newUsage) {
        usage = newUsage;
    }

    public void printUsageAndExit() {
        System.out.println("Usage:   mixer_tools " + usage);
        System.exit(0);
    }

    public void printUsageAndExit(int exitcode) {
        System.out.println("Usage:   mixer_tools " + usage);
        System.exit(exitcode);
    }

    protected void setDatasetAndNorm(String files, String normType, boolean allowPrinting) {
        dataset = HiCFileTools.extractDatasetForCLT(Arrays.asList(files.split("\\+")), allowPrinting);

        norm = dataset.getNormalizationHandler().getNormTypeFromString(normType);
        if (norm == null) {
            System.err.println("Normalization type " + norm + " unrecognized.  Normalization type must be one of \n" +
                    "\"NONE\", \"VC\", \"VC_SQRT\", \"KR\", \"GW_KR\"," +
                    " \"GW_VC\", \"INTER_KR\", \"INTER_VC\", or a custom added normalization.");
            System.exit(16);
        }
    }
}
