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

package mixer.algos;

import javastraw.feature1D.GenomeWide1DList;
import javastraw.reader.Dataset;
import javastraw.reader.basics.ChromosomeHandler;
import javastraw.reader.type.NormalizationHandler;
import javastraw.reader.type.NormalizationType;
import javastraw.tools.HiCFileTools;
import javastraw.tools.UNIXTools;
import mixer.clt.CommandLineParserForMixer;
import mixer.clt.MixerCLT;
import mixer.utils.BedTools;
import mixer.utils.shuffle.Shuffle;
import mixer.utils.tracks.SubcompartmentInterval;

import java.io.File;
import java.util.List;
import java.util.Random;

public class Chic extends MixerCLT {

    private final Random generator = new Random(22871L);
    protected NormalizationType norm = NormalizationHandler.INTER_KR;
    private Dataset ds;
    private int resolution = 100000;
    private File outputDirectory;
    private String[] prefix;
    private String[] referenceBedFiles;
    private int numSplits = 2;

    public Chic() {
        super("chic [-r resolution] [-k NONE/INTER_KR/INTER_SCALE] [--verbose] " +
                "<file.hic> <outfolder> <file1.bed,file2.bed,...> <name1,name2,...>");
    }

    @Override
    protected void readMixerArguments(String[] args, CommandLineParserForMixer mixerParser) {
        if (args.length != 5) {
            printUsageAndExit(51);
        }

        ds = HiCFileTools.extractDatasetForCLT(args[1], false, false, true);
        outputDirectory = HiCFileTools.createValidDirectory(args[2]);
        referenceBedFiles = args[3].split(",");
        prefix = args[4].split(",");
        if (referenceBedFiles.length != prefix.length) {
            System.err.println("Number of bed files should match number of prefixes");
            printUsageAndExit(53);
        }

        NormalizationType preferredNorm = mixerParser.getNormalizationTypeOption(ds.getNormalizationHandler());
        if (preferredNorm != null) norm = preferredNorm;

        List<Integer> possibleResolutions = mixerParser.getMultipleResolutionOptions();
        if (possibleResolutions != null) {
            if (possibleResolutions.size() > 1)
                System.err.println("Only one resolution can be specified\nUsing " + possibleResolutions.get(0));
            resolution = possibleResolutions.get(0);
        }
        numSplits = mixerParser.getWindowSizeOption(2);
        updateGeneratorSeed(mixerParser, generator);
    }

    @Override
    public void run() {

        ChromosomeHandler chromosomeHandler = ds.getChromosomeHandler();
        UNIXTools.makeDir(outputDirectory);

        for (int i = 0; i < referenceBedFiles.length; i++) {
            GenomeWide1DList<SubcompartmentInterval> subcompartments =
                    BedTools.loadBedFile(chromosomeHandler, referenceBedFiles[i]);
            System.out.println("Processing " + prefix[i]);
            Shuffle matrix = new Shuffle(ds, norm, resolution);
            matrix.runGWStats(subcompartments, outputDirectory, prefix[i], numSplits);
            matrix = null;
            subcompartments = null;
        }
        System.out.println("Shuffle complete");
    }
}
