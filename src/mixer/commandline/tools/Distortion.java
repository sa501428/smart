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

package mixer.commandline.tools;

import mixer.commandline.handling.CommandLineParserForMixer;
import mixer.commandline.handling.MixerCLT;
import mixer.commandline.utils.mixer.DistortionGenerator;
import mixer.commandline.utils.mixer.DistortionMixedGenerator;
import mixer.data.ChromosomeHandler;
import mixer.data.Dataset;
import mixer.data.HiCFileTools;
import mixer.windowui.NormalizationType;
import org.broad.igv.feature.Chromosome;
import org.broad.igv.util.Pair;

import java.io.File;
import java.util.List;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;

/**
 * Generating Regions of Interest for Network Discovery
 */

public class Distortion extends MixerCLT {

    private Dataset ds, ds2;
    private int x, y, z, stride, resolution;
    private File outputDirectory;
    private boolean onlyOneFileBeingUsed = true;

    public Distortion() {
        super("distort [-k NONE/KR/VC/VC_SQRT] [-r resolution] [--stride increment] " +
                "<hic file> <size,manipulations,examples> <directory>");
    }

    @Override
    protected void readMixerArguments(String[] args, CommandLineParserForMixer mixerParser) {
        if (args.length != 4) {
            printUsageAndExit();
        }

        String[] files = args[1].split("\\+");
        onlyOneFileBeingUsed = files.length == 1;
        ds = HiCFileTools.extractDatasetForCLT(files[0], true);

        if (!onlyOneFileBeingUsed) {
            ds2 = HiCFileTools.extractDatasetForCLT(files[1], true);
        }

        // split on commas
        // save the dimensions
        String[] dimensions = args[2].split(",");
        x = Integer.parseInt(dimensions[0]);
        y = Integer.parseInt(dimensions[1]);
        z = Integer.parseInt(dimensions[2]);

        stride = mixerParser.getStride();
        outputDirectory = HiCFileTools.createValidDirectory(args[3]);

        NormalizationType preferredNorm = mixerParser.getNormalizationTypeOption(ds.getNormalizationHandler());
        if (preferredNorm != null) norm = preferredNorm;

        List<String> possibleResolutions = mixerParser.getMultipleResolutionOptions();
        if (possibleResolutions != null) {
            for (String num : possibleResolutions) {
                resolution = Integer.parseInt(num);
            }
        }
    }

    @Override
    public void run() {

        ChromosomeHandler chromosomeHandler = ds.getChromosomeHandler();

        if (givenChromosomes != null)
            chromosomeHandler = HiCFileTools.stringToChromosomes(givenChromosomes, chromosomeHandler);

        ExecutorService executor = Executors.newFixedThreadPool(Runtime.getRuntime().availableProcessors());

        for (Pair<Chromosome, Chromosome> chromosomePair : chromosomeHandler.getAutosomalPairs()) {
            Runnable worker = new Runnable() {
                @Override
                public void run() {
                    if (onlyOneFileBeingUsed) {
                        DistortionGenerator generator = new DistortionGenerator(ds, chromosomePair, x, y, z, resolution, norm, stride, outputDirectory);
                        generator.makeExamples();
                    } else {
                        DistortionMixedGenerator generator = new DistortionMixedGenerator(ds, ds2, chromosomePair, x, y, z, resolution, norm, stride, outputDirectory);
                        generator.makeExamples();
                    }
                }
            };
            executor.execute(worker);
        }
        executor.shutdown();

        // Wait until all threads finish
        while (!executor.isTerminated()) {
        }
    }

}
