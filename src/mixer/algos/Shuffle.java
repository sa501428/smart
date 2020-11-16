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

package mixer.algos;

import javastraw.featurelist.GenomeWideList;
import javastraw.reader.ChromosomeHandler;
import javastraw.reader.Dataset;
import javastraw.reader.HiCFileTools;
import javastraw.tools.UNIXTools;
import javastraw.type.NormalizationType;
import mixer.clt.CommandLineParserForMixer;
import mixer.clt.MixerCLT;
import mixer.utils.custommetrics.RobustCosineMetric;
import mixer.utils.shuffle.ShuffleMatrix;
import mixer.utils.slice.structures.SliceUtils;
import mixer.utils.slice.structures.SubcompartmentInterval;
import tagbio.umap.metric.Metric;

import java.io.File;
import java.util.ArrayList;
import java.util.List;
import java.util.Random;

/**
 * experimental code
 * <p>
 * Created by muhammadsaadshamim on 9/14/15.
 */
public class Shuffle extends MixerCLT {

    public static Metric metric = RobustCosineMetric.SINGLETON;
    private final List<String> inputHicFilePaths = new ArrayList<>();
    private final Random generator = new Random(22871L);
    private Dataset ds;
    private int resolution = 100000;
    private int compressionFactor = 8;
    private File outputDirectory;
    private String[] prefix;
    private String[] referenceBedFiles;

    // subcompartment lanscape identification via clustering enrichment
    public Shuffle() {
        super("shuffle [-r resolution] [-k NONE/VC/VC_SQRT/KR/SCALE] [-w window] [--verbose] " +
                "<file.hic> <subcompartment.bed> <outfolder> <prefix>");
    }

    @Override
    protected void readMixerArguments(String[] args, CommandLineParserForMixer mixerParser) {
        if (args.length != 5) {
            printUsageAndExit(51);
        }

        ds = HiCFileTools.extractDatasetForCLT(args[1], true, false);

        referenceBedFiles = args[2].split(",");
        outputDirectory = HiCFileTools.createValidDirectory(args[3]);
        prefix = args[4].split(",");

        NormalizationType preferredNorm = mixerParser.getNormalizationTypeOption(ds.getNormalizationHandler());
        if (preferredNorm != null) norm = preferredNorm;

        List<Integer> possibleResolutions = mixerParser.getMultipleResolutionOptions();
        if (possibleResolutions != null) {
            if (possibleResolutions.size() > 1)
                System.err.println("Only one resolution can be specified\nUsing " + possibleResolutions.get(0));
            resolution = possibleResolutions.get(0);
        }

        long[] possibleSeeds = mixerParser.getMultipleSeedsOption();
        if (possibleSeeds != null && possibleSeeds.length > 0) {
            for (long seed : possibleSeeds) {
                generator.setSeed(seed);
            }
        }

        int minSize = mixerParser.getAPAWindowSizeOption();
        if (minSize > 1) {
            compressionFactor = minSize;
        }
    }

    @Override
    public void run() {

        ChromosomeHandler chromosomeHandler = ds.getChromosomeHandler();
        if (givenChromosomes != null)
            chromosomeHandler = HiCFileTools.stringToChromosomes(givenChromosomes, chromosomeHandler);

        for (int i = 0; i < referenceBedFiles.length; i++) {
            GenomeWideList<SubcompartmentInterval> subcompartments =
                    SliceUtils.loadFromSubcompartmentBEDFile(chromosomeHandler, referenceBedFiles[i]);
            System.out.println("Processing " + prefix[i]);
            File newFolder = new File(outputDirectory, prefix[i]);
            UNIXTools.makeDir(newFolder);
            ShuffleMatrix matrix = new ShuffleMatrix(ds, norm, resolution, compressionFactor);
            matrix.runAnalysis(subcompartments, newFolder);
            matrix.savePlotsAndResults(newFolder);
        }
        System.out.println("Shuffle complete");
    }
}
