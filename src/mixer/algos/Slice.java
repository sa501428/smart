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

import javastraw.reader.Dataset;
import javastraw.reader.basics.Chromosome;
import javastraw.reader.basics.ChromosomeHandler;
import javastraw.reader.norm.NormalizationPicker;
import javastraw.reader.type.NormalizationType;
import javastraw.tools.HiCFileTools;
import mixer.clt.CommandLineParserForMixer;
import mixer.clt.MixerCLT;
import mixer.utils.cleaning.BadIndexFinder;
import mixer.utils.cleaning.MatrixPreprocessor;
import mixer.utils.drive.BinMappings;
import mixer.utils.drive.MatrixAndWeight;
import mixer.utils.drive.MatrixBuilder;
import mixer.utils.intra.IndexOrderer;
import mixer.utils.kmeans.ClusteringMagic;
import mixer.utils.translocations.SimpleTranslocationFinder;

import java.io.File;
import java.util.Map;
import java.util.Random;
import java.util.Set;

/**
 * experimental code
 * <p>
 * Created by muhammadsaadshamim on 9/14/15.
 */
public class Slice extends MixerCLT {

    public static final int INTRA_SCALE_INDEX = 0;
    public static final int INTER_SCALE_INDEX = 1;
    private boolean useScale = false;
    private final Random generator = new Random(22871L);
    private int resolution = 100000;
    private Dataset ds;
    private File outputDirectory;
    private NormalizationType[] norms;
    private String prefix = "";

    // subcompartment landscape identification via compressing enrichments
    public Slice() {
        super("slice [-r resolution] [--verbose] [--scale] " +
                //"<-k NONE/VC/VC_SQRT/KR/SCALE> [--compare reference.bed] [--has-translocation] " +
                "<file.hic> <K0,KF> <outfolder> <prefix_>\n" +
                "   K0 - minimum number of clusters\n" +
                "   KF - maximum number of clusters");
    }

    @Override
    protected void readMixerArguments(String[] args, CommandLineParserForMixer mixerParser) {
        if (args.length != 5) {
            printUsageAndExit(5);
        }

        resolution = updateResolution(mixerParser, resolution);
        ds = HiCFileTools.extractDatasetForCLT(args[1], true, false, resolution > 100);

        try {
            String[] valString = args[2].split(",");
            ClusteringMagic.startingClusterSizeK = Integer.parseInt(valString[0]);
            ClusteringMagic.numClusterSizeKValsUsed = Integer.parseInt(valString[1]) - ClusteringMagic.startingClusterSizeK;
        } catch (Exception e) {
            printUsageAndExit(5);
        }

        useScale = mixerParser.getScaleOption();
        outputDirectory = HiCFileTools.createValidDirectory(args[3]);
        prefix = args[4];
        norms = populateNormalizations(ds);

        updateGeneratorSeed(mixerParser, generator);
    }

    private NormalizationType[] populateNormalizations(Dataset ds) {
        NormalizationType[] norms = new NormalizationType[2];
        norms[INTRA_SCALE_INDEX] = NormalizationPicker.getFirstValidNormInThisOrder(ds, new String[]{"SCALE", "KR"});
        norms[INTER_SCALE_INDEX] = NormalizationPicker.getFirstValidNormInThisOrder(ds, new String[]{"INTER_SCALE", "INTER_KR"});
        return norms;
    }

    @Override
    public void run() {
        ChromosomeHandler handler = ds.getChromosomeHandler();
        Chromosome[] chromosomes = handler.getAutosomalChromosomesArray();

        Map<Integer, Set<Integer>> badIndices = BadIndexFinder.getBadIndices(ds, chromosomes, resolution);

        // todo should be at lower res
        SimpleTranslocationFinder translocations = new SimpleTranslocationFinder(ds, norms, outputDirectory,
                badIndices, resolution);

        // todo should be at lower res
        BinMappings mappings = IndexOrderer.getInitialMappings(ds, chromosomes, resolution,
                badIndices, norms[INTRA_SCALE_INDEX], generator.nextLong(), outputDirectory);

        MatrixAndWeight slice0 = MatrixBuilder.populateMatrix(ds, chromosomes, resolution,
                norms[INTER_SCALE_INDEX], mappings, translocations, outputDirectory, useScale);

        slice0.export(outputDirectory, "pre-clean");
        runWithSettings2(slice0, false, handler, chromosomes);
        runWithSettings2(slice0, true, handler, chromosomes);
        System.out.println("\nSLICE complete");
    }

    private void runWithSettings(MatrixAndWeight slice0, int cutoff, boolean useExp, boolean useZscore,
                                 ChromosomeHandler handler, Chromosome[] chromosomes) {
        MatrixAndWeight slice = MatrixPreprocessor.clean(slice0.deepCopy(), chromosomes, cutoff, useExp, useZscore);
        //slice.export(outputDirectory, "slice");
        if (slice.notEmpty()) {
            ClusteringMagic clustering = new ClusteringMagic(slice, outputDirectory, handler, generator.nextLong());
            clustering.extractFinalGWSubcompartments(getNewPrefix(prefix, cutoff, useExp, useZscore));
            System.out.println("*");
        }
    }

    private void runWithSettings2(MatrixAndWeight slice0, boolean useLog,
                                  ChromosomeHandler handler, Chromosome[] chromosomes) {
        String stem = getNewPrefix2(prefix, useLog);
        MatrixAndWeight slice = MatrixPreprocessor.clean2(slice0.deepCopy(), chromosomes, useLog, stem, outputDirectory);
        if (slice.notEmpty()) {
            ClusteringMagic clustering = new ClusteringMagic(slice, outputDirectory, handler, generator.nextLong());
            clustering.extractFinalGWSubcompartments(stem);
            System.out.println("*");
        }
    }

    private String getNewPrefix(String prefix, int cutoff, boolean useExp, boolean useZscore) {
        String newPrefix = prefix + "_c" + cutoff + "_";
        if (useExp) {
            newPrefix += "exp_";
        } else {
            newPrefix += "log_";
        }
        if (useZscore) {
            newPrefix += "Zscore_";
        } else {
            newPrefix += "noZ_";
        }
        return newPrefix;
    }

    private String getNewPrefix2(String prefix, boolean useLog) {
        String newPrefix = prefix + "_";
        if (useLog) {
            newPrefix += "log_";
        } else {
            newPrefix += "exp_";
        }
        return newPrefix;
    }
}
