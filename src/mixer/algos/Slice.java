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
import javastraw.reader.basics.Chromosome;
import javastraw.reader.basics.ChromosomeHandler;
import javastraw.reader.norm.NormalizationPicker;
import javastraw.reader.type.NormalizationType;
import javastraw.tools.HiCFileTools;
import mixer.SmartTools;
import mixer.clt.CommandLineParserForMixer;
import mixer.clt.MixerCLT;
import mixer.utils.cleaning.BadIndexFinder;
import mixer.utils.cleaning.MatrixPreprocessor;
import mixer.utils.drive.BinMappings;
import mixer.utils.drive.FinalMatrix;
import mixer.utils.drive.MatrixAndWeight;
import mixer.utils.drive.MatrixBuilder;
import mixer.utils.intra.IndexOrderer;
import mixer.utils.kmeans.ClusteringMagic;
import mixer.utils.refinement.InternalShuffle;
import mixer.utils.tracks.SubcompartmentInterval;
import mixer.utils.translocations.SimpleTranslocationFinder;

import java.io.File;
import java.util.List;
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
    private final Random generator = new Random(22871L);
    private int resolution = 100000;
    private Dataset ds;
    private File outputDirectory;
    private NormalizationType[] norms;
    boolean useExpandedIntraOE = true;

    // subcompartment landscape identification via compressing enrichments
    public Slice() {
        super("slice [-r resolution] [--verbose] [--both-norms] [--append] " +
                "<file.hic> <K0,KF> <outfolder>\n" +
                "   K0 - minimum number of clusters\n" +
                "   KF - maximum number of clusters");
        //"<-k NONE/VC/VC_SQRT/KR/SCALE> [--compare reference.bed] [--has-translocation] " +
    }

    @Override
    protected void readMixerArguments(String[] args, CommandLineParserForMixer mixerParser) {
        if (args.length != 4) {
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

        outputDirectory = HiCFileTools.createValidDirectory(args[3]);
        norms = populateNormalizations(ds);
        updateGeneratorSeed(mixerParser, generator);
    }

    private NormalizationType[] populateNormalizations(Dataset ds) {
        NormalizationType[] norms = new NormalizationType[2];
        norms[INTRA_SCALE_INDEX] = NormalizationPicker.getFirstValidNormInThisOrder(ds, new String[]{"SCALE", "KR", "VC"});
        norms[INTER_SCALE_INDEX] = NormalizationPicker.getFirstValidNormInThisOrder(ds, new String[]{"INTER_SCALE", "INTER_KR", "INTER_VC"});
        System.out.println("Using normalizations: " + norms[INTRA_SCALE_INDEX].getLabel() + " and " + norms[INTER_SCALE_INDEX].getLabel());
        return norms;
    }

    @Override
    public void run() {
        ChromosomeHandler handler = ds.getChromosomeHandler();
        Chromosome[] chromosomes = handler.getAutosomalChromosomesArray();

        Map<Integer, Set<Integer>> badIndices = BadIndexFinder.getBadIndices(ds, chromosomes, resolution);
        SimpleTranslocationFinder translocations = new SimpleTranslocationFinder(ds, norms, outputDirectory,
                badIndices, resolution);
        BinMappings mappings = IndexOrderer.getInitialMappings(ds, chromosomes, resolution,
                badIndices, norms[INTRA_SCALE_INDEX], generator.nextLong(), outputDirectory,
                useExpandedIntraOE);

        MatrixAndWeight slice0 = MatrixBuilder.populateMatrix(ds, chromosomes, resolution,
                norms[INTER_SCALE_INDEX], norms[INTRA_SCALE_INDEX], mappings, translocations, outputDirectory);

        runWithSettings(slice0, handler, chromosomes,
                true, false, false, true, false, false);
        System.out.println("\nSLICE complete");
    }

    private void runWithSettings(MatrixAndWeight slice0, ChromosomeHandler handler, Chromosome[] chromosomes,
                                 boolean includeIntra, boolean useLog, boolean useBothNorms,
                                 boolean appendIntra, boolean useRowZ, boolean shouldRegularize) {
        String stem = getNewPrefix(includeIntra, useLog, useBothNorms, appendIntra, useRowZ, shouldRegularize);
        FinalMatrix slice = MatrixPreprocessor.preprocess(slice0.deepCopy(), chromosomes, includeIntra, useLog,
                useBothNorms, appendIntra, useRowZ, shouldRegularize);

        if (slice.notEmpty()) {
            if (SmartTools.printVerboseComments) {
                slice.export(outputDirectory, stem);
            }
            ClusteringMagic clustering = new ClusteringMagic(slice, outputDirectory, handler, generator.nextLong());
            Map<Integer, List<String>> bedFiles = clustering.extractFinalGWSubcompartments(stem);
            System.out.print("*");
            Map<Integer, GenomeWide1DList<SubcompartmentInterval>> bestClusterings =
                    InternalShuffle.determineBest(bedFiles, resolution,
                            handler, ds, norms[INTER_SCALE_INDEX]);
            for (int k : bestClusterings.keySet()) {
                File outBedFile = new File(outputDirectory, stem + "_k" + k + "_best_clusters.bed"); // "_wcss" + wcss +
                bestClusterings.get(k).simpleExport(outBedFile);
            }
            System.out.println(">>");
        }
    }

    private String getNewPrefix(boolean includeIntra, boolean useLog, boolean useBothNorms, boolean appendIntra,
                                boolean useRowZ, boolean shouldRegularize) {
        String stem = "SLICE";
        /*
        if (includeIntra) {
            stem += "_GW";
        } else {
            stem += "_INTER";
        }
        if (useLog) {
            stem += "_LOG";
        } else {
            stem += "_EXP";
        }

        if (shouldRegularize) {
            stem += "_REG";
        } else {
            stem += "_IRR";
        }
        if (useBothNorms) {
            stem += "_2NORM";
        } else {
            stem += "_1NORM";
        }
        if (appendIntra) {
            stem += "_APPEND";
        } else {
            stem += "_INSERT";
        }
        if (useRowZ) {
            stem += "_RowZ";
        } else {
            stem += "_ColZ";
        }
         */
        return stem;
    }
}
