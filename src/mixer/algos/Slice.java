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
import javastraw.reader.type.NormalizationHandler;
import javastraw.reader.type.NormalizationType;
import javastraw.tools.HiCFileTools;
import javastraw.tools.UNIXTools;
import mixer.SmartTools;
import mixer.clt.CommandLineParserForMixer;
import mixer.clt.MixerCLT;
import mixer.utils.cleaning.BadIndexFinder;
import mixer.utils.cleaning.MatrixPreprocessor;
import mixer.utils.drive.FinalMatrix;
import mixer.utils.drive.Mappings;
import mixer.utils.drive.MatrixAndWeight;
import mixer.utils.drive.MatrixBuilder;
import mixer.utils.intra.IndexOrderer;
import mixer.utils.kmeans.ClusteringMagic;
import mixer.utils.refinement.InternalShuffle;
import mixer.utils.tracks.SubcompartmentInterval;
import mixer.utils.translocations.SimpleTranslocationFinder;
import mixer.utils.translocations.TranslocationSet;

import java.io.File;
import java.util.*;

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
    private File parentDirectory;
    private NormalizationType[] norms;
    boolean useExpandedIntraOE = true, appendIntra = true;
    private boolean doPostNorm = false;
    private NormalizationType internalShuffleNorm = NormalizationHandler.INTER_VC;
    private boolean skipCheck = false;

    // subcompartment landscape identification via compressing enrichments
    public Slice() {
        super("slice [-r resolution] [--post-norm] [--skip-check] [--verbose] [-k INTRA_NORM,INTER_NORM] " +
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

        if (resolution < 25000) {
            SmartTools.NUM_ENTRIES_TO_SKIP_MEDIAN = 30000 / resolution;
        }

        try {
            String[] valString = args[2].split(",");
            ClusteringMagic.startingClusterSizeK = Integer.parseInt(valString[0]);
            ClusteringMagic.numClusterSizeKValsUsed = Integer.parseInt(valString[1]) - ClusteringMagic.startingClusterSizeK;
        } catch (Exception e) {
            printUsageAndExit(5);
        }

        parentDirectory = HiCFileTools.createValidDirectory(args[3]);
        norms = populateNormalizations(ds, mixerParser);
        internalShuffleNorm = norms[INTER_SCALE_INDEX];
        skipCheck = mixerParser.getSkipCheckOption();
        if (mixerParser.getPostNormOption()) {
            doPostNorm = true;
            norms[INTER_SCALE_INDEX] = NormalizationHandler.NONE;
        }
        System.out.println("Using normalizations: " + norms[INTRA_SCALE_INDEX].getLabel() + " and " + norms[INTER_SCALE_INDEX].getLabel());

        updateGeneratorSeed(mixerParser, generator);
    }

    public static NormalizationType[] populateNormalizations(Dataset ds, CommandLineParserForMixer parser) {

        NormalizationType[] potentialNorms = parser.getTwoNormsTypeOption(ds.getNormalizationHandler());
        if (potentialNorms != null && potentialNorms.length == 2) {
            if (potentialNorms[0] != null && potentialNorms[1] != null) {
                return potentialNorms;
            }
        }

        NormalizationType[] norms = new NormalizationType[2];
        norms[INTRA_SCALE_INDEX] = NormalizationPicker.getFirstValidNormInThisOrder(ds, new String[]{"SCALE", "KR", "VC"});
        norms[INTER_SCALE_INDEX] = NormalizationPicker.getFirstValidNormInThisOrder(ds, new String[]{"INTER_SCALE", "INTER_KR", "INTER_VC"});
        return norms;
    }

    @Override
    public void run() {
        ChromosomeHandler handler = ds.getChromosomeHandler();
        Chromosome[] chromosomes = handler.getAutosomalChromosomesArray();

        Map<Integer, Set<Integer>> badIndices = BadIndexFinder.getBadIndices(ds, chromosomes, resolution, norms[INTRA_SCALE_INDEX]);
        File tempOutputDirectory = new File(parentDirectory, "work");
        UNIXTools.makeDir(tempOutputDirectory);
        TranslocationSet translocations = SimpleTranslocationFinder.find(ds, norms, tempOutputDirectory,
                badIndices, resolution);
        Mappings mappings = IndexOrderer.getInitialMappings(ds, chromosomes, resolution,
                badIndices, norms[INTRA_SCALE_INDEX], generator.nextLong(), tempOutputDirectory,
                useExpandedIntraOE);

        MatrixAndWeight slice0 = MatrixBuilder.populateMatrix(ds, chromosomes, resolution,
                norms[INTER_SCALE_INDEX], norms[INTRA_SCALE_INDEX], mappings, translocations, true);

        Map<Integer, List<String>> bedFiles = new HashMap<>();
        runWithSettings(slice0, handler, chromosomes, tempOutputDirectory, bedFiles);

        Map<Integer, GenomeWide1DList<SubcompartmentInterval>> bestClusterings;
        if (skipCheck) {
            bestClusterings = InternalShuffle.getDefault(bedFiles, resolution,
                    handler);
        } else {
            bestClusterings = InternalShuffle.determineBest(bedFiles, resolution,
                    handler, ds, internalShuffleNorm);
        }

        for (int k : bestClusterings.keySet()) {
            File outBedFile = new File(parentDirectory, "SLICE_k" + k + "_best_clusters.bed"); // "_wcss" + wcss +
            bestClusterings.get(k).simpleExport(outBedFile);
        }
        System.out.println("\nSLICE complete");
    }

    private void runWithSettings(MatrixAndWeight slice0, ChromosomeHandler handler, Chromosome[] chromosomes,
                                 File tempOutputDirectory, Map<Integer, List<String>> bedFiles) {
        String stem = "SLICE";
        FinalMatrix slice = MatrixPreprocessor.preprocess(slice0, chromosomes, doPostNorm,
                appendIntra); // slice0.deepCopy()

        if (slice.notEmpty()) {
            if (SmartTools.printVerboseComments) {
                slice.export(tempOutputDirectory, stem);
            }
            ClusteringMagic clustering = new ClusteringMagic(slice, tempOutputDirectory,
                    handler, generator.nextLong(), skipCheck);
            clustering.extractFinalGWSubcompartments(stem, bedFiles);
        }
    }
}
