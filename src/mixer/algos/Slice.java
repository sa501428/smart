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
import mixer.utils.cleaning.MatrixPreprocessor;
import mixer.utils.drive.BinMappings;
import mixer.utils.drive.MatrixAndWeight;
import mixer.utils.drive.MatrixBuilder;
import mixer.utils.intra.IndexOrderer;
import mixer.utils.magic.ClusteringMagic;
import mixer.utils.nv.BadIndexFinder;
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
    private final int maxIters = 200;
    public static final int INTRA_SCALE_INDEX = 0;
    public static final int INTER_SCALE_INDEX = 1;
    public static final int GW_SCALE_INDEX = 2;
    public static final boolean PROJECT_TO_UMAP = true;
    public static final boolean USE_WEIGHTED_MEAN = false;
    public static boolean FILTER_OUTLIERS = true;
    private final Random generator = new Random(22871L);
    private int resolution = 100000;
    private Dataset ds;
    private File outputDirectory;
    private NormalizationType[] norms;
    private String prefix = "";
    public static boolean USE_KMEANS = false, USE_KMEDIANS = true;
    public static boolean USE_ENCODE_MODE = false;

    // subcompartment lanscape identification via clustering enrichment
    public Slice(String command) {
        super("slice [-r resolution] [--verbose] " +
                //"<-k NONE/VC/VC_SQRT/KR/SCALE> [--compare reference.bed] [--has-translocation] " +
                "<file.hic> <K0,KF,nK> <outfolder> <prefix_>\n" +
                "   K0 - minimum number of clusters\n" +
                "   KF - maximum number of clusters\n" +
                "   nK - number of times to rerun kmeans");
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
            ClusteringMagic.numAttemptsForKMeans = Integer.parseInt(valString[2]);
        } catch (Exception e) {
            printUsageAndExit(5);
        }

        outputDirectory = HiCFileTools.createValidDirectory(args[3]);
        prefix = args[4];
        norms = populateNormalizations(ds);

        updateGeneratorSeed(mixerParser, generator);
    }

    private NormalizationType[] populateNormalizations(Dataset ds) {
        NormalizationType[] norms = new NormalizationType[3];
        norms[INTRA_SCALE_INDEX] = NormalizationPicker.getFirstValidNormInThisOrder(ds, new String[]{"SCALE", "KR"});
        norms[INTER_SCALE_INDEX] = NormalizationPicker.getFirstValidNormInThisOrder(ds, new String[]{"INTER_SCALE", "INTER_KR"});
        norms[GW_SCALE_INDEX] = NormalizationPicker.getFirstValidNormInThisOrder(ds, new String[]{"GW_SCALE", "GW_KR"});
        return norms;
    }

    @Override
    public void run() {

        ChromosomeHandler handler = ds.getChromosomeHandler();
        Chromosome[] chromosomes = handler.getAutosomalChromosomesArray();

        Map<Integer, Set<Integer>> badIndices = BadIndexFinder.getBadIndices(ds, chromosomes, resolution);

        // todo should be at lower res
        SimpleTranslocationFinder translocations = new SimpleTranslocationFinder(ds, norms, outputDirectory, resolution, badIndices);

        // todo should be at lower res
        BinMappings mappings = IndexOrderer.getInitialMappings(ds, chromosomes, resolution,
                badIndices, norms[INTRA_SCALE_INDEX], generator.nextLong(), outputDirectory);

        MatrixAndWeight slice = MatrixBuilder.populateMatrix(ds, chromosomes, resolution,
                norms[INTER_SCALE_INDEX], mappings, translocations, outputDirectory);

        MatrixPreprocessor.clean(slice, mappings, chromosomes, false, true);

        slice.export(outputDirectory, "magic");

        ClusteringMagic clustering = new ClusteringMagic(slice, outputDirectory, handler, generator.nextLong());
        clustering.extractFinalGWSubcompartments(prefix);

        System.out.println("\nSLICE complete");
    }
}
