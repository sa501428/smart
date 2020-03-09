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

package mixer.tools.dev;

import mixer.MixerGlobals;
import mixer.data.ChromosomeHandler;
import mixer.data.Dataset;
import mixer.data.HiCFileTools;
import mixer.data.feature.GenomeWideList;
import mixer.tools.clt.CommandLineParserForMixer;
import mixer.tools.clt.MixerCLT;
import mixer.tools.utils.common.UNIXTools;
import mixer.tools.utils.dev.drink.*;
import mixer.windowui.NormalizationType;
import org.broad.igv.util.Pair;

import java.io.File;
import java.util.*;

/**
 * experimental code
 *
 * Created by muhammadsaadshamim on 9/14/15.
 */
public class Drink extends MixerCLT {

    public static final int USE_ONLY_DERIVATIVE = 1;
    public static final int IGNORE_DERIVATIVE = 2;

    private int resolution = 100000;
    private Dataset ds;
    private File outputDirectory;
    private final int numIntraIters = 3;
    private int numIntraClusters = 10;
    private final int whichApproachtoUse = 0;
    private int numInterClusters = 8;
    private final List<Dataset> datasetList = new ArrayList<>();
    private List<String> inputHicFilePaths = new ArrayList<>();
    private final boolean compareOnlyNotSubcompartment;
    private final float oeThreshold = 3f;
    private double[] convolution1d = null;
    private Random generator = new Random(22871L);
    private int derivativeStatus = 0;
    private boolean useNormalizationOfRows = false;

    public Drink(boolean compareOnlyNotSubcompartment) {
        super("drink [-r resolution] [-k NONE/VC/VC_SQRT/KR] [-m num_clusters] <input1.hic+input2.hic+input3.hic...> <output_file>");
        MixerGlobals.useCache = false;
        this.compareOnlyNotSubcompartment = compareOnlyNotSubcompartment;
    }

    @Override
    protected void readMixerArguments(String[] args, CommandLineParserForMixer mixerParser) {
        if (args.length != 3) {
            printUsageAndExit();
        }

        determineNumClusters(mixerParser);

        if (whichApproachtoUse == 0) {
            for (String path : args[1].split("\\+")) {
                System.out.println("Extracting " + path);
                inputHicFilePaths.add(path);
                List<String> tempList = new ArrayList<>();
                tempList.add(path);
                datasetList.add(HiCFileTools.extractDatasetForCLT(tempList, true));
            }
            ds = datasetList.get(0);
        } else {
            ds = HiCFileTools.extractDatasetForCLT(Arrays.asList(args[1].split("\\+")), true);
        }
        outputDirectory = HiCFileTools.createValidDirectory(args[2]);

        NormalizationType preferredNorm = mixerParser.getNormalizationTypeOption(ds.getNormalizationHandler());
        if (preferredNorm != null) norm = preferredNorm;

        List<String> possibleResolutions = mixerParser.getMultipleResolutionOptions();
        if (possibleResolutions != null) {
            if (possibleResolutions.size() > 1)
                System.err.println("Only one resolution can be specified for Drink\nUsing " + possibleResolutions.get(0));
            resolution = Integer.parseInt(possibleResolutions.get(0));
        }

        long[] possibleSeeds = mixerParser.getMultipleSeedsOption();
        if (possibleSeeds != null && possibleSeeds.length > 0) {
            for (long seed : possibleSeeds) {
                generator.setSeed(seed);
            }
        }

        convolution1d = mixerParser.getConvolutionOption();
        derivativeStatus = mixerParser.getUsingDerivativeStatus();
        useNormalizationOfRows = mixerParser.getUsingRowNomalizationStatus();
    }

    private void determineNumClusters(CommandLineParserForMixer mixerParser) {
        int n = mixerParser.getMatrixSizeOption();
        if (n > 1) {
            numInterClusters = n;
            numIntraClusters = n + 5;
        }
    }

    @Override
    public void run() {

        ChromosomeHandler chromosomeHandler = ds.getChromosomeHandler();
        if (givenChromosomes != null)
            chromosomeHandler = HiCFileTools.stringToChromosomes(givenChromosomes, chromosomeHandler);

        if (datasetList.size() < 1) return;

        InitialClusterer clusterer = new InitialClusterer(datasetList, chromosomeHandler, resolution, norm, numIntraClusters, generator, oeThreshold, convolution1d, numIntraIters);

        File initialClusteringOut = new File(outputDirectory, "initial_clustering");
        UNIXTools.makeDir(initialClusteringOut);
        Pair<List<GenomeWideList<SubcompartmentInterval>>, Map<Integer, float[]>> initialClustering = clusterer.extractAllComparativeIntraSubcompartmentsTo(initialClusteringOut);
        System.out.println("\nInitial clustering done");

        for (int i = 0; i < datasetList.size(); i++) {
            initialClustering.getFirst().get(i).simpleExport(new File(initialClusteringOut, DrinkUtils.cleanUpPath(inputHicFilePaths.get(i)) + "." + i + ".init.bed"));
        }

        if (compareOnlyNotSubcompartment) {
            ComparativeSubcompartmentsProcessor processor = new ComparativeSubcompartmentsProcessor(initialClustering,
                    chromosomeHandler, resolution);

            // process differences for diff vector
            processor.writeDiffVectorsRelativeToBaselineToFiles(outputDirectory, inputHicFilePaths,
                    "drink_r_" + resolution + "_k_" + numIntraClusters + "_diffs");

            processor.writeConsensusSubcompartmentsToFile(outputDirectory);

            processor.writeFinalSubcompartmentsToFiles(outputDirectory, inputHicFilePaths);
        } else {

            for (int i = 0; i < datasetList.size(); i++) {
                FullGenomeOEWithinClusters withinClusters = new FullGenomeOEWithinClusters(datasetList.get(i),
                        chromosomeHandler, resolution, norm, initialClustering.getFirst().get(i), oeThreshold, derivativeStatus, useNormalizationOfRows);

                Map<Integer, GenomeWideList<SubcompartmentInterval>> gwListMap = withinClusters.extractFinalGWSubcompartments(outputDirectory, generator);

                for (Integer key : gwListMap.keySet()) {
                    GenomeWideList<SubcompartmentInterval> gwList = gwListMap.get(key);
                    DrinkUtils.collapseGWList(gwList);
                    gwList.simpleExport(new File(outputDirectory, "gw_full_" + key + "_clusters_" + DrinkUtils.cleanUpPath(inputHicFilePaths.get(i)) + ".subcompartment.bed"));
                }
                System.out.println("\nClustering complete");
            }
        }
    }
}
