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

import mixer.MixerGlobals;
import mixer.commandline.handling.CommandLineParserForMixer;
import mixer.commandline.handling.MixerCLT;
import mixer.commandline.utils.common.UNIXTools;
import mixer.commandline.utils.drink.*;
import mixer.data.ChromosomeHandler;
import mixer.data.Dataset;
import mixer.data.HiCFileTools;
import mixer.data.feature.GenomeWideList;
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
    private final float oeThreshold = 3f;
    private final int whichApproachtoUse = 0;
    private final List<Dataset> datasetList = new ArrayList<>();
    private List<String> inputHicFilePaths = new ArrayList<>();
    private final boolean compareOnlyNotSubcompartment;
    private int[] numIntraClusters = new int[]{10, 10, 10};
    private double[] convolution1d = null;
    private Random generator = new Random(22871L);
    private int derivativeStatus = 0;
    private boolean useNormalizationOfRows = false;
    private boolean useStackingAlongRow = false;
    private int minIntervalSizeAllowed = 3; // 1
    private String prefix = "";

    public Drink(String command) {
        super("drink [-r resolution] [-k NONE/VC/VC_SQRT/KR] [-m num_clusters] <input1.hic+input2.hic+input3.hic...> <output_file>");
        MixerGlobals.useCache = false;
        this.compareOnlyNotSubcompartment = command.equalsIgnoreCase("drink");
        useStackingAlongRow = command.contains("2");
    }

    @Override
    protected void readMixerArguments(String[] args, CommandLineParserForMixer mixerParser) {
        if (args.length != 4) {
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
        prefix = args[3];

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

        int minSize = mixerParser.getAPAWindowSizeOption();
        if (minSize > 0) {
            minIntervalSizeAllowed = minSize;
        }

        convolution1d = mixerParser.getConvolutionOption();
        derivativeStatus = mixerParser.getUsingDerivativeStatus();
        useNormalizationOfRows = mixerParser.getUsingRowNomalizationStatus();
    }

    private void determineNumClusters(CommandLineParserForMixer mixerParser) {
        int n = mixerParser.getMatrixSizeOption();
        if (n > 1) {
            numIntraClusters = new int[]{n, n, n};
        }
    }

    @Override
    public void run() {

        ChromosomeHandler chromosomeHandler = ds.getChromosomeHandler();
        if (givenChromosomes != null)
            chromosomeHandler = HiCFileTools.stringToChromosomes(givenChromosomes, chromosomeHandler);

        if (datasetList.size() < 1) return;

        File initialClusteringOutDir = new File(outputDirectory, "initial_clustering");
        UNIXTools.makeDir(initialClusteringOutDir);


        if (compareOnlyNotSubcompartment) {
            InitialClusterer clusterer = new InitialClusterer(datasetList, chromosomeHandler, resolution, norm,
                    numIntraClusters, generator, oeThreshold, convolution1d, numIntraIters, useStackingAlongRow);

            Pair<List<GenomeWideList<SubcompartmentInterval>>, Map<Integer, float[]>> initialClustering =
                    clusterer.extractIntraSubcompartmentsTo(initialClusteringOutDir);
            for (int i = 0; i < datasetList.size(); i++) {
                initialClustering.getFirst().get(i).simpleExport(new File(initialClusteringOutDir, DrinkUtils.cleanUpPath(inputHicFilePaths.get(i)) + "." + i + ".init.bed"));
            }

            ComparativeSubcompartmentsProcessor processor = new ComparativeSubcompartmentsProcessor(initialClustering,
                    chromosomeHandler, resolution);

            // process differences for diff vector
            processor.writeDiffVectorsRelativeToBaselineToFiles(outputDirectory, inputHicFilePaths,
                    "drink_r_" + resolution + "_k_" + numIntraClusters + "_diffs");

            processor.writeConsensusSubcompartmentsToFile(outputDirectory);

            processor.writeFinalSubcompartmentsToFiles(outputDirectory, inputHicFilePaths);
        } else {

            /*
            DistanceSplitter splitter = new DistanceSplitter(datasetList, chromosomeHandler, resolution, norm, oeThreshold, useStackingAlongRow);
            List<GenomeWideList<SubcompartmentInterval>> initialSplit = splitter.extractIntraSubcompartmentsTo(initialClusteringOutDir);
             */

            InitialClusterer clusterer = new InitialClusterer(datasetList, chromosomeHandler, resolution, norm,
                    numIntraClusters, generator, oeThreshold, convolution1d, numIntraIters, useStackingAlongRow);

            List<GenomeWideList<SubcompartmentInterval>> initialSplit =
                    clusterer.extractIntraSubcompartmentsTo(initialClusteringOutDir).getFirst();
            for (int i = 0; i < datasetList.size(); i++) {
                GenomeWideList<SubcompartmentInterval> temp = initialSplit.get(i);
                DrinkUtils.collapseGWList(temp);
                temp.simpleExport(new File(initialClusteringOutDir, DrinkUtils.cleanUpPath(inputHicFilePaths.get(i)) + "." + i + ".init.bed"));
            }

            if (useStackingAlongRow) {

                initialSplit.get(0).simpleExport(new File(initialClusteringOutDir, DrinkUtils.cleanUpPath(inputHicFilePaths.get(0)) + ".init.split.bed"));
                /*
                FullGenomeOEWithinClusters withinClusters = new StackedFullGenomeOEWithinClusters(datasetList,
                        chromosomeHandler, resolution, norm, initialClustering.getFirst().get(0), oeThreshold, derivativeStatus, useNormalizationOfRows);

                Map<Integer, GenomeWideList<SubcompartmentInterval>> gwListMap = withinClusters.extractFinalGWSubcompartments(outputDirectory, generator);

                for (Integer key : gwListMap.keySet()) {
                    GenomeWideList<SubcompartmentInterval> gwList = gwListMap.get(key);
                    DrinkUtils.collapseGWList(gwList);
                    gwList.simpleExport(new File(outputDirectory, "gw_full_" + key + "_clusters_" + DrinkUtils.cleanUpPath(inputHicFilePaths.get(i)) + ".subcompartment.bed"));
                }

                 */
            } else {

                for (int i = 0; i < datasetList.size(); i++) {
                    initialSplit.get(i).simpleExport(new File(initialClusteringOutDir, DrinkUtils.cleanUpPath(inputHicFilePaths.get(i)) + "." + i + ".init.split.bed"));
                }

                for (int i = 0; i < datasetList.size(); i++) {
                    FullGenomeOEWithinClusters withinClusters = new FullGenomeOEWithinClusters(datasetList.get(i),
                            chromosomeHandler, resolution, norm, initialSplit.get(i), oeThreshold, minIntervalSizeAllowed);

                    withinClusters.extractFinalGWSubcompartments(outputDirectory, generator, inputHicFilePaths, prefix, i);
                }
                System.out.println("\nClustering complete");
            }
        }
    }
}
