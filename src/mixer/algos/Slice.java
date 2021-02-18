/*
 * The MIT License (MIT)
 *
 * Copyright (c) 2011-2021 Rice University, Baylor College of Medicine, Aiden Lab
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

import javastraw.reader.ChromosomeHandler;
import javastraw.reader.Dataset;
import javastraw.reader.HiCFileTools;
import javastraw.type.NormalizationType;
import mixer.clt.CommandLineParserForMixer;
import mixer.clt.MixerCLT;
import mixer.utils.similaritymeasures.RobustMedianAbsoluteError;
import mixer.utils.similaritymeasures.SimilarityMetric;
import mixer.utils.slice.cleaning.MatrixCleanerAndProjector;
import mixer.utils.slice.kmeans.FullGenomeOEWithinClusters;
import mixer.utils.slice.matrices.SliceMatrix;
import mixer.utils.slice.structures.HiCInterTools;

import java.io.File;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;
import java.util.Random;

/**
 * experimental code
 * <p>
 * Created by muhammadsaadshamim on 9/14/15.
 */
public class Slice extends MixerCLT {

    private final List<Dataset> datasetList = new ArrayList<>();
    private final List<String> inputHicFilePaths = new ArrayList<>();
    private final Random generator = new Random(22871L);
    private int resolution = 100000;
    private Dataset ds;
    private File outputDirectory;
    private NormalizationType[] normArray;
    public static SimilarityMetric metric = RobustMedianAbsoluteError.SINGLETON; //RobustManhattanDistance
    private String prefix = "";
    private String[] referenceBedFiles;
    private final boolean compareMaps;

    // subcompartment lanscape identification via clustering enrichment
    public Slice(String command) {
        super("slice [-r resolution] <-k NONE/VC/VC_SQRT/KR/SCALE>  [-w window] " +
                "[--compare reference.bed] [--verbose] " + //"[--has-translocation] " +
                "<input1.hic+input2.hic...> <K0,KF,nK> <outfolder> <prefix_>");
        compareMaps = command.contains("dice");
    }

    @Override
    protected void readMixerArguments(String[] args, CommandLineParserForMixer mixerParser) {
        if (args.length != 5) {
            printUsageAndExit(5);
        }

        for (String path : args[1].split(",")) {
            System.out.println("Extracting " + path);
            inputHicFilePaths.add(path);
            datasetList.add(HiCFileTools.extractDatasetForCLT(path, true, false));
        }

        try {
            String[] valString = args[2].split(",");
            FullGenomeOEWithinClusters.startingClusterSizeK = Integer.parseInt(valString[0]);
            FullGenomeOEWithinClusters.numClusterSizeKValsUsed = Integer.parseInt(valString[1]) -
                    FullGenomeOEWithinClusters.startingClusterSizeK;
            FullGenomeOEWithinClusters.numAttemptsForKMeans = Integer.parseInt(valString[2]);
        } catch (Exception e) {
            printUsageAndExit(5);
        }

        ds = datasetList.get(0);
        outputDirectory = HiCFileTools.createValidDirectory(args[3]);
        prefix = args[4];

        normArray = new NormalizationType[datasetList.size()];
        NormalizationType[] preferredNorm = mixerParser.getNormalizationTypeOption(datasetList);
        if (preferredNorm != null && preferredNorm.length > 0) {
            if (preferredNorm.length == normArray.length) {
                normArray = preferredNorm;
            } else {
                Arrays.fill(normArray, preferredNorm[0]);
            }
        } else {
            System.err.println("Normalization must be specified");
            System.exit(9);
        }


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
        if (minSize > 0) {
            SliceMatrix.numColumnsToPutTogether = minSize;
        } else {
            int numColumns = HiCInterTools.calculateIdealWidth(ds, resolution);
            for (int i = 1; i < datasetList.size(); i++) {
                int numColumns2 = HiCInterTools.calculateIdealWidth(datasetList.get(i), resolution);
                if (numColumns2 > numColumns) {
                    numColumns = numColumns2;
                }
            }
            SliceMatrix.numColumnsToPutTogether = numColumns;
            System.out.println("Using compression width: " + numColumns);
        }

        int subsampling = mixerParser.getSubsamplingOption();
        if (subsampling > 0) {
            MatrixCleanerAndProjector.NUM_PER_CENTROID = subsampling;
        }

        String bedFiles = mixerParser.getCompareReferenceOption();
        if (bedFiles != null && bedFiles.length() > 1) {
            referenceBedFiles = bedFiles.split(",");
        }
    }

    @Override
    public void run() {

        ChromosomeHandler chromosomeHandler = ds.getChromosomeHandler();
        if (givenChromosomes != null)
            chromosomeHandler = HiCFileTools.stringToChromosomes(givenChromosomes, chromosomeHandler);

        if (datasetList.size() < 1) return;

        FullGenomeOEWithinClusters withinClusters = new FullGenomeOEWithinClusters(datasetList,
                chromosomeHandler, resolution, normArray, outputDirectory, generator.nextLong(), referenceBedFiles, metric);
        withinClusters.extractFinalGWSubcompartments(inputHicFilePaths, prefix, 0, compareMaps);

        System.out.println("\nClustering complete");
    }
}
