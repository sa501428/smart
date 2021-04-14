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

import javastraw.reader.Dataset;
import javastraw.reader.basics.ChromosomeHandler;
import javastraw.reader.type.NormalizationType;
import javastraw.tools.HiCFileTools;
import mixer.clt.CommandLineParserForMixer;
import mixer.clt.MixerCLT;
import mixer.utils.similaritymeasures.RobustManhattanDistance;
import mixer.utils.similaritymeasures.SimilarityMetric;
import mixer.utils.slice.cleaning.SliceMatrixCleaner;
import mixer.utils.slice.kmeans.FullGenomeOEWithinClusters;

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
    private NormalizationType[] intraNormArray;
    private NormalizationType[] interNormArray;
    public static SimilarityMetric metric = RobustManhattanDistance.SINGLETON; //  RobustMedianAbsoluteError
    private String prefix = "";
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

        intraNormArray = new NormalizationType[datasetList.size()];
        interNormArray = new NormalizationType[datasetList.size()];

        NormalizationType[][] preferredNorm = mixerParser.getNormalizationTypeOption(datasetList);
        if (preferredNorm != null && preferredNorm.length > 0 && preferredNorm[0].length > 0) {
            if (preferredNorm[0].length == intraNormArray.length) {
                intraNormArray = preferredNorm[0];
                interNormArray = preferredNorm[1];
            } else {
                Arrays.fill(intraNormArray, preferredNorm[0]);
                Arrays.fill(interNormArray, preferredNorm[0]);
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

        int subsampling = mixerParser.getSubsamplingOption();
        if (subsampling > 0) {
            SliceMatrixCleaner.NUM_PER_CENTROID = subsampling;
        }
    }

    @Override
    public void run() {

        ChromosomeHandler chromosomeHandler = ds.getChromosomeHandler();
        if (givenChromosomes != null)
            chromosomeHandler = HiCFileTools.stringToChromosomes(givenChromosomes, chromosomeHandler);

        if (datasetList.size() < 1) return;

        FullGenomeOEWithinClusters withinClusters = new FullGenomeOEWithinClusters(datasetList,
                chromosomeHandler, resolution, intraNormArray, interNormArray, outputDirectory, generator.nextLong(), metric);
        withinClusters.extractFinalGWSubcompartments(prefix);

        System.out.println("\nClustering complete");
    }
}
