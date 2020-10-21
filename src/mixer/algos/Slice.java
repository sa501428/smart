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

import javastraw.reader.ChromosomeHandler;
import javastraw.reader.Dataset;
import javastraw.reader.HiCFileTools;
import javastraw.tools.UNIXTools;
import javastraw.type.NormalizationType;
import mixer.clt.CommandLineParserForMixer;
import mixer.clt.MixerCLT;
import mixer.utils.slice.FullGenomeOEWithinClusters;

import java.io.File;
import java.util.ArrayList;
import java.util.List;
import java.util.Random;

/**
 * experimental code
 * <p>
 * Created by muhammadsaadshamim on 9/14/15.
 */
public class Slice extends MixerCLT {
    
    private int resolution = 100000;
    private Dataset ds;
    private File outputDirectory;
    private final List<Dataset> datasetList = new ArrayList<>();
    private final List<String> inputHicFilePaths = new ArrayList<>();
    private final Random generator = new Random(22871L);
    private boolean useStackingAlongRow = false;
    private int minIntervalSizeAllowed = 1; // 1
    private String prefix = "";
    private String[] referenceBedFiles;

    // subcompartment lanscape identification via clustering enrichment
    public Slice(String command) {
        super("slice [-r resolution] [-k NONE/VC/VC_SQRT/KR] [-m num_clusters] " +
                "<input1.hic+input2.hic+input3.hic...> <output_file> <exp_name> [reference.bed]");
        useStackingAlongRow = command.contains("2");
    }
    
    @Override
    protected void readMixerArguments(String[] args, CommandLineParserForMixer mixerParser) {
        if (args.length != 4 && args.length != 5) {
            printUsageAndExit(5);
        }

        for (String path : args[1].split("\\+")) {
            System.out.println("Extracting " + path);
            inputHicFilePaths.add(path);
            datasetList.add(HiCFileTools.extractDatasetForCLT(path, true, false));
        }
        ds = datasetList.get(0);
        outputDirectory = HiCFileTools.createValidDirectory(args[2]);
        prefix = args[3];

        NormalizationType preferredNorm = mixerParser.getNormalizationTypeOption(ds.getNormalizationHandler());
        if (preferredNorm != null) norm = preferredNorm;

        List<String> possibleResolutions = mixerParser.getMultipleResolutionOptions();
        if (possibleResolutions != null) {
            if (possibleResolutions.size() > 1)
                System.err.println("Only one resolution can be specified\nUsing " + possibleResolutions.get(0));
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

        if (args.length == 5) {
            referenceBedFiles = args[4].split("\\+");
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
        
        if (useStackingAlongRow) {
            FullGenomeOEWithinClusters withinClusters = new FullGenomeOEWithinClusters(datasetList.get(0),
                    chromosomeHandler, resolution, norm, minIntervalSizeAllowed, outputDirectory, generator, referenceBedFiles);
            for (int i = 1; i < datasetList.size(); i++) {
                withinClusters.appendGWDataFromAdditionalDataset(datasetList.get(i));
            }
            withinClusters.extractFinalGWSubcompartments(generator, inputHicFilePaths, prefix, 0);
        } else {
            for (int i = 0; i < datasetList.size(); i++) {
                FullGenomeOEWithinClusters withinClusters = new FullGenomeOEWithinClusters(datasetList.get(i),
                        chromosomeHandler, resolution, norm, minIntervalSizeAllowed, outputDirectory, generator, referenceBedFiles);
                withinClusters.extractFinalGWSubcompartments(generator, inputHicFilePaths, prefix, i);
            }
            System.out.println("\nClustering complete");
        }
        
    }
}
