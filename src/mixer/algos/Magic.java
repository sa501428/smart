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
import javastraw.reader.type.NormalizationHandler;
import javastraw.reader.type.NormalizationType;
import javastraw.tools.HiCFileTools;
import mixer.clt.CommandLineParserForMixer;
import mixer.clt.MixerCLT;
import mixer.utils.SubcompartmentBedFile;
import mixer.utils.magic.ClusteringMagic;
import mixer.utils.magic.MagicMatrix;
import mixer.utils.slice.structures.SubcompartmentInterval;

import java.io.File;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.Random;

/**
 * experimental code
 * <p>
 * Created by muhammadsaadshamim on 9/14/15.
 */
public class Magic extends MixerCLT {

    public static final boolean PROJECT_TO_UMAP = true;
    private final Random generator = new Random(22871L);
    private Dataset ds;
    private int resolution = 1000;
    private File outputDirectory;
    private String prefix, bedpath;
    private NormalizationType norm = NormalizationHandler.NONE;
    private int numCols = -1;
    private int numRows = -1;

    // subcompartment lanscape identification via clustering enrichment
    public Magic(String command) {
        super("magic [-r resolution] [--verbose] [-k norm] " +
                //"<-k NONE/VC/VC_SQRT/KR/SCALE> [--compare reference.bed] [--has-translocation] " +
                "<file.hic> <K0,KF,nK> <initial_clusters.bed> <outfolder> <prefix_>\n" +
                "   K0 - minimum number of clusters\n" +
                "   KF - maximum number of clusters\n" +
                "   nK - number of times to rerun k-medians");
    }

    @Override
    protected void readMixerArguments(String[] args, CommandLineParserForMixer mixerParser) {
        if (args.length != 6) {
            printUsageAndExit(5);
        }

        ds = HiCFileTools.extractDatasetForCLT(args[1], true, false);
        try {
            String[] valString = args[2].split(",");
            ClusteringMagic.startingClusterSizeK = Integer.parseInt(valString[0]);
            ClusteringMagic.numClusterSizeKValsUsed = Integer.parseInt(valString[1]) -
                    ClusteringMagic.startingClusterSizeK;
            ClusteringMagic.numAttemptsForKMeans = Integer.parseInt(valString[2]);
        } catch (Exception e) {
            printUsageAndExit(5);
        }

        bedpath = args[3];
        outputDirectory = HiCFileTools.createValidDirectory(args[4]);
        prefix = args[5];

        NormalizationType potentialNorm = mixerParser.getNormalizationTypeOption(ds.getNormalizationHandler());
        if (potentialNorm != null) {
            norm = potentialNorm;
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
    }

    @Override
    public void run() {

        ChromosomeHandler chromosomeHandler = ds.getChromosomeHandler();
        if (givenChromosomes != null)
            chromosomeHandler = HiCFileTools.stringToChromosomes(givenChromosomes, chromosomeHandler);

        int[] offset = generateFromChromosomes(chromosomeHandler);

        Map<Integer, Integer> binToColumn = populateFromBedFile(bedpath, chromosomeHandler, offset, resolution);

        MagicMatrix matrix = new MagicMatrix(ds, chromosomeHandler, resolution, norm,
                outputDirectory, generator.nextLong(), offset, binToColumn, numRows, numCols);
        matrix.export(new File(outputDirectory, "magic.npy").getAbsolutePath());

        ClusteringMagic clustering = new ClusteringMagic(matrix, outputDirectory, chromosomeHandler, 10L);
        clustering.extractFinalGWSubcompartments(prefix);
        System.out.println("\nClustering complete");
    }

    private int[] generateFromChromosomes(ChromosomeHandler handler) {
        Chromosome[] chromosomes = handler.getAutosomalChromosomesArray();
        int[] offsets = new int[chromosomes.length];
        int offset = 0;
        for (int z = 0; z < chromosomes.length; z++) {
            Chromosome chromosome = chromosomes[z];
            offsets[z] = offset;
            offset += (int) (chromosome.getLength() / resolution) + 1;
        }
        numRows = offset;
        return offsets;
    }

    private Map<Integer, Integer> populateFromBedFile(String bedpath, ChromosomeHandler handler,
                                                      int[] offsets, int resolution) {
        SubcompartmentBedFile bedFile = new SubcompartmentBedFile(bedpath, handler);
        if (bedFile.getNumTotalFeatures() < 1) {
            System.err.println("bed file is empty or incorrect path provided.");
            System.exit(3);
        }

        Map<Integer, Integer> binToColumnNumber = new HashMap<>();
        Map<Integer, Integer> clusterNumToColIndex = new HashMap<>();

        Chromosome[] chromosomes = handler.getAutosomalChromosomesArray();
        for (int z = 0; z < chromosomes.length; z++) {
            Chromosome chromosome = chromosomes[z];
            List<SubcompartmentInterval> intervals = bedFile.get(chromosome.getIndex());
            int offset = offsets[z];
            for (SubcompartmentInterval si : intervals) {
                for (int x = si.getX1() / resolution; x < si.getX2() / resolution; x++) {
                    int currIndex = offset + x;
                    int cluster = si.getClusterID();
                    if (!clusterNumToColIndex.containsKey(cluster)) {
                        clusterNumToColIndex.put(cluster, numCols);
                        numCols++;
                    }
                    int col = clusterNumToColIndex.get(cluster);
                    binToColumnNumber.put(currIndex, col);
                }
            }
        }
        return binToColumnNumber;
    }
}
