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
import javastraw.reader.type.NormalizationType;
import javastraw.tools.HiCFileTools;
import javastraw.tools.UNIXTools;
import mixer.clt.CommandLineParserForMixer;
import mixer.clt.MixerCLT;
import mixer.utils.BedTools;
import mixer.utils.shuffle.Partition;
import mixer.utils.shuffle.ShuffleAction;
import mixer.utils.tracks.SubcompartmentInterval;

import java.io.File;
import java.util.*;

/**
 * experimental code
 * <p>
 * Created by muhammadsaadshamim on 9/14/15.
 */
public class ChicScore extends MixerCLT {

    //private final Random generator = new Random(22871L);
    private final Random generator = new Random(22871L);
    private final boolean useOriginal;
    private Dataset ds;
    private int resolution = 100000;
    private int compressionFactor = 8;
    private File outputDirectory;
    private String[] prefix;
    private String[] referenceBedFiles;
    private NormalizationType norm;

    // subcompartment lanscape identification via clustering enrichment
    public ChicScore(String name) {
        super("shuffle [-r resolution] [-k NONE/VC/VC_SQRT/KR/SCALE] [-w window] [--verbose] " +
                "<file.hic> <subcompartment.bed(s)> <outfolder> <prefix>");
        useOriginal = name.contains("original");
    }

    @Override
    protected void readMixerArguments(String[] args, CommandLineParserForMixer mixerParser) {
        if (args.length != 5) {
            printUsageAndExit(51);
        }

        ds = HiCFileTools.extractDatasetForCLT(args[1], false, false, true);

        referenceBedFiles = args[2].split(",");
        outputDirectory = HiCFileTools.createValidDirectory(args[3]);
        prefix = args[4].split(",");

        NormalizationType preferredNorm = mixerParser.getNormalizationTypeOption(ds.getNormalizationHandler());
        if (preferredNorm != null) norm = preferredNorm;

        List<Integer> possibleResolutions = mixerParser.getMultipleResolutionOptions();
        if (possibleResolutions != null) {
            if (possibleResolutions.size() > 1)
                System.err.println("Only one resolution can be specified\nUsing " + possibleResolutions.get(0));
            resolution = possibleResolutions.get(0);
        }

        updateGeneratorSeed(mixerParser, generator);

        int minSize = mixerParser.getWindowSizeOption(1600000 / resolution);
        if (minSize > 1) {
            compressionFactor = minSize;
        }
    }

    public static void ensureSameLoci(List<GenomeWide1DList<SubcompartmentInterval>> subcompartments,
                                      ChromosomeHandler handler, int resolution) {

        int numFiles = subcompartments.size();

        Map<String, int[]> chromToLociToCounts = new HashMap<>();
        for (Chromosome chromosome : handler.getAutosomalChromosomesArray()) {
            int n = (int) ((chromosome.getLength() / resolution) + 1);
            chromToLociToCounts.put("" + chromosome.getIndex(), new int[n]);
        }

        for (GenomeWide1DList<SubcompartmentInterval> list : subcompartments) {
            list.processLists((s, list1) -> {
                int[] counts = chromToLociToCounts.get(s);
                for (SubcompartmentInterval interval : list1) {
                    counts[interval.getX1() / resolution]++;
                }
            });
        }

        for (GenomeWide1DList<SubcompartmentInterval> list : subcompartments) {
            list.filterLists((s, list1) -> {
                int[] counts = chromToLociToCounts.get(s);
                List<SubcompartmentInterval> toKeep = new ArrayList<>(list1.size());
                for (SubcompartmentInterval interval : list1) {
                    if (counts[interval.getX1() / resolution] == numFiles) {
                        toKeep.add(interval);
                    }
                }
                list1.clear();
                return toKeep;
            });
        }

    }

    @Override
    public void run() {

        ChromosomeHandler handler = ds.getChromosomeHandler();

        UNIXTools.makeDir(outputDirectory);

        Partition.Type[] mapTypes = {Partition.Type.SPLIT1, Partition.Type.SPLIT2,
                Partition.Type.SPLIT3, Partition.Type.SPLIT4, Partition.Type.SPLIT5};
        if (useOriginal) {
            mapTypes = new Partition.Type[]{Partition.Type.ODDS_VS_EVENS,
                    Partition.Type.SKIP_BY_TWOS, Partition.Type.FIRST_HALF_VS_SECOND_HALF};
        }

        List<GenomeWide1DList<SubcompartmentInterval>> subcompartments = new ArrayList<>(referenceBedFiles.length);
        for (String referenceBedFile : referenceBedFiles) {
            subcompartments.add(BedTools.loadBedFileAtResolution(handler, referenceBedFile, resolution));
        }

        ensureSameLoci(subcompartments, handler, resolution);

        /*
        for (int i = 0; i < referenceBedFiles.length; i++) {
            System.out.println("Filtering for common loci " + prefix[i]);
            File outBedFile = new File(referenceBedFiles[i]+ "_filtered_for_common_loci.bed"); // "_wcss" + wcss +
            subcompartments.get(i).simpleExport(outBedFile);
        }
        */

        for (int i = 0; i < referenceBedFiles.length; i++) {
            System.out.println("Processing " + prefix[i]);
            File newFolder = new File(outputDirectory, "shuffle_" + prefix[i]);
            UNIXTools.makeDir(newFolder);
            ShuffleAction matrix = new ShuffleAction(ds, norm, resolution, compressionFactor, mapTypes);
            matrix.runInterAnalysis(subcompartments.get(i), newFolder, generator);
            matrix.savePlotsAndResults(newFolder, prefix[i]);
        }
        System.out.println("Shuffle complete");
    }
}
