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
import javastraw.reader.basics.ChromosomeHandler;
import javastraw.reader.basics.ChromosomeTools;
import javastraw.tools.MatrixTools;
import mixer.clt.CommandLineParserForMixer;
import mixer.clt.MixerCLT;
import mixer.utils.BedTools;
import mixer.utils.tracks.ARITools;
import mixer.utils.tracks.Concensus2DTools;
import mixer.utils.tracks.SubcompartmentInterval;

import java.util.ArrayList;
import java.util.List;

public class Compare extends MixerCLT {

    private List<GenomeWide1DList<SubcompartmentInterval>> files;
    private int resolution = 1000;
    private final boolean perChromosome, doARI;
    private String outputName;

    public Compare(String name) {
        super("compare[-per-chrom] [-r resolution] <genomeID> <output.npy> <file1.bed> <file2.bed> ... <fileN.bed>\n" +
                "ari [-r resolution] <genomeID> <output.npy> <file1.bed> <file2.bed> ... <fileN.bed>");
        perChromosome = name.contains("per") && name.contains("chrom");
        doARI = name.contains("ari");
    }

    @Override
    protected void readMixerArguments(String[] args, CommandLineParserForMixer mixerParser) {
        if (args.length < 5) {
            printUsageAndExit(51);
        }

        List<Integer> possibleResolutions = mixerParser.getMultipleResolutionOptions();
        if (possibleResolutions != null) {
            if (possibleResolutions.size() > 1)
                System.err.println("Only one resolution can be specified\nUsing " + possibleResolutions.get(0));
            resolution = possibleResolutions.get(0);
        }

        ChromosomeHandler handler = ChromosomeTools.loadChromosomes(args[1]);
        outputName = args[2];
        files = new ArrayList<>();
        for (int k = 3; k < args.length; k++) {
            files.add(BedTools.loadBedFileAtResolution(handler, args[k], resolution));
        }
    }

    @Override
    public void run() {
        if (files.size() == 2) {
            if (doARI) {
                System.out.println("ARI = " + ARITools.getARI(files.get(0), files.get(1)));
            } else {
                if (perChromosome) {
                    Concensus2DTools.checkOverlapPerChrom(files.get(0), files.get(1));
                }
                double accuracy = Concensus2DTools.checkOverlap(files.get(0), files.get(1));
                System.out.println("Accuracy " + (100 * accuracy));
            }
        } else {
            double[][] result = new double[files.size()][files.size()];
            for (int i = 0; i < result.length; i++) {
                for (int j = i; j < result.length; j++) {
                    if (doARI) {
                        result[i][j] = ARITools.getARI(files.get(i), files.get(j));
                    } else {
                        result[i][j] = Concensus2DTools.checkOverlap(files.get(i), files.get(j));
                    }
                    result[j][i] = result[i][j];
                }
            }
            MatrixTools.saveMatrixTextNumpy(outputName, result);
        }
    }
}
