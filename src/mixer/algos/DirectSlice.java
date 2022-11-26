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
import javastraw.reader.type.NormalizationHandler;
import javastraw.tools.HiCFileTools;
import javastraw.tools.MatrixTools;
import mixer.clt.CommandLineParserForMixer;
import mixer.clt.MixerCLT;
import mixer.utils.BedTools;
import mixer.utils.drive.*;
import mixer.utils.tracks.SubcompartmentInterval;
import mixer.utils.translocations.TranslocationSet;

import java.io.File;

/**
 * experimental code
 * <p>
 * Created by muhammadsaadshamim on 9/14/15.
 */
public class DirectSlice extends MixerCLT {

    public static final int INTRA_SCALE_INDEX = 0;
    public static final int INTER_SCALE_INDEX = 1;
    private int resolution = 100000;
    private Dataset ds;
    private ChromosomeHandler handler;
    private File parentDirectory;
    private GenomeWide1DList<SubcompartmentInterval> clusters;

    // subcompartment landscape identification via compressing enrichments
    public DirectSlice(String command) {
        super("direct-slice [-r resolution] [--verbose] [-k INTER_NORM] " +
                "<file.hic> <input.bed> <outfolder>");
        //"<-k NONE/VC/VC_SQRT/KR/SCALE> [--compare reference.bed] [--has-translocation] " +
    }

    @Override
    protected void readMixerArguments(String[] args, CommandLineParserForMixer mixerParser) {
        if (args.length != 4) {
            printUsageAndExit(7);
        }

        resolution = updateResolution(mixerParser, resolution);
        ds = HiCFileTools.extractDatasetForCLT(args[1], true, false, resolution > 100);
        handler = ds.getChromosomeHandler();

        clusters = BedTools.loadBedFileAtResolution(handler, args[2], resolution);

        parentDirectory = HiCFileTools.createValidDirectory(args[3]);
    }

    @Override
    public void run() {
        Chromosome[] chromosomes = handler.getAutosomalChromosomesArray();

        Mappings mappings = new BedFileMappings(chromosomes, resolution, clusters);

        MatrixAndWeight slice = MatrixBuilder.populateMatrix(ds, chromosomes, resolution,
                NormalizationHandler.NONE, NormalizationHandler.NONE, mappings,
                new TranslocationSet(), false);

        slice.doSimpleVCNorm();

        FinalMatrix result = slice.getFinalMatrix(false);
        result.removeAllNanRows();
        result.removeAllNanCols();

        String path1 = new File(parentDirectory, "slice.matrix.npy").getAbsolutePath();
        MatrixTools.saveMatrixTextNumpy(path1, result.matrix);

        String path2 = new File(parentDirectory, "genome.indices.npy").getAbsolutePath();
        MatrixTools.saveMatrixTextNumpy(path2, result.getGenomeIndices());

        System.out.println("\nDirect SLICE Compression complete");
    }
}
