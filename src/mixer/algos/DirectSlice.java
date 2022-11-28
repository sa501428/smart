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
import javastraw.reader.basics.ChromosomeTools;
import javastraw.reader.type.NormalizationHandler;
import javastraw.tools.HiCFileTools;
import javastraw.tools.MatrixTools;
import mixer.clt.CommandLineParserForMixer;
import mixer.clt.MixerCLT;
import mixer.utils.BedTools;
import mixer.utils.cleaning.SimilarityMatrixTools;
import mixer.utils.common.FloatMatrixTools;
import mixer.utils.drive.*;
import mixer.utils.similaritymeasures.RobustManhattanDistance;
import mixer.utils.tracks.SubcompartmentInterval;
import mixer.utils.translocations.TranslocationSet;

import java.io.*;
import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.nio.file.Paths;

/**
 * experimental code
 * <p>
 * Created by muhammadsaadshamim on 9/14/15.
 */
public class DirectSlice extends MixerCLT {
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
        String genomeID = mixerParser.getGenomeOption();
        if (genomeID != null && genomeID.length() > 2) {
            handler = ChromosomeTools.loadChromosomes(genomeID);
        }

        clusters = BedTools.loadBedFileAtResolution(handler, args[2], resolution);

        parentDirectory = HiCFileTools.createValidDirectory(args[3]);
    }

    public static void saveMatrixText(String filename, float[][] matrix) {
        Writer writer = null;

        try {
            writer = new BufferedWriter(new OutputStreamWriter(Files.newOutputStream(Paths.get(filename)), StandardCharsets.UTF_8));
            for (float[] row : matrix) {
                for (int k = 0; k < row.length; k++) {
                    writer.write("" + row[k]);
                    if (k == row.length - 1) {
                        writer.write("\n");
                    } else {
                        writer.write(" ");
                    }
                }
            }
        } catch (IOException e) {
            e.printStackTrace();
        } finally {
            try {
                if (writer != null) {
                    writer.close();
                }
            } catch (Exception e) {
                e.printStackTrace();
            }
        }
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
        FloatMatrixTools.log(result.matrix, 1);

        String path = new File(parentDirectory, "slice.matrix.npy").getAbsolutePath();
        MatrixTools.saveMatrixTextNumpy(path, result.matrix);

        System.out.print("matrix size " + result.matrix.length + " x ");
        System.out.println(result.matrix[0].length);

        path = new File(parentDirectory, "genome.indices.npy").getAbsolutePath();
        MatrixTools.saveMatrixTextNumpy(path, result.getGenomeIndices());

        System.out.print("indices size " + result.getGenomeIndices().length + " x ");
        System.out.println(result.getGenomeIndices()[0].length);

        float[][] dist = SimilarityMatrixTools.getSymmetricDistanceMatrix(result.matrix,
                RobustManhattanDistance.SINGLETON);
        export("slice.l1.distance.matrix", dist);

        /*
        dist = SimilarityMatrixTools.getSymmetricDistanceMatrix(result.matrix,
                RobustEuclideanDistance.SINGLETON);
        export("slice.l2.distance.matrix", dist);

        dist = SimilarityMatrixTools.getSymmetricDistanceMatrix(result.matrix,
                RobustCosineSimilarity.SINGLETON);
        export("slice.cosine.distance.matrix", dist);

        dist = SimilarityMatrixTools.getSymmetricDistanceMatrix(result.matrix,
                RobustCorrelationSimilarity.SINGLETON);
        export("slice.corr.distance.matrix", dist);
        */

        System.out.println("\nDirect SLICE Compression complete");
    }

    private void export(String name, float[][] dist) {
        long matrixSize = dist.length;
        matrixSize *= dist.length;
        if (matrixSize < Integer.MAX_VALUE && matrixSize > 0) {
            String path = new File(parentDirectory, name + ".npy").getAbsolutePath();
            MatrixTools.saveMatrixTextNumpy(path, dist);
        } else {
            String path = new File(parentDirectory, name + ".txt").getAbsolutePath();
            saveMatrixText(path, dist);
        }
    }
}
