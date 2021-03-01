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

package mixer.utils.shuffle;

import javastraw.featurelist.GenomeWideList;
import javastraw.reader.ChromosomeHandler;
import javastraw.reader.Dataset;
import javastraw.reader.basics.Chromosome;
import javastraw.tools.MatrixTools;
import javastraw.tools.UNIXTools;
import javastraw.type.NormalizationType;
import mixer.utils.common.FloatMatrixTools;
import mixer.utils.similaritymeasures.SimilarityMetric;
import mixer.utils.slice.cleaning.IntraMatrixCleaner;
import mixer.utils.slice.structures.SliceUtils;
import mixer.utils.slice.structures.SubcompartmentInterval;
import tagbio.umap.Umap;

import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.List;
import java.util.Random;

public class UMAPIntraMatrices {

    private final Random generator = new Random(0);
    private final Dataset ds;
    private final NormalizationType norm;
    private final int compressionFactor;
    private final int resolution;
    private final InterOnlyMatrix.INTRA_TYPE intra_type;
    private final SimilarityMetric metric;

    public UMAPIntraMatrices(Dataset ds, NormalizationType norm, int resolution, int compressionFactor,
                             InterOnlyMatrix.INTRA_TYPE intra_type, SimilarityMetric metric) {
        this.resolution = resolution;
        this.compressionFactor = compressionFactor;
        this.ds = ds;
        this.norm = norm;
        this.intra_type = intra_type;
        this.metric = metric;
    }

    public void runAnalysis(String[] bedFiles, File outputDirectory,
                            ChromosomeHandler chromosomeHandler) {

        Chromosome[] chromosomes = chromosomeHandler.getAutosomalChromosomesArray();
        List<GenomeWideList<SubcompartmentInterval>> allSubcompartments = new ArrayList<>();
        for (int i = 0; i < bedFiles.length; i++) {
            GenomeWideList<SubcompartmentInterval> subcompartments =
                    SliceUtils.loadFromSubcompartmentBEDFile(chromosomeHandler, bedFiles[i]);
            SliceUtils.collapseGWList(subcompartments);
            allSubcompartments.add(subcompartments);
        }

        try {
            writeToFile(outputDirectory, bedFiles);
        } catch (Exception e) {
            e.printStackTrace();
        }

        int numBedFiles = allSubcompartments.size();
        // todo, using only big size? todo sorting picture

        for (int y = 0; y < chromosomes.length; y++) {
            final InterOnlyMatrix interMatrix = new InterOnlyMatrix(ds, norm, resolution, chromosomes[y],
                    intra_type, metric);
            interMatrix.applySimpleLog();

            float[][] matrix = interMatrix.getMatrix();
            IntraMatrixCleaner.rollingAverage(matrix, compressionFactor);

            int[][] rowIndicesToClusterIDs = new int[matrix.length][numBedFiles];

            for (int z = 0; z < numBedFiles; z++) {
                populateIndexToClusterIDMap(interMatrix.getRowChromosomes(), interMatrix.getRowOffsets(),
                        allSubcompartments.get(z), rowIndicesToClusterIDs, z);
            }

            File outfolder = new File(outputDirectory, chromosomes[y].getName());
            UNIXTools.makeDir(outfolder);
            runUmapAndSaveMatrices(matrix, outfolder, "rows", rowIndicesToClusterIDs);
        }
    }

    private void runUmapAndSaveMatrices(float[][] data, File outputDirectory, String outstem,
                                        int[][] indexToIDs) {

        int numThreads = Math.max(1, Runtime.getRuntime().availableProcessors());
        System.out.println("Running UMAP with " + outstem);
        final Umap umap = new Umap();
        umap.setNumberComponents(2);
        umap.setNumberNearestNeighbours(50); // 15 -> 50 for more global picture
        umap.setThreads(numThreads);
        umap.setMinDist(0.2f);  //0.1f -> 0.2f for more general features
        umap.setVerbose(false);
        umap.setSeed(generator.nextLong());
        final float[][] result = umap.fitTransform(data);
        File temp = new File(outputDirectory, "umap_" + outstem + "_embedding.npy");
        FloatMatrixTools.saveMatrixTextNumpy(temp.getAbsolutePath(), result);
        System.gc();

        temp = new File(outputDirectory, outstem + "genome_indices.npy");
        MatrixTools.saveMatrixTextNumpy(temp.getAbsolutePath(), indexToIDs);
    }

    private void populateIndexToClusterIDMap(Chromosome[] chromosomes, int[] offsets,
                                             GenomeWideList<SubcompartmentInterval> subcompartments,
                                             int[][] indicesToClusterIDs, int colIndex) {
        for (int x = 0; x < chromosomes.length; x++) {
            Chromosome chrom = chromosomes[x];
            List<SubcompartmentInterval> intervalList = subcompartments.getFeatures("" + chrom.getIndex());
            for (SubcompartmentInterval interval : intervalList) {
                int xStart = interval.getX1() / resolution;
                int xEnd = interval.getX2() / resolution;
                int clusterID = interval.getClusterID();

                List<Integer> tempList = new ArrayList<>();
                for (int k = xStart; k < xEnd; k++) {
                    final int actualPosition = k + offsets[x];
                    tempList.add(actualPosition);
                    indicesToClusterIDs[actualPosition][colIndex] = clusterID;
                }
            }
        }
    }

    private void writeToFile(File outfolder, String[] filenames) throws IOException {
        FileWriter myWriter = new FileWriter(new File(outfolder, "bedfile_order.txt"));
        for (String filename : filenames) {
            myWriter.write(filename + "\n");
        }
        myWriter.close();
    }
}
