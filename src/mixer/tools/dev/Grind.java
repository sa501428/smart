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

package mixer.tools.dev;

import mixer.data.HiCFileTools;
import mixer.tools.clt.CommandLineParserForMixer;
import mixer.tools.clt.MixerCLT;
import mixer.tools.utils.mixer.grind.*;
import mixer.track.feature.Feature2DParser;
import mixer.windowui.NormalizationType;

import java.util.Arrays;
import java.util.List;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;

/**
 * Generating Regions of Interest for Network Discovery
 */

public class Grind extends MixerCLT {

    public static final int LIST_ITERATION_OPTION = 1;
    public static final int DOMAIN_OPTION = 2;
    public static final int DOWN_DIAGONAL_OPTION = 3;
    public static final int DISTORTION_OPTION = 4;
    private ParameterConfigurationContainer container = new ParameterConfigurationContainer();

    public Grind() {
        super("grind [-k NONE/KR/VC/VC_SQRT] [-r resolution] [--stride increment] " +
                "[--off-from-diagonal max-dist-from-diag] " +
                "--observed-over-expected --dense-labels --ignore-feature-orientation --only-make-positives " +
                "<mode> <hic file> <bedpe positions> <x,y,z> <directory>" +
                "     \n" +
                "     mode: --iterate-down-diagonal --iterate-on-list --iterate-distortions --iterate-domains");
    }

    @Override
    protected void readMixerArguments(String[] args, CommandLineParserForMixer mixerParser) {
        if (args.length != 5) {
            printUsageAndExit();
        }

        container.ds = HiCFileTools.
                extractDatasetForCLT(Arrays.asList(args[1].split("\\+")), true);

        container.featureListPath = args[2];

        // split on commas
        // save the dimensions
        String[] dimensions = args[3].split(",");
        container.x = Integer.parseInt(dimensions[0]);
        container.y = Integer.parseInt(dimensions[1]);
        container.z = Integer.parseInt(dimensions[2]);

        container.useObservedOverExpected = mixerParser.getUseObservedOverExpectedOption();
        container.featureDirectionOrientationIsImportant = mixerParser.getDontIgnoreDirectionOrientationOption();
        container.useAmorphicPixelLabeling = mixerParser.getUseAmorphicLabelingOption();
        container.onlyMakePositiveExamples = mixerParser.getUseOnlyMakePositiveExamplesOption();
        container.useDenseLabelsNotBinary = mixerParser.getDenseLabelsOption();
        container.wholeGenome = mixerParser.getUseWholeGenome();
        container.offsetOfCornerFromDiagonal = mixerParser.getCornerOffBy();
        container.stride = mixerParser.getStride();
        container.outputDirectory = HiCFileTools.createValidDirectory(args[4]);
        container.useDiagonal = mixerParser.getUseGenomeDiagonal();
        container.useTxtInsteadOfNPY = mixerParser.getUseTxtInsteadOfNPY();

        NormalizationType preferredNorm = mixerParser.getNormalizationTypeOption(container.ds.getNormalizationHandler());
        if (preferredNorm != null) norm = preferredNorm;
        container.norm = norm;

        List<String> possibleResolutions = mixerParser.getMultipleResolutionOptions();
        if (possibleResolutions != null) {
            for (String num : possibleResolutions) {
                container.resolutions.add(Integer.parseInt(num));
            }
        } else {
            container.resolutions.add(10000);
        }

        container.grindIterationTypeOption = mixerParser.getGrindDataSliceOption();
        container.imgFileType = mixerParser.getGenerateImageFormatPicturesOption();
    }

    @Override
    public void run() {

        container.chromosomeHandler = container.ds.getChromosomeHandler();
        container.feature2DList = null;
        try {
            container.feature2DList = Feature2DParser.loadFeatures(container.featureListPath, container.chromosomeHandler, false, null, false);
        } catch (Exception e) {
            if (container.grindIterationTypeOption != 4) {
                System.err.println("Feature list failed to load");
                e.printStackTrace();
                System.exit(-9);
            }
        }

        if (givenChromosomes != null)
            container.chromosomeHandler = HiCFileTools.stringToChromosomes(givenChromosomes, container.chromosomeHandler);

        RegionFinder finder = null;
        if (container.grindIterationTypeOption == LIST_ITERATION_OPTION) {
            finder = new IterateOnFeatureListFinder(container);
        } else if (container.grindIterationTypeOption == DOMAIN_OPTION) {
            finder = new DomainFinder(container);
        } else if (container.grindIterationTypeOption == DOWN_DIAGONAL_OPTION) {
            finder = new IterateDownDiagonalFinder(container);
        } else {
            runDistortionTypeOfIteration();
        }

        if (finder != null) {
            finder.makeExamples();
        }
    }

    private void runDistortionTypeOfIteration() {
        ExecutorService executor = Executors.newFixedThreadPool(container.resolutions.size());
        for (final int resolution : container.resolutions) {
            Runnable worker = new Runnable() {
                @Override
                public void run() {
                    RegionFinder finder = new DistortionFinder(resolution, container);
                    finder.makeExamples();
                }
            };
            executor.execute(worker);
        }
        executor.shutdown();

        // Wait until all threads finish
        while (!executor.isTerminated()) {
        }
    }

}
