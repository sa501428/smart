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

package mixer.clt;

import jargs.gnu.CmdLineParser;
import javastraw.type.NormalizationHandler;
import javastraw.type.NormalizationType;
import mixer.algos.Slice;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

/**
 * Command Line Parser for Mixer commands (hiccups, arrowhead, apa)
 * @author Muhammad Shamim
 */
public class CommandLineParserForMixer extends CmdLineParser {

    // used flags
    // wmnxcrplafdptkqbvuhgjyz

    // available flags
    // oes
    // universal
    private final Option verboseOption = addBooleanOption('v', "verbose");
    private final Option helpOption = addBooleanOption('h', "help");
    private final Option versionOption = addBooleanOption('V', "version");
    private final Option normalizationTypeOption = addStringOption('k', "normalization");


    // General
    private final Option matrixSizeOption = addIntegerOption('m', "matrix-window-width");
    private final Option multipleChromosomesOption = addStringOption('c', "chromosomes");
    private final Option multipleResolutionsOption = addStringOption('r', "resolutions");
    private final Option threadNumOption = addIntegerOption('z', "threads");
    private final Option randomSeedsOption = addStringOption("random-seeds");
    private final Option convolutionOption = addStringOption("conv1d");

    // APA
    private final Option apaWindowOption = addIntegerOption('w', "window");
    private final Option apaMinValOption = addDoubleOption('n', "min_dist");
    private final Option apaMaxValOption = addDoubleOption('x', "max_dist");
    private final Option multipleCornerRegionDimensionsOption = addStringOption('q', "corner-width");

    // HICCUPS
    private final Option fdrOption = addStringOption('f', "fdr-thresholds");
    private final Option windowOption = addStringOption('i', "window-width");
    private final Option peakOption = addStringOption('p', "peak-width");
    private final Option clusterRadiusOption = addStringOption('d', "centroid-radii");
    private final Option thresholdOption = addStringOption('t', "postprocessing-thresholds");
    private final Option cpuVersionHiCCUPSOption = addBooleanOption('j', "cpu");
    private final Option restrictSearchRegionsOption = addBooleanOption('y', "restrict");

    // previously for AFA
    private final Option multipleAttributesOption = addStringOption('a', "attributes");

    // for GRIND
    private final Option useObservedOverExpectedOption = addBooleanOption("observed-over-expected");
    private final Option useDenseLabelsOption = addBooleanOption("dense-labels");
    private final Option useWholeGenome = addBooleanOption("whole-genome");
    private final Option useDiagonalOption = addBooleanOption("diagonal");
    private final Option cornerOffBy = addIntegerOption("off-from-diagonal");
    private final Option stride = addIntegerOption("stride");
    private final Option useDontIgnoreDirectionOrientationOption = addBooleanOption("use-feature-orientation");
    private final Option useOnlyMakePositiveExamplesOption = addBooleanOption("only-make-positives");
    private final Option generateImageFormatPicturesOption = addStringOption("img");
    private final Option useAmorphicLabelingOption = addBooleanOption("amorphic-labeling");
    private final Option useTxtInsteadOfNPYOption = addBooleanOption("text-output");

    //iterate-down-diagonal, iterate-on-list, iterate-distortions, iterate-domains
    private final Option useListIterationOption = addBooleanOption("iterate-on-list");
    private final Option useDomainOption = addBooleanOption("iterate-domains");
    private final Option useIterationDownDiagonalOption = addBooleanOption("iterate-down-diagonal");
    private final Option useDistortionOption = addBooleanOption("iterate-distortions");

    // DRINKS
    private final Option useDerivativeOption = addBooleanOption("derivative");
    private final Option ignoreDerivativeOption = addBooleanOption("ignore-derivative");
    private final Option usingRowNormalizationOption = addBooleanOption("normalize-rows");


    public CommandLineParserForMixer() {
    }

    public int getUsingDerivativeStatus() {
        Object opt = getOptionValue(useDerivativeOption);
        if (opt != null) return Slice.USE_ONLY_DERIVATIVE;
        opt = getOptionValue(ignoreDerivativeOption);
        if (opt != null) return Slice.IGNORE_DERIVATIVE;
        return 0;
    }

    public boolean getUsingRowNomalizationStatus() {
        return optionToBoolean(usingRowNormalizationOption);
    }


    // for GRIND
    public boolean getUseObservedOverExpectedOption() {
        return optionToBoolean(useObservedOverExpectedOption);
    }

    public boolean getUseAmorphicLabelingOption() {
        return optionToBoolean(useAmorphicLabelingOption);
    }

    public boolean getUseWholeGenome() {
        return optionToBoolean(useWholeGenome);
    }

    public boolean getUseGenomeDiagonal() {
        return optionToBoolean(useDiagonalOption);
    }

    public boolean getDenseLabelsOption() {
        return optionToBoolean(useDenseLabelsOption);
    }

    public boolean getDontIgnoreDirectionOrientationOption() {
        return optionToBoolean(useDontIgnoreDirectionOrientationOption);
    }

    public boolean getUseOnlyMakePositiveExamplesOption() {
        return optionToBoolean(useOnlyMakePositiveExamplesOption);
    }

    /**
     * String flags
     */

    public NormalizationType getNormalizationTypeOption(NormalizationHandler normalizationHandler) {
        return retrieveNormalization(optionToString(normalizationTypeOption), normalizationHandler);
    }

    /**
     * int flags
     */
    public int getAPAWindowSizeOption() {
        return optionToInt(apaWindowOption);
    }

    public int getCornerOffBy() {
        return optionToInt(cornerOffBy);
    }

    public int getStride() {
        return optionToInt(stride);
    }

    public int getMatrixSizeOption() {
        return optionToInt(matrixSizeOption);
    }

    public int getNumThreads() {
        return optionToInt(threadNumOption);
    }

    /**
     * String Set flags
     */

    List<String> getChromosomeListOption() {
        return optionToStringList(multipleChromosomesOption);
    }

    // todo fix to return list of ints
    public List<String> getMultipleResolutionOptions() {
        return optionToStringList(multipleResolutionsOption);
    }

    public String getGenerateImageFormatPicturesOption() {
        return optionToString(generateImageFormatPicturesOption);
    }

    public boolean getUseTxtInsteadOfNPY() {
        return optionToBoolean(useTxtInsteadOfNPYOption);
    }

    public long[] getMultipleSeedsOption() {
        List<String> possibleSeeds = optionToStringList(randomSeedsOption);
        if (possibleSeeds != null) {
            long[] seeds = new long[possibleSeeds.size()];
            for (int i = 0; i < seeds.length; i++) {
                seeds[i] = Long.parseLong(possibleSeeds.get(i));
            }
            return seeds;
        }
        return null;
    }

    public double[] getConvolutionOption() {
        List<String> conv1d = optionToStringList(convolutionOption);
        if (conv1d != null) {
            double[] values = new double[conv1d.size()];
            for (int i = 0; i < values.length; i++) {
                values[i] = Double.parseDouble(conv1d.get(i));
            }
            return values;
        }
        return null;
    }

    /**
     * boolean flags
     */
    private boolean optionToBoolean(Option option) {
        Object opt = getOptionValue(option);
        return opt != null && (Boolean) opt;
    }

    public boolean getHelpOption() {
        return optionToBoolean(helpOption);
    }

    public boolean getVerboseOption() {
        return optionToBoolean(verboseOption);
    }

    public boolean getVersionOption() {
        return optionToBoolean(versionOption);
    }

    /**
     * String flags
     */
    private String optionToString(Option option) {
        Object opt = getOptionValue(option);
        return opt == null ? null : opt.toString();
    }

    /**
     * int flags
     */
    private int optionToInt(Option option) {
        Object opt = getOptionValue(option);
        return opt == null ? 0 : ((Number) opt).intValue();
    }

    /**
     * double flags
     */
    private double optionToDouble(Option option) {
        Object opt = getOptionValue(option);
        return opt == null ? 0 : ((Number) opt).doubleValue();
    }

    private List<String> optionToStringList(Option option) {
        Object opt = getOptionValue(option);
        return opt == null ? null : new ArrayList<>(Arrays.asList(opt.toString().split(",")));
    }

    private NormalizationType retrieveNormalization(String norm, NormalizationHandler normalizationHandler) {
        if (norm == null || norm.length() < 1)
            return null;

        try {
            return normalizationHandler.getNormTypeFromString(norm);
        } catch (IllegalArgumentException error) {
            System.err.println("Normalization must be one of \"NONE\", \"VC\", \"VC_SQRT\", \"KR\", \"GW_KR\", \"GW_VC\", \"INTER_KR\", or \"INTER_VC\".");
            System.exit(7);
        }
        return null;
    }
}