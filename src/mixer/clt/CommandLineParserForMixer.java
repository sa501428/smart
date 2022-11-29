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

package mixer.clt;

import jargs.gnu.CmdLineParser;
import javastraw.reader.type.NormalizationHandler;
import javastraw.reader.type.NormalizationType;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

/**
 * Command Line Parser for Mixer commands
 *
 * @author Muhammad Shamim
 */
public class CommandLineParserForMixer extends CmdLineParser {

    private final Option verboseOption = addBooleanOption('v', "verbose");
    private final Option helpOption = addBooleanOption('h', "help");
    private final Option versionOption = addBooleanOption('V', "version");
    private final Option normalizationTypeOption = addStringOption('k', "normalization");
    private final Option multipleResolutionsOption = addStringOption('r', "resolutions");
    private final Option randomSeedsOption = addStringOption("seed");
    private final Option genomeOption = addStringOption("genome");
    private final Option zScoreOption = addBooleanOption("zscore");
    private final Option bothNormsOption = addBooleanOption("both-norms");
    private final Option appendOption = addBooleanOption("append");
    private final Option scaleOption = addBooleanOption("scale");
    private final Option postNormOption = addBooleanOption("post-norm");
    private final Option skipCheckOption = addBooleanOption("skip-check");
    private final Option windowOption = addIntegerOption('w', "window");
    private final Option skipIntraOption = addBooleanOption("skip-intra");


    public CommandLineParserForMixer() {
    }

    /**
     * String flags
     */
    public String getGenomeOption() {
        return optionToString(genomeOption);
    }

    public NormalizationType getNormalizationTypeOption(NormalizationHandler normalizationHandler) {
        return retrieveNormalization(optionToString(normalizationTypeOption), normalizationHandler);
    }

    public NormalizationType[] getTwoNormsTypeOption(NormalizationHandler normalizationHandler) {
        String normString = optionToString(normalizationTypeOption);
        if (normString != null && normString.length() > 3 && normString.contains(",")) {
            String[] nStrings = normString.split(",");
            if (nStrings.length == 2) {
                return new NormalizationType[]{
                        retrieveNormalization(nStrings[0], normalizationHandler),
                        retrieveNormalization(nStrings[1], normalizationHandler)
                };
            }
        }
        return null;
    }

    public List<Integer> getMultipleResolutionOptions() {
        return optionToIntegerList(multipleResolutionsOption);
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

    private List<String> optionToStringList(Option option) {
        Object opt = getOptionValue(option);
        return opt == null ? null : new ArrayList<>(Arrays.asList(opt.toString().split(",")));
    }

    private List<Integer> optionToIntegerList(Option option) {
        Object opt = getOptionValue(option);
        if (opt == null) return null;
        List<String> tempList = new ArrayList<>(Arrays.asList(opt.toString().split(",")));
        List<Integer> intList = new ArrayList<>();
        for (String temp : tempList) {
            intList.add(Integer.parseInt(temp));
        }
        return intList;
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

    public boolean getZScoreOption() {
        return optionToBoolean(zScoreOption);
    }

    public boolean getScaleOption() {
        return optionToBoolean(scaleOption);
    }

    public boolean getUseBothNormsOption() {
        return optionToBoolean(bothNormsOption);
    }

    public boolean getAppendOption() {
        return optionToBoolean(appendOption);
    }

    public boolean getSkipCheckOption() {
        return optionToBoolean(skipCheckOption);
    }

    public boolean getPostNormOption() {
        return optionToBoolean(postNormOption);
    }

    private int optionToInt(Option option, int defaultNum) {
        Object opt = getOptionValue(option);
        return opt == null ? defaultNum : ((Number) opt).intValue();
    }

    public int getWindowSizeOption(int defaultNum) {
        return optionToInt(windowOption, defaultNum);
    }

    public boolean getSkipIntraOption() {
        return optionToBoolean(skipIntraOption);
    }
}