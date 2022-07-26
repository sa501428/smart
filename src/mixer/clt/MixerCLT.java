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

import java.util.List;
import java.util.Random;

/**
 * Created by muhammadsaadshamim on 9/21/15.
 */
public abstract class MixerCLT {
    private static String usage;

    protected MixerCLT(String usage) {
        setUsage(usage);
    }

    public void readArguments(String[] args, CommandLineParserForMixer parser) {
        readMixerArguments(args, parser);
    }

    protected abstract void readMixerArguments(String[] args, CommandLineParserForMixer mixerParser);

    public abstract void run();

    private void setUsage(String newUsage) {
        usage = newUsage;
    }

    public void printUsageAndExit(int exitcode) {
        System.out.println("Usage:   mixer_tools " + usage);
        System.exit(exitcode);
    }

    protected int updateResolution(CommandLineParserForMixer mixerParser, int r0) {
        List<Integer> possibleResolutions = mixerParser.getMultipleResolutionOptions();
        if (possibleResolutions != null) {
            if (possibleResolutions.size() > 1)
                System.err.println("Only one resolution can be specified\nUsing " + possibleResolutions.get(0));
            return possibleResolutions.get(0);
        }
        return r0;
    }

    protected void updateGeneratorSeed(CommandLineParserForMixer mixerParser, Random generator) {
        long[] possibleSeeds = mixerParser.getMultipleSeedsOption();
        if (possibleSeeds != null && possibleSeeds.length > 0) {
            for (long seed : possibleSeeds) {
                generator.setSeed(seed);
            }
        }
    }
}
