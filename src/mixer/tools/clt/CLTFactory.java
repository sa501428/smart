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

package mixer.tools.clt;

import mixer.MixerGlobals;
import mixer.tools.clt.mixer.CompareLists;
import mixer.tools.clt.mixer.LoopDomains;
import mixer.tools.dev.Distortion;
import mixer.tools.dev.Drink;
import mixer.tools.dev.Grind;
import mixer.tools.dev.Shuffle;


/**
 * Factory for command line tools to call different functions
 *
 * @author Muhammad Shamim
 * @since 1/30/2015
 */
public class CLTFactory {

    // Commenting some out because we're not going to release all these when we release CLT
    private final static String[] commandLineToolUsages = {
            "grind---"
    };

    public static void generalUsage() {
        System.out.println("Mixer Tools Version " + MixerGlobals.versionNum);
        System.out.println("Usage:");
        for (String usage : commandLineToolUsages) {
            System.out.println("\t" + usage);
        }
        System.out.println("\t" + "-h, --help print help");
        System.out.println("\t" + "-v, --verbose verbose mode");
        System.out.println("\t" + "-V, --version print version");
        System.out.println("Type mixer_tools <commandName> for more detailed usage instructions");
    }

    public static MixerCLT getCLTCommand(String cmd) {

        cmd = cmd.toLowerCase();
        if (cmd.equals("compare")) {
            return new CompareLists();
        } else if (cmd.equals("shuffle")) {
            return new Shuffle();
        } else if (cmd.equals("loop_domains")) {
            return new LoopDomains();
        } else if (cmd.startsWith("drink")) {
            return new Drink(cmd);
        } else if (cmd.equals("grind")) {
            return new Grind();
        } else if (cmd.equals("distort")) {
            return new Distortion();
        }

        return null;
    }
}
