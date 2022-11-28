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

package mixer;

import jargs.gnu.CmdLineParser;
import mixer.clt.CLTFactory;
import mixer.clt.CommandLineParserForMixer;
import mixer.clt.MixerCLT;


public class SmartTools {

    public static final String versionNum = "0.30.2";
    public static final int bufferSize = 2097152;
    public static boolean printVerboseComments = false;
    public static int NUM_ENTRIES_TO_SKIP_MEDIAN = 0;

    public static void main(String[] argv) throws CmdLineParser.UnknownOptionException, CmdLineParser.IllegalOptionValueException {

        if (argv.length == 0 || argv[0].equals("-h") || argv[0].equals("--help") || argv[0].equals("-V") || argv[0].equals("--version")) {
            CLTFactory.generalUsage();
            System.exit(0);
        }

        CommandLineParserForMixer parser = new CommandLineParserForMixer();
        boolean help;
        boolean version;
        parser.parse(argv);

        help = parser.getHelpOption();
        version = parser.getVersionOption();
        printVerboseComments = printVerboseComments || parser.getVerboseOption();

        String[] args = parser.getRemainingArgs();

        MixerCLT instanceOfCLT;
        String cmd = "";
        if (args.length == 0) {
            instanceOfCLT = null;
        } else {
            cmd = args[0];
            instanceOfCLT = CLTFactory.getCLTCommand(cmd);
        }
        if (instanceOfCLT != null) {
            if (version) {
                System.out.println("Smart tools version " + versionNum);
            }
            if (args.length == 1 || help) {
                instanceOfCLT.printUsageAndExit(1);
            }

            instanceOfCLT.readArguments(args, parser);
            instanceOfCLT.run();
        } else {
            throw new RuntimeException("Unknown command: " + cmd);
        }
    }
}
