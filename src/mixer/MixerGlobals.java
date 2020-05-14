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

package mixer;

/**
 * @author Muhammad Shamim
 * @since 11/25/14
 */
public class MixerGlobals {

    public static final String versionNum = "2.06.05";
    public static final int minVersion = 6;
    public static final int bufferSize = 2097152;

    // whether MatrixZoomData should cache or not
    public static boolean useCache = false;
    public static boolean printVerboseComments = false;

    // whether instance was linked before mouse press or not
    public static boolean isLegacyOutputPrintingEnabled = false;

    public static void verifySupportedHiCFileVersion(int version) throws RuntimeException {
        if (version < minVersion) {
            throw new RuntimeException("This file is version " + version +
                    ". Only versions " + minVersion + " and greater are supported at this time.");
        }
    }
}