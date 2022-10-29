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


import mixer.SmartTools;

/**
 * Created for testing multiple CLTs at once
 * Basically scratch space
 */
@SuppressWarnings({"UnusedAssignment", "unused"})
public class AggregateProcessing {

    public static void main(String[] argv) throws Exception {
        String[] args = new String[]{"slice", "-r", "25000", "/Users/muhammad/Desktop/hicfiles/gm12878_rh14_30.hic",
                "5,7", "/Users/muhammad/Desktop/slice2/final_hammer/result_25k_v2", "gm14_"};

        args = new String[]{"slice", "-r", "50000", "/Volumes/AidenLabWD6/hicfiles/hg38/HCT116_RAD21_30.hic",
                "5,10", "/Users/muhammad/Desktop/slice2/final_hammer/result_degron_50k_v2", "degron_"};

        SmartTools.main(args);
    }
}
