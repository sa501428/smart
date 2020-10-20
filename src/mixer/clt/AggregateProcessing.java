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

import mixer.utils.slice.SliceMatrix;


/**
 * Created for testing multiple CLTs at once
 * Basically scratch space
 */
public class AggregateProcessing {

    public static boolean useDerivative = false;
    public static boolean useL1Norm = false;

    public static void main(String[] argv) throws Exception {


        String[] strings;
        String refs = "/Users/muhammad/Desktop/research/drinks/existingmethods/GSE63525_GM12878_subcompartments.bed" + "+" +
                "/Users/muhammad/Desktop/research/drinks/existingmethods/ultra_res_100k_default_clean_outliers.bed" + "+" +
                "/Users/muhammad/Desktop/research/drinks/existingmethods/New_STRICT_gold_standard.bed";
        String file42 = "/Volumes/AidenLabWD7/Backup/AidenLab/LocalFiles/gm12878/GM12878_ultra_42B_1k_30.hic";
        String file14 = "/Users/muhammad/Desktop/insitumboi/combined_GM12878_insitu_combined_30.hic";
        String file15B = "/Volumes/AidenLabWD7/Backup/AidenLab/LocalFiles/gm12878/GM12878_insitu_combined_15B_30.hic";
        String folder;
        String primary14 = "/Volumes/AidenLabWD7/Backup/AidenLab/LocalFiles/gm12878/CELL2014/GSE63525_GM12878_insitu_primary_30.hic";
        String replicate14 = "/Volumes/AidenLabWD7/Backup/AidenLab/LocalFiles/gm12878/CELL2014/GSE63525_GM12878_insitu_replicate_30.hic";

        file14 = "/Users/mss/Desktop/gm12878_rh14_30.hic";

        SliceMatrix.USE_CORRELATION = true;

        folder = "alpha_r100k_w2_gwkr";
        strings = new String[]{"slice", "-r", "100000", "-k", "GW_KR", "-w", "2", //"--verbose",
                file14, "/Users/mss/Desktop/slice/" + folder, folder + "_", refs};
        System.out.println("-----------------------------------------------------");
        //MixerTools.main(strings);
        System.gc();
    
        String[] prefixes = new String[]{"Primary", "Replicate", "RH14"};
        String[] files = new String[]{primary14, replicate14, file14};
    
        for (int i = 0; i < prefixes.length; i++) {
            String myfile = files[i];
            String myprefix = prefixes[i];
            for (int width : new int[]{4, 8}) {
                folder = "Delphi_" + myprefix + "_50k_KR_w" + width;
                strings = new String[]{"slice", "-r", "100000", "-k", "KR", "-w", "" + width, //"--verbose",
                        myfile, "/Users/muhammad/Desktop/research/slice/" + folder, folder + "_", refs};
                System.out.println("-----------------------------------------------------");
                //MixerTools.main(strings);
                System.gc();
            }
        }
    
    }
}
