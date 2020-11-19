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

import mixer.MixerTools;


/**
 * Created for testing multiple CLTs at once
 * Basically scratch space
 */
public class AggregateProcessing {

    public static boolean useL1Norm = false;

    public static void main(String[] argv) throws Exception {

        String[] strings;

        String refs = "/Users/mss/Desktop/SLICE.work/subcompartment_analysis/slice/existing/GSE63525_GM12878_subcompartments.bed," +
                "/Users/mss/Desktop/SLICE.work/subcompartment_analysis/slice/existing/GM12878_SCI_sub_compartments.bed," +
                "/Users/mss/Desktop/SLICE.work/subcompartment_analysis/slice/existing/GM12878_track_hg19.bed";
        String file14 = "/Users/mss/Desktop/hic_files/gm12878_rh14_30.hic";

        String folder = "collins_dice_test";
        strings = new String[]{"dice", "-r", "100000", "-k", "KR",
                "--type", "zscore-cosine", "--subsample", "50",
                "/Volumes/AidenLabWD7/Backup/AidenLab/LocalFiles/collins/map1_30.hic,/Volumes/AidenLabWD7/Backup/AidenLab/LocalFiles/collins/map2_30.hic",
                "2,11,10",
                "/Users/mss/Desktop/SLICE.work/tempslice/" + folder, folder + "_"
        };
        System.out.println("-----------------------------------------------------");
        //MixerTools.main(strings);
        System.gc();


        strings = new String[]{"shuffle", "-r", "100000", "-k", "KR", "-w", "16",
                "/Volumes/AidenLabWD7/Backup/AidenLab/LocalFiles/collins/map1_30.hic",
                "/Users/mss/Desktop/SLICE.work/tempslice/collins_dice_test/collins_dice_test_10_clusters_map1_30.subcompartment.bed",
                "/Users/mss/Desktop/SLICE.work/tempslice/shuffle4",
                "collins1_10_again"
        };
        System.out.println("-----------------------------------------------------");
        MixerTools.main(strings);
        System.gc();

        strings = new String[]{"shuffle", "-r", "100000", "-k", "KR", "-w", "16",
                "/Volumes/AidenLabWD7/Backup/AidenLab/LocalFiles/collins/map2_30.hic",
                "/Users/mss/Desktop/SLICE.work/tempslice/collins_dice_test/collins_dice_test_10_clusters_map1_30.subcompartment.bed",
                "/Users/mss/Desktop/SLICE.work/tempslice/shuffle4",
                "collins2_10_again"
        };
        System.out.println("-----------------------------------------------------");
        MixerTools.main(strings);
        System.gc();

    }
}
