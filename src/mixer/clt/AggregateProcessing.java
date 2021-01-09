/*
 * The MIT License (MIT)
 *
 * Copyright (c) 2011-2021 Rice University, Baylor College of Medicine, Aiden Lab
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
@SuppressWarnings({"UnusedAssignment", "unused"})
public class AggregateProcessing {

    public static void main(String[] argv) throws Exception {

        String pro = "/Users/mss/Desktop/hic_files/mss_pro_30.hic";
        String sen = "/Users/mss/Desktop/hic_files/mss_sen_30.hic";

        String folder = "prosen_dice";
        String[] strings = new String[]{"dice", "-r", "100000", "-k", "KR", "--subsample", "30",
                pro + "," + sen, "2,11,10",
                "/Users/mss/Desktop/hic_files/prosen_100K_V11_s30", folder + "_"
        };
        System.out.println("-----------------------------------------------------");
        //MixerTools.main(strings);
        System.gc();

        folder = "covid_dice";
        strings = new String[]{"dice", "-r", "100000", "-k", "KR", "--subsample", "30",
                "/Users/mss/Desktop/hic_files/marianna/HIC68_30.hic," +
                        "/Users/mss/Desktop/hic_files/marianna/HIC80_30.hic",
                "2,11,10", // 2,11,10
                "/Users/mss/Desktop/hic_files/marianna/new_dice_V11_100K_s30", folder + "_"
        };
        System.out.println("-----------------------------------------------------");
        //MixerTools.main(strings);
        System.gc();


        String stem = "/Users/mss/Desktop/hic_files/prosen_100K_V11_s30/prosen_dice_";
        String dice100K = stem + "10_clusters_mss_pro_30.subcompartment.bed," +
                stem + "2_clusters_mss_pro_30.subcompartment.bed," +
                stem + "3_clusters_mss_pro_30.subcompartment.bed," +
                stem + "4_clusters_mss_pro_30.subcompartment.bed," +
                stem + "5_clusters_mss_pro_30.subcompartment.bed," +
                stem + "6_clusters_mss_pro_30.subcompartment.bed," +
                stem + "7_clusters_mss_pro_30.subcompartment.bed," +
                stem + "8_clusters_mss_pro_30.subcompartment.bed," +
                stem + "9_clusters_mss_pro_30.subcompartment.bed";

        strings = new String[]{"shuffle", "-r", "100000", "-k", "KR", "-w", "30",
                pro,
                dice100K,
                "/Users/mss/Desktop/hic_files/prosen_100K_V11_s30/shuffle_dice_PRO",
                "C10,C2,C3,C4,C5,C6,C7,C8,C9"
        };
        System.out.println("-----------------------------------------------------");
        MixerTools.main(strings);
        System.gc();

        strings = new String[]{"shuffle", "-r", "100000", "-k", "KR", "-w", "30",
                sen,
                dice100K,
                "/Users/mss/Desktop/hic_files/prosen_100K_V11_s30/shuffle_dice_SEN",
                "C10,C2,C3,C4,C5,C6,C7,C8,C9"
        };
        System.out.println("-----------------------------------------------------");
        MixerTools.main(strings);
        System.gc();


        stem = "/Users/mss/Desktop/hic_files/marianna/new_dice_V11_100K_s30/";
        dice100K = stem + "covid_dice_10_clusters_HIC68_30.subcompartment.bed," +
                stem + "covid_dice_2_clusters_HIC68_30.subcompartment.bed," +
                stem + "covid_dice_3_clusters_HIC68_30.subcompartment.bed," +
                stem + "covid_dice_4_clusters_HIC68_30.subcompartment.bed," +
                stem + "covid_dice_5_clusters_HIC68_30.subcompartment.bed," +
                stem + "covid_dice_6_clusters_HIC68_30.subcompartment.bed," +
                stem + "covid_dice_7_clusters_HIC68_30.subcompartment.bed," +
                stem + "covid_dice_8_clusters_HIC68_30.subcompartment.bed," +
                stem + "covid_dice_9_clusters_HIC68_30.subcompartment.bed";

        strings = new String[]{"shuffle", "-r", "100000", "-k", "KR", "-w", "30",
                "/Users/mss/Desktop/hic_files/marianna/HIC68_30.hic",
                dice100K,
                "/Users/mss/Desktop/hic_files/marianna/new_dice_V11_100K_s30/shuffle_dice_100K_68.hic",
                "C10,C2,C3,C4,C5,C6,C7,C8,C9"
        };
        System.out.println("-----------------------------------------------------");
        MixerTools.main(strings);
        System.gc();

        strings = new String[]{"shuffle", "-r", "100000", "-k", "KR", "-w", "30",
                "/Users/mss/Desktop/hic_files/marianna/HIC80_30.hic",
                dice100K,
                "/Users/mss/Desktop/hic_files/marianna/new_dice_V11_100K_s30/shuffle_dice_100K_80.hic",
                "C10,C2,C3,C4,C5,C6,C7,C8,C9"
        };
        System.out.println("-----------------------------------------------------");
        MixerTools.main(strings);
        System.gc();


    }
}
