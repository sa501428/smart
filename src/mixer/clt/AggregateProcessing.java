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

        String folder = "phoenix_similarity_reboot_cosine_sub50";
        strings = new String[]{"slice", "-r", "100000", "-k", "KR", "-w", "2",
                "--type", "cosine", "--zscore", "--subsample", "50",
                "--compare", refs,
                file14, "2,11,10",
                "/Users/mss/Desktop/SLICE.work/tempslice/" + folder, folder + "_"
        };
        System.out.println("-----------------------------------------------------");
        //MixerTools.main(strings);
        System.gc();

        //====================

        strings = new String[]{"shuffle", "-r", "100000", "-k", "KR", "-w", "16",
                file14, refs,
                "/Users/mss/Desktop/SLICE.work/tempslice/shuffle1",
                "RH2014,SCI,SNIPER"
        };

        strings = new String[]{"shuffle", "-r", "100000", "-k", "KR", "-w", "16",
                file14,
                "/Users/mss/Desktop/SLICE.work/tempslice/metric4umap_v1_cosine/metric4umap_v1_cosine_5_clusters_gm12878_rh14_30.subcompartment.bed,/Users/mss/Desktop/SLICE.work/tempslice/metric4umap_v1_cosine_zscore/metric4umap_v1_cosine_5_clusters_gm12878_rh14_30.subcompartment.bed,/Users/mss/Desktop/SLICE.work/tempslice/metric4umap_v1_l2_zscore/metric4umap_v1_l2_zscore_5_clusters_gm12878_rh14_30.subcompartment.bed,/Users/mss/Desktop/SLICE.work/tempslice/metric4umap_v1_js/metric4umap_v1_js_5_clusters_gm12878_rh14_30.subcompartment.bed,/Users/mss/Desktop/SLICE.work/tempslice/final_baseline_lots_to_umap/final_baseline_lots_to_umap_5_clusters_gm12878_rh14_30.subcompartment.bed,/Users/mss/Desktop/SLICE.work/tempslice/final_baseline_JSDist/final_baseline_JSDist_5_clusters_gm12878_rh14_30.subcompartment.bed",
                "/Users/mss/Desktop/SLICE.work/tempslice/shuffle1",
                "newcosine,newzscorecosine,newzscoreL2,newJS,prevSlice,distMatrixJSPrev"
        };


        System.out.println("-----------------------------------------------------");
        //MixerTools.main(strings);
        System.gc();

        strings = new String[]{"shuffle", "-r", "25000", "-k", "KR", "-w", "32",
                file14,
                "/Users/mss/Desktop/SLICE.work/tempslice/metric4umap_v1_cosine/metric4umap_v1_cosine_5_clusters_gm12878_rh14_30.subcompartment.bed,/Users/mss/Desktop/SLICE.work/tempslice/metric4umap_v1_cosine_zscore/metric4umap_v1_cosine_5_clusters_gm12878_rh14_30.subcompartment.bed,/Users/mss/Desktop/SLICE.work/tempslice/metric4umap_v1_l2_zscore/metric4umap_v1_l2_zscore_5_clusters_gm12878_rh14_30.subcompartment.bed,/Users/mss/Desktop/SLICE.work/tempslice/metric4umap_v1_js/metric4umap_v1_js_5_clusters_gm12878_rh14_30.subcompartment.bed,/Users/mss/Desktop/SLICE.work/tempslice/final_baseline_lots_to_umap/final_baseline_lots_to_umap_5_clusters_gm12878_rh14_30.subcompartment.bed,/Users/mss/Desktop/SLICE.work/tempslice/final_baseline_JSDist/final_baseline_JSDist_5_clusters_gm12878_rh14_30.subcompartment.bed",
                "/Users/mss/Desktop/SLICE.work/tempslice/" + folder,
                "newcosine,newzscorecosine,newzscoreL2,newJS,prevSlice,distMatrixJSPrev"
        };


        System.out.println("-----------------------------------------------------");
        //MixerTools.main(strings);
        System.gc();

        String degronbeds = "/Users/mss/Desktop/slice_degron/slice_wt_degron_30.hic_50000_16_9_clusters_wt_degron_30.subcompartment.bed," +
                "/Users/mss/Desktop/slice_degron/slice_treated_degron_30.hic_50000_16_9_clusters_treated_degron_30.subcompartment.bed";

        String degronfolder = "/Users/mss/Desktop/degron_shuffle";
        String wthic = "/Volumes/AidenLabWD7/Backup/AidenLab/LocalFiles/degron/notreat_nosync_combined.hic";
        String treatedhic = "/Volumes/AidenLabWD7/Backup/AidenLab/LocalFiles/degron/6hrtreat_nosync_combined.hic";


        strings = new String[]{"shuffle", "-r", "50000", "-k", "KR", "-w", "50",
                wthic, degronbeds, degronfolder,
                "wtbed_wtnic,trbed_wthic"
        };
        System.out.println("-----------------------------------------------------");
        //MixerTools.main(strings);
        System.gc();

        strings = new String[]{"shuffle", "-r", "50000", "-k", "KR", "-w", "50",
                treatedhic, degronbeds, degronfolder,
                "wtbed_trnic,trbed_trhic"
        };
        System.out.println("-----------------------------------------------------");
        //MixerTools.main(strings);
        System.gc();


        strings = new String[]{"shuffle", "-r", "100000", "-k", "KR", "-w", "16",
                file14,
                "/Users/mss/Desktop/SLICE.work/tempslice/phoenix_similarity_reboot_cosine/phoenix_similarity_reboot_cosine_5_clusters_gm12878_rh14_30.subcompartment.bed," +
                        "/Users/mss/Desktop/SLICE.work/subcompartment_analysis/slice/existing/GSE63525_GM12878_subcompartments.bed," +
                        "/Users/mss/Desktop/SLICE.work/tempslice/phoenix_similarity_reboot_cosine_sub50/phoenix_similarity_reboot_cosine_sub50_5_clusters_gm12878_rh14_30.subcompartment.bed",
                "/Users/mss/Desktop/SLICE.work/tempslice/shuffle2",
                "phoenix_cosine_zscore_again,phoenix_rh14,phoenix_cosine_zscore_sub50"
        };
        System.out.println("-----------------------------------------------------");
        MixerTools.main(strings);
        System.gc();

    }
}
