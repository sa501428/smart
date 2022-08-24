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


import mixer.MixerTools;

/**
 * Created for testing multiple CLTs at once
 * Basically scratch space
 */
@SuppressWarnings({"UnusedAssignment", "unused"})
public class AggregateProcessing {

    public static void main(String[] argv) throws Exception {

        String path = "/Users/mshamim/Desktop/hg38_files/";
        String[] stems = new String[]{"gm", "hepg2", "imr", "hct", "k562"};
        String[] kStems = new String[]{"GM12878", "HepG2", "IMR90", "HCT116", "K562"};

        for (int f : new int[]{0, 3}) {//= 0; f < stems.length; f++
            String stem = stems[f];
            String kStem = kStems[f];
            String file = "/Volumes/AidenLabWD6/hicfiles/hg38/" + kStem + "_30.hic";
            for (String k : new String[]{"INTER_SCALE"}) {// "GW_SCALE", "KR" ,normtype2[f]
                for (int r : new int[]{100}) {
                    String beds = path + "kyle/hmm_clustering/" + kStem + "_imputed_slice_hg38_6_clusters.bed,"
                            + path + "kyle/hmm_clustering2/" + kStem + "_imputed.bed,"
                            + path + stem + "_100k/c2/slice_GHMM_impute_1.bed,"
                            + path + stem + "_100k/c2/slice_GHMM_impute_3.bed,"
                            + path + stem + "_100k/c2/slice_GHMM_impute_4.bed,"
                            //+path+stem+"_100k/kmedians/"+stem+"_5_kmedians_clusters.bed,"
                            + path + stem + "_100k/kmedians/" + stem + "_6_kmedians_clusters.bed";
                    String labels = "GS,GS2,I1,I3,I4,S6";

                    String[] strings = new String[]{"shuffle",
                            "-r", r + "000", "-k", k, "-w", "16",
                            file,
                            beds,
                            path + stem + "_100k/shuffle_w16_redo_" + r + "_" + k,
                            labels
                    };
                    System.out.println("-----------------------------------------------------");
                    MixerTools.main(strings);
                    System.gc();
                }
            }
        }

        /*

        String[] files = new String[]{
                //"/Users/mshamim/Desktop/hicfiles/gm12878_rh14_30.hic" //,
                //"/Users/mshamim/Desktop/hicfiles/GM12878_intact_18.7B_8.15.20_30.hic",
                //"/Volumes/AidenLabWD7/Backup/AidenLab/LocalFiles/gm12878/GM12878_insitu_combined_15B_30.hic"
        };

        files = new String[]{
                "/Users/mshamim/Desktop/hicfiles/SCALE/hap1_SCALE_30.hic",
                "/Users/mshamim/Desktop/hicfiles/SCALE/imr90_rh14_SCALE_30.hic",
                "/Users/mshamim/Desktop/hicfiles/SCALE/K562_2014_SCALE_30.hic",
                "/Users/mshamim/Desktop/hicfiles/gm12878_rh14_30.hic"
        };

        String[] stems = new String[]{
                "hap1",
                "imr",
                "k562",
                "gm"
        };

        files = new String[]{
                "/Users/mshamim/Desktop/subsampling_experiment/primary_15.hic",
                "/Users/mshamim/Desktop/subsampling_experiment/primary_29.hic",
                "/Users/mshamim/Desktop/subsampling_experiment/primary_43.hic",
                "/Users/mshamim/Desktop/subsampling_experiment/primary_58.hic",
                "/Users/mshamim/Desktop/hicfiles/GM12878_primary14_30.hic",
                "/Users/mshamim/Desktop/hicfiles/gm12878_rh14_30.hic"
        };

        stems = new String[]{"p15", "p29", "p43", "p58", "primary", "gmMega"};



        for (int f = 5; f < files.length; f++) {// files.length
            String file = files[f];
            String stem = stems[f];
            for (int res : new int[]{100}) { //  ,100000,   50000,25000,10000 100000 100000 50000
                String folder = stem;
                String[] strings = new String[]{"slice", "-r", res + "000", "--encode-mode",
                        file, "2,12,4",
                        "/Users/mshamim/Desktop/reSLICE/encode_z8_" + res + "000_" + folder,
                        folder + "_"
                };
                System.out.println("-----------------------------------------------------");
                //MixerTools.main(strings);
                System.gc();
            }
        }
        System.gc();

        // 100 exp +/-
        // 101 ^3
        // 102 only corr of standard
        // 103 fix mask
        // 104 clean corr
        // 105 clean corr 10
        // 107 fifty mb split
        // 108 redid post analysis
        // restructure
        int id = 123;

        String[] strings = new String[]{"slice", "-r", "100000", "--encode-mode",
                "/Users/mshamim/Desktop/various_hic_files/HCT116_Degron/treated_degron_30_25k.hic",
                "2,14,4",
                "/Users/mshamim/Desktop/reSLICE/encode_HCT_Z" + id + "_100K",
                "hct"
        };
        System.out.println("-----------------------------------------------------");
        MixerTools.main(strings);

        System.gc();

        strings = new String[]{"slice", "-r", "100000", "--encode-mode",
                "/Users/mshamim/Desktop/hicfiles/gm12878_rh14_30.hic",
                "2,14,4",
                "/Users/mshamim/Desktop/reSLICE/encode_GM_Z" + id + "_100K",
                "hct"
        };
        System.out.println("-----------------------------------------------------");
        MixerTools.main(strings);

        /*

        String beds = "/Users/mshamim/Desktop/reSLICE/80X_beds/gmMega_SLICE_800__5_kmeans_clusters.bed,/Users/mshamim/Desktop/reSLICE/80X_beds/gmMega_SLICE_803__5_kmeans_clusters.bed,/Users/mshamim/Desktop/reSLICE/80X_beds/gmMega_SLICE_804__5_kmeans_clusters.bed,/Users/mshamim/Desktop/reSLICE/80X_beds/p15_SLICE_803__5_kmeans_clusters.bed,/Users/mshamim/Desktop/reSLICE/80X_beds/p15_SLICE_804__5_kmeans_clusters.bed,/Users/mshamim/Desktop/reSLICE/existing/GSE63525_GM12878_subcompartments.bed";
        String labels = "gmMega_800,gmMega_803,gmMega_804,p15_803,p15_804,rh2014";

        beds = "/Users/mshamim/Desktop/reSLICE/SLICE804/GM12878_SCI_sub_compartments.bed,/Users/mshamim/Desktop/reSLICE/SLICE804/GM12878_track_hg19.bed,/Users/mshamim/Desktop/reSLICE/SLICE804/GSE63525_GM12878_subcompartments.bed,/Users/mshamim/Desktop/reSLICE/SLICE804/SLICE_GM_MEGA_100K.bed,/Users/mshamim/Desktop/reSLICE/SLICE804/SLICE_P15_100K.bed,/Users/mshamim/Desktop/reSLICE/SLICE804/SLICE_P29_100K.bed,/Users/mshamim/Desktop/reSLICE/SLICE804/SLICE_P43_100K.bed,/Users/mshamim/Desktop/reSLICE/SLICE804/SLICE_P58_100K.bed,/Users/mshamim/Desktop/reSLICE/SLICE804/SLICE_PRIMARY_100K.bed";
        labels = "SCI,SNIPER,RH2014,SLICE_MEGA,SLICE_P15,SLICE_P29,SLICE_P43,SLICE_P58,SLICE_PRIMARY";

        beds = "/Users/mshamim/Desktop/reSLICE/SLICE804/GM12878_SCI_sub_compartments.bed," +
                "/Users/mshamim/Desktop/reSLICE/SLICE804/GM12878_track_hg19.bed," +
                "/Users/mshamim/Desktop/reSLICE/SLICE804/GSE63525_GM12878_subcompartments.bed";
        labels = "SCI,SNIPER,RH2014";


        beds = "/Users/mshamim/Desktop/reSLICE/SLICE804/GM12878_SCI_sub_compartments.bed," +
                "/Users/mshamim/Desktop/reSLICE/SLICE804/GM12878_track_hg19.bed," +
                "/Users/mshamim/Desktop/reSLICE/SLICE804/GSE63525_GM12878_subcompartments.bed," +
                "/Users/mshamim/Desktop/various_hic_files/gm_451/gm_5_kmedians_clusters.bed," +
                "/Users/mshamim/Desktop/various_hic_files/gm_452/gm_5_kmeans_clusters.bed";
        labels = "SCI,SNIPER,RH2014,SLICE_KMEDIANS,SLICE_KMEANS";


        for (int f = 5; f < files.length; f++) {//
            String file = files[f];
            String stem = stems[f];
            for (String k : new String[]{"INTER_KR"}) {// "GW_SCALE", "KR" ,normtype2[f]
                for (int r : new int[]{100}) { // 100 50, 25
                    String[] strings = new String[]{"shuffle",
                            "-r", r + "000", "-k", k, "-w", "" + 16 * (100 / r),
                            file,
                            beds,
                            "/Users/mshamim/Desktop/shuffle_vault/shuffle_905_" + stem + "_" + r + "_" + k,
                            labels
                    };
                    System.out.println("-----------------------------------------------------");
                    //MixerTools.main(strings);
                    System.gc();
                }
            }
        }


        String[] strings;
        String changes = "2,A1,34,139,34;" +
                "3,A2,152,251,152;" +
                "5,B1,220,20,60;" +
                "4,B2,255,255,0;" +
                "1,B3,112,128,144";

        String[] bedfiles = new String[]{
                "/Users/mshamim/Desktop/reSLICE/80X_beds/p15_SLICE_804__5_kmeans_clusters.bed",
                "/Users/mshamim/Desktop/reSLICE/80X_beds/p29_SLICE_804__5_kmeans_clusters.bed",
                "/Users/mshamim/Desktop/reSLICE/80X_beds/p43_SLICE_804__5_kmeans_clusters.bed",
                "/Users/mshamim/Desktop/reSLICE/80X_beds/p58_SLICE_804__5_kmeans_clusters.bed",
                "/Users/mshamim/Desktop/reSLICE/80X_beds/primary_SLICE_804__5_kmeans_clusters.bed"
        };

        String[] bedstem = new String[]{
                "P15",
                "P29",
                "P43",
                "P58",
                "PRIMARY"
        };

        for (int z = 0; z < bedfiles.length; z++) {
            strings = new String[]{"rename",
                    changes,
                    bedfiles[z],
                    "/Users/mshamim/Desktop/reSLICE/fin_slice/SLICE_" + bedstem[z] + "_100K.bed"
            };
            //MixerTools.main(strings);
            System.gc();
        }

        /*
        String[] labels = new String[4];
        Arrays.fill(labels, "SLICE,SNIPER");
        labels[3] = "SLICE,SNIPER,RH2014,SCI";

        String[] normtype1 = new String[]{
                "INTER_SCALE", "INTER_SCALE", "INTER_SCALE", "INTER_KR"
        };
        String[] normtype2 = new String[]{
                "GW_SCALE", "GW_SCALE", "GW_SCALE", "GW_KR"
        };

        for (int f = 1; f < files.length; f++) {//
            String file = files[f];
            String stem = stems[f];
            for (String k : new String[]{normtype1[f]}) {// "GW_SCALE", "KR" ,normtype2[f]
                for (int r : new int[]{100}) { // 100 50, 25

                    String[] strings = new String[]{"shuffle",
                            "-r", r + "000", "-k", k, "-w", "" + 16 * (100 / r),
                            file,
                            //beds[f],
                            "/Users/mshamim/Desktop/reSLICE/shuffle2_" + "_slice_vs_sniper_" + stem + "_" + r + "_" + k,
                            labels[f]
                    };
                    System.out.println("-----------------------------------------------------");
                    //MixerTools.main(strings);
                    System.gc();
                }
            }
        }

        /*

        String[] types = new String[]{"RAW", "LOG", "O_E", "LOG_E(O)"};
        String[] corrs = new String[]{"IDENTITY", "COSINE", "MedianAbsDev", "PEARSON", "MANHATTAN",
                "EUCLIDEAN", "COSINE_ZSCORE", "PEARSON_ZSCORE", "COSINE_ARC", "PEARSON_ARC"};
        String[] filePrefixes = new String[]{"GM_RH_14_", "GM_Intact_18B_", "GM_15B_"};
        String[] norms = new String[]{"INTER_KR", "SCALE", "KR"};

        int[] metrics = new int[]{1, 3, 6, 7, 8, 9};//1, 3, 4, 5 , 7

        HiCMatrix.NUM_CENTROIDS = 10;

        for (int q = 0; q < 1; q++) {
            for (int x : metrics) {
                for (int y = 2; y < 3; y++) {
                    String[] strings = new String[]{"gw-umap", //umap shuffle
                            "-r", 100000 + "", // 50000
                            "-k", norms[q],
                            "-w", "" + 3,
                            "--type", "" + y,
                            "--corr", "" + x,
                            //"-c", "1,2,10,14,17,18,19",
                            files[q],
                            //"/Users/mshamim/Desktop/SLICE.Reboot/existing/SLICE_PRIMARY_100K.bed," +
                            "/Users/mshamim/Desktop/SLICE.Reboot/INTER-KR_GM_res_100000v3.18.07/INTER-KR_GM_res_100000v3.18.07__5_clusters_gm12878_rh14_30.subcompartment.bed," +
                                    "/Users/mshamim/Desktop/SLICE.Reboot/existing/SLICE_GM_100K_no2Pass.bed," +
                                    "/Users/mshamim/Desktop/SLICE.Reboot/existing/GSE63525_GM12878_subcompartments.bed",
                            "/Users/mshamim/Desktop/SLICE5/intra_umaps/NEW_SUBSET_CENT10_REALCOS_UMAP_" + filePrefixes[q] + "100K_" + types[y] + "_" + corrs[x]
                            //"/Users/mshamim/Desktop/SLICE.Reboot/intra_shuffles/REBOOT" + filePrefixes[q] + "50K_" + types[y] + "_" + corrs[x]
                            , "InterSlice100,SLICE100,RH2014"
                    };

                    System.out.println("-----------------------------------------------------");
                    MixerTools.main(strings);
                    System.gc();
                }
            }
        }

        String file14 = "/Users/mshamim/Desktop/hicfiles/gm12878_rh14_30.hic";
        for (int res : new int[]{100000}) { //100000
            String folder = "GM12878_DualNorm2_ZscoreNoWidth_" + res + "v" + MixerTools.versionNum;
            String[] strings = new String[]{"slice", "-r", res + "",
                    "-k", "KR:INTER_KR", //"--verbose",
                    file14, "2,13,10",
                    "/Users/mshamim/Desktop/SLICE.Reboot/" + folder, folder + "_"
            };
            System.out.println("-----------------------------------------------------");
            //MixerTools.main(strings);
            System.gc();
        }

        /*






        for (int q = 0; q < files.length - 1; q++) {
            for (int x = 0; x < corrs.length; x++) {
                for (int y = 0; y < types.length; y++) {
                    String[] strings = new String[]{"intra-umap", //umap shuffle
                            "-r", 100000 + "", // 50000
                            "-k", norms[q],
                            "-w", "" + 2,
                            "--type", "" + y,
                            "--corr", "" + x,
                            "-c", "1,14,19",
                            files[q],
                            "/Users/mshamim/Desktop/SLICE.Reboot/existing/SLICE_PRIMARY_100K.bed," +
                                    "/Users/mshamim/Desktop/SLICE.Reboot/existing/SLICE_GM_100K_no2Pass.bed," +
                                    "/Users/mshamim/Desktop/SLICE.Reboot/existing/GSE63525_GM12878_subcompartments.bed",
                            "/Users/mshamim/Desktop/SLICE.Reboot/intra_umaps/NEW_UMAP_PARAMS_" + filePrefixes[q] + "100K_" + types[y] + "_" + corrs[x]
                            //"/Users/mshamim/Desktop/SLICE.Reboot/intra_shuffles/REBOOT" + filePrefixes[q] + "50K_" + types[y] + "_" + corrs[x]
                            //                             , "PRIMARY100,SLICE100,RH2014"
                    };

                    System.out.println("-----------------------------------------------------");
                    //MixerTools.main(strings);
                    System.gc();
                }
            }
        }

        String file14 = "/Users/mshamim/Desktop/hicfiles/gm12878_rh14_30.hic";

        for (String norm : new String[]{"KR:INTER_KR"}) { //"GW_KR",
            String folder = "SHUFFLE_again6_DualNorm_W16_GM12878_" + norm + "_100K_v" + MixerTools.versionNum;
            String[] strings = new String[]{"shuffle-umap", //umap
                    "-r", 100000 + "",
                    "-k", norm, "-w", "" + 16,
                    file14,
                    //refs + "," +
                    "/Users/mshamim/Desktop/SLICE.Reboot/existing/GSE63525_GM12878_subcompartments.bed," +
                            "/Users/mshamim/Desktop/SLICE.Reboot/existing/SLICE_GM_100K_no2Pass.bed," +
                            "/Users/mshamim/Desktop/SLICE.Reboot/existing/SLICE_PRIMARY_100K.bed," +
                            "/Users/mshamim/Desktop/SLICE.Reboot/INTER-KR_GM_res_100000v3.18.07/INTER-KR_GM_res_100000v3.18.07__5_clusters_gm12878_rh14_30.subcompartment.bed," +
                            "/Users/mshamim/Desktop/SLICE.Reboot/existing/GM12878_track_hg19.bed," +
                            "/Users/mshamim/Desktop/SLICE.Reboot/existing/GM12878_SCI_sub_compartments.bed",
                    "/Users/mshamim/Desktop/SLICE.Reboot/shuffles/" +
                            folder
                    , "RH2014,SLICE100,PRIMARY100,INTER_SLICE,SNIPER,SCI"
            };
            System.out.println("-----------------------------------------------------");
            //MixerTools.main(strings);
            System.gc();
        }

        for (int res : new int[]{100000}) { //100000
            String folder = "INTER-KR_GM_res_" + res + "v" + MixerTools.versionNum;
            String[] strings = new String[]{"slice", "-r", res + "",
                    "-k", "INTER_KR", //"--verbose",
                    file14, "2,8,10",
                    "/Users/mshamim/Desktop/SLICE.Reboot/" + folder, folder + "_"
            };
            System.out.println("-----------------------------------------------------");
            //MixerTools.main(strings);
            System.gc();
        }

        file14 = "/Users/mshamim/Desktop/hicfiles/gm12878_rh14_30.hic";
        for (int res : new int[]{100000, 50000}) { //100000
            String folder = "GM12878_DualNorm2_" + res + "v" + MixerTools.versionNum;
            String[] strings = new String[]{"slice", "-r", res + "",
                    "-k", "KR:INTER_KR", //"--verbose",
                    file14, "2,20,10",
                    "/Users/mshamim/Desktop/SLICE.Reboot/" + folder, folder + "_"
            };
            System.out.println("-----------------------------------------------------");
            //MixerTools.main(strings);
            System.gc();

            strings = new String[]{"shuffle", "-r", res + "", "-k", "KR", "-w", "" + (16 * (100000 / res)), file14,
                    //refs + "," +
                    "/Users/mshamim/Desktop/SLICE.Reboot/" + folder +
                            "/" + folder + "__5_clusters_K562_2014_SCALE_30.subcompartment.bed",
                    "/Users/mshamim/Desktop/SLICE.Reboot/shuffle_K562_Final_res_" + res + "v" + MixerTools.versionNum,
                    "SLICE3"
            };
            System.out.println("-----------------------------------------------------");
            //MixerTools.main(strings);
            System.gc();
        }


        files = new String[]{
                "/Users/mshamim/Desktop/hicfiles/SCALE/K562_2014_SCALE_30.hic",
                "/Users/mshamim/Desktop/hicfiles/SCALE/imr90_rh14_SCALE_30.hic",
                "/Users/mshamim/Desktop/hicfiles/SCALE/hap1_SCALE_30.hic"
        };

        String[] prefixes = new String[]{
                "k562", "imr90", "hap1"
        };

        for (int res : new int[]{50000}) { //100000
            for (int q = 0; q < files.length; q++) {
                String folder = prefixes[q] + "_" + res + "v" + MixerTools.versionNum;
                String[] strings = new String[]{"slice", "-r", res + "",
                        "-k", "SCALE:INTER_SCALE", //"--verbose",
                        files[q], "2,11,10",
                        "/Users/mshamim/Desktop/SLICE.Reboot/DualNorm2/" + folder, folder + "_"
                };
                System.out.println("-----------------------------------------------------");
                //MixerTools.main(strings);
                System.gc();
            }
        }

        for (int res : new int[]{50000}) { //100000
            String folder = "GM12878_" + res + "v" + MixerTools.versionNum;
            String[] strings = new String[]{"slice", "-r", res + "",
                    "-k", "KR:INTER_KR", //"--verbose",
                    file14, "2,16,10",
                    "/Users/mshamim/Desktop/SLICE.Reboot/DualNorm2/" + folder, folder + "_"
            };
            System.out.println("-----------------------------------------------------");
            //MixerTools.main(strings);
            System.gc();
        }

        /*

        String file14 = "/Users/mshamim/Desktop/hicfiles/K562_2014_SCALE_30.hic";

        for (int res : new int[]{50000}) { //100000
            String folder = "K562_res_" + res + "v" + MixerTools.versionNum;
            String[] strings = new String[]{"slice", "-r", res + "", "-k", "SCALE", //"--verbose",
                    file14, "2,8,10",
                    "/Users/mshamim/Desktop/SLICE.Reboot/" + folder, folder + "_"
            };
            System.out.println("-----------------------------------------------------");
            //MixerTools.main(strings);
            System.gc();

            strings = new String[]{"shuffle", "-r", res + "", "-k", "KR", "-w", "" + (16 * (100000 / res)), file14,
                    //refs + "," +
                    "/Users/mshamim/Desktop/SLICE.Reboot/" + folder +
                            "/" + folder + "__5_clusters_K562_2014_SCALE_30.subcompartment.bed",
                    "/Users/mshamim/Desktop/SLICE.Reboot/shuffle_K562_Final_res_" + res + "v" + MixerTools.versionNum,
                    "SLICE3"
            };
            System.out.println("-----------------------------------------------------");
            MixerTools.main(strings);
            System.gc();
        }

        /*

        String[] strings = new String[]{"intra-shuffle", "-r", 100000 + "", "-k", "SCALE", "-w", "" + 4,
                file14,
                "/Users/mshamim/Desktop/SLICE.Reboot/existing/SLICE_PRIMARY_100K.bed," +
                        "/Users/mshamim/Desktop/SLICE.Reboot/existing/GSE63525_GM12878_subcompartments.bed",
                "/Users/mshamim/Desktop/SLICE.Reboot/intra_shuffles/INTACT18B_W4_REAL_OE_" + 100 + "K",
                "PRIMARY,RH2014"
        };

        System.out.println("-----------------------------------------------------");
        MixerTools.main(strings);
        System.gc();


        strings = new String[]{"intra-shuffle", "-r", 100000 + "", "-k", "SCALE", "-w", "" + 4,
                file14,
                "/Users/mshamim/Desktop/SLICE.Reboot/existing/SLICE_PRIMARY_100K.bed," +
                        "/Users/mshamim/Desktop/SLICE.Reboot/existing/GSE63525_GM12878_subcompartments.bed",
                "/Users/mshamim/Desktop/SLICE.Reboot/intra_shuffles/INTACT18B_W4_REAL_LOGBASE_OE_" + 100 + "K",
                "PRIMARY,RH2014"
        };

        System.out.println("-----------------------------------------------------");
        MixerTools.main(strings);
        System.gc();


        strings = new String[]{"intra-shuffle", "-r", 100000 + "", "-k", "SCALE", "-w", "" + 4,
                file14,
                "/Users/mshamim/Desktop/SLICE.Reboot/existing/SLICE_PRIMARY_100K.bed," +
                        "/Users/mshamim/Desktop/SLICE.Reboot/existing/GSE63525_GM12878_subcompartments.bed",
                "/Users/mshamim/Desktop/SLICE.Reboot/intra_shuffles/INTACT18B_W4_SIMPLE_" + 100 + "K",
                "PRIMARY,RH2014"
        };

        System.out.println("-----------------------------------------------------");
        MixerTools.main(strings);
        System.gc();


        String[] strings;
        String changes = "3,A1,34,139,34;" +
                "2,A2,152,251,152;" +
                "5,B1,220,20,60;" +
                "4,B2,255,255,0;" +
                "1,B3,112,128,144";

                strings = new String[]{"rename",
                changes,
                "/Users/mshamim/Desktop/SLICE.Reboot/GMPrimary_Final_100000_v3.17.05/GMPrimary_Final_100000_v3.17.05_5_clusters_GM12878_primary14_30.subcompartment.bed",
                "/Users/mshamim/Desktop/SLICE.Reboot/existing/SLICE_PRIMARY_100K.bed"
        };
        //MixerTools.main(strings);
        System.gc();

        String outfolder = "/Users/mshamim/Desktop/SLICE.Reboot/deepcalls";

        String[] bigFiles = new String[]{
                "/Users/mshamim/Desktop/hicfiles/GM12878_intact_18.7B_8.15.20_30.hic",
                "/Volumes/AidenLabWD7/Backup/AidenLab/LocalFiles/gm12878/GM_2019_mega_30.hic",
                "/Volumes/AidenLabWD7/Backup/AidenLab/LocalFiles/gm12878/GM12878_insitu_combined_15B_30.hic",
                "/Volumes/AidenLabWD7/Backup/AidenLab/LocalFiles/gm12878/GM12878_ultra_42B_1k_30.hic",
                "/Users/mshamim/Desktop/hicfiles/gm12878_rh14_30.hic",
                "/Users/mshamim/Desktop/hicfiles/GM12878_primary14_30.hic",
                "/Volumes/AidenLabWD7/Backup/AidenLab/LocalFiles/gm12878/GM12878_mnase_9.5B_2.11.20_30.hic"
        };

        String[] norms = new String[]{
                "SCALE", "KR", "KR", "KR", "KR", "KR", "KR"
        };

        String[] fileNames = new String[]{
                "GM_Intact_18B", "GM_2019_Mega", "GM12878_15B",
                "Ultra_42B", "RH14_Mega", "RH14_Primary", "GM_MNASE_9.5B"
        };

        for(int q = 5; q < fileNames.length; q++) {
            String folder = fileNames[q] + "_50K_v" + MixerTools.versionNum;
            strings = new String[]{"slice", "-r", "50000", "-k", norms[q],
                    bigFiles[q], "2,13,4",
                    "/Users/mshamim/Desktop/SLICE.Reboot/ULTRAS/" + folder, fileNames[q]
            };
            System.out.println("-----------------------------------------------------");
            //MixerTools.main(strings);
            System.gc();
        }




        /*

        String file14 = "/Users/mshamim/Desktop/hicfiles/gm12878_rh14_30.hic";
        strings = new String[]{"shuffle", "-r", 100000 + "", "-k", "KR", "-w", "" + 40,
                file14,
                "/Users/mshamim/Desktop/SLICE.Reboot/existing/SLICE_PRIMARY_100K.bed",
                "/Users/mshamim/Desktop/SLICE.Reboot/shuffles/PRIMARY_Balanced7_BED_Mega14_HIC_" + 100 + "K_v",
                "SLICE_PRIMARY"
        };
        System.out.println("-----------------------------------------------------");
        MixerTools.main(strings);
        System.gc();

        /*

        String file14 = "/Users/mshamim/Desktop/hicfiles/imr90_rh14_30.hic";
        String bedfiles = "/Users/mshamim/Desktop/SLICE.Reboot/existing/GSE63525_GM12878_subcompartments.bed," +
                "/Users/mshamim/Desktop/SLICE.Reboot/existing/GM12878_track_hg19.bed," +
                "/Users/mshamim/Desktop/SLICE.Reboot/existing/GM12878_SCI_sub_compartments.bed," +
                "/Users/mshamim/Desktop/SLICE.Reboot/existing/SLICE_GM_100K.bed," +
                "/Users/mshamim/Desktop/SLICE.Reboot/existing/SLICE_GM_100K_no2Pass.bed," +
                "/Users/mshamim/Desktop/SLICE.Reboot/existing/SLICE_PRIMARY_100K.bed";

        file14 = "/Users/mshamim/Desktop/hicfiles/imr90_rh14_30.hic";
        bedfiles = "/Users/mshamim/Desktop/SLICE.Reboot/existing/sniper/IMR90_track_hg19.bed," +
                "/Users/mshamim/Desktop/SLICE.Reboot/existing/SLICE_IMR90_100K.bed";

        strings = new String[]{"umap",
                "-r", "100000", "-k", "KR", "-w", "8",
                file14,
                bedfiles,
                "/Users/mshamim/Desktop/SLICE.Reboot/UMAP_IMR90"
        };
        MixerTools.main(strings);
        System.gc();

        file14 = "/Users/mshamim/Desktop/hicfiles/hap1_30.hic";
        bedfiles = "/Users/mshamim/Desktop/SLICE.Reboot/existing/sniper/HAP1_track_hg19.bed," +
                "/Users/mshamim/Desktop/SLICE.Reboot/existing/SLICE_HAP1_100K.bed";

        strings = new String[]{"umap",
                "-r", "100000", "-k", "KR", "-w", "8",
                file14,
                bedfiles,
                "/Users/mshamim/Desktop/SLICE.Reboot/UMAP_HAP1"
        };
        MixerTools.main(strings);
        System.gc();


        /*

        strings = new String[]{"rename",
                changes,
                "/Users/mshamim/Desktop/SLICE.Reboot/GMPrimary_Final_100000_v3.17.05/GMPrimary_Final_100000_v3.17.05_5_clusters_GM12878_primary14_30.subcompartment.bed",
                "/Users/mshamim/Desktop/SLICE.Reboot/existing/SLICE_PRIMARY_100K.bed"
        };
        //MixerTools.main(strings);
        System.gc();

        strings = new String[]{"rename",
                changes,
                "/Users/mshamim/Desktop/SLICE.Reboot/GMPrimary_Final_50000_v3.17.05/GMPrimary_Final_50000_v3.17.05_5_clusters_GM12878_primary14_30.subcompartment.bed",
                "/Users/mshamim/Desktop/SLICE.Reboot/existing/SLICE_PRIMARY_50K.bed"
        };
        //MixerTools.main(strings);
        System.gc();

        String file14 = "/Users/mshamim/Desktop/hicfiles/gm12878_rh14_30.hic";
        strings = new String[]{"shuffle", "-r", 100000 + "", "-k", "KR", "-w", "" + 16,
                file14,
                "/Users/mshamim/Desktop/SLICE.Reboot/existing/SLICE_PRIMARY_100K.bed",
                "/Users/mshamim/Desktop/SLICE.Reboot/shuffles/PRIMARY_BED_Mega14_HIC_" + 100 + "K_v" + MixerTools.versionNum,
                "SLICE_PRIMARY"
        };
        System.out.println("-----------------------------------------------------");
        MixerTools.main(strings);
        System.gc();



/*

        String primary = "/Users/mshamim/Desktop/hicfiles/GM12878_primary14_30.hic";

        for (int res : new int[]{100000,50000}) { //
            String folder = "GMPrimary_Final_" + res + "_v" + MixerTools.versionNum;
            strings = new String[]{"slice", "-r", res + "", "-k", "KR", //"--verbose",
                    primary, "2,8,10",
                    "/Users/mshamim/Desktop/SLICE.Reboot/" + folder, folder + "_"
            };
            System.out.println("-----------------------------------------------------");
            MixerTools.main(strings);
            System.gc();

            strings = new String[]{"shuffle", "-r", res + "", "-k", "KR", "-w", "" + (16 * (100000 / res)),
                    primary,
                    //refs + "," +
                    "/Users/mshamim/Desktop/SLICE.Reboot/" + folder +
                            "/" + folder + "_5_clusters_GM12878_primary14_30.subcompartment.bed",
                    "/Users/mshamim/Desktop/SLICE.Reboot/shuffle_GM_Primary_" + res + "v" + MixerTools.versionNum,
                    "SLICE3"
            };
            System.out.println("-----------------------------------------------------");
            MixerTools.main(strings);
            System.gc();
        }


        /*
        changes = "4,A1,34,139,34;" +
                "2,A2,152,251,152;" +
                "3,B1,220,20,60;" +
                "5,B2,255,255,0;" +
                "1,B3,112,128,144";
        strings = new String[]{"rename",
                changes,
                "/Users/mshamim/Desktop/SLICE.Reboot/calls/HAP1_res_100000v3.17.03/HAP1_res_100000v3.17.03_5_clusters_hap1_30.subcompartment.bed",
                "/Users/mshamim/Desktop/SLICE.Reboot/existing/SLICE_HAP1_100K.bed"
        };
        System.out.println("-----------------------------------------------------");
        MixerTools.main(strings);
        System.gc();

        String file14 = "/Users/mshamim/Desktop/hicfiles/hap1_30.hic";
        strings = new String[]{"shuffle", "-r", 100000 + "", "-k", "KR", "-w", "" + 16,
                file14,
                "/Users/mshamim/Desktop/SLICE.Reboot/existing/SLICE_HAP1_100K.bed",
                "/Users/mshamim/Desktop/SLICE.Reboot/shuffles/SLICE_HAP1_100K_FIN",
                "SLICE3"
        };
        MixerTools.main(strings);
        System.gc();

        /*

        String file14 = "/Users/mshamim/Desktop/hicfiles/imr90_rh14_30.hic";
        strings = new String[]{"shuffle", "-r", 100000 + "", "-k", "KR", "-w", "" + 16,
                file14,
                "/Users/mshamim/Desktop/SLICE.Reboot/existing/sniper/IMR90_track_hg19.bed",
                "/Users/mshamim/Desktop/SLICE.Reboot/existing/sniper/IMR90_100K_Shuffle",
                "SNIPER"
        };
        MixerTools.main(strings);
        System.gc();

        file14 = "/Users/mshamim/Desktop/hicfiles/hap1_30.hic";
        strings = new String[]{"shuffle", "-r", 100000 + "", "-k", "KR", "-w", "" + 16,
                file14,
                "/Users/mshamim/Desktop/SLICE.Reboot/existing/sniper/HAP1_track_hg19.bed",
                "/Users/mshamim/Desktop/SLICE.Reboot/existing/sniper/HAP1_100K_Shuffle",
                "SNIPER"
        };
        MixerTools.main(strings);
        System.gc();

/*

        strings = new String[]{"rename",
                changes,
                "/Users/mshamim/Desktop/SLICE.Reboot/calls/GM_100000v3.17.01/new_test_res_100000v3.17.01_5_clusters_gm12878_rh14_30.subcompartment.bed",
                "/Users/mshamim/Desktop/SLICE.Reboot/existing/archive/SLICE_GM_100K.bed"
        };
        System.out.println("-----------------------------------------------------");
        MixerTools.main(strings);
        System.gc();

        String file14 = "/Users/mshamim/Desktop/hicfiles/gm12878_rh14_30.hic";
        strings = new String[]{"shuffle", "-r", 100000 + "", "-k", "KR", "-w", "" + 16,
                file14,
                "/Users/mshamim/Desktop/SLICE.Reboot/existing/archive/SLICE_GM_100K.bed",
                "/Users/mshamim/Desktop/SLICE.Reboot/existing/archive/SLICE_GM_100K_FIN",
                "SLICE3"
        };
        MixerTools.main(strings);
        System.gc();

/*



        strings = new String[]{"rename",
                changes,
                "/Users/mshamim/Desktop/SLICE.Reboot/GM_Final_res_100000v3.17.04/GM_Final_res_100000v3.17.04_5_clusters_gm12878_rh14_30.subcompartment.bed",
                "/Users/mshamim/Desktop/SLICE.Reboot/existing/SLICE_GM_100K.bed"
        };
        System.out.println("-----------------------------------------------------");
        MixerTools.main(strings);
        System.gc();

        String file14 = "/Users/mshamim/Desktop/hicfiles/gm12878_rh14_30.hic";
        strings = new String[]{"shuffle", "-r", 100000 + "", "-k", "KR", "-w", "" + 16,
                file14,
                "/Users/mshamim/Desktop/SLICE.Reboot/existing/SLICE_GM_100K.bed",
                "/Users/mshamim/Desktop/SLICE.Reboot/shuffles/SLICE_GM_100K_FIN2",
                "SLICE3"
        };
        MixerTools.main(strings);
        System.gc();



        strings = new String[]{"rename",
                changes,
                "/Users/mshamim/Desktop/SLICE.Reboot/GM_Final_res_50000v3.17.04/GM_Final_res_50000v3.17.04_5_clusters_gm12878_rh14_30.subcompartment.bed",
                "/Users/mshamim/Desktop/SLICE.Reboot/existing/SLICE_GM_50K.bed"
        };
        System.out.println("-----------------------------------------------------");
        MixerTools.main(strings);
        System.gc();

        strings = new String[]{"shuffle", "-r", 50000 + "", "-k", "KR", "-w", "" + 32,
                file14,
                "/Users/mshamim/Desktop/SLICE.Reboot/existing/SLICE_GM_50K.bed," +
                        "/Users/mshamim/Desktop/SLICE.Reboot/existing/SLICE_GM_100K.bed," +
                        "/Users/mshamim/Desktop/SLICE.Reboot/existing/GSE63525_GM12878_subcompartments.bed",
                "/Users/mshamim/Desktop/SLICE.Reboot/shuffles/SLICE_GM_50K_FIN2",
                "SLICE50,SLICE100,RH2014"
        };
        MixerTools.main(strings);
        System.gc();

        /*

        String refs = "/Users/mshamim/Desktop/SLICE.Reboot/existing/GSE63525_GM12878_subcompartments.bed," +
                "/Users/mshamim/Desktop/SLICE.Reboot/existing/GM12878_SCI_sub_compartments.bed," +
                "/Users/mshamim/Desktop/SLICE.Reboot/existing/GM12878_track_hg19.bed";
        file14 = "/Users/mshamim/Desktop/hicfiles/gm12878_rh14_30.hic";

        //file14 = "/Users/mshamim/Desktop/hicfiles/hap1_30.hic";

        for (int res : new int[]{50000}) { //
            String folder = "GM_Final_res_" + res + "v" + MixerTools.versionNum;
            strings = new String[]{"slice", "-r", res + "", "-k", "KR", //"--verbose",
                    file14, "2,8,10",
                    "/Users/mshamim/Desktop/SLICE.Reboot/" + folder, folder + "_"
            };
            System.out.println("-----------------------------------------------------");
            MixerTools.main(strings);
            System.gc();

            strings = new String[]{"shuffle", "-r", res + "", "-k", "KR", "-w", "" + (16 * (100000 / res)), file14,
                    //refs + "," +
                    "/Users/mshamim/Desktop/SLICE.Reboot/" + folder +
                            "/" + folder + "_5_clusters_gm12878_rh14_30.subcompartment.bed",
                    "/Users/mshamim/Desktop/SLICE.Reboot/shuffle_GM_Final_res_" + res + "v" + MixerTools.versionNum,
                    "SLICE3"
            };
            System.out.println("-----------------------------------------------------");
            MixerTools.main(strings);
            System.gc();
        }

        //ExtractionB4.extract();

       */
    }
}
