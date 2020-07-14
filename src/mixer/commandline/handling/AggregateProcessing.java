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

package mixer.commandline.handling;

import mixer.commandline.MixerTools;
import mixer.commandline.utils.drink.ExtractingOEDataUtils;
import mixer.commandline.utils.drink.GenomewideBadIndexFinder;
import mixer.commandline.utils.drink.SliceMatrix;


/**
 * Created for testing multiple CLTs at once
 * Basically scratch space
 */
public class AggregateProcessing {


    public static ExtractingOEDataUtils.ThresholdType beforeThresholdType = ExtractingOEDataUtils.ThresholdType.TRUE_OE;
    public static ExtractingOEDataUtils.ThresholdType afterThresholdType = ExtractingOEDataUtils.ThresholdType.TRUE_OE;
    public static boolean useDerivative = false;
    public static boolean useL1Norm = false;

    public static void main(String[] argv) throws Exception {
    
    
        String[] strings;
        String refs = "/Users/muhammad/Desktop/drinks/existingmethods/GSE63525_GM12878_subcompartments.bed" + "+" +
                "/Users/muhammad/Desktop/drinks/existingmethods/ultra_res_100k_default_clean_outliers.bed" + "+" +
                "/Users/muhammad/Desktop/drinks/existingmethods/New_STRICT_gold_standard.bed";
        String file42 = "/Volumes/AidenLabWD7/Backup/AidenLab/LocalFiles/gm12878/GM12878_ultra_42B_1k_30.hic";
        String file14 = "/Users/muhammad/Desktop/insitumboi/combined_GM12878_insitu_combined_30.hic";
        String folder;
    
        String file1B = "/Users/muhammad/Desktop/local_hic_files/sub1bil.hic";
    
        for (int w = 3; w < 3.1; w++) {
            SliceMatrix.USE_DERIV = true;
            folder = "redo_sparseV_slice_14w" + w + "_deriv";
            strings = new String[]{"slice", "-r", "100000", "-w", "" + w, //"--verbose",
                    file14, "/Users/muhammad/Desktop/slice/" + folder, folder + "_", refs};
            System.out.println("-----------------------------------------------------");
            //MixerTools.main(strings);
            System.gc();
        
            SliceMatrix.USE_DERIV = false;
            folder = "redo_sparseV_slice_14w" + w + "";
            strings = new String[]{"slice", "-r", "100000", "-w", "" + w, //"--verbose",
                    file14, "/Users/muhammad/Desktop/slice/" + folder, folder + "_", refs};
            System.out.println("-----------------------------------------------------");
            //MixerTools.main(strings);
            System.gc();
        }
    
        for (float dd = -3f; dd < -2.9f; dd += 0.5f) {
            SliceMatrix.USE_DERIV = false;
            GenomewideBadIndexFinder.ZSCORE_THRESHOLD_HIGHER = dd;
            folder = "redo_slice_14_d_" + dd + "_SPARSE";
            strings = new String[]{"slice", "-r", "100000", "-w", "4", //"--verbose",
                    file14, "/Users/muhammad/Desktop/slice/" + folder, folder + "_", refs};
            System.out.println("----------------------------------------------------");
            MixerTools.main(strings);
            System.gc();
        }
        
        /*
        
        SliceMatrix.USE_DERIV = true;
        folder = "redo_slice_14_july12_deriv";
        strings = new String[]{"slice", "-r", "100000", "-w", "2", "--verbose",
                file14, "/Users/muhammad/Desktop/slice/" + folder, folder + "_", refs};
        System.out.println("-----------------------------------------------------");
        MixerTools.main(strings);
        System.gc();
        
        SliceMatrix.USE_DERIV = false;
        folder = "redo_slice_14_july12";
        strings = new String[]{"slice", "-r", "100000", "-w", "2", "--verbose",
                file14, "/Users/muhammad/Desktop/slice/" + folder, folder + "_", refs};
        System.out.println("-----------------------------------------------------");
        MixerTools.main(strings);
        System.gc();
        
        
        SliceMatrix.USE_DERIV = false;
        folder = "redo_slice_42_july12";
        strings = new String[]{"slice", "-r", "100000", "-w", "2", "--verbose",
                file42, "/Users/muhammad/Desktop/slice/" + folder, folder + "_", refs};
        System.out.println("-----------------------------------------------------");
        MixerTools.main(strings);
        System.gc();
    
        
        SliceMatrix.USE_DERIV = true;
        folder = "redo_slice_42_july12_deriv";
        strings = new String[]{"slice", "-r", "100000", "-w", "2", "--verbose",
                file42, "/Users/muhammad/Desktop/slice/" + folder, folder + "_", refs};
        System.out.println("-----------------------------------------------------");
        MixerTools.main(strings);
        System.gc();
        
/*
        int num = 10;//50;

        String folder = "digest_" + num + "_42_july7";
        String prefix = folder + "_";
        MatrixCleanupReduction.USE_DERIV = false;
        MatrixCleanupReduction.DERIV_THEN_ZSCORE = false;
        MatrixCleanupReduction.USE_RANDOM3_TYPE = true;
        MatrixCleanupReduction.reductionScalar = num;

        strings = new String[]{"links", "-r", "100000", "-w", "1",
                "--verbose",
                "/Volumes/AidenLabWD7/Backup/AidenLab/LocalFiles/gm12878/GM12878_ultra_42B_1k_30.hic",
                //"/Users/muhammad/Desktop/insitumboi/combined_GM12878_insitu_combined_30.hic",
                "/Users/muhammad/Desktop/drinks/hundred/" + folder, prefix,
                "/Users/muhammad/Desktop/drinks/existingmethods/GSE63525_GM12878_subcompartments.bed" + "+" +
                        "/Users/muhammad/Desktop/drinks/existingmethods/ultra_res_100k_default_clean_outliers.bed" + "+" +
                        "/Users/muhammad/Desktop/drinks/existingmethods/New_STRICT_gold_standard.bed"};

        System.out.println("-----------------------------------------------------");
        MixerTools.main(strings);

        folder = "subc50k_v_deriv_" + num + "_42";
        prefix = folder + "_";
        MatrixCleanupReduction.USE_DERIV = true;
        MatrixCleanupReduction.DERIV_THEN_ZSCORE = false;
        MatrixCleanupReduction.USE_RANDOM3_TYPE = true;
        MatrixCleanupReduction.reductionScalar = num;

        strings = new String[]{"links", "-r", "50000", "-w", "2",
                "--verbose",
                "/Volumes/AidenLabWD7/Backup/AidenLab/LocalFiles/gm12878/GM12878_ultra_42B_1k_30.hic",
                //"/Users/muhammad/Desktop/insitumboi/combined_GM12878_insitu_combined_30.hic",
                "/Users/muhammad/Desktop/drinks/fifty/" + folder, prefix};

        System.out.println("-----------------------------------------------------");
        //MixerTools.main(strings);

        /*


            folder = "rebootderiv_"+num+"_42";
            prefix = folder + "_";
            MatrixCleanupReduction.USE_DERIV = true;
            MatrixCleanupReduction.DERIV_ON_ZSCORE = false;
        MatrixCleanupReduction.DERIV_THEN_ZSCORE = false;
            MatrixCleanupReduction.USE_RANDOM3_TYPE = true;
            MatrixCleanupReduction.reductionScalar = num;

            strings = new String[]{"links", "-r", "100000", "-w", "2",
                    //"--verbose",
                    "/Volumes/AidenLabWD7/Backup/AidenLab/LocalFiles/gm12878/GM12878_ultra_42B_1k_30.hic",
                    //"/Users/muhammad/Desktop/insitumboi/combined_GM12878_insitu_combined_30.hic",
                    "/Users/muhammad/Desktop/drinks/pfol/" + folder, prefix};

            System.out.println("-----------------------------------------------------");
            //MixerTools.main(strings);


            folder = "rebootderivONZ_"+num+"_42";
            prefix = folder + "_";
            MatrixCleanupReduction.USE_DERIV = true;
            MatrixCleanupReduction.DERIV_ON_ZSCORE = true;
            MatrixCleanupReduction.USE_RANDOM3_TYPE = true;
        MatrixCleanupReduction.DERIV_THEN_ZSCORE = false;
            MatrixCleanupReduction.reductionScalar = num;

            strings = new String[]{"links", "-r", "100000", "-w", "2",
                    //"--verbose",
                    "/Volumes/AidenLabWD7/Backup/AidenLab/LocalFiles/gm12878/GM12878_ultra_42B_1k_30.hic",
                    //"/Users/muhammad/Desktop/insitumboi/combined_GM12878_insitu_combined_30.hic",
                    "/Users/muhammad/Desktop/drinks/pfol/" + folder, prefix};

            System.out.println("-----------------------------------------------------");
            //MixerTools.main(strings);


        folder = "reboot_reset_w3_"+num+"_14";
        prefix = folder + "_";
        MatrixCleanupReduction.USE_DERIV = false;
        MatrixCleanupReduction.DERIV_ON_ZSCORE = false;
        MatrixCleanupReduction.USE_RANDOM3_TYPE = true;
        MatrixCleanupReduction.DERIV_THEN_ZSCORE = true;
        MatrixCleanupReduction.reductionScalar = num;

        strings = new String[]{"links", "-r", "100000", "-w", "3",
                //"--verbose",
                //"/Volumes/AidenLabWD7/Backup/AidenLab/LocalFiles/gm12878/GM12878_ultra_42B_1k_30.hic",
                "/Users/muhammad/Desktop/insitumboi/combined_GM12878_insitu_combined_30.hic",
                "/Users/muhammad/Desktop/drinks/pfol/" + folder, prefix};

        System.out.println("-----------------------------------------------------");
        //MixerTools.main(strings);

        folder = "reboot_smooth_deriv_"+num+"_42";
        prefix = folder + "_";
        MatrixCleanupReduction.USE_DERIV = true;
        MatrixCleanupReduction.DERIV_ON_ZSCORE = false;
        MatrixCleanupReduction.USE_RANDOM3_TYPE = true;
        MatrixCleanupReduction.DERIV_THEN_ZSCORE = true;
        MatrixCleanupReduction.reductionScalar = num;

        strings = new String[]{"links", "-r", "100000", "-w", "2",
                //"--verbose",
                "/Volumes/AidenLabWD7/Backup/AidenLab/LocalFiles/gm12878/GM12878_ultra_42B_1k_30.hic",
                //"/Users/muhammad/Desktop/insitumboi/combined_GM12878_insitu_combined_30.hic",
                "/Users/muhammad/Desktop/drinks/pfol/" + folder, prefix};

        System.out.println("-----------------------------------------------------");
        //MixerTools.main(strings);


        /*

        String folder = "plinks3_42";
        String prefix = folder + "_";
        LinksMatrix.USE_DERIV = false;
        LinksMatrix.USE_RANDOM3_TYPE = false;

        strings = new String[]{"links", "-r", "100000", "-w", "2",
                "--verbose",
                "/Volumes/AidenLabWD7/Backup/AidenLab/LocalFiles/gm12878/GM12878_ultra_42B_1k_30.hic",
                //"/Users/muhammad/Desktop/insitumboi/combined_GM12878_insitu_combined_30.hic",
                "/Users/muhammad/Desktop/drinks/" + folder, prefix};

        System.out.println("-----------------------------------------------------");
        MixerTools.main(strings);



        System.gc();
        LinksMatrix.USE_DERIV = true;
        LinksMatrix.USE_RANDOM3_TYPE = false;

        folder = "plinks3_42_deriv";
        prefix = folder + "_";

        strings = new String[]{"links", "-r", "100000", "-w", "2",
                "--verbose",
                "/Volumes/AidenLabWD7/Backup/AidenLab/LocalFiles/gm12878/GM12878_ultra_42B_1k_30.hic",
                //"/Users/muhammad/Desktop/insitumboi/combined_GM12878_insitu_combined_30.hic",
                "/Users/muhammad/Desktop/drinks/" + folder, prefix};

        System.out.println("-----------------------------------------------------");
        //MixerTools.main(strings);

        LinksMatrix.USE_DERIV = false;
        LinksMatrix.USE_RANDOM3_TYPE = true;

        folder = "plinks3_42_T3";
        prefix = folder + "_";

        strings = new String[]{"links", "-r", "100000", "-w", "2",
                "--verbose",
                //"/Users/muhammad/Desktop/insitumboi/combined_GM12878_insitu_combined_30.hic",
                "/Volumes/AidenLabWD7/Backup/AidenLab/LocalFiles/gm12878/GM12878_ultra_42B_1k_30.hic",
                "/Users/muhammad/Desktop/drinks/" + folder, prefix};

        System.out.println("-----------------------------------------------------");
        MixerTools.main(strings);

        LinksMatrix.USE_DERIV = true;
        LinksMatrix.USE_RANDOM3_TYPE = true;

        folder = "plinks3_42_deriv_T3";
        prefix = folder + "_";

        strings = new String[]{"links", "-r", "100000", "-w", "2",
                "--verbose",
                //"/Users/muhammad/Desktop/insitumboi/combined_GM12878_insitu_combined_30.hic",
                "/Volumes/AidenLabWD7/Backup/AidenLab/LocalFiles/gm12878/GM12878_ultra_42B_1k_30.hic",
                "/Users/muhammad/Desktop/drinks/" + folder, prefix};

        System.out.println("-----------------------------------------------------");
        //MixerTools.main(strings);

/*
        String folder = "emerrald42_w2_all";
        String prefix = folder+"_";

         strings = new String[]{"links", "-r", "100000", "-w", "2", // 1
         //"--verbose",
         //"-c", "2,3,4,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22",
         "/Volumes/AidenLabWD7/Backup/AidenLab/LocalFiles/gm12878/GM12878_ultra_42B_1k_30.hic",
         //    "/Volumes/AidenLabWD7/Backup/AidenLab/LocalFiles/gm12878/GM12878_ultra_42B_1k.hic",
         "/Users/muhammad/Desktop/drinks/" + folder, prefix};

         System.out.println("-----------------------------------------------------");
         MixerTools.main(strings);

        System.gc();


        folder = "emerrald42_w2_50k_all";
        prefix = folder+"_";

        strings = new String[]{"links", "-r", "50000", "-w", "2", // 1
                //"--verbose",
                //"-c", "2,3,4,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22",
                "/Volumes/AidenLabWD7/Backup/AidenLab/LocalFiles/gm12878/GM12878_ultra_42B_1k_30.hic",
                //    "/Volumes/AidenLabWD7/Backup/AidenLab/LocalFiles/gm12878/GM12878_ultra_42B_1k.hic",
                "/Users/muhammad/Desktop/drinks/" + folder, prefix};

        System.out.println("-----------------------------------------------------");
        MixerTools.main(strings);

        System.gc();



/*

        folder = "ruby42_for_1";
        prefix = folder + "_";

        strings = new String[]{"links", "-r", "100000", "-w", "1",
                "--verbose",
                "-c", "1,2,4,6,7,8,9,10,11,12,13,15,16,17,18,19,20,21,22",
                "/Volumes/AidenLabWD7/Backup/AidenLab/LocalFiles/gm12878/GM12878_ultra_42B_1k_30.hic",
                "/Users/muhammad/Desktop/drinks/ultra42/" + folder, prefix};

        System.out.println("-----------------------------------------------------");
        MixerTools.main(strings);


        System.gc();

        folder = "ruby42_for_5";
        prefix = folder + "_";

        strings = new String[]{"links", "-r", "100000", "-w", "1",
                "--verbose",
                "-c", "2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22",
                "/Volumes/AidenLabWD7/Backup/AidenLab/LocalFiles/gm12878/GM12878_ultra_42B_1k_30.hic",
                "/Users/muhammad/Desktop/drinks/ultra42/" + folder, prefix};

        System.out.println("-----------------------------------------------------");
        MixerTools.main(strings);

        System.gc();


        folder = "ruby42_all";
        prefix = folder + "_";

        strings = new String[]{"links", "-r", "100000", "-w", "1",
                "--verbose",
                "/Volumes/AidenLabWD7/Backup/AidenLab/LocalFiles/gm12878/GM12878_ultra_42B_1k_30.hic",
                "/Users/muhammad/Desktop/drinks/ultra42/" + folder, prefix};

        System.out.println("-----------------------------------------------------");
        MixerTools.main(strings);

        System.gc();

        */

        /*
        DrinkUtils.createUnanimousList(
                "/Users/muhammad/Desktop/drinks/existingmethods/GSE63525_GM12878_subcompartments.bed",
                "/Users/muhammad/Desktop/drinks/existingmethods/SNIPER_GM12878_track_hg19.bed",
                        "/Users/muhammad/Desktop/drinks/existingmethods/SCI_GM12878_SCI_sub_compartments.bed");
 */

    }
}
