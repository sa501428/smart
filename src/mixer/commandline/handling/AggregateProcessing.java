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


/**
 * Created for testing multiple CLTs at once
 * Basically scratch space
 */
public class AggregateProcessing {


    public static ExtractingOEDataUtils.ThresholdType beforeThresholdType = ExtractingOEDataUtils.ThresholdType.TRUE_OE;
    public static ExtractingOEDataUtils.ThresholdType afterThresholdType = ExtractingOEDataUtils.ThresholdType.TRUE_OE;
    public static boolean useDerivative = false;
    public static boolean useL1Norm = false;

    public static double scalar = 0;


    public static void main(String[] argv) throws Exception {


        String[] strings = new String[]{"distort", "-r", "250000", "--stride", "50",
                "-c", "1,2,3,4,5,7,8,9,10,12,13,14,15,16,17,18,19,20,21,22,X", //6,11
                "/Volumes/AidenLabWD7/Backup/AidenLab/LocalFiles/gm12878/GSE63525_GM12878_insitu_primary_30.hic+" +
                        "/Volumes/AidenLabWD7/Backup/AidenLab/LocalFiles/gm12878/GSE63525_GM12878_insitu_replicate_30.hic",
                "200,20,100", "/Users/muhammad/Desktop/findsv/train_set_2_250kb_m200_full_minus_6_11"};


        /*
        +/- OE vs logeo   - rreal oe
        +/- deriv              no
        +/- L1                 no
        log vs linear vs real before  real before
        log vs linear vs real after   real after


        +/- clean diagonal
        new tests
        more initial clusters for intra
        try log(o+1/e+1) for deterministic splitter





        Lessons learned

don't use L1 norm before/after
don't use LOGEO before/after (REAL_OE+1 woks better)
don't use derivative before/after
log vs linear vs real before  - real works
log vs linear vs real after  - rela works
more initial clusters for intra - unequal verrsion doesn't work

***** zebrra offset helped a little +

try log(o+1/e+1) for deterministic splitter





(length/1,000,000  -- failed


post cleanuper?




try log rround?

         */


        String folder = "canyon_0_reboot_mega42_again_with_augment3";//""one_round_log_2r";
        String prefix = folder + "_";

        strings = new String[]{"drinks", "-r", "100000", "-w", "3", "--verbose",
                //"/Users/muhammad/Desktop/insitumboi/combined_GM12878_insitu_combined_30.hic"
                "/Volumes/AidenLabWD7/Backup/AidenLab/LocalFiles/gm12878/GM12878_ultra_42B_1k_30.hic"
                , "/Users/muhammad/Desktop/drinks/" + folder, prefix};
        System.out.println("-----------------------------------------------------");
        MixerTools.main(strings);

        System.gc();


        /*


        DrinkUtils.createUnanimousList(
                "/Users/muhammad/Desktop/drinks/existingmethods/GSE63525_GM12878_subcompartments.bed",
                "/Users/muhammad/Desktop/drinks/existingmethods/SNIPER_GM12878_track_hg19.bed",
                        "/Users/muhammad/Desktop/drinks/existingmethods/SCI_GM12878_SCI_sub_compartments.bed");

 */

    }
}
