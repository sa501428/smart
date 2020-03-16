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



/**
 * Created for testing multiple CLTs at once
 * Basically scratch space
 */
class AggregateProcessing {


    public static void main(String[] argv) throws Exception {


        String[] strings = new String[]{"distort", "-r", "250000", "--stride", "50",
                "-c", "1,2,3,4,5,7,8,9,10,12,13,14,15,16,17,18,19,20,21,22,X", //6,11
                "/Volumes/AidenLabWD7/Backup/AidenLab/LocalFiles/gm12878/GSE63525_GM12878_insitu_primary_30.hic+" +
                        "/Volumes/AidenLabWD7/Backup/AidenLab/LocalFiles/gm12878/GSE63525_GM12878_insitu_replicate_30.hic",
                "200,20,100", "/Users/muhammad/Desktop/findsv/train_set_2_250kb_m200_full_minus_6_11"};

        MixerTools.main(strings);

        // load the model

        /*
        String simpleMlp = "/Users/muhammad/Desktop/deeplearning/models/Clean64DistortionDiffHalfLocalizerV0BinCross.h5";
        MultiLayerNetwork model = KerasModelImport.importKerasSequentialModelAndWeights(simpleMlp);


        // make a random sample
        int inputs = 10;
        INDArray features = Nd4j.zeros(inputs);
        for (int i=0; i<inputs; i++) {
            features.putScalar(new int[]{i}, Math.random() < 0.5 ? 0 : 1);
        }
// get the prediction
        //double prediction = model.output(features).getDouble(0);

         */



    }
}
