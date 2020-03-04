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

package mixer.track.feature;

import java.awt.*;
import java.util.List;

/**
 * Created by muhammadsaadshamim on 10/27/15.
 */
public class Feature2DTools {

    public static Feature2DList subtract(final Feature2DList listA, final Feature2DList listB) {
        Feature2DList result = new Feature2DList(listA);
        result.filterLists(new FeatureFilter() {
            @Override
            public List<Feature2D> filter(String chr, List<Feature2D> feature2DList) {
                if (listB.containsKey(chr)) {
                    feature2DList.removeAll(listB.getFeatureList(chr));
                }
                return feature2DList;
            }
        });
        result.removeDuplicates();
        return result;
    }


    public static boolean loopIsUpstreamOfDomain(Feature2D loop, Feature2D domain, int threshold) {
        return loop.getEnd1() < domain.getStart1() - threshold &&
                loop.getEnd2() < domain.getStart2() - threshold;
    }

    public static boolean loopIsDownstreamOfDomain(Feature2D loop, Feature2D domain, int threshold) {
        return loop.getStart1() > domain.getEnd1() + threshold &&
                loop.getStart2() > domain.getEnd2() + threshold;
    }

    public static boolean domainContainsLoopWithinExpandedTolerance(Feature2D loop, Feature2D domain, int threshold) {

        Rectangle bounds = new Rectangle(domain.getStart1() - threshold, domain.getStart2() - threshold,
                domain.getWidth1() + 2 * threshold, domain.getWidth2() + 2 * threshold);
        Point point = new Point(loop.getMidPt1(), loop.getMidPt2());

        return bounds.contains(point);
    }

    /**
     * Compares a feature against all other features in list
     *
     * @param feature
     * @param existingFeatures
     * @return
     */
    public static boolean doesOverlap(Feature2D feature, List<Feature2D> existingFeatures) {
        boolean repeat = false;
        for (Feature2D existingFeature : existingFeatures) {
            if (existingFeature.overlapsWith(feature)) {
                repeat = true;
            }
        }
        return repeat;
    }

    public static boolean isResolutionPresent(final Feature2DList feature2DList, final int resolution) {
        final boolean[] returnValue = new boolean[1];
        returnValue[0] = false;
        feature2DList.processLists(new FeatureFunction() {
            @Override
            public void process(String chr, List<Feature2D> feature2DList) {
                for (Feature2D feature : feature2DList) {
                    if (feature.getWidth1() == resolution) {
                        returnValue[0] = true;
                        return;
                    }
                }
            }
        });
        return returnValue[0];
    }
}
