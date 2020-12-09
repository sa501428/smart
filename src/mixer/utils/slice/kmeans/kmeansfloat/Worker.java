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

package mixer.utils.slice.kmeans.kmeansfloat;

import mixer.utils.similaritymeasures.RobustEuclideanDistance;

import java.util.concurrent.BrokenBarrierException;

/**
 * The class which does the hard work of the subtasks.
 */
public class Worker implements Runnable {

    static final int DOING_NOTHING = 0;
    static final int COMPUTING_DISTANCES = 1;
    static final int MAKING_ASSIGNMENTS = 2;
    // Codes used to identify what step is being done.
    // What the object is currently doing
    public static int mDoing = Worker.DOING_NOTHING;
    // Defines range of coordinates to cover.
    private final int mStartCoord;
    private final int mNumCoords;
    // Number of moves made by this worker in the last call
    // to workerMakeAssignments().  The SubtaskManager totals up
    // this value from all the workers in numberOfMoves().
    private int mMoves;

    /**
     * Constructor
     *
     * @param startCoord index of the first coordinate covered by
     *                   this Worker.
     * @param numCoords  the number of coordinates covered.
     */
    Worker(int startCoord, int numCoords) {
        mStartCoord = startCoord;
        mNumCoords = numCoords;
    }

    /**
     * Returns the number of moves this worker made in the last
     * execution of workerMakeAssignments()
     */
    int numberOfMoves() {
        return mMoves;
    }

    /**
     * The run method.  It accesses the SubtaskManager field mDoing
     * to determine what subtask to perform.
     */
    public void run() {
        try {
            switch (mDoing) {
                case COMPUTING_DISTANCES:
                    workerComputeDistances(ConcurrentKMeans.mProtoClusters);
                    break;
                case MAKING_ASSIGNMENTS:
                    workerMakeAssignments();
                    break;
            }
        } finally {
            // If there's a barrier, call its await() method.  To ensure it
            // gets done, it's placed in the finally clause.
            if (SubtaskManager.mBarrier != null) {
                try {
                    SubtaskManager.mBarrier.await();
                    // barrier.isBroken() will return true if either of these
                    // exceptions happens, so the SubtaskManager will detect
                    // the problem.
                } catch (InterruptedException | BrokenBarrierException ignored) {
                }
            }
        }

    }

    /**
     * Compute the distances for the covered coordinates
     * to the updated centers.
     */
    private void workerComputeDistances(ProtoCluster[] mProtoClusters) {
        int lim = mStartCoord + mNumCoords;
        for (int i = mStartCoord; i < lim; i++) {
            int numClusters = mProtoClusters.length;
            for (int c = 0; c < numClusters; c++) {
                ProtoCluster cluster = mProtoClusters[c];
                if (cluster.getConsiderForAssignment() && cluster.needsUpdate()) {
                    ConcurrentKMeans.mDistanceCache[i][c] = distanceL2Norm(ConcurrentKMeans.mCoordinates[i], cluster.getCenter());
                }
            }
        }
    }

    /**
     * Assign each covered coordinate to the nearest cluster.
     */
    private void workerMakeAssignments() {
        mMoves = 0;
        int lim = mStartCoord + mNumCoords;
        for (int i = mStartCoord; i < lim; i++) {
            int c = nearestCluster(i);
            ConcurrentKMeans.mProtoClusters[c].add(i);
            if (ConcurrentKMeans.mClusterAssignments[i] != c) {
                ConcurrentKMeans.mClusterAssignments[i] = c;
                mMoves++;
            }
        }
    }

    /**
     * Compute the euclidean distance between the two arguments.
     */
    private float distanceL2Norm(float[] coord, float[] center) {
        if (ConcurrentKMeans.useNonNanVersion) {
            return RobustEuclideanDistance.SINGLETON.distance(coord, center);
        }
        double sumSquared = 0.0;
        for (int i = 0; i < coord.length; i++) {
            float v = coord[i] - center[i];
            sumSquared += (v * v);
        }
        return (float) Math.sqrt(sumSquared);
    }

    /**
     * Find the nearest cluster to the coordinate identified by
     * the specified index.
     */
    private int nearestCluster(int ndx) {
        int nearest = -1;
        double min = Double.MAX_VALUE;
        int numClusters = ConcurrentKMeans.mProtoClusters.length;
        for (int c = 0; c < numClusters; c++) {
            if (ConcurrentKMeans.mProtoClusters[c].getConsiderForAssignment()) {
                double d = ConcurrentKMeans.mDistanceCache[ndx][c];
                if (d < min) {
                    min = d;
                    nearest = c;
                }
            }
        }
        return nearest;
    }

}