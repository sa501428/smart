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

package mixer.utils.slice.kmeans.kmeansfloat;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;
import java.util.Random;

/**
 * The version of K-means clustering adapted for true concurrency
 * or simultaneous multithreading (SMT).  The subtasks of
 * computing distances and making assignments are delegate to
 * a subtask manager which oversees a thread pool.
 */
@SuppressWarnings("ForLoopReplaceableByForEach")
public class ConcurrentKMeans implements KMeans {

    // 2D array holding the coordinates to be clustered.
    public static float[][] mCoordinates;
    public static boolean useNonNanVersion = false;
    // The desired number of clusters and maximum number
    // of iterations.
    private final int mK;
    private final int mMaxIterations;
    // Seed for the random number generator used to select
    // coordinates for the initial cluster centers.
    private final long mRandomSeed;
    // The number of threads used to perform the subtasks.
    private final int mThreadCount;
    // Listeners to be notified of significant happenings.
    private final List<KMeansListener> mListeners = new ArrayList<>(1);
    // Temporary clusters used during the clustering process.  Converted to
    // an array of the simpler class Cluster at the conclusion.
    public static ProtoCluster[] mProtoClusters;
    // Cache of coordinate-to-cluster distances. Number of entries =
    // number of clusters X number of coordinates.
    public static float[][] mDistanceCache;
    // Used in makeAssignments() to figure out how many moves are made
    // during each iteration -- the cluster assignment for coordinate n is
    // found in mClusterAssignments[n] where the N coordinates are numbered
    // 0 ... (N-1)
    public static int[] mClusterAssignments;
    // Subtask manager that handles the thread pool to which
    // time-consuming tasks are delegated.
    private SubtaskManager mSubtaskManager;
    // An array of Cluster objects: the output of k-means.
    private Cluster[] mClusters;

    /**
     * Constructor
     *
     * @param coordinates   two-dimensional array containing the coordinates to be clustered.
     * @param k             the number of desired clusters.
     * @param maxIterations the maximum number of clustering iterations.
     * @param randomSeed    seed used with the random number generator.
     * @param threadCount   the number of threads to be used for computing time-consuming steps.
     */
    private ConcurrentKMeans(float[][] coordinates, int k, int maxIterations,
                             long randomSeed, int threadCount) {
        mCoordinates = coordinates;
        // Can't have more clusters than coordinates.
        mK = Math.min(k, mCoordinates.length);
        mMaxIterations = maxIterations;
        mRandomSeed = randomSeed;
        mThreadCount = threadCount;
    }

    /**
     * Constructor that uses the return from
     * <tt>Runtime.getRuntime().availableProcessors()</tt> as the number
     * of threads for time-consuming steps.
     *
     * @param coordinates   two-dimensional array containing the coordinates to be clustered.
     * @param k             the number of desired clusters.
     * @param maxIterations the maximum number of clustering iterations.
     * @param randomSeed    seed used with the random number generator.
     */
    public ConcurrentKMeans(float[][] coordinates, int k, int maxIterations,
                            long randomSeed) {
        this(coordinates, k, maxIterations, randomSeed,
                Runtime.getRuntime().availableProcessors());
    }

    /**
     * Adds a KMeansListener to be notified of significant happenings.
     *
     * @param l the listener to be added.
     */
    public void addKMeansListener(KMeansListener l) {
        synchronized (mListeners) {
            if (!mListeners.contains(l)) {
                mListeners.add(l);
            }
        }
    }

    /**
     * Removes a KMeansListener
     *
     * @param l the listener to be removed.
     */
    public void removeKMeansListener(KMeansListener l) {
        synchronized (mListeners) {
            mListeners.remove(l);
        }
    }

    /**
     * Posts a message to registered KMeansListeners.
     *
     * @param message
     */
    private void postKMeansMessage(String message) {
        if (mListeners.size() > 0) {
            synchronized (mListeners) {
                for (KMeansListener mListener : mListeners) {
                    mListener.kmeansMessage(message);
                }
            }
        }
    }

    /**
     * Notifies registered listeners that k-means is complete.
     *
     * @param clusters      the output of clustering.
     * @param executionTime the number of milliseconds taken to cluster.
     */
    private void postKMeansComplete(Cluster[] clusters, long executionTime) {
        if (mListeners.size() > 0) {
            synchronized (mListeners) {
                for (KMeansListener mListener : mListeners) {
                    mListener.kmeansComplete(clusters, executionTime);
                }
            }
        }
    }

    /**
     * Notifies registered listeners that k-means has failed because of
     * a Throwable caught in the run method.
     *
     * @param err
     */
    private void postKMeansError(Throwable err) {
        if (mListeners.size() > 0) {
            synchronized (mListeners) {
                for (KMeansListener mListener : mListeners) {
                    mListener.kmeansError(err);
                }
            }
        }
    }

    /**
     * Get the clusters computed by the algorithm.  This method should
     * not be called until clustering has completed successfully.
     *
     * @return an array of Cluster objects.
     */
    public Cluster[] getClusters() {
        return mClusters;
    }

    /**
     * Run the clustering algorithm.
     */
    public void run() {

        try {

            // Note the start time.
            long startTime = System.currentTimeMillis();

            postKMeansMessage("K-Means clustering started");

            // Randomly initialize the cluster centers creating the
            // array mProtoClusters.
            initCenters();
            postKMeansMessage("... centers initialized");

            // Instantiate the subtask manager.
            mSubtaskManager = new SubtaskManager(mThreadCount);

            // Post a message about the state of concurrent subprocessing.
            if (mThreadCount > 1) {
                postKMeansMessage("... concurrent processing mode with "
                        + mThreadCount + " subtask threads");
            } else {
                postKMeansMessage("... non-concurrent processing mode");
            }

            // Perform the initial computation of distances.
            computeDistances();

            // Make the initial cluster assignments.
            makeAssignments();

            // Number of moves in the iteration and the iteration counter.
            int moves, it = 0;

            // Main Loop:
            //
            // Two stopping criteria:
            // - no moves in makeAssignments
            //   (moves == 0)
            // OR
            // - the maximum number of iterations has been reached
            //   (it == mMaxIterations)
            //
            do {

                // Compute the centers of the clusters that need updating.
                computeCenters();

                // Compute the stored distances between the updated clusters and the
                // coordinates.
                computeDistances();

                // Make this iteration's assignments.
                moves = makeAssignments();

                it++;

                postKMeansMessage("... iteration " + it + " moves = " + moves);

            } while (moves > 0 && it < mMaxIterations);

            // Transform the array of ProtoClusters to an array
            // of the simpler class Cluster.
            mClusters = generateFinalClusters();

            long executionTime = System.currentTimeMillis() - startTime;

            postKMeansComplete(mClusters, executionTime);

        } catch (Throwable t) {

            postKMeansError(t);

        } finally {

            // Clean up temporary data structures used during the algorithm.
            cleanup();

        }
    }

    /**
     * Randomly select coordinates to be the initial cluster centers.
     */
    private void initCenters() {

        Random random = new Random(mRandomSeed);

        int coordCount = mCoordinates.length;

        // The array mClusterAssignments is used only to keep track of the cluster
        // membership for each coordinate.  The method makeAssignments() uses it
        // to keep track of the number of moves.
        if (mClusterAssignments == null) {
            mClusterAssignments = new int[coordCount];
            // Initialize to -1 to indicate that they haven't been assigned yet.
            Arrays.fill(mClusterAssignments, -1);
        }


        int[] indices = new SmartInitialization(mCoordinates, mK, random.nextInt(coordCount)).getSmartKmeansInitialization();
        mProtoClusters = new ProtoCluster[mK];
        for (int i = 0; i < mK; i++) {
            int coordIndex = indices[i];
            mProtoClusters[i] = new ProtoCluster(mCoordinates[coordIndex], coordIndex);
            mClusterAssignments[indices[i]] = i;
        }
    }

    /**
     * Recompute the centers of the protoclusters with
     * update flags set to true.
     */
    private void computeCenters() {

        // Sets the update flags of the protoclusters that haven't been deleted and
        // whose memberships have changed in the iteration just completed.
        //
        for (int q = 0; q < mProtoClusters.length; q++) {
            ProtoCluster cluster = mProtoClusters[q];
            if (cluster.getConsiderForAssignment()) {
                if (cluster.isNotEmpty()) {
                    // This sets the protocluster's update flag to
                    // true only if its membership changed in last call
                    // to makeAssignments().
                    cluster.setUpdateFlag();
                    // If the update flag was set, update the center.
                    if (cluster.needsUpdate()) {
                        cluster.updateCenter(mCoordinates);
                    }
                } else {
                    // When a cluster loses all of its members, it
                    // falls out of contention.  So it is possible for
                    // k-means to return fewer than k clusters.
                    cluster.setConsiderForAssignment(false);
                }
            }
        }
    }

    /**
     * Compute distances between coodinates and cluster centers,
     * storing them in the distanceChi2 cache.  Only distances that
     * need to be computed are computed.  This is determined by
     * distanceChi2 update flags in the protocluster objects.
     */
    private void computeDistances() throws InsufficientMemoryException {

        if (mDistanceCache == null) {
            int numCoords = mCoordinates.length;
            int numClusters = mProtoClusters.length;
            // Explicit garbage collection to reduce likelihood of insufficient
            // memory.
            System.gc();
            // Ensure there is enough memory available for the distances.
            // Throw an exception if not.
            long memRequired = 8L * numCoords * numClusters;
            if (Runtime.getRuntime().freeMemory() < memRequired) {
                throw new InsufficientMemoryException("Not enough memory for compute distances");
            }
            // Instantiate an array to hold the distances between coordinates
            // and cluster centers
            mDistanceCache = new float[numCoords][numClusters];
        }

        // Bulk of the work is delegated to the
        // SubtaskManager.
        mSubtaskManager.computeDistances();
    }

    /**
     * Assign each coordinate to the nearest cluster.  Called once
     * per iteration.  Returns the number of coordinates that have
     * changed their cluster membership.
     */
    private int makeAssignments() {

        // Checkpoint the clusters, so we'll be able to tell
        // which one have changed after all the assignments have been
        // made.
        for (int q = 0; q < mProtoClusters.length; q++) {
            ProtoCluster mProtoCluster = mProtoClusters[q];
            if (mProtoCluster.getConsiderForAssignment()) {
                mProtoCluster.checkPoint();
            }
        }

        // Bulk of the work is delegated to the SubtaskManager.
        mSubtaskManager.makeAssignments();
        // Get the number of moves from the SubtaskManager.
        return mSubtaskManager.numberOfMoves();
    }

    /**
     * Generate an array of Cluster objects from mProtoClusters.
     *
     * @return array of Cluster object references.
     */
    private Cluster[] generateFinalClusters() {

        int numClusters = mProtoClusters.length;

        // Convert the proto-clusters to the final Clusters.
        //
        // - accumulate in a list.
        List<Cluster> clusterList = new ArrayList<>(numClusters);
        for (ProtoCluster pcluster : mProtoClusters) {
            if (pcluster.isNotEmpty()) {
                Cluster cluster = new Cluster(pcluster.getMembership(), pcluster.getCenter());
                clusterList.add(cluster);
            }
        }

        // - convert list to an array.
        Cluster[] clusters = new Cluster[clusterList.size()];
        clusterList.toArray(clusters);

        return clusters;
    }

    /**
     * Clean up items used by the clustering algorithm that are no longer needed.
     */
    private void cleanup() {
        mProtoClusters = null;
        mDistanceCache = null;
        mClusterAssignments = null;
        if (mSubtaskManager != null) {
            mSubtaskManager.shutdown();
            mSubtaskManager = null;
        }
    }
}
