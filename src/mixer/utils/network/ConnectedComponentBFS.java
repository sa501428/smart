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

package mixer.utils.network;

import java.util.*;

public class ConnectedComponentBFS {
    private final List<Set<Integer>> largeNetworks = new ArrayList<>();
    private final int BIG_SIZE = 25;

    private final int[] nodeStatus;
    private final List<Integer> clusterSizes = new ArrayList<>();

    public ConnectedComponentBFS(Map<Integer, Set<Integer>> adjacencyMatrix, int numNodes) {
        nodeStatus = new int[numNodes];
        Arrays.fill(nodeStatus, -1);


        int clusterCounter = 0;
        for (int k = 0; k < nodeStatus.length; k++) {
            if (nodeIsNotAssigned(k)) {
                int size = BFS(adjacencyMatrix, k, clusterCounter);
                clusterSizes.add(size);
                clusterCounter++;
            }
        }
    }

    private boolean nodeIsNotAssigned(int k) {
        return nodeStatus[k] < 0;
    }

    private int BFS(Map<Integer, Set<Integer>> adjacencyMatrix, int startingVertex, int clusterID) {
        Set<Integer> nodesThisRound = new HashSet<>();
        Queue<Integer> queue = new LinkedList<>();
        queue.add(startingVertex);
        nodesThisRound.add(startingVertex);
        nodeStatus[startingVertex] = clusterID;
        int size = 1;

        while (!queue.isEmpty()) {
            Integer current = queue.remove();
            for (int i = 0; i < nodeStatus.length; i++)
                if (isConnected(adjacencyMatrix, current, i) && nodeIsNotAssigned(i)) {
                    nodeStatus[i] = clusterID;
                    queue.add(i);
                    nodesThisRound.add(i);
                    size++;
                }
        }

        if (size > BIG_SIZE) {
            largeNetworks.add(nodesThisRound);
        }

        return size;
    }

    private boolean isConnected(Map<Integer, Set<Integer>> adjacencyMatrix, Integer node1, Integer node2) {
        if (adjacencyMatrix.containsKey(node1)) {
            return adjacencyMatrix.get(node1).contains(node2);
        }
        return false;
    }

    public List<Integer> getClusterSizes() {
        return clusterSizes;
    }

    public List<Set<Integer>> getLargeNetworks() {
        return largeNetworks;
    }
}
