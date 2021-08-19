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

public class LoopAnchor {
    private final List<Integer> connectedLoops = new ArrayList<>();
    private final int buffer;
    private long position;

    public LoopAnchor(long midBin, int id, int buffer) {
        connectedLoops.add(id);
        position = midBin;
        this.buffer = buffer;
    }


    public boolean overlaps(long midBin) {
        return Math.abs(midBin - position) < buffer;
    }

    public void add(long midBin, int id) {
        long currentTotal = position * connectedLoops.size();
        currentTotal += midBin;
        connectedLoops.add(id);
        position = currentTotal / connectedLoops.size();
    }

    public int size() {
        return connectedLoops.size();
    }

    public void populateAdjacencyMatrix(Map<Integer, Set<Integer>> adjacencyMatrix) {
        for (int i : connectedLoops) {
            if (adjacencyMatrix.containsKey(i)) {
                adjacencyMatrix.get(i).addAll(connectedLoops);
            } else {
                adjacencyMatrix.put(i, new HashSet<>(connectedLoops));
            }
        }
    }
}
