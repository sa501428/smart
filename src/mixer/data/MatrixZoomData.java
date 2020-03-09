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


package mixer.data;

import mixer.HiC;
import mixer.MixerGlobals;
import mixer.matrix.BasicMatrix;
import mixer.track.HiCFixedGridAxis;
import mixer.track.HiCFragmentAxis;
import mixer.track.HiCGridAxis;
import mixer.windowui.HiCZoom;
import mixer.windowui.NormalizationHandler;
import mixer.windowui.NormalizationType;
import org.broad.igv.feature.Chromosome;
import org.broad.igv.util.collections.LRUCache;

import java.io.IOException;
import java.util.*;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;
import java.util.concurrent.TimeUnit;
import java.util.concurrent.atomic.AtomicInteger;



/**
 * @author jrobinso
 * @since Aug 10, 2010
 */
public class MatrixZoomData {

    final Chromosome chr1;  // Chromosome on the X axis
    final Chromosome chr2;  // Chromosome on the Y axis
    final HiCZoom zoom;    // Unit and bin size
    private final HiCGridAxis xGridAxis;
    private final HiCGridAxis yGridAxis;
    // Observed values are organized into sub-matrices ("blocks")
    private final int blockBinCount;   // block size in bins
    private final int blockColumnCount;     // number of block columns
    // Cache the last 20 blocks loaded
    private final LRUCache<String, Block> blockCache = new LRUCache<>(500);
    private final HashMap<NormalizationType, BasicMatrix> pearsonsMap;
    private final HashMap<NormalizationType, BasicMatrix> normSquaredMaps;
    private final HashSet<NormalizationType> missingPearsonFiles;
    DatasetReader reader;
    private double averageCount = -1;

    /**
     * Constructor, sets the grid axes.  Called when read from file.
     *
     * @param chr1             Chromosome 1
     * @param chr2             Chromosome 2
     * @param zoom             Zoom (bin size and BP or FRAG)
     * @param blockBinCount    Number of bins divided by number of columns (around 1000)
     * @param blockColumnCount Number of bins divided by 1000 (BLOCK_SIZE)
     * @param chr1Sites        Used for looking up fragment
     * @param chr2Sites        Used for looking up fragment
     * @param reader           Pointer to file reader
     */
    public MatrixZoomData(Chromosome chr1, Chromosome chr2, HiCZoom zoom, int blockBinCount, int blockColumnCount,
                          int[] chr1Sites, int[] chr2Sites, DatasetReader reader) {
        this.chr1 = chr1;
        this.chr2 = chr2;
        this.zoom = zoom;
        this.reader = reader;
        this.blockBinCount = blockBinCount;
        this.blockColumnCount = blockColumnCount;

        int correctedBinCount = blockBinCount;
        if (reader.getVersion() < 8 && chr1.getLength() < chr2.getLength()) {
            boolean isFrag = zoom.getUnit() == HiC.Unit.FRAG;
            int len1 = chr1.getLength();
            int len2 = chr2.getLength();
            if (chr1Sites != null && chr2Sites != null && isFrag) {
                len1 = chr1Sites.length + 1;
                len2 = chr2Sites.length + 1;
            }
            int nBinsX = Math.max(len1, len2) / zoom.getBinSize() + 1;
            correctedBinCount = nBinsX / blockColumnCount + 1;
        }

        if (this instanceof CustomMatrixZoomData) {
            this.xGridAxis = new HiCFixedGridAxis(chr1.getLength() / zoom.getBinSize() + 1, zoom.getBinSize(), null);
            this.yGridAxis = new HiCFixedGridAxis(chr2.getLength() / zoom.getBinSize() + 1, zoom.getBinSize(), null);
        } else if (zoom.getUnit() == HiC.Unit.BP) {
            this.xGridAxis = new HiCFixedGridAxis(correctedBinCount * blockColumnCount, zoom.getBinSize(), chr1Sites);
            this.yGridAxis = new HiCFixedGridAxis(correctedBinCount * blockColumnCount, zoom.getBinSize(), chr2Sites);
        } else {
            this.xGridAxis = new HiCFragmentAxis(zoom.getBinSize(), chr1Sites, chr1.getLength());
            this.yGridAxis = new HiCFragmentAxis(zoom.getBinSize(), chr2Sites, chr2.getLength());
        }

        pearsonsMap = new HashMap<>();
        normSquaredMaps = new HashMap<>();
        missingPearsonFiles = new HashSet<>();
    }

    public Chromosome getChr1() {
        return chr1;
    }

    public Chromosome getChr2() {
        return chr2;
    }

    public HiCGridAxis getXGridAxis() {
        return xGridAxis;
    }

    public HiCGridAxis getYGridAxis() {
        return yGridAxis;
    }

    public int getBinSize() {
        return zoom.getBinSize();
    }

    public int getChr1Idx() {
        return chr1.getIndex();
    }

    public int getChr2Idx() {
        return chr2.getIndex();
    }

    public HiCZoom getZoom() {
        return zoom;
    }

    private int getBlockColumnCount() {
        return blockColumnCount;
    }

    public String getKey() {
        return chr1.getName() + "_" + chr2.getName() + "_" + zoom.getKey();
    }

    // i think this is how it should be? todo sxxgrc please confirm use case
    private String getKey(int chr1, int chr2) {
        return chr1 + "_" + chr2 + "_" + zoom.getKey();
    }

    public String getBlockKey(int blockNumber, NormalizationType no) {
        return getKey() + "_" + blockNumber + "_" + no;
    }

    public String getNormLessBlockKey(Block block) {
        return getKey() + "_" + block.getNumber() + "_";
    }

    private String getBlockKey(int blockNumber, NormalizationType no, int chr1, int chr2) {
        return getKey(chr1, chr2) + "_" + blockNumber + "_" + no;
    }

    /**
     * Return the blocks of normalized, observed values overlapping the rectangular region specified.
     * The units are "bins"
     *
     * @param binY1 leftmost position in "bins"
     * @param binX2 rightmost position in "bins"
     * @param binY2 bottom position in "bins"
     * @param no    normalization type
     * @param isImportant used for debugging
     * @return List of overlapping blocks, normalized
     */
    public List<Block> getNormalizedBlocksOverlapping(int binX1, int binY1, int binX2, int binY2, final NormalizationType no,
                                                      boolean isImportant, boolean fillUnderDiagonal) {
        final List<Block> blockList = new ArrayList<>();
        Block b = new Block(1, getBlockKey(1, no));
        return addNormalizedBlocksToList(blockList, binX1, binY1, binX2, binY2, no, fillUnderDiagonal);
    }

    private void populateBlocksToLoad(int r, int c, NormalizationType no, List<Block> blockList, Set<Integer> blocksToLoad) {
        int blockNumber = r * getBlockColumnCount() + c;
        String key = getBlockKey(blockNumber, no);
        Block b;
        if (MixerGlobals.useCache && blockCache.containsKey(key)) {
            b = blockCache.get(key);
            blockList.add(b);
        } else {
            blocksToLoad.add(blockNumber);
        }
    }

    /**
     * Return the blocks of normalized, observed values overlapping the rectangular region specified.
     *
     * @param binY1 leftmost position in "bins"
     * @param binX2 rightmost position in "bins"
     * @param binY2 bottom position in "bins"
     * @param no    normalization type
     * @return List of overlapping blocks, normalized
     */
    private List<Block> addNormalizedBlocksToList(final List<Block> blockList, int binX1, int binY1, int binX2, int binY2,
                                                  final NormalizationType no, boolean getBelowDiagonal) {

        Set<Integer> blocksToLoad = new HashSet<>();
      
        // have to do this regardless (just in case)
        int col1 = binX1 / blockBinCount;
        int row1 = binY1 / blockBinCount;
        int col2 = binX2 / blockBinCount;
        int row2 = binY2 / blockBinCount;

        for (int r = row1; r <= row2; r++) {
            for (int c = col1; c <= col2; c++) {
                populateBlocksToLoad(r, c, no, blockList, blocksToLoad);
            }
        }

        if (getBelowDiagonal && binY1 < binX2) {
            for (int r = row1; r <= row2; r++) {
                for (int c = col1; c <= col2; c++) {
                    populateBlocksToLoad(c, r, no, blockList, blocksToLoad);
                }
            }
        }

        actuallyLoadGivenBlocks(blockList, blocksToLoad, no);

        return new ArrayList<>(new HashSet<>(blockList));
    }

    private List<Block> addNormalizedBlocksToList(final List<Block> blockList, int binX1, int binY1, int binX2, int binY2,
                                                  final NormalizationType no, int chr1, int chr2) {

        Set<Integer> blocksToLoad = new HashSet<>();

        // have to do this regardless (just in case)
        int col1 = binX1 / blockBinCount;
        int row1 = binY1 / blockBinCount;
        int col2 = binX2 / blockBinCount;
        int row2 = binY2 / blockBinCount;

        for (int r = row1; r <= row2; r++) {
            for (int c = col1; c <= col2; c++) {
                populateBlocksToLoad(r, c, no, blockList, blocksToLoad);
            }
        }

        actuallyLoadGivenBlocks(blockList, blocksToLoad, no, chr1, chr2);
        System.out.println("I am block size: " + blockList.size());
        System.out.println("I am first block: " + blockList.get(0).getNumber());
        return new ArrayList<>(new HashSet<>(blockList));
    }


    private void actuallyLoadGivenBlocks(final List<Block> blockList, Set<Integer> blocksToLoad,
                                         final NormalizationType no) {
        final AtomicInteger errorCounter = new AtomicInteger();

        ExecutorService service = Executors.newFixedThreadPool(200);

        final int binSize = getBinSize();
        final int chr1Index = chr1.getIndex();
        final int chr2Index = chr2.getIndex();

        for (final int blockNumber : blocksToLoad) {
            Runnable loader = new Runnable() {
                @Override
                public void run() {
                    try {
                        String key = getBlockKey(blockNumber, no);
                        Block b = reader.readNormalizedBlock(blockNumber, MatrixZoomData.this, no);
                        if (b == null) {
                            b = new Block(blockNumber, key);   // An empty block
                        }
                        //Run out of memory if do it here
                        if (MixerGlobals.useCache) {
                            blockCache.put(key, b);
                        }
                        blockList.add(b);
                    } catch (IOException e) {
                        errorCounter.incrementAndGet();
                    }
                }
            };

            service.submit(loader);
        }

        // done submitting all jobs
        service.shutdown();

        // wait for all to finish
        try {
            service.awaitTermination(Long.MAX_VALUE, TimeUnit.NANOSECONDS);
        } catch (InterruptedException e) {
            System.err.println("Error loading mzd data " + e.getLocalizedMessage());
            if (MixerGlobals.printVerboseComments) {
                e.printStackTrace();
            }
        }

        // error printing
        if (errorCounter.get() > 0) {
            System.err.println(errorCounter.get() + " errors while reading blocks");
        }
    }

    private void actuallyLoadGivenBlocks(final List<Block> blockList, Set<Integer> blocksToLoad,
                                         final NormalizationType no, final int chr1Id, final int chr2Id) {
        final AtomicInteger errorCounter = new AtomicInteger();

        ExecutorService service = Executors.newFixedThreadPool(200);

        final int binSize = getBinSize();

        for (final int blockNumber : blocksToLoad) {
            Runnable loader = new Runnable() {
                @Override
                public void run() {
                    try {
                        String key = getBlockKey(blockNumber, no, chr1Id, chr2Id);
                        Block b = reader.readNormalizedBlock(blockNumber, MatrixZoomData.this, no);
                        if (b == null) {
                            b = new Block(blockNumber, key);   // An empty block
                        }
                        //Run out of memory if do it here
                        if (MixerGlobals.useCache) {
                            blockCache.put(key, b);
                        }
                        blockList.add(b);
                    } catch (IOException e) {
                        errorCounter.incrementAndGet();
                    }
                }
            };

            service.submit(loader);
        }

        // done submitting all jobs
        service.shutdown();

        // wait for all to finish
        try {
            service.awaitTermination(Long.MAX_VALUE, TimeUnit.NANOSECONDS);
        } catch (InterruptedException e) {
            System.err.println("Error loading mzd data " + e.getLocalizedMessage());
            if (MixerGlobals.printVerboseComments) {
                e.printStackTrace();
            }
        }

        // error printing
        if (errorCounter.get() > 0) {
            System.err.println(errorCounter.get() + " errors while reading blocks");
        }
    }


    /**
     * Return the observed value at the specified location. Supports tooltip text
     * This implementation is naive, but might get away with it for tooltip.
     *
     * @param binX              X bin
     * @param binY              Y bin
     * @param normalizationType Normalization type
     */
    public float getObservedValue(int binX, int binY, NormalizationType normalizationType) {

        // Intra stores only lower diagonal
        if (chr1 == chr2) {
            if (binX > binY) {
                int tmp = binX;
                //noinspection SuspiciousNameCombination
                binX = binY;
                binY = tmp;
            }
        }

        List<Block> blocks = getNormalizedBlocksOverlapping(binX, binY, binX, binY, normalizationType, false, false);
        if (blocks == null) return 0;
        for (Block b : blocks) {
            for (ContactRecord rec : b.getContactRecords()) {
                if (rec.getBinX() == binX && rec.getBinY() == binY) {
                    return rec.getCounts();
                }
            }
        }
        // No record found for this bin
        return 0;
    }

//    /**
//     * Return a slice of the matrix at the specified Y been as a list of wiggle scores
//     *
//     * @param binY
//     */
//    public List<BasicScore> getSlice(int startBinX, int endBinX, int binY, NormalizationType normalizationType) {
//
//        // Intra stores only lower diagonal
//        if (chr1 == chr2) {
//            if (binX > binY) {
//                int tmp = binX;
//                binX = binY;
//                binY = tmp;
//
//            }
//        }
//
//        List<Block> blocks = getNormalizedBlocksOverlapping(binX, binY, binX, binY, normalizationType);
//        if (blocks == null) return 0;
//        for (Block b : blocks) {
//            for (ContactRecord rec : b.getContactRecords()) {
//                if (rec.getBinX() == binX && rec.getBinY() == binY) {
//                    return rec.getCounts();
//                }
//            }
//        }
//        // No record found for this bin
//        return 0;
//    }



    /**
     * Utility for printing description of this matrix.
     */
    public String getDescription() {
        return chr1.getName() + " - " + chr2.getName() + " - " + getZoom();
    }



    /**
     * For a specified region, select the block numbers corresponding to it
     * @param regionIndices
     * @return
     */
    List<Integer> getBlockNumbersForRegionFromGenomePosition(int[] regionIndices) {
        int resolution = zoom.getBinSize();
        int[] regionBinIndices = new int[4];
        for (int i = 0; i < regionBinIndices.length; i++) {
            regionBinIndices[i] = regionIndices[i] / resolution;
        }
        return getBlockNumbersForRegionFromBinPosition(regionBinIndices);
    }

    private List<Integer> getBlockNumbersForRegionFromBinPosition(int[] regionIndices) {

        int col1 = regionIndices[0] / blockBinCount;
        int col2 = (regionIndices[1] + 1) / blockBinCount;
        int row1 = regionIndices[2] / blockBinCount;
        int row2 = (regionIndices[3] + 1) / blockBinCount;

        // first check the upper triangular matrix
        Set<Integer> blocksSet = new HashSet<>();
        for (int r = row1; r <= row2; r++) {
            for (int c = col1; c <= col2; c++) {
                int blockNumber = r * getBlockColumnCount() + c;
                blocksSet.add(blockNumber);
            }
        }
        // check region part that overlaps with lower left triangle
        // but only if intrachromosomal
        if (chr1.getIndex() == chr2.getIndex()) {
            for (int r = col1; r <= col2; r++) {
                for (int c = row1; c <= row2; c++) {
                    int blockNumber = r * getBlockColumnCount() + c;
                    blocksSet.add(blockNumber);
                }
            }
        }

        List<Integer> blocksToIterateOver = new ArrayList<>(blocksSet);
        Collections.sort(blocksToIterateOver);
        return blocksToIterateOver;
    }





    /**
     * Returns the average count
     *
     * @return Average count
     */
    public double getAverageCount() {
        return averageCount;
    }

    /**
     * Sets the average count
     *
     * @param averageCount Average count to set
     */
    public void setAverageCount(double averageCount) {
        this.averageCount = averageCount;
    }

    /**
     * Returns iterator for contact records
     *
     * @return iterator for contact records
     */
    public Iterator<ContactRecord> getNewContactRecordIterator() {
        return new ContactRecordIterator();
    }

    public List<ContactRecord> getContactRecordList() {
        List<ContactRecord> records = new ArrayList<>();
        Iterator<ContactRecord> iterator = getNewContactRecordIterator();
        while (iterator.hasNext()) {
            ContactRecord cr = iterator.next();
            records.add(cr);
        }
        return records;
    }

    public void clearCache() {
        blockCache.clear();
    }


    /**
     * Class for iterating over the contact records
     */
    class ContactRecordIterator implements Iterator<ContactRecord> {

        final List<Integer> blockNumbers;
        int blockIdx;
        Iterator<ContactRecord> currentBlockIterator;

        /**
         * Initializes the iterator
         */
        ContactRecordIterator() {
            this.blockIdx = -1;
            this.blockNumbers = reader.getBlockNumbers(MatrixZoomData.this);
        }

        /**
         * Indicates whether or not there is another block waiting; checks current block
         * iterator and creates a new one if need be
         *
         * @return true if there is another block to be read
         */
        @Override
        public boolean hasNext() {

            if (currentBlockIterator != null && currentBlockIterator.hasNext()) {
                return true;
            } else {
                blockIdx++;
                if (blockIdx < blockNumbers.size()) {
                    try {
                        int blockNumber = blockNumbers.get(blockIdx);

                        // Optionally check the cache
                        // TODO why is this always NONE, should trace to ensure hard coding doesn't cause bug?
                        String key = getBlockKey(blockNumber, NormalizationHandler.NONE);
                        Block nextBlock;
                        if (MixerGlobals.useCache && blockCache.containsKey(key)) {
                            nextBlock = blockCache.get(key);
                        } else {
                            nextBlock = reader.readNormalizedBlock(blockNumber, MatrixZoomData.this, NormalizationHandler.NONE);
                        }
                        currentBlockIterator = nextBlock.getContactRecords().iterator();
                        return true;
                    } catch (IOException e) {
                        System.err.println("Error fetching block " + e.getMessage());
                        return false;
                    }
                }
            }

            return false;
        }

        /**
         * Returns the next contact record
         *
         * @return The next contact record
         */
        @Override
        public ContactRecord next() {
            return currentBlockIterator == null ? null : currentBlockIterator.next();
        }

        /**
         * Not supported
         */
        @Override
        public void remove() {
            //Not supported
            throw new RuntimeException("remove() is not supported");
        }
    }
}
