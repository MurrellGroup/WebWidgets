/**
 * Needleman-Wunsch Alignment Library
 * Hosted at: https://murrellgroup.github.io/WebWidgets/nw.js
 */

function affineNWAlign(s1, s2, gapOpen = -10.0, gapExtend = -0.2, matchCost = 1.0, mismatchCost = -0.7, boundaryGapFactor = 10) {
    const s1arr = s1.split('');
    const s1len = s1arr.length;
    const s2arr = s2.split('');
    const s2len = s2arr.length;
    const M = Array(s1len + 1).fill().map(() => Array(s2len + 1).fill(0));
    const IX = Array(s1len + 1).fill().map(() => Array(s2len + 1).fill(0));
    const IY = Array(s1len + 1).fill().map(() => Array(s2len + 1).fill(0));
    for (let i = 0; i <= s1len; i++) {
        IX[i][0] = gapOpen + gapExtend * i;
        IY[i][0] = -Infinity;
        M[i][0] = gapOpen + gapExtend * i;
    }
    for (let j = 0; j <= s2len; j++) {
        IX[0][j] = -Infinity;
        IY[0][j] = gapOpen + gapExtend * j;
        M[0][j] = gapOpen + gapExtend * j;
    }
    const traceM = Array(s1len + 1).fill().map(() => Array(s2len + 1).fill(0));
    const traceIX = Array(s1len + 1).fill().map(() => Array(s2len + 1).fill(0));
    const traceIY = Array(s1len + 1).fill().map(() => Array(s2len + 1).fill(0));
    traceM[0].fill(3);
    traceIX[0].fill(3);
    traceIY[0].fill(3);
    for (let i = 1; i <= s1len; i++) {
        traceM[i][0] = 2;
        traceIX[i][0] = 2;
        traceIY[i][0] = 2;
    }
    for (let i = 1; i <= s1len; i++) {
        for (let j = 1; j <= s2len; j++) {
            const diagCost = s1arr[i - 1] === s2arr[j - 1] ? matchCost : mismatchCost;
            const diagM = M[i - 1][j - 1] + diagCost;
            const IX2M = IX[i - 1][j - 1] + diagCost;
            const IY2M = IY[i - 1][j - 1] + diagCost;
            [M[i][j], traceM[i][j]] = findMax([diagM, IX2M, IY2M]);

            const boundaryGapExtendX = (i === s1len) ? gapExtend / boundaryGapFactor : gapExtend;
            const boundaryGapExtendY = (j === s2len) ? gapExtend / boundaryGapFactor : gapExtend;
            
            const M2IX = M[i - 1][j] + gapOpen;
            const IXextend = IX[i - 1][j] + boundaryGapExtendY;
            [IX[i][j], traceIX[i][j]] = findMax([M2IX, IXextend]);

            const M2IY = M[i][j - 1] + gapOpen;
            const IYextend = IY[i][j - 1] + boundaryGapExtendX;
            [IY[i][j], traceIY[i][j]] = findMax([M2IY, -Infinity, IYextend]);
        }
    }
    const revArr1 = [];
    const revArr2 = [];
    const mats = [traceM, traceIX, traceIY];
    let xI = s1len;
    let yI = s2len;
    let mI = findMax([M[xI][yI], IX[xI][yI], IY[xI][yI]])[1];
    while (xI > 0 && yI > 0) {
        const nextMI = mats[mI - 1][xI][yI];
        if (mI === 1) {
            revArr1.push(s1arr[xI - 1]);
            revArr2.push(s2arr[yI - 1]);
            xI--;
            yI--;
        } else if (mI === 2) {
            revArr1.push(s1arr[xI - 1]);
            revArr2.push('-');
            xI--;
        } else if (mI === 3) {
            revArr1.push('-');
            revArr2.push(s2arr[yI - 1]);
            yI--;
        }
        mI = nextMI;
    }
    while (xI > 0) {
        revArr1.push(s1arr[xI - 1]);
        revArr2.push('-');
        xI--;
    }
    while (yI > 0) {
        revArr1.push('-');
        revArr2.push(s2arr[yI - 1]);
        yI--;
    }
    return [revArr1.reverse().join(''), revArr2.reverse().join('')];
}

function findMax(arr) {
    let maxVal = arr[0];
    let maxIndex = 1;
    for (let i = 1; i < arr.length; i++) {
        if (arr[i] > maxVal) {
            maxVal = arr[i];
            maxIndex = i + 1;
        }
    }
    return [maxVal, maxIndex];
}









/**
 * Constrained NW Aligner
 * Used for aligning gaps between fixed blocks.
 * 
 * isStart: true if this segment is the very beginning of the sequences (Left Flank).
 * isEnd:   true if this segment is the very end of the sequences (Right Flank).
 * 
 * If isEnd is false, boundaryGapFactor is ignored (treated as 1.0) to prevent free gaps 
 * where we should be connecting to a match.
 * M[0][0] is explicitly 0 to ensure we start from a "Match" state unless isStart implies otherwise.
 */
function constrainedNWAlign(s1, s2, gapOpen, gapExtend, matchCost, mismatchCost, boundaryGapFactor, isStart, isEnd) {
    const s1arr = s1.split('');
    const s1len = s1arr.length;
    const s2arr = s2.split('');
    const s2len = s2arr.length;

    // Use boundary factor ONLY if we are at the true end of the alignment
    const useBoundaryFactor = isEnd ? boundaryGapFactor : 1;

    const M = Array(s1len + 1).fill().map(() => Array(s2len + 1).fill(0));
    const IX = Array(s1len + 1).fill().map(() => Array(s2len + 1).fill(0));
    const IY = Array(s1len + 1).fill().map(() => Array(s2len + 1).fill(0));

    // Initialization
    // M[0][0] = 0 implies we are in a match state (or start of seq)
    M[0][0] = 0; 
    IX[0][0] = -Infinity;
    IY[0][0] = -Infinity;

    for (let i = 1; i <= s1len; i++) {
        IX[i][0] = gapOpen + gapExtend * i;
        IY[i][0] = -Infinity;
        M[i][0] = gapOpen + gapExtend * i;
    }
    for (let j = 1; j <= s2len; j++) {
        IX[0][j] = -Infinity;
        IY[0][j] = gapOpen + gapExtend * j;
        M[0][j] = gapOpen + gapExtend * j;
    }

    const traceM = Array(s1len + 1).fill().map(() => Array(s2len + 1).fill(0));
    const traceIX = Array(s1len + 1).fill().map(() => Array(s2len + 1).fill(0));
    const traceIY = Array(s1len + 1).fill().map(() => Array(s2len + 1).fill(0));
    
    // Traceback initialization
    traceM[0].fill(3);
    traceIX[0].fill(3);
    traceIY[0].fill(3);
    for (let i = 1; i <= s1len; i++) {
        traceM[i][0] = 2;
        traceIX[i][0] = 2;
        traceIY[i][0] = 2;
    }

    for (let i = 1; i <= s1len; i++) {
        for (let j = 1; j <= s2len; j++) {
            const diagCost = s1arr[i - 1] === s2arr[j - 1] ? matchCost : mismatchCost;
            const diagM = M[i - 1][j - 1] + diagCost;
            const IX2M = IX[i - 1][j - 1] + diagCost;
            const IY2M = IY[i - 1][j - 1] + diagCost;
            [M[i][j], traceM[i][j]] = findMax([diagM, IX2M, IY2M]);

            // Apply factor only if isEnd is true (via useBoundaryFactor)
            const boundaryGapExtendX = (i === s1len) ? gapExtend / useBoundaryFactor : gapExtend;
            const boundaryGapExtendY = (j === s2len) ? gapExtend / useBoundaryFactor : gapExtend;
            
            const M2IX = M[i - 1][j] + gapOpen;
            const IXextend = IX[i - 1][j] + boundaryGapExtendY;
            [IX[i][j], traceIX[i][j]] = findMax([M2IX, IXextend]);

            const M2IY = M[i][j - 1] + gapOpen;
            const IYextend = IY[i][j - 1] + boundaryGapExtendX;
            [IY[i][j], traceIY[i][j]] = findMax([M2IY, -Infinity, IYextend]);
        }
    }

    const revArr1 = [];
    const revArr2 = [];
    const mats = [traceM, traceIX, traceIY];
    let xI = s1len;
    let yI = s2len;
    let mI = findMax([M[xI][yI], IX[xI][yI], IY[xI][yI]])[1];

    while (xI > 0 && yI > 0) {
        const nextMI = mats[mI - 1][xI][yI];
        if (mI === 1) {
            revArr1.push(s1arr[xI - 1]);
            revArr2.push(s2arr[yI - 1]);
            xI--;
            yI--;
        } else if (mI === 2) {
            revArr1.push(s1arr[xI - 1]);
            revArr2.push('-');
            xI--;
        } else if (mI === 3) {
            revArr1.push('-');
            revArr2.push(s2arr[yI - 1]);
            yI--;
        }
        mI = nextMI;
    }
    while (xI > 0) {
        revArr1.push(s1arr[xI - 1]);
        revArr2.push('-');
        xI--;
    }
    while (yI > 0) {
        revArr1.push('-');
        revArr2.push(s2arr[yI - 1]);
        yI--;
    }
    return [revArr1.reverse().join(''), revArr2.reverse().join('')];
}

/**
 * Task 1: K-mer Seeded Alignment
 * Uses LIS to find consistent anchors and Erodes edges to prevent boundary issues.
 */
function kmer_seeded_nwalign(s1, s2, gapOpen = -10.0, gapExtend = -0.2, matchCost = 1.0, mismatchCost = -0.7, boundaryGapFactor = 10) {
    const K = 25; 
    const minLen = Math.min(s1.length, s2.length);
    
    // Fallback if sequences are too short to support seeds + erosion
    if (minLen < 3 * K) {
        return affineNWAlign(s1, s2, gapOpen, gapExtend, matchCost, mismatchCost, boundaryGapFactor);
    }

    // 1. Identify Unique K-mers
    const getUniqueKmers = (str) => {
        const counts = new Map();
        const indices = new Map();
        for (let i = 0; i <= str.length - K; i++) {
            const sub = str.substring(i, i + K);
            counts.set(sub, (counts.get(sub) || 0) + 1);
            indices.set(sub, i);
        }
        return { counts, indices };
    };

    const map1 = getUniqueKmers(s1);
    const map2 = getUniqueKmers(s2);

    let matches = [];
    for (const [kmer, count] of map1.counts) {
        if (count === 1 && map2.counts.get(kmer) === 1) {
            matches.push({
                i: map1.indices.get(kmer),
                j: map2.indices.get(kmer),
                len: K
            });
        }
    }

    if (matches.length === 0) {
        return affineNWAlign(s1, s2, gapOpen, gapExtend, matchCost, mismatchCost, boundaryGapFactor);
    }

    // 2. LIS for consistency (Sort by I, LIS on J)
    matches.sort((a, b) => a.i - b.i);
    const lisMatches = getLIS(matches);
    
    // 3. Merge and Erode to create safe anchors
    const anchors = mergeAndErodeAnchors(lisMatches, K);
    
    // 4. Stitch using constrained NW
    return stitchAlignments(s1, s2, anchors, [gapOpen, gapExtend, matchCost, mismatchCost, boundaryGapFactor]);
}

/**
 * Task 2: Double DP Alignment
 * Uses Sparse DP to chain small seeds, then erodes edges.
 */
function doubleDP_nwalign(s1, s2, gapOpen = -10.0, gapExtend = -0.2, matchCost = 1.0, mismatchCost = -0.7, boundaryGapFactor = 10) {
    const K = 11;
    const minLen = Math.min(s1.length, s2.length);

    if (minLen < 3 * K) {
        return affineNWAlign(s1, s2, gapOpen, gapExtend, matchCost, mismatchCost, boundaryGapFactor);
    }

    // 1. Index s1
    const s1Map = new Map();
    for (let i = 0; i <= s1.length - K; i++) {
        const sub = s1.substring(i, i + K);
        if (!s1Map.has(sub)) s1Map.set(sub, []);
        s1Map.get(sub).push(i);
    }

    // 2. Find Seeds in s2
    let matches = [];
    const maxHits = 25; 
    
    for (let j = 0; j <= s2.length - K; j++) {
        const sub = s2.substring(j, j + K);
        const hits = s1Map.get(sub);
        if (hits && hits.length <= maxHits) {
            for (let i of hits) {
                matches.push({ i, j, len: K });
            }
        }
    }

    if (matches.length === 0) {
        return affineNWAlign(s1, s2, gapOpen, gapExtend, matchCost, mismatchCost, boundaryGapFactor);
    }

    matches.sort((a, b) => (a.i - b.i) || (a.j - b.j));

    // 3. Sparse DP Chain
    const scores = new Float64Array(matches.length);
    const parents = new Int32Array(matches.length).fill(-1);
    const lookback = 80;

    for (let cur = 0; cur < matches.length; cur++) {
        const mCurr = matches[cur];
        let maxScore = mCurr.len; 
        
        const startSearch = Math.max(0, cur - lookback);
        for (let prev = cur - 1; prev >= startSearch; prev--) {
            const mPrev = matches[prev];

            // Strict time-forward check
            if (mPrev.i < mCurr.i && mPrev.j < mCurr.j) {
                
                const diagDiff = Math.abs((mCurr.j - mCurr.i) - (mPrev.j - mPrev.i));
                const gapI = mCurr.i - (mPrev.i + mPrev.len);
                const gapJ = mCurr.j - (mPrev.j + mPrev.len);
                const dist = Math.max(0, gapI) + Math.max(0, gapJ);
                
                // Heuristic scoring
                const penalty = (diagDiff * 3.0) + (dist * 0.1); 
                
                const newScore = scores[prev] + mCurr.len - penalty;
                if (newScore > maxScore) {
                    maxScore = newScore;
                    parents[cur] = prev;
                }
            }
        }
        scores[cur] = maxScore;
    }

    let bestIdx = 0;
    let maxS = scores[0];
    for(let i=1; i<matches.length; i++){
        if(scores[i] > maxS){
            maxS = scores[i];
            bestIdx = i;
        }
    }

    const chain = [];
    let curr = bestIdx;
    while(curr !== -1) {
        chain.push(matches[curr]);
        curr = parents[curr];
    }
    chain.reverse();

    // 4. Merge and Erode
    const anchors = mergeAndErodeAnchors(chain, K);

    return stitchAlignments(s1, s2, anchors, [gapOpen, gapExtend, matchCost, mismatchCost, boundaryGapFactor]);
}

// --- Helpers ---

function getLIS(matches) {
    if (matches.length === 0) return [];
    const tails = []; 
    const parent = new Int32Array(matches.length).fill(-1);
    for (let i = 0; i < matches.length; i++) {
        const val = matches[i].j;
        let left = 0, right = tails.length;
        while (left < right) {
            const mid = (left + right) >>> 1;
            if (matches[tails[mid]].j < val) left = mid + 1;
            else right = mid;
        }
        if (left < tails.length) {
            tails[left] = i;
            parent[i] = (left > 0) ? tails[left - 1] : -1;
        } else {
            tails.push(i);
            parent[i] = (tails.length > 1) ? tails[tails.length - 2] : -1;
        }
    }
    const result = [];
    let curr = tails[tails.length - 1];
    while (curr !== -1) {
        result.push(matches[curr]);
        curr = parent[curr];
    }
    return result.reverse();
}

/**
 * Merges continuous/overlapping blocks on same diagonal, 
 * then erodes ends by 'margin' to ensure safe stitching.
 */
function mergeAndErodeAnchors(matches, margin) {
    if (matches.length === 0) return [];
    
    // Phase 1: Merge
    const merged = [];
    let curr = { ...matches[0] };

    for (let k = 1; k < matches.length; k++) {
        const next = matches[k];
        
        const diagCurr = curr.j - curr.i;
        const diagNext = next.j - next.i;
        
        // Check connectivity/overlap
        const isConnected = next.i <= (curr.i + curr.len);

        if (diagCurr === diagNext && isConnected) {
            const newEnd = Math.max(curr.i + curr.len, next.i + next.len);
            curr.len = newEnd - curr.i;
        } else {
            merged.push(curr);
            curr = { ...next };
        }
    }
    merged.push(curr);

    // Phase 2: Erode
    const eroded = [];
    for (const m of merged) {
        // "Chew back" margin from both sides
        const newLen = m.len - (2 * margin);
        
        if (newLen > 0) {
            eroded.push({
                i: m.i + margin,
                j: m.j + margin,
                len: newLen
            });
        }
    }
    return eroded;
}

function stitchAlignments(s1, s2, anchors, nwParams) {
    let finalS1 = "";
    let finalS2 = "";
    let idx1 = 0;
    let idx2 = 0;

    // nwParams: [gapOpen, gapExtend, matchCost, mismatchCost, boundaryGapFactor]

    for (let k = 0; k < anchors.length; k++) {
        const anchor = anchors[k];

        // --- GAP REGION ---
        const sub1 = s1.substring(idx1, anchor.i);
        const sub2 = s2.substring(idx2, anchor.j);
        
        // Check if we are at Start or End of entire sequence
        // Left Flank (Start): k === 0 && idx1 === 0
        const isStart = (k === 0 && idx1 === 0);
        // Internal Gap: isStart = false, isEnd = false (since we are connecting to an anchor)
        const isEnd = false;

        // Perform constrained alignment on the gap
        // Even if strings are empty, run it if we need to ensure state transition logic, 
        // though typically empty aligns to empty.
        if (sub1.length > 0 || sub2.length > 0) {
            const [aln1, aln2] = constrainedNWAlign(sub1, sub2, ...nwParams, isStart, isEnd);
            finalS1 += aln1;
            finalS2 += aln2;
        }

        // --- ANCHOR REGION (Perfect Match) ---
        const anchorSeq = s1.substring(anchor.i, anchor.i + anchor.len);
        finalS1 += anchorSeq;
        finalS2 += anchorSeq;

        idx1 = anchor.i + anchor.len;
        idx2 = anchor.j + anchor.len;
    }

    // --- TAIL REGION (Right Flank) ---
    if (idx1 < s1.length || idx2 < s2.length) {
        const sub1 = s1.substring(idx1);
        const sub2 = s2.substring(idx2);
        
        // This is the Right Flank: isStart=false, isEnd=true
        const [aln1, aln2] = constrainedNWAlign(sub1, sub2, ...nwParams, false, true);
        finalS1 += aln1;
        finalS2 += aln2;
    }

    return [finalS1, finalS2];
}
