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
 * K-mer Seeded NW Alignment (Task 1)
 *
 * Strategy:
 * 1. Find all 'exact unique' k-mers (default k=25) that appear exactly once in s1 and exactly once in s2.
 * 2. Sort these anchors.
 * 3. Filter for strictly increasing consistency (Longest Increasing Subsequence logic, though usually consistent by nature of uniqueness).
 * 4. Merge adjacent/overlapping anchors on the same diagonal into larger blocks.
 * 5. Run standard affineNWAlign on the gaps between blocks.
 */
function kmer_seeded_nwalign(s1, s2, gapOpen = -10.0, gapExtend = -0.2, matchCost = 1.0, mismatchCost = -0.7, boundaryGapFactor = 10) {
    // 1. Configuration
    const K = 25; // Large K for uniqueness
    const minLen = Math.min(s1.length, s2.length);
    
    // Fallback to standard NW if sequences are too short for seeding
    if (minLen < K) {
        return affineNWAlign(s1, s2, gapOpen, gapExtend, matchCost, mismatchCost, boundaryGapFactor);
    }

    // 2. Identify Unique K-mers in O(N)
    const getKmerCounts = (str) => {
        const counts = new Map();
        const indices = new Map();
        for (let i = 0; i <= str.length - K; i++) {
            const sub = str.substring(i, i + K);
            counts.set(sub, (counts.get(sub) || 0) + 1);
            indices.set(sub, i);
        }
        return { counts, indices };
    };

    const map1 = getKmerCounts(s1);
    const map2 = getKmerCounts(s2);

    // Collect matches that appear exactly once in both sequences
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

    // 3. Sort by position in s1
    matches.sort((a, b) => a.i - b.i);

    // 4. LIS / Filter for consistency
    // Since k-mers are unique, we just need to ensure s2 indices (j) are increasing
    const validMatches = [];
    let lastJ = -1;
    for (const m of matches) {
        if (m.j > lastJ) {
            validMatches.push(m);
            lastJ = m.j;
        }
    }

    // 5. Merge overlapping/adjacent anchors on the same diagonal
    // This turns individual k-mer seeds into long exact match blocks
    const anchors = mergeAnchors(validMatches);

    // 6. Stitch alignments
    return stitchAlignments(s1, s2, anchors, [gapOpen, gapExtend, matchCost, mismatchCost, boundaryGapFactor]);
}


/**
 * Double DP NW Alignment (Task 2)
 *
 * Strategy:
 * 1. Find all matches for small k-mers (k=9). Allow multiple matches.
 * 2. Prune high-frequency k-mers to keep complexity linear-ish.
 * 3. Run a Sparse DP (Chain of Seeds) to find the lowest cost path through the matches.
 *    - Cost heuristic combines match reward vs gap penalty.
 * 4. Convert the optimal chain into anchors.
 * 5. Run standard affineNWAlign on the gaps.
 */
function doubleDP_nwalign(s1, s2, gapOpen = -10.0, gapExtend = -0.2, matchCost = 1.0, mismatchCost = -0.7, boundaryGapFactor = 10) {
    const K = 9; // Small K for sensitivity
    const minLen = Math.min(s1.length, s2.length);

    if (minLen < K) {
        return affineNWAlign(s1, s2, gapOpen, gapExtend, matchCost, mismatchCost, boundaryGapFactor);
    }

    // 1. Index s1 k-mers
    const s1Map = new Map();
    for (let i = 0; i <= s1.length - K; i++) {
        const sub = s1.substring(i, i + K);
        if (!s1Map.has(sub)) s1Map.set(sub, []);
        s1Map.get(sub).push(i);
    }

    // 2. Scan s2 and collect matches
    let matches = [];
    const maxHits = 10; // Pruning: Ignore repetitive regions (e.g. poly-A) to avoid O(N^2)
    
    for (let j = 0; j <= s2.length - K; j++) {
        const sub = s2.substring(j, j + K);
        const hits = s1Map.get(sub);
        if (hits && hits.length <= maxHits) {
            for (let i of hits) {
                matches.push({ i, j, len: K });
            }
        }
    }

    // If no matches, fallback
    if (matches.length === 0) {
        return affineNWAlign(s1, s2, gapOpen, gapExtend, matchCost, mismatchCost, boundaryGapFactor);
    }

    // Sort by S1 index
    matches.sort((a, b) => a.i - b.i);

    // 3. Sparse DP (Chain of Seeds)
    // dp[x] = max score ending at match x
    // To keep it fast, we look back at a limited window of previous matches
    const lookback = 50; 
    const scores = new Float64Array(matches.length);
    const parents = new Int32Array(matches.length).fill(-1);

    // Approximate scoring for the chaining phase
    // We want to maximize coverage (len) while minimizing gap distance
    // We treat the gap penalty roughly linear here for speed in seed selection
    for (let cur = 0; cur < matches.length; cur++) {
        const mCurr = matches[cur];
        let maxScore = mCurr.len; // Base score is just the match length
        
        // Look at previous matches to extend
        // We iterate backwards. Because it's sorted by i, i_prev <= i_curr
        let startSearch = Math.max(0, cur - lookback);
        
        for (let prev = cur - 1; prev >= startSearch; prev--) {
            const mPrev = matches[prev];

            // Must be strictly previous in both dimensions to chain
            // Note: We allow slight overlap in the chain logic, but will flatten later.
            // Strict ordering: i and j must increase
            if (mPrev.i < mCurr.i && mPrev.j < mCurr.j) {
                
                // Calculate Gap
                const di = mCurr.i - (mPrev.i + mPrev.len);
                const dj = mCurr.j - (mPrev.j + mPrev.len);
                
                // We only link if gaps are non-negative (no overlapping blocks in time)
                // or small overlaps (which we handle by merging). 
                // For simplicity in this heuristic, we prefer non-overlapping starts.
                
                // Distance in base pairs
                const dist = Math.abs(di) + Math.abs(dj); 
                
                // Heuristic score: Previous Score + New Match Length - Penalty
                // Penalty is roughly 1.0 per base gap distance
                const penalty = (dist > 0) ? (dist * 0.5) : 0; 
                
                const newScore = scores[prev] + mCurr.len - penalty;
                
                if (newScore > maxScore) {
                    maxScore = newScore;
                    parents[cur] = prev;
                }
            }
        }
        scores[cur] = maxScore;
    }

    // Find best end point
    let bestIdx = 0;
    let maxS = scores[0];
    for (let i = 1; i < matches.length; i++) {
        if (scores[i] > maxS) {
            maxS = scores[i];
            bestIdx = i;
        }
    }

    // Backtrack to build the chain
    let chain = [];
    let curr = bestIdx;
    while (curr !== -1) {
        chain.push(matches[curr]);
        curr = parents[curr];
    }
    chain.reverse();

    // 4. Merge Anchors
    const anchors = mergeAnchors(chain);

    // 5. Stitch
    return stitchAlignments(s1, s2, anchors, [gapOpen, gapExtend, matchCost, mismatchCost, boundaryGapFactor]);
}

// --- Helpers ---

/**
 * Merges a list of match objects {i, j, len} if they are adjacent
 * or overlapping on the same diagonal.
 */
function mergeAnchors(matches) {
    if (matches.length === 0) return [];
    
    const anchors = [];
    let curr = { ...matches[0] };

    for (let k = 1; k < matches.length; k++) {
        const next = matches[k];
        
        // Check if on same diagonal: (j - i) must be constant
        const diagCurr = curr.j - curr.i;
        const diagNext = next.j - next.i;

        // Check continuity: next starts before or exactly where curr ends
        // We allow overlaps because our DP might select (0,0,10) and (5,5,10).
        // Since they are on the same diagonal, this is just one long match (0,0,15).
        const overlapOrTouch = next.i <= (curr.i + curr.len);

        if (diagCurr === diagNext && overlapOrTouch) {
            // Extend current
            const newEnd = Math.max(curr.i + curr.len, next.i + next.len);
            curr.len = newEnd - curr.i;
        } else {
            anchors.push(curr);
            curr = { ...next };
        }
    }
    anchors.push(curr);
    return anchors;
}

/**
 * Takes the original sequences and a list of "perfect match" anchors.
 * Runs affineNWAlign on the gaps between anchors and constructs the final strings.
 */
function stitchAlignments(s1, s2, anchors, nwParams) {
    let finalS1 = "";
    let finalS2 = "";
    let idx1 = 0;
    let idx2 = 0;

    for (const anchor of anchors) {
        // 1. Align the gap before this anchor
        const sub1 = s1.substring(idx1, anchor.i);
        const sub2 = s2.substring(idx2, anchor.j);

        if (sub1.length > 0 || sub2.length > 0) {
            const [aln1, aln2] = affineNWAlign(sub1, sub2, ...nwParams);
            finalS1 += aln1;
            finalS2 += aln2;
        }

        // 2. Add the anchor (perfect match)
        // Note: s1 and s2 matches at anchor are identical by definition of the anchor
        const anchorSeq = s1.substring(anchor.i, anchor.i + anchor.len);
        finalS1 += anchorSeq;
        finalS2 += anchorSeq;

        // 3. Advance indices
        idx1 = anchor.i + anchor.len;
        idx2 = anchor.j + anchor.len;
    }

    // 4. Align any remaining tail
    if (idx1 < s1.length || idx2 < s2.length) {
        const sub1 = s1.substring(idx1);
        const sub2 = s2.substring(idx2);
        const [aln1, aln2] = affineNWAlign(sub1, sub2, ...nwParams);
        finalS1 += aln1;
        finalS2 += aln2;
    }

    return [finalS1, finalS2];
}
