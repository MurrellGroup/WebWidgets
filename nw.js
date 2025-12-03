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
 * ROBUST K-mer Seeded NW Alignment (Task 1)
 * 
 * Improvements:
 * - Uses O(N log N) Longest Increasing Subsequence (LIS) to find the maximal consistent set of unique anchors.
 * - Handles overlapping anchors gracefully by trimming.
 * - Robust index management for the gap filling.
 */
function kmer_seeded_nwalign(s1, s2, gapOpen = -10.0, gapExtend = -0.2, matchCost = 1.0, mismatchCost = -0.7, boundaryGapFactor = 10) {
    const K = 25; // High specificity for unique matching
    const minLen = Math.min(s1.length, s2.length);
    
    // Fallback if too short
    if (minLen < K) {
        return affineNWAlign(s1, s2, gapOpen, gapExtend, matchCost, mismatchCost, boundaryGapFactor);
    }

    // 1. Identify Unique K-mers
    // We utilize a single pass map strategy
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

    // Collect matches that appear exactly once in both
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

    // 2. Sort by S1 index to prepare for LIS
    matches.sort((a, b) => a.i - b.i);

    // 3. Longest Increasing Subsequence (LIS) based on S2 index (j)
    // Since we sorted by i, finding the LIS on j gives us the longest chain of non-crossing matches.
    const lis = getLIS(matches);
    
    // 4. Resolve Overlaps and Stitch
    // We treat the LIS anchors as "Truth". We trim them if they overlap to maintain strict monotonicity.
    const safeAnchors = resolveOverlaps(lis);
    
    return stitchAlignments(s1, s2, safeAnchors, [gapOpen, gapExtend, matchCost, mismatchCost, boundaryGapFactor]);
}

/**
 * ROBUST Double DP NW Alignment (Task 2)
 * 
 * Improvements:
 * - Uses a sparse DP with a "Diagonal Penalty" to prefer matches on the same alignment trajectory.
 * - Prunes high-frequency k-mers to prevent performance degradation.
 * - Trims overlapping seeds in the final chain.
 */
function doubleDP_nwalign(s1, s2, gapOpen = -10.0, gapExtend = -0.2, matchCost = 1.0, mismatchCost = -0.7, boundaryGapFactor = 10) {
    const K = 11; // Smaller K for sensitivity, but > 9 to reduce noise
    const minLen = Math.min(s1.length, s2.length);

    if (minLen < K) {
        return affineNWAlign(s1, s2, gapOpen, gapExtend, matchCost, mismatchCost, boundaryGapFactor);
    }

    // 1. Index s1 with a Map
    const s1Map = new Map();
    // Optimization: Step by K/2 to reduce index density while ensuring coverage
    const step = 1; 
    for (let i = 0; i <= s1.length - K; i += step) {
        const sub = s1.substring(i, i + K);
        if (!s1Map.has(sub)) s1Map.set(sub, []);
        s1Map.get(sub).push(i);
    }

    // 2. Find Seeds
    let matches = [];
    const maxHits = 20; // Ignore highly repetitive kmers (e.g. polyA)
    
    for (let j = 0; j <= s2.length - K; j += step) {
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

    // Sort by i (primary) then j (secondary)
    matches.sort((a, b) => (a.i - b.i) || (a.j - b.j));

    // 3. Sparse Chain DP
    // We want to find a path matches[p1] -> matches[p2] ...
    // Penalty is based heavily on *diagonal deviation*.
    const scores = new Float64Array(matches.length);
    const parents = new Int32Array(matches.length).fill(-1);
    const lookback = 100; // How many previous seeds to check

    for (let cur = 0; cur < matches.length; cur++) {
        const mCurr = matches[cur];
        let maxScore = mCurr.len; // Base score
        
        // Search backwards
        const startSearch = Math.max(0, cur - lookback);
        
        for (let prev = cur - 1; prev >= startSearch; prev--) {
            const mPrev = matches[prev];

            // Strict ordering required for time flow
            // Allow slight overlap in indices (will fix later), but start points must differ
            if (mPrev.i < mCurr.i && mPrev.j < mCurr.j) {
                
                // Diagonal Difference:
                // If we stay on the same diagonal, (j - i) is constant.
                // Deviation = |(j2 - i2) - (j1 - i1)|
                const diagDiff = Math.abs((mCurr.j - mCurr.i) - (mPrev.j - mPrev.i));
                
                // Gap distance (Chebyshev distance approximation for gap length)
                const dist = Math.max(mCurr.i - (mPrev.i + mPrev.len), mCurr.j - (mPrev.j + mPrev.len));
                const gapPenalty = (dist > 0) ? (dist * 0.05) : 0; // Small penalty for distance
                
                // Heavy penalty for changing diagonals (indel approximation)
                const diagPenalty = diagDiff * 2.0; 

                const chainScore = scores[prev] + mCurr.len - diagPenalty - gapPenalty;

                if (chainScore > maxScore) {
                    maxScore = chainScore;
                    parents[cur] = prev;
                }
            }
        }
        scores[cur] = maxScore;
    }

    // 4. Backtrack
    let bestIdx = 0;
    let maxVal = -Infinity;
    for(let i=0; i<matches.length; i++){
        if(scores[i] > maxVal){
            maxVal = scores[i];
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

    // 5. Clean and Stitch
    const safeAnchors = resolveOverlaps(chain);
    return stitchAlignments(s1, s2, safeAnchors, [gapOpen, gapExtend, matchCost, mismatchCost, boundaryGapFactor]);
}

// --- Helpers ---

/**
 * Standard O(N log N) Longest Increasing Subsequence algorithm
 * Adapted to work on the 'j' property of match objects (since 'i' is already sorted).
 */
function getLIS(matches) {
    if (matches.length === 0) return [];
    
    // tails[k] stores the index in 'matches' of the smallest tail of all increasing subsequences of length k+1.
    const tails = []; 
    // parent[x] stores the index of the predecessor of matches[x] in the LIS
    const parent = new Int32Array(matches.length).fill(-1);

    for (let i = 0; i < matches.length; i++) {
        const val = matches[i].j;
        
        // Binary search for val in tails
        let left = 0, right = tails.length;
        while (left < right) {
            const mid = (left + right) >>> 1;
            if (matches[tails[mid]].j < val) {
                left = mid + 1;
            } else {
                right = mid;
            }
        }
        
        // If we found a valid extension
        if (left < tails.length) {
            tails[left] = i;
            parent[i] = (left > 0) ? tails[left - 1] : -1;
        } else {
            tails.push(i);
            parent[i] = (tails.length > 1) ? tails[tails.length - 2] : -1;
        }
    }

    // Reconstruct
    const result = [];
    let curr = tails[tails.length - 1];
    while (curr !== -1) {
        result.push(matches[curr]);
        curr = parent[curr];
    }
    return result.reverse();
}

/**
 * Ensures anchors do not physically overlap in a way that reverses time.
 * If Anchor B starts before Anchor A ends, we trim Anchor B's start.
 * If B is completely swallowed, it is removed.
 */
function resolveOverlaps(anchors) {
    if (anchors.length === 0) return [];
    
    const valid = [];
    let prev = anchors[0];
    valid.push({ ...prev }); // Copy to avoid mutation issues

    for (let k = 1; k < anchors.length; k++) {
        let curr = { ...anchors[k] };
        const last = valid[valid.length - 1];

        // Check s1 overlap
        let overlapI = (last.i + last.len) - curr.i;
        // Check s2 overlap
        let overlapJ = (last.j + last.len) - curr.j;
        
        let overlap = Math.max(overlapI, overlapJ);

        if (overlap > 0) {
            // Trim start of current
            curr.i += overlap;
            curr.j += overlap;
            curr.len -= overlap;
        }

        // Only add if valid length remains
        if (curr.len > 0) {
            valid.push(curr);
        }
    }
    return valid;
}

/**
 * Runs the expensive affineNWAlign only on the gaps between anchors.
 * Concatenates the results.
 */
function stitchAlignments(s1, s2, anchors, nwParams) {
    let finalS1 = "";
    let finalS2 = "";
    let idx1 = 0;
    let idx2 = 0;

    for (const anchor of anchors) {
        // Gap Region
        // Ensure we don't pass negative length strings (should be handled by resolveOverlaps, but safety first)
        const gapLen1 = anchor.i - idx1;
        const gapLen2 = anchor.j - idx2;

        if (gapLen1 > 0 || gapLen2 > 0) {
            const sub1 = (gapLen1 > 0) ? s1.substring(idx1, anchor.i) : "";
            const sub2 = (gapLen2 > 0) ? s2.substring(idx2, anchor.j) : "";
            const [aln1, aln2] = affineNWAlign(sub1, sub2, ...nwParams);
            finalS1 += aln1;
            finalS2 += aln2;
        }

        // Anchor Region (Perfect Match)
        const anchorSeq = s1.substring(anchor.i, anchor.i + anchor.len);
        finalS1 += anchorSeq;
        finalS2 += anchorSeq;

        idx1 = anchor.i + anchor.len;
        idx2 = anchor.j + anchor.len;
    }

    // Tail Region
    if (idx1 < s1.length || idx2 < s2.length) {
        const sub1 = s1.substring(idx1);
        const sub2 = s2.substring(idx2);
        const [aln1, aln2] = affineNWAlign(sub1, sub2, ...nwParams);
        finalS1 += aln1;
        finalS2 += aln2;
    }

    return [finalS1, finalS2];
}
