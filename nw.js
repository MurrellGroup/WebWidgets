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
