/*!
FrameClean — assign a robust coding frame to an existing nucleotide alignment, despite within-codon indels and rare 1–2 bp “noise insertions”.

Core idea
- Treat alignment columns as either:
  (B) backbone columns that advance a shared coding coordinate, or
  (I) insertion-only columns relative to that shared coordinate (they do not advance the global frame).
- Infer a B/I mask with a small Viterbi DP driven by:
  - a column model (defaults to occupancy-based backbone-vs-insertion likelihood),
  - an indel prior that strongly prefers insertion runs whose lengths are multiples of 3,
  - optional hard/soft constraints from a reference sequence with trusted frame.
- Then assign a phase (0/1/2) to each backbone column; insertion columns get null (⊥).
  - If reference provided: choose the offset that best matches the reference’s phase.
  - Else: choose the offset that maximizes a codon-coherence score (defaults: stop-codon penalty + amino-acid consensus reward).

Outputs
- frameVec: length L array with entries {0,1,2} for backbone columns, null for insertion columns.
- backboneMask: length L boolean array, true for backbone columns.
- Convenience: buildCorrectedAlignment() to produce a new alignment where “noise insertions” do not break the frame:
  - mode="pad" (default): keep insertion columns and optionally add all-gap padding columns so every insertion-run length is a multiple of 3.
  - mode="remove-insertions": drop columns classified as insertions (output sequences will differ from input).

API (high level)
- inferFrame(alignment, options) -> FrameResult
- buildCorrectedAlignment(alignment, frameResult, options) -> { alignment, frameVec, columnMap }

Alignment format
- alignment: array of strings OR array of { id, seq } objects.
  All seq strings must be same length and contain A/C/G/T/- (case-insensitive).

No external deps. ES module + CommonJS compatible.
*/

(function (root, factory) {
  if (typeof module === "object" && typeof module.exports === "object") {
    module.exports = factory();
  } else {
    root.FrameClean = factory();
  }
})(typeof self !== "undefined" ? self : this, function () {
  "use strict";

  // ---------------------------
  // Utilities
  // ---------------------------

  function assert(cond, msg) {
    if (!cond) throw new Error(msg);
  }

  function normalizeAlignment(alignment) {
    assert(Array.isArray(alignment) && alignment.length > 0, "alignment must be a non-empty array");
    let seqs = alignment.map((x, idx) => {
      if (typeof x === "string") return { id: `seq${idx + 1}`, seq: x };
      assert(x && typeof x.seq === "string", "alignment entries must be strings or {id, seq} objects");
      return { id: x.id ?? `seq${idx + 1}`, seq: x.seq };
    });

    const L = seqs[0].seq.length;
    assert(L > 0, "alignment sequences must be non-empty");
    for (const s of seqs) {
      assert(s.seq.length === L, "all alignment sequences must have the same length");
      s.seq = s.seq.toUpperCase();
      assert(/^[ACGTN\-]+$/.test(s.seq), "sequences may only contain A,C,G,T,N,-");
    }
    return { seqs, L, N: seqs.length };
  }

  function getColumnChars(aln, j) {
    const col = new Array(aln.N);
    for (let i = 0; i < aln.N; i++) col[i] = aln.seqs[i].seq[j];
    return col;
  }

  function columnOccupancy(col) {
    let m = 0;
    for (const c of col) if (c !== "-") m++;
    return m;
  }

  function argMax3(a0, a1, a2) {
    if (a0 >= a1 && a0 >= a2) return 0;
    if (a1 >= a0 && a1 >= a2) return 1;
    return 2;
  }

  // ---------------------------
  // Genetic code (standard)
  // ---------------------------

  const STANDARD_CODE = (() => {
    const m = Object.create(null);
    // U -> T. Only canonical codons. N handling occurs elsewhere.
    const pairs = [
      ["TTT","F"],["TTC","F"],["TTA","L"],["TTG","L"],
      ["TCT","S"],["TCC","S"],["TCA","S"],["TCG","S"],
      ["TAT","Y"],["TAC","Y"],["TAA","*"],["TAG","*"],
      ["TGT","C"],["TGC","C"],["TGA","*"],["TGG","W"],

      ["CTT","L"],["CTC","L"],["CTA","L"],["CTG","L"],
      ["CCT","P"],["CCC","P"],["CCA","P"],["CCG","P"],
      ["CAT","H"],["CAC","H"],["CAA","Q"],["CAG","Q"],
      ["CGT","R"],["CGC","R"],["CGA","R"],["CGG","R"],

      ["ATT","I"],["ATC","I"],["ATA","I"],["ATG","M"],
      ["ACT","T"],["ACC","T"],["ACA","T"],["ACG","T"],
      ["AAT","N"],["AAC","N"],["AAA","K"],["AAG","K"],
      ["AGT","S"],["AGC","S"],["AGA","R"],["AGG","R"],

      ["GTT","V"],["GTC","V"],["GTA","V"],["GTG","V"],
      ["GCT","A"],["GCC","A"],["GCA","A"],["GCG","A"],
      ["GAT","D"],["GAC","D"],["GAA","E"],["GAG","E"],
      ["GGT","G"],["GGC","G"],["GGA","G"],["GGG","G"],
    ];
    for (const [c, aa] of pairs) m[c] = aa;
    return m;
  })();

  function translateCodon(codon, codeMap = STANDARD_CODE) {
    // codon: string length 3, containing A/C/G/T. If any N or -, treat as unknown.
    if (codon.length !== 3) return "X";
    if (codon.indexOf("-") !== -1) return "X";
    if (codon.indexOf("N") !== -1) return "X";
    return codeMap[codon] ?? "X";
  }

  // ---------------------------
  // Models (modular components)
  // ---------------------------

  /**
   * Default occupancy-based column model.
   * Assumes:
   *  - Backbone columns have higher non-gap occupancy (pBackbone),
   *  - Insertion columns have lower non-gap occupancy (pInsertion).
   * Uses binomial log-likelihood of occupancy; binomial coefficient cancels in comparisons.
   */
  class OccupancyColumnModel {
    constructor(opts = {}) {
      this.pBackbone = clamp01(opts.pBackbone ?? 0.90);
      this.pInsertion = clamp01(opts.pInsertion ?? 0.05);
      this.extraBackboneBonus = opts.extraBackboneBonus ?? 0.0; // e.g. encourage backbone slightly
      this.extraInsertionBonus = opts.extraInsertionBonus ?? 0.0;
    }

    logLikeBackbone(col) {
      const N = col.length;
      const m = columnOccupancy(col);
      return (
        m * Math.log(this.pBackbone) +
        (N - m) * Math.log(1 - this.pBackbone) +
        this.extraBackboneBonus
      );
    }

    logLikeInsertion(col) {
      const N = col.length;
      const m = columnOccupancy(col);
      return (
        m * Math.log(this.pInsertion) +
        (N - m) * Math.log(1 - this.pInsertion) +
        this.extraInsertionBonus
      );
    }
  }

  function clamp01(x) {
    if (x <= 0) return 1e-12;
    if (x >= 1) return 1 - 1e-12;
    return x;
  }

  /**
   * Default indel prior for backbone/insertion segmentation.
   *
   * - insertionStartPenalty: penalty for entering insertion run (B->I)
   * - insertionExtendPenalty: penalty per insertion column while staying in insertion (I->I)
   * - insertionEndResidPenalty[r]: penalty paid when exiting insertion run (I->B) with residue r = runLength mod 3
   *     r=0 => length multiple of 3 (preferred)
   *     r=1 or 2 => strongly penalized (noise/seq error)
   */
  class Mod3InsertionPrior {
    constructor(opts = {}) {
      this.insertionStartPenalty = opts.insertionStartPenalty ?? -2.0;
      this.insertionExtendPenalty = opts.insertionExtendPenalty ?? -0.2;
      // Defaults: strongly disfavor residue != 0 at exit.
      this.insertionEndResidPenalty = opts.insertionEndResidPenalty ?? [0.0, -6.0, -6.0];
      this.backboneStayBonus = opts.backboneStayBonus ?? 0.0;
      this.backboneAfterInsertionBonus = opts.backboneAfterInsertionBonus ?? 0.0;
      this.endOfAlignmentResidPenalty = opts.endOfAlignmentResidPenalty ?? [0.0, -6.0, -6.0];
    }
  }

  /**
   * Default codon coherence model.
   * - stopCodonPenalty: penalty per stop codon in fully observed codons (per-sequence, per-codon)
   * - aaConsensusWeight: reward proportional to most common amino acid count among sequences (excluding stop and X)
   * - ignoreTerminalCodonStops: if true, do not penalize stop codons in the final codon (useful for genuine terminal stop)
   */
  class SimpleCodonCoherenceModel {
    constructor(opts = {}) {
      this.stopCodonPenalty = opts.stopCodonPenalty ?? 8.0;
      this.aaConsensusWeight = opts.aaConsensusWeight ?? 0.2;
      this.ignoreTerminalCodonStops = opts.ignoreTerminalCodonStops ?? true;
      this.codeMap = opts.codeMap ?? STANDARD_CODE;
    }

    codonScore(codonsPerSeq, codonIndex, codonCount) {
      // codonsPerSeq: array length N containing codon strings (len 3) or null (missing / gapped)
      const aas = [];
      let stops = 0;
      for (const c of codonsPerSeq) {
        if (!c) continue;
        const aa = translateCodon(c, this.codeMap);
        if (aa === "*") stops++;
        aas.push(aa);
      }
      let score = 0.0;
      const isTerminal = codonIndex === codonCount - 1;
      if (!(this.ignoreTerminalCodonStops && isTerminal)) {
        score -= this.stopCodonPenalty * stops;
      }

      // Consensus reward (exclude stop, X)
      const counts = Object.create(null);
      for (const aa of aas) {
        if (aa === "*" || aa === "X") continue;
        counts[aa] = (counts[aa] ?? 0) + 1;
      }
      let best = 0;
      for (const k in counts) if (counts[k] > best) best = counts[k];
      score += this.aaConsensusWeight * best;
      return score;
    }
  }

  // ---------------------------
  // Reference handling
  // ---------------------------

  function resolveReference(aln, refOpt) {
    // refOpt: { id, index, mode, frameOffset }
    if (!refOpt) return null;

    let idx = null;
    if (typeof refOpt.index === "number") idx = refOpt.index;
    if (idx == null && typeof refOpt.id === "string") {
      idx = aln.seqs.findIndex(s => s.id === refOpt.id);
      assert(idx >= 0, `reference id not found: ${refOpt.id}`);
    }
    assert(idx != null && idx >= 0 && idx < aln.N, "invalid reference specification");

    const mode = refOpt.mode ?? "hard"; // "hard" or "soft"
    assert(mode === "hard" || mode === "soft", "reference.mode must be 'hard' or 'soft'");

    // frameOffset: phase at reference ungapped position 1 (0/1/2).
    const frameOffset = refOpt.frameOffset ?? 0;
    assert(frameOffset === 0 || frameOffset === 1 || frameOffset === 2, "reference.frameOffset must be 0,1,2");

    // soft penalties (applied if a column violates the ref-based forced type)
    const softPenalty = refOpt.softPenalty ?? -20.0;

    return { index: idx, mode, frameOffset, softPenalty };
  }

  function forcedColumnTypeFromReference(ref, colChar) {
    // If reference base present => backbone, if gap => insertion.
    // This is a strong heuristic that matches typical “reference defines coordinate” behavior.
    return colChar === "-" ? "I" : "B";
  }

  // ---------------------------
  // Viterbi for backboneMask (B vs I)
  // ---------------------------

  /**
   * Infer backboneMask via Viterbi DP with small state:
   *  - B state: not in insertion run
   *  - I[r] states: in insertion run with residue r = runLength mod 3 (r in {0,1,2})
   *
   * Emissions: columnModel log-likelihoods for B vs I
   * Transitions: prior penalties and mod-3 exit penalties
   * Reference: optional hard/soft constraints per column
   */
  function inferBackboneMask(aln, options) {
    const columnModel = options.columnModel ?? new OccupancyColumnModel(options.columnModelOptions ?? {});
    const prior = options.priorModel ?? new Mod3InsertionPrior(options.priorModelOptions ?? {});
    const ref = resolveReference(aln, options.reference);

    const L = aln.L;

    // State indices:
    // 0: B
    // 1: I0 (insertion, residue 0)
    // 2: I1 (residue 1)
    // 3: I2 (residue 2)
    const S = 4;

    const NEG_INF = -1e300;
    let dpPrev = new Float64Array(S);
    let dpCur = new Float64Array(S);
    for (let s = 0; s < S; s++) dpPrev[s] = NEG_INF;
    dpPrev[0] = 0.0; // start in backbone

    // Backpointers
    // prevState[j][s] and action[j][s] (action: 1 backbone, 0 insertion)
    const prevState = Array.from({ length: L }, () => new Int16Array(S));
    const action = Array.from({ length: L }, () => new Uint8Array(S));

    // Helper: update dpCur for a candidate transition
    function relax(j, sFrom, sTo, act, score) {
      if (score > dpCur[sTo]) {
        dpCur[sTo] = score;
        prevState[j][sTo] = sFrom;
        action[j][sTo] = act;
      }
    }

    for (let j = 0; j < L; j++) {
      for (let s = 0; s < S; s++) dpCur[s] = NEG_INF;

      const col = getColumnChars(aln, j);
      const eB = columnModel.logLikeBackbone(col);
      const eI = columnModel.logLikeInsertion(col);

      // Reference constraint for this column, if any.
      let forced = null; // "B" or "I" or null
      if (ref) forced = forcedColumnTypeFromReference(ref, col[ref.index]);

      for (let sFrom = 0; sFrom < S; sFrom++) {
        const baseScore = dpPrev[sFrom];
        if (baseScore <= NEG_INF / 10) continue;

        const fromIsB = (sFrom === 0);
        const fromIsI = !fromIsB;
        const fromResid = fromIsI ? (sFrom === 1 ? 0 : (sFrom === 2 ? 1 : 2)) : 0;

        // Option 1: choose backbone for this column
        if (!forced || forced === "B") {
          let trans = 0.0;
          let refAdj = 0.0;
          if (ref && forced === "B" && ref.mode === "soft" && col[ref.index] === "-") {
            // contradicts expected I, but forced would be I in that case; this block won't run
            // kept for completeness
            refAdj += ref.softPenalty;
          }

          if (fromIsB) {
            trans += prior.backboneStayBonus;
          } else {
            // exiting insertion run, pay residue penalty
            trans += prior.insertionEndResidPenalty[fromResid] + prior.backboneAfterInsertionBonus;
          }
          relax(j, sFrom, 0, 1, baseScore + trans + eB + refAdj);
        } else if (ref && forced === "B" && ref.mode === "soft") {
          // allow violating forced backbone via insertion, penalize below (handled in insertion option)
        }

        // Option 2: choose insertion for this column
        if (!forced || forced === "I") {
          let trans = 0.0;
          let sTo = 0;

          if (fromIsB) {
            trans += prior.insertionStartPenalty;
            // entering insertion with length 1 => residue 1
            sTo = 2; // I1
          } else {
            trans += prior.insertionExtendPenalty;
            // advance residue
            const nextResid = (fromResid + 1) % 3;
            sTo = (nextResid === 0 ? 1 : (nextResid === 1 ? 2 : 3));
          }

          relax(j, sFrom, sTo, 0, baseScore + trans + eI);
        } else if (ref && forced === "I" && ref.mode === "soft") {
          // allow violating forced insertion via backbone, penalize in backbone option instead
        }

        // Soft reference penalties: if mode=soft, we allow both actions but penalize violations
        // The above "forced" logic blocks the other action when forced != null.
        // For "soft" we want to allow both. Implement by re-running both with penalties.
      }

      // If soft reference, redo transitions allowing both actions with soft penalties.
      if (ref && ref.mode === "soft" && forced) {
        // Recompute dpCur from dpPrev but without forcing, and add penalty if action mismatches forced.
        // Take max between existing forced version and this soft version.
        const dpSoft = new Float64Array(S);
        for (let s = 0; s < S; s++) dpSoft[s] = NEG_INF;
        const prevSoft = new Int16Array(S);
        const actSoft = new Uint8Array(S);

        function relaxSoft(sFrom, sTo, act, score) {
          if (score > dpSoft[sTo]) {
            dpSoft[sTo] = score;
            prevSoft[sTo] = sFrom;
            actSoft[sTo] = act;
          }
        }

        for (let sFrom = 0; sFrom < S; sFrom++) {
          const baseScore = dpPrev[sFrom];
          if (baseScore <= NEG_INF / 10) continue;

          const fromIsB = (sFrom === 0);
          const fromIsI = !fromIsB;
          const fromResid = fromIsI ? (sFrom === 1 ? 0 : (sFrom === 2 ? 1 : 2)) : 0;

          // backbone action, with penalty if forced === "I"
          {
            let trans = 0.0;
            if (fromIsB) trans += prior.backboneStayBonus;
            else trans += prior.insertionEndResidPenalty[fromResid] + prior.backboneAfterInsertionBonus;
            const viol = (forced === "I") ? ref.softPenalty : 0.0;
            relaxSoft(sFrom, 0, 1, baseScore + trans + eB + viol);
          }

          // insertion action, with penalty if forced === "B"
          {
            let trans = 0.0;
            let sTo = 0;
            if (fromIsB) { trans += prior.insertionStartPenalty; sTo = 2; }
            else { trans += prior.insertionExtendPenalty; const nextResid = (fromResid + 1) % 3; sTo = (nextResid === 0 ? 1 : (nextResid === 1 ? 2 : 3)); }
            const viol = (forced === "B") ? ref.softPenalty : 0.0;
            relaxSoft(sFrom, sTo, 0, baseScore + trans + eI + viol);
          }
        }

        // Merge: for each state, if soft version better, overwrite dpCur + backpointers.
        for (let sTo = 0; sTo < S; sTo++) {
          if (dpSoft[sTo] > dpCur[sTo]) {
            dpCur[sTo] = dpSoft[sTo];
            prevState[j][sTo] = prevSoft[sTo];
            action[j][sTo] = actSoft[sTo];
          }
        }
      }

      // Swap
      const tmp = dpPrev; dpPrev = dpCur; dpCur = tmp;
    }

    // Apply end-of-alignment residue penalty if we end in insertion.
    let bestFinal = -1e300;
    let bestState = 0;
    for (let s = 0; s < S; s++) {
      let sc = dpPrev[s];
      if (s !== 0) {
        const resid = (s === 1 ? 0 : (s === 2 ? 1 : 2));
        sc += prior.endOfAlignmentResidPenalty[resid];
      }
      if (sc > bestFinal) { bestFinal = sc; bestState = s; }
    }

    // Backtrack to recover action path (s_j)
    const backboneMask = new Array(L);
    let s = bestState;
    for (let j = L - 1; j >= 0; j--) {
      const act = action[j][s]; // 1 backbone, 0 insertion
      backboneMask[j] = (act === 1);
      s = prevState[j][s];
    }

    return { backboneMask, viterbiScore: bestFinal };
  }

  // ---------------------------
  // Assign phase vector (0/1/2 or null) given backboneMask
  // ---------------------------

  function assignFrames(aln, backboneMask, options) {
    const ref = resolveReference(aln, options.reference);
    const codonModel = options.codonModel ?? new SimpleCodonCoherenceModel(options.codonModelOptions ?? {});
    const L = aln.L;

    // backboneIndex[j] = index among backbone columns (0-based), or -1 if insertion
    const backboneIndex = new Int32Array(L);
    let k = 0;
    for (let j = 0; j < L; j++) {
      if (backboneMask[j]) {
        backboneIndex[j] = k;
        k++;
      } else {
        backboneIndex[j] = -1;
      }
    }
    const backboneCount = k;
    const codonCount = Math.floor(backboneCount / 3);

    // If no codons, trivial
    if (backboneCount === 0) {
      const frameVec = new Array(L).fill(null);
      return { frameVec, offset: 0, backboneIndex, backboneCount, scoresByOffset: [0,0,0] };
    }

    // Score each offset delta in {0,1,2}.
    function scoreOffset(delta) {
      let score = 0.0;

      // Reference agreement (if present)
      if (ref) {
        // Compute reference ungapped position for each column where backboneMask && ref has base.
        let refPos = 0; // 1-based coordinate in effect; we'll increment then compute phase
        // First pass: precompute refPosAtCol[j] for speed if desired; keep simple.
        for (let j = 0; j < L; j++) {
          if (aln.seqs[ref.index].seq[j] !== "-") refPos++;
          if (!backboneMask[j]) continue;
          if (aln.seqs[ref.index].seq[j] === "-") continue;
          const refPhase = ((refPos - 1 + ref.frameOffset) % 3);
          const inferred = ((backboneIndex[j] + delta) % 3);
          score += (refPhase === inferred) ? 1.0 : -1.0; // simple agreement score
        }
      }

      // Codon coherence score (backbone columns grouped by 3)
      for (let t = 0; t < codonCount; t++) {
        const b0 = 3 * t + 0;
        const b1 = 3 * t + 1;
        const b2 = 3 * t + 2;
        // We only score if these are actual backbone positions (they are by construction)
        const cols = backboneColsForBackboneIndices(aln, backboneIndex, b0, b1, b2);
        // Phase offset delta defines where backboneIndex 0 lands, but codon grouping by 3 is independent of delta.
        // delta instead tells you the phase label, not the codon partition, in this design.
        // If you want codon partition to shift with delta, you can shift b0/b1/b2 by delta.
        // Here we *do* shift codon partition with delta: codon boundaries are defined by (backboneIndex + delta) mod 3.
        // That means codon t contains the three backbone indices whose phase labels are 0,1,2.
      }

      // Recompute codon partition that depends on delta:
      // Build list of backbone column indices in alignment order:
      const backboneCols = [];
      backboneCols.length = backboneCount;
      for (let j = 0; j < L; j++) {
        const bi = backboneIndex[j];
        if (bi >= 0) backboneCols[bi] = j;
      }

      // Find codon starts: backboneIndex bi where (bi + delta) mod 3 === 0
      // Then codon bases are bi, bi+1, bi+2.
      const codonStarts = [];
      for (let bi = 0; bi + 2 < backboneCount; bi++) {
        if (((bi + delta) % 3) === 0) codonStarts.push(bi);
      }
      const codonTotal = codonStarts.length;

      for (let ci = 0; ci < codonTotal; ci++) {
        const bi0 = codonStarts[ci];
        const bi1 = bi0 + 1;
        const bi2 = bi0 + 2;
        const j0 = backboneCols[bi0], j1 = backboneCols[bi1], j2 = backboneCols[bi2];

        const codonsPerSeq = new Array(aln.N);
        for (let i = 0; i < aln.N; i++) {
          const c0 = aln.seqs[i].seq[j0];
          const c1 = aln.seqs[i].seq[j1];
          const c2 = aln.seqs[i].seq[j2];
          if (c0 === "-" || c1 === "-" || c2 === "-") { codonsPerSeq[i] = null; continue; }
          if (c0 === "N" || c1 === "N" || c2 === "N") { codonsPerSeq[i] = null; continue; }
          codonsPerSeq[i] = c0 + c1 + c2;
        }
        score += codonModel.codonScore(codonsPerSeq, ci, codonTotal);
      }

      return score;
    }

    const s0 = scoreOffset(0);
    const s1 = scoreOffset(1);
    const s2 = scoreOffset(2);
    const best = argMax3(s0, s1, s2);

    // If reference is provided, allow forcing offset from reference if requested.
    let offset = best;
    if (options.reference && options.reference.forceOffset === true) {
      // Choose offset that maximizes only reference agreement.
      // This is implicit in scoreOffset(); but if codonModel dominates, user may want strict reference.
      // Implement by temporarily zeroing codon component: easiest is a separate function.
      offset = bestReferenceOffset(aln, backboneMask, backboneIndex, resolveReference(aln, options.reference));
    }

    // Build frameVec
    const frameVec = new Array(L);
    for (let j = 0; j < L; j++) {
      const bi = backboneIndex[j];
      frameVec[j] = (bi >= 0) ? ((bi + offset) % 3) : null;
    }

    return {
      frameVec,
      offset,
      backboneIndex,
      backboneCount,
      scoresByOffset: [s0, s1, s2],
    };
  }

  function bestReferenceOffset(aln, backboneMask, backboneIndex, ref) {
    if (!ref) return 0;
    const L = aln.L;
    const scores = [0, 0, 0];
    for (let delta = 0; delta < 3; delta++) {
      let score = 0;
      let refPos = 0;
      for (let j = 0; j < L; j++) {
        if (aln.seqs[ref.index].seq[j] !== "-") refPos++;
        if (!backboneMask[j]) continue;
        if (aln.seqs[ref.index].seq[j] === "-") continue;
        const refPhase = ((refPos - 1 + ref.frameOffset) % 3);
        const inferred = ((backboneIndex[j] + delta) % 3);
        score += (refPhase === inferred) ? 1.0 : -1.0;
      }
      scores[delta] = score;
    }
    return argMax3(scores[0], scores[1], scores[2]);
  }

  function backboneColsForBackboneIndices(aln, backboneIndex, b0, b1, b2) {
    // Not used in final; kept for potential custom models.
    return null;
  }

  // ---------------------------
  // Corrected alignment builder
  // ---------------------------

  /**
   * Build a “corrected-frame” alignment from an input alignment and a FrameResult.
   *
   * Options:
   * - mode: "pad" (default) | "remove-insertions"
   * - padToMultipleOf3: (default true) add all-gap columns at end of each insertion-run so run length % 3 == 0.
   * - noiseMaxSupport: (default 1) insertion columns with occupancy <= noiseMaxSupport are considered "noise".
   * - noisePolicy (pad mode only): "keep" (default) | "mask"
   *      keep: leave inserted bases as-is
   *      mask: replace bases in noise columns with gaps (i.e. delete the noise insertion content but keep columns/padding)
   *
   * Returns:
   * - alignment: array of {id, seq} with new sequences
   * - frameVec: phase labels (0/1/2) for every output column (no nulls). This is the “corrected frame”.
   * - columnMap: array mapping output columns -> input column index (or -1 for padding columns)
   */
  function buildCorrectedAlignment(alignment, frameResult, options = {}) {
    const aln = normalizeAlignment(alignment);
    const { backboneMask, frameVec: inputFrameVec } = frameResult;
    assert(backboneMask && backboneMask.length === aln.L, "frameResult.backboneMask must match alignment length");

    const mode = options.mode ?? "pad";
    assert(mode === "pad" || mode === "remove-insertions", "mode must be 'pad' or 'remove-insertions'");

    const padToMultipleOf3 = options.padToMultipleOf3 ?? true;
    const noiseMaxSupport = options.noiseMaxSupport ?? 1;
    const noisePolicy = options.noisePolicy ?? "keep";
    assert(noisePolicy === "keep" || noisePolicy === "mask", "noisePolicy must be 'keep' or 'mask'");

    // Determine insertion runs in input alignment (contiguous columns where backboneMask[j] == false)
    const runs = [];
    let j = 0;
    while (j < aln.L) {
      const isIns = !backboneMask[j];
      const start = j;
      while (j < aln.L && (!backboneMask[j]) === isIns) j++;
      const end = j; // exclusive
      runs.push({ isIns, start, end });
    }

    // Build output columns as a list of actions:
    // each item: { kind: "orig", j } or { kind: "pad" }
    const outCols = [];
    for (const r of runs) {
      if (!r.isIns) {
        for (let jj = r.start; jj < r.end; jj++) {
          // backbone columns always included
          outCols.push({ kind: "orig", j: jj });
        }
      } else {
        if (mode === "remove-insertions") {
          // skip all insertion columns
          continue;
        }
        // pad mode: include insertion columns, optionally mask noise
        for (let jj = r.start; jj < r.end; jj++) outCols.push({ kind: "orig", j: jj });

        if (padToMultipleOf3) {
          const len = r.end - r.start;
          const rem = len % 3;
          if (rem !== 0) {
            const toAdd = 3 - rem;
            for (let t = 0; t < toAdd; t++) outCols.push({ kind: "pad" });
          }
        }
      }
    }

    const outL = outCols.length;

    // Build output sequences
    const outSeqs = aln.seqs.map(s => ({ id: s.id, seqArr: new Array(outL) }));
    const columnMap = new Int32Array(outL); // -1 for pad
    for (let outJ = 0; outJ < outL; outJ++) {
      const item = outCols[outJ];
      if (item.kind === "pad") {
        columnMap[outJ] = -1;
        for (let i = 0; i < aln.N; i++) outSeqs[i].seqArr[outJ] = "-";
        continue;
      }

      const inJ = item.j;
      columnMap[outJ] = inJ;

      // In pad mode, optionally mask noise insertions:
      // noise defined only for insertion columns (inputFrameVec[inJ] == null OR backboneMask false)
      let maskThisCol = false;
      if (mode === "pad" && !backboneMask[inJ] && noisePolicy === "mask") {
        const col = getColumnChars(aln, inJ);
        const m = columnOccupancy(col);
        if (m <= noiseMaxSupport) maskThisCol = true;
      }

      for (let i = 0; i < aln.N; i++) {
        const c = aln.seqs[i].seq[inJ];
        outSeqs[i].seqArr[outJ] = maskThisCol ? "-" : c;
      }
    }

    // Finalize sequences
    const outAlignment = outSeqs.map(s => ({ id: s.id, seq: s.seqArr.join("") }));

    // Build corrected frame vector for output columns (no nulls):
    // Walk through outCols, assigning phases by:
    // - For columns mapping to backbone input columns: use inputFrameVec (0/1/2)
    // - For insertion/pad columns: assign phases sequentially, but ensure they don't disrupt backbone phases.
    //
    // Approach:
    // - Maintain current phase "p" that always equals the last assigned output phase.
    // - When we hit a backbone-mapped column, we set its phase to the inferred inputFrameVec.
    //   If insertion/pad columns exist before it, and padToMultipleOf3=true, phases should align automatically.
    // - For insertion/pad columns, increment p each column (cyclic), but if this causes mismatch at next backbone,
    //   the mismatch reflects inconsistent inputs; we do not silently “fix” beyond padding.
    const outFrameVec = new Int8Array(outL);
    let lastPhase = 0;
    let haveLast = false;

    for (let outJ = 0; outJ < outL; outJ++) {
      const inJ = columnMap[outJ];
      if (inJ >= 0 && backboneMask[inJ]) {
        const ph = inputFrameVec[inJ];
        // ph should be 0/1/2
        outFrameVec[outJ] = ph;
        lastPhase = ph;
        haveLast = true;
      } else {
        // insertion or pad
        if (!haveLast) {
          // If alignment begins with insertion/pad, start at phase 0 by convention.
          outFrameVec[outJ] = lastPhase;
          lastPhase = (lastPhase + 1) % 3;
          haveLast = true;
        } else {
          lastPhase = (lastPhase + 1) % 3;
          outFrameVec[outJ] = lastPhase;
        }
      }
    }

    return { alignment: outAlignment, frameVec: Array.from(outFrameVec), columnMap: Array.from(columnMap) };
  }

  // ---------------------------
  // High-level API
  // ---------------------------

  /**
   * Infer robust backbone mask + frame assignment for an existing nucleotide alignment.
   *
   * options:
   * - reference: {
   *     id?: string,
   *     index?: number,
   *     mode?: "hard"|"soft",
   *     softPenalty?: number,
   *     frameOffset?: 0|1|2,          // phase at ref ungapped position 1
   *     forceOffset?: boolean         // if true, choose offset purely by reference agreement
   *   }
   *
   * - columnModel: custom object with methods { logLikeBackbone(colChars), logLikeInsertion(colChars) }
   * - columnModelOptions: opts passed to default OccupancyColumnModel if columnModel not supplied
   *
   * - priorModel: custom prior object, or use priorModelOptions for default Mod3InsertionPrior
   * - priorModelOptions: { insertionStartPenalty, insertionExtendPenalty, insertionEndResidPenalty, ... }
   *
   * - codonModel: custom object with method { codonScore(codonsPerSeq, codonIndex, codonCount) }
   * - codonModelOptions: opts passed to default SimpleCodonCoherenceModel if codonModel not supplied
   *
   * returns FrameResult:
   * {
   *   backboneMask: boolean[],
   *   frameVec: (0|1|2|null)[],
   *   offset: 0|1|2,
   *   diagnostics: { viterbiScore, scoresByOffset, backboneCount, insertionCount }
   * }
   */
  function inferFrame(alignment, options = {}) {
    const aln = normalizeAlignment(alignment);

    const { backboneMask, viterbiScore } = inferBackboneMask(aln, options);
    const { frameVec, offset, scoresByOffset, backboneCount } = assignFrames(aln, backboneMask, options);

    const insertionCount = aln.L - backboneCount;

    return {
      backboneMask,
      frameVec,
      offset,
      diagnostics: {
        viterbiScore,
        scoresByOffset,
        backboneCount,
        insertionCount,
      },
    };
  }

  // ---------------------------
  // Exports
  // ---------------------------

  return {
    inferFrame,
    buildCorrectedAlignment,

    // Export default models for ablation/customization
    models: {
      OccupancyColumnModel,
      Mod3InsertionPrior,
      SimpleCodonCoherenceModel,
      STANDARD_CODE,
      translateCodon,
    },

    // Low-level utilities (useful for testing)
    _internal: {
      normalizeAlignment,
      inferBackboneMask,
      assignFrames,
    },
  };
});
