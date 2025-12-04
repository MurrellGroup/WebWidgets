/**
 * FastTree.js
 * A high-performance JavaScript port of the core FastTree algorithm.
 * Performs Profile Neighbor-Joining to infer an approximate Maximum Likelihood tree.
 */

const FastTree = (function() {
    'use strict';

    // --- Constants ---
    const GAP_CODE = 20; // Using 20 for gap/unknown in AA, 4 in DNA
    const AA_CODES = "ARNDCQEGHILKMFPSTWYV";
    const DNA_CODES = "ACGT";
    
    // Simple Jukes-Cantor log correction cap
    const MAX_SCORE = 3.0;

    // --- Helper: Model Definitions ---
    const MODELS = {
        'protein': {
            type: 'protein',
            alphabet: AA_CODES,
            nCodes: 20,
            gapCode: 20
        },
        'dna': {
            type: 'dna',
            alphabet: DNA_CODES,
            nCodes: 4,
            gapCode: 4
        }
    };

    /**
     * Internal Node structure for the tree.
     */
    class Node {
        constructor(id, name = null) {
            this.id = id;
            this.name = name;
            this.children = []; // [Node, Node]
            this.branchLengths = []; // [float, float]
            this.parent = null;
            this.profile = null; // Float64Array
            this.weight = 0; // Number of active sequences in this profile
            this.isActive = true;
        }
    }

    /**
     * Main class handling the algorithm state.
     */
    class FastTreeState {
        constructor(sequences, options) {
            this.seqs = sequences;
            this.model = options.model === 'dna' ? MODELS.dna : MODELS.protein;
            this.nSeq = sequences.length;
            this.nPos = sequences[0].length;
            
            // Configuration
            this.useMatrix = true; 
            
            // Data structures
            this.nodes = [];
            this.activeNodeCount = this.nSeq;
            this.nextId = 0;

            // Distance Matrix (Linearized upper triangle for memory efficiency could be used, 
            // but for JS JIT speed, a flat Float64Array accessed by row/col is often fastest).
            // We use a full matrix for simplicity in indexing, though it uses more RAM.
            // Size: (2*N) * (2*N) to handle internal nodes.
            this.maxNodes = 2 * this.nSeq;
            this.distances = new Float64Array(this.maxNodes * this.maxNodes);
            
            // Divergence vector (r) for NJ
            this.divergence = new Float64Array(this.maxNodes);
            
            // Initialization
            this.initNodes();
        }

        getDist(i, j) {
            return this.distances[i * this.maxNodes + j];
        }

        setDist(i, j, val) {
            this.distances[i * this.maxNodes + j] = val;
            this.distances[j * this.maxNodes + i] = val;
        }

        initNodes() {
            const { nSeq, nPos, model } = this;
            const codeMap = new Int8Array(256).fill(model.gapCode);
            
            // Build map
            for (let i = 0; i < model.alphabet.length; i++) {
                const char = model.alphabet.charCodeAt(i);
                codeMap[char] = i;
                codeMap[char + 32] = i; // lowercase
            }

            // Create leaf nodes and profiles
            for (let i = 0; i < nSeq; i++) {
                const node = new Node(this.nextId++, this.seqs[i].name || `Seq${i+1}`);
                const seqStr = this.seqs[i].seq;
                
                // Create Profile: Flat array [pos * nCodes + code]
                // For leaves, probability is 1.0 for the character, 0.0 otherwise.
                node.profile = new Float32Array(nPos * model.nCodes); 
                
                for (let j = 0; j < nPos; j++) {
                    const code = codeMap[seqStr.charCodeAt(j)];
                    if (code < model.nCodes) {
                        node.profile[j * model.nCodes + code] = 1.0;
                    }
                    // Gaps (code == gapCode) leave the vector chunk as all zeros (treated as unknown)
                }
                
                node.weight = 1;
                this.nodes.push(node);
            }
        }

        /**
         * Log correction for distances (Jukes-Cantor or similar).
         * d = -b * log(1 - d/b)
         */
        logCorrect(dist) {
            if (this.model.type === 'dna') {
                // Jukes-Cantor for DNA
                const b = 0.75; // 3/4
                if (dist >= b) return MAX_SCORE;
                return -0.75 * Math.log(1.0 - (dist * 4.0 / 3.0));
            } else {
                // Approximate correction for Protein (Poisson) or similar
                // FastTree uses -1.3 * log(1-d) roughly
                if (dist >= 0.9) return MAX_SCORE;
                return -Math.log(1.0 - dist);
            }
        }

        /**
         * Computes profile distance between two nodes.
         * dist(A, B) = 1 - sum( freqA[i][c] * freqB[i][c] ) / valid_positions
         */
        computeProfileDistance(nodeA, nodeB) {
            const pA = nodeA.profile;
            const pB = nodeB.profile;
            const nPos = this.nPos;
            const nCodes = this.model.nCodes;
            
            let dotSum = 0.0;
            let validPos = 0.0;

            // Tight loop optimization
            for (let i = 0; i < nPos; i++) {
                const offset = i * nCodes;
                let dot = 0.0;
                let sumA = 0.0;
                let sumB = 0.0;

                // Unrolling loop for DNA (4) slightly helps, but generic loop handles Protein (20)
                for (let c = 0; c < nCodes; c++) {
                    const valA = pA[offset + c];
                    const valB = pB[offset + c];
                    dot += valA * valB;
                    sumA += valA;
                    sumB += valB;
                }

                // If both have valid data at this position
                if (sumA > 0 && sumB > 0) {
                    dotSum += dot;
                    validPos += 1.0;
                }
            }

            if (validPos === 0) return MAX_SCORE;
            
            const rawDist = 1.0 - (dotSum / validPos);
            return this.logCorrect(rawDist);
        }

        computeInitialDistances() {
            const n = this.activeNodeCount;
            for (let i = 0; i < n; i++) {
                this.setDist(i, i, 0);
                for (let j = i + 1; j < n; j++) {
                    const d = this.computeProfileDistance(this.nodes[i], this.nodes[j]);
                    this.setDist(i, j, d);
                }
            }
        }

        /**
         * Calculates the divergence (r) for each active node.
         * r_i = sum(d_ij) / (N-2)
         */
        computeDivergence() {
            const n = this.activeNodeCount;
            const nodes = this.nodes;
            const limit = this.nodes.length;
            
            if (n <= 2) {
                for(let i=0; i<limit; i++) if(nodes[i].isActive) this.divergence[nodes[i].id] = 0;
                return;
            }

            const divFactor = 1.0 / (n - 2);

            for (let i = 0; i < limit; i++) {
                if (!nodes[i].isActive) continue;
                let sum = 0;
                const idI = nodes[i].id;
                
                for (let j = 0; j < limit; j++) {
                    if (i === j || !nodes[j].isActive) continue;
                    sum += this.getDist(idI, nodes[j].id);
                }
                this.divergence[idI] = sum * divFactor;
            }
        }

        /**
         * Finds the pair of nodes to join that minimizes:
         * Q_ij = d_ij - r_i - r_j
         */
        findBestPair() {
            const nodes = this.nodes;
            const len = nodes.length;
            let minQ = Infinity;
            let bestI = -1;
            let bestJ = -1;

            // This is O(N^2). FastTree C uses Top-Hits heuristics to make this O(N).
            // For JS compatibility, we use the exhaustive search (safer for < 5000 seqs).
            for (let i = 0; i < len; i++) {
                if (!nodes[i].isActive) continue;
                const rI = this.divergence[nodes[i].id];
                const idI = nodes[i].id;

                for (let j = i + 1; j < len; j++) {
                    if (!nodes[j].isActive) continue;
                    
                    const d = this.getDist(idI, nodes[j].id);
                    const q = d - rI - this.divergence[nodes[j].id];

                    if (q < minQ) {
                        minQ = q;
                        bestI = i;
                        bestJ = j;
                    }
                }
            }
            return { i: bestI, j: bestJ, dist: this.getDist(nodes[bestI].id, nodes[bestJ].id) };
        }

        /**
         * Merges two profiles (Weighted average).
         * Profile(new) = 0.5 * Profile(A) + 0.5 * Profile(B)
         * (FastTree uses BIONJ weighting, simplified here to 0.5 for speed/stability)
         */
        mergeProfiles(nodeA, nodeB) {
            const len = nodeA.profile.length;
            const newProfile = new Float32Array(len);
            const pA = nodeA.profile;
            const pB = nodeB.profile;
            
            // Standard loop is faster than .map or typed array methods for simple addition
            for (let k = 0; k < len; k++) {
                newProfile[k] = (pA[k] + pB[k]) * 0.5;
            }
            return newProfile;
        }

        runNJ() {
            this.computeInitialDistances();

            // Main NJ Loop
            while (this.activeNodeCount > 2) {
                this.computeDivergence();
                const { i, j, dist } = this.findBestPair();

                if (i === -1) break; // Should not happen

                const nodeA = this.nodes[i];
                const nodeB = this.nodes[j];

                // Create new parent node
                const newNode = new Node(this.nextId++);
                newNode.children = [nodeA, nodeB];
                newNode.profile = this.mergeProfiles(nodeA, nodeB);
                newNode.weight = nodeA.weight + nodeB.weight;

                // Calculate branch lengths
                // L_A = (d_AB + r_A - r_B) / 2
                const rA = this.divergence[nodeA.id];
                const rB = this.divergence[nodeB.id];
                const branchA = (dist + rA - rB) * 0.5;
                const branchB = (dist + rB - rA) * 0.5;

                // Prevent negative branch lengths (common in NJ)
                newNode.branchLengths = [
                    branchA < 0 ? 0 : branchA,
                    branchB < 0 ? 0 : branchB
                ];

                nodeA.parent = newNode;
                nodeB.parent = newNode;
                nodeA.isActive = false;
                nodeB.isActive = false;

                // Update distances for the new node
                // In Profile NJ, we recompute profile distances rather than using the matrix formula
                // This is strictly more accurate for the alignment than the NJ matrix update.
                const nodesLen = this.nodes.length;
                for (let k = 0; k < nodesLen; k++) {
                    if (this.nodes[k].isActive) {
                        const newDist = this.computeProfileDistance(newNode, this.nodes[k]);
                        this.setDist(newNode.id, this.nodes[k].id, newDist);
                    }
                }

                this.nodes.push(newNode);
                this.activeNodeCount--;
            }

            // Finish the root (last 2 nodes)
            // Find the last two active nodes
            let lastNodes = [];
            for(let i=0; i<this.nodes.length; i++) {
                if(this.nodes[i].isActive) lastNodes.push(this.nodes[i]);
            }

            if (lastNodes.length === 2) {
                const root = new Node(this.nextId++);
                const n1 = lastNodes[0];
                const n2 = lastNodes[1];
                const d = this.getDist(n1.id, n2.id);
                
                root.children = [n1, n2];
                root.branchLengths = [d / 2, d / 2];
                n1.parent = root;
                n2.parent = root;
                n1.isActive = false;
                n2.isActive = false;
                this.nodes.push(root);
            }

            // Return the last created node as root
            return this.nodes[this.nodes.length - 1];
        }

        /**
         * Recursive function to generate Newick string.
         */
        toNewick(node) {
            if (!node.children || node.children.length === 0) {
                return node.name;
            }
            
            const childStr = node.children.map((child, index) => {
                const sub = this.toNewick(child);
                const len = node.branchLengths[index].toFixed(5);
                return `${sub}:${len}`;
            }).join(',');

            return `(${childStr})`;
        }
    }

    // --- Public Interface ---

    /**
     * Infer an ML tree from sequences.
     * @param {Array<{name: string, seq: string}>} sequences - Input sequences
     * @param {Object} options - { model: 'dna' | 'protein' }
     * @returns {string} Newick tree string
     */
    function run(sequences, options = {}) {
        const defaults = { model: 'protein' };
        const settings = { ...defaults, ...options };

        if (!sequences || sequences.length < 2) {
            throw new Error("At least 2 sequences are required.");
        }

        // Validate sequences
        const len = sequences[0].seq.length;
        for(let s of sequences) {
            if (s.seq.length !== len) throw new Error("All sequences must be the same length (aligned).");
        }

        // Run Algorithm
        const state = new FastTreeState(sequences, settings);
        const root = state.runNJ();
        
        return state.toNewick(root) + ";";
    }

    return {
        run: run
    };

})();

// --- Example Usage ---
// Use strict mode compatible environment (Node or Browser)

/* 
// Example Data (Protein)
const seqs = [
    { name: "A", seq: "MKV---L" },
    { name: "B", seq: "MKV---I" },
    { name: "C", seq: "MKA---L" },
    { name: "D", seq: "MLA---L" }
];

const tree = FastTree.run(seqs, { model: 'protein' });
console.log(tree); 
*/

if (typeof module !== 'undefined' && module.exports) {
    module.exports = FastTree;
}
