/**
 * DNA Thermodynamics & Tm Calculation Module
 * Implements NN (Nearest Neighbor) model with salt corrections.
 */
const TMCalculator = (function() {
    // NN Delta H (kcal/mol) and Delta S (cal/K/mol) parameters (SantaLucia 1998)
    const NN_PARAMS = {
        'AA': {dH: -7.9, dS: -22.2}, 'TT': {dH: -7.9, dS: -22.2},
        'AT': {dH: -7.2, dS: -20.4}, 'TA': {dH: -7.2, dS: -21.3},
        'CA': {dH: -8.5, dS: -22.7}, 'TG': {dH: -8.5, dS: -22.7},
        'GT': {dH: -8.4, dS: -22.4}, 'AC': {dH: -8.4, dS: -22.4},
        'CT': {dH: -7.8, dS: -21.0}, 'AG': {dH: -7.8, dS: -21.0},
        'GA': {dH: -8.2, dS: -22.2}, 'TC': {dH: -8.2, dS: -22.2},
        'CG': {dH: -10.6, dS: -27.2}, 'GC': {dH: -9.8, dS: -24.4},
        'GG': {dH: -8.0, dS: -19.9}, 'CC': {dH: -8.0, dS: -19.9},
        'sym': {dH: 0, dS: -1.4}
    };

    /**
     * Standard default configurations for different applications.
     */
    const PRESETS = {
        gibson: {
            dnaConc: 5,     // nM (Oligo concentration)
            naConc: 50,     // mM (Monovalent salt)
            mgConc: 10,     // mM (Divalent salt)
            dntpConc: 0.8   // mM (dNTP concentration)
        },
        primer: {
            dnaConc: 200,   // nM (0.2 uM)
            naConc: 50,     // mM
            mgConc: 3,      // mM
            dntpConc: 0.8   // mM
        }
    };

    /**
     * Internal utility to check for perfect self-complementarity.
     */
    function isSelfComplementary(seq) {
        const comp = { A: 'T', T: 'A', C: 'G', G: 'C' };
        let rc = "";
        for (let i = seq.length - 1; i >= 0; i--) {
            const b = seq[i];
            const c = comp[b];
            if (!c) return false;
            rc += c;
        }
        return rc === seq;
    }

    /**
     * Internal utility to compute free Mg2+ after dNTP binding.
     */
    function computeFreeMg(totalMg, totalDntp, Ka) {
        if (totalMg <= 0 || Ka <= 0) return 0;
        const A = Ka * totalDntp - Ka * totalMg + 1;
        const D = A * A + 4 * Ka * totalMg;
        const mg = (-A + Math.sqrt(Math.max(0, D))) / (2 * Ka);
        return Math.min(totalMg, Math.max(0, mg));
    }

    /**
     * Calculates the melting temperature (Tm) of a DNA duplex.
     * 
     * @param {string} seq - The DNA sequence (A, C, G, T only).
     * @param {Object|string} options - Configuration object with dnaConc, naConc, mgConc, dntpConc, 
     *                                   or a preset name ('gibson', 'primer').
     * @returns {number} The calculated Tm in degrees Celsius.
     */
    function calculateTm(seq, options = 'gibson') {
        let config;
        if (typeof options === 'string') {
            config = PRESETS[options] || PRESETS.gibson;
        } else {
            // Start with gibson defaults and override with provided options
            config = Object.assign({}, PRESETS.gibson, options);
        }

        seq = (seq || "").toUpperCase();
        if (seq.length < 2) return 0;
        if (/[^ACGT]/.test(seq)) return 0;

        const R_gas = 1.987; // cal·K^-1·mol^-1
        const K0 = 273.15;
        const Ka_Mg_dNTP = 3.0e4; // M^-1

        // Convert concentrations to Molar for calculation
        const oligoEff = Math.max(1e-30, (config.dnaConc ?? 5) * 1e-9);
        const naConc = Math.max(0, (config.naConc ?? 50) * 1e-3);
        const tMgConc = Math.max(0, (config.mgConc ?? 10) * 1e-3);
        const dntpConc = Math.max(0, (config.dntpConc ?? 0.8) * 1e-3);

        let H_kcal = 0.0;
        let S_calK = 0.0;

        // NN transitions
        for (let i = 0; i < seq.length - 1; i++) {
            const dimer = seq.substring(i, i + 2);
            const p = NN_PARAMS[dimer];
            if (!p) return 0;
            H_kcal += p.dH;
            S_calK += p.dS;
        }

        // Initiation terms
        const end5 = seq[0];
        const end3 = seq[seq.length - 1];
        [end5, end3].forEach(end => {
            if (end === 'G' || end === 'C') {
                H_kcal += 0.1;
                S_calK += -2.8;
            } else if (end === 'A' || end === 'T') {
                H_kcal += 2.3;
                S_calK += 4.1;
            }
        });

        // Symmetry correction
        if (isSelfComplementary(seq)) {
            S_calK += NN_PARAMS.sym.dS;
        }

        const H_cal = 1000.0 * H_kcal;
        
        // Base Tm in 1M NaCl
        const Tm_1M_K = H_cal / (S_calK + R_gas * Math.log(oligoEff));

        const gcCount = (seq.match(/[GC]/g) || []).length;
        const f_GC = gcCount / seq.length;

        const mgFree = computeFreeMg(tMgConc, dntpConc, Ka_Mg_dNTP);

        let useMonovalent = false;
        let useDivalent = false;

        if (naConc > 0) {
            if (mgFree > 0) {
                const R_ratio = Math.sqrt(mgFree) / naConc;
                if (R_ratio < 0.22) useMonovalent = true;
                else useDivalent = true;
            } else {
                useMonovalent = true;
            }
        } else {
            if (mgFree > 0) useDivalent = true;
        }

        // Salt corrections (IDT-style)
        if (useMonovalent) {
            const L = Math.log(naConc);
            const invTm_Na = (1 / Tm_1M_K) + (((4.29 * f_GC - 3.95) * L + 0.940 * (L * L)) * 1e-5);
            return (1 / invTm_Na) - K0;
        }

        if (useDivalent) {
            const N_bp = seq.length;
            const M = Math.log(mgFree);
            let a, d, g;

            if (naConc > 0) {
                const R_ratio = Math.sqrt(mgFree) / naConc;
                if (R_ratio >= 0.22 && R_ratio <= 6.0) {
                    const N = Math.log(naConc);
                    a = 3.92 * (0.843 - 0.352 * Math.sqrt(naConc) * N);
                    d = 1.42 * (1.279 - 4.03e-3 * N - 8.03e-3 * (N * N));
                    g = 8.31 * (0.486 - 0.258 * N + 5.25e-3 * (N * N * N));
                } else {
                    a = 3.92; d = 1.42; g = 8.31;
                }
            } else {
                a = 3.92; d = 1.42; g = 8.31;
            }

            const term = a - 0.911 * M + f_GC * (6.26 + d * M) + (1 / (2 * (N_bp - 1))) * (-48.2 + 52.5 * M + g * (M * M));
            const invTm_Mg = (1 / Tm_1M_K) + term * 1e-5;
            return (1 / invTm_Mg) - K0;
        }

        return Tm_1M_K - K0;
    }

    // Public API
    return {
        calculateTm,
        PRESETS
    };
})();

// Export for various environments
if (typeof module !== 'undefined' && module.exports) {
    module.exports = TMCalculator;
}
