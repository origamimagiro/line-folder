import { M } from "./flatfolder/math.js";
import { X } from "./flatfolder/conversion.js";
import { SVG } from "./flatfolder/svg.js";
import { NOTE } from "./flatfolder/note.js";
import { SOLVER } from "./flatfolder/solver.js";
import { CON } from "./flatfolder/constraints.js";

export const TYPE_LABEL = [
    "Pureland",
    "Inside Reverse",   "Outside Reverse",  "Mixed Reverse",
    "Open Sink",        "Closed Sink",      "Mixed Sink",
    "Complex",
];
export const TYPE = Object.fromEntries(TYPE_LABEL
    .map(t => t.toUpperCase().replaceAll(" ", "_"))
    .map((t, i) => [t, i]));

export const COMP = {
    augment_FOLD_FO: (FOLD) => {
        const {V, FV, EV, EF, FE, Ff, FO, FOO, FM} = FOLD;
        [FOLD.H, FOLD.EA] = COMP.FO_Ff_EF_2_H_EA(FO, Ff, EF);
        FOLD.Vf = X.V_FV_EV_EA_2_Vf_Ff(V, FV, EV, FOLD.EA)[0];
        if (M.polygon_area2(M.expand(FV[0], FOLD.Vf)) < 0) {
            FOLD.Vf = FOLD.Vf.map(v => M.refY(v));
        }
        const v0 = FOLD.Vf[0];
        FOLD.Vf = FOLD.Vf.map(p => M.sub(p, v0));
        const [c1, s1] = FOLD.Vf[1];
        FOLD.Vf = FOLD.Vf.map(p => M.rotate_cos_sin(p, c1, -s1));
        FOLD.Vf = FOLD.Vf.map(p => M.rotate_cos_sin(p, 0, 1));
        FOLD.Vf = COMP.V_P_transform(FOLD.Vf, FOLD.Vf, false);
        [FOLD.L, FOLD.LL] = (() => { // linearize state
            const {H} = FOLD;
            const n = Ff.length;
            const Adj = Array(n).fill(0).map(() => []);
            for (const [k, o] of H) {
                if (o != -1) { continue; }
                const [f1, f2] = M.decode(k);
                Adj[f1].push(f2);
            }
            const L = [];
            const seen = Array(n).fill(false);
            const dfs = (i) => {
                if (seen[i]) { return; }
                seen[i] = true;
                for (const j of Adj[i]) { dfs(j); }
                L.push(i);
            };
            for (let i = 0; i < n; ++i) { dfs(i); }
            L.reverse();
            console.assert(L.length == n);
            const idx_map = Array(n).fill(undefined);
            for (let i = 0; i < n; ++i) { idx_map[L[i]] = i; }
            for (const [s, o] of H) {
                if (o != -1) { continue; }
                const [f1, f2] = M.decode(s);
                if (idx_map[f1] <= idx_map[f2]) { continue; }
                return [undefined, undefined]; // cycle
            }
            for (let i = 0; i < n; ++i) { seen[i] = false; }
            const layers = [];
            for (let i = 0; i < n; ++i) {
                const fi = L[i];
                if (seen[fi]) { continue; }
                seen[fi] = true;
                const layer = [fi];
                const Adj_set = new Set();
                for (const fj of Adj[fi]) { Adj_set.add(fj); }
                for (let j = i + 1; j < L.length; ++j) {
                    const fj = L[j];
                    if (seen[fj]) { continue; }
                    if (!Adj_set.has(fj)) {
                        seen[fj] = true;
                        layer.push(fj);
                    }
                    for (const fk of Adj[fj]) { Adj_set.add(fk); }
                }
                layers.push(layer);
            }
            return [L, layers];
        })();
        if (FOO == undefined) { return; }
        const [type, RF, FR, EC, V_sink, V_border] = COMP.classify(
            V, EV, EF, FE, Ff, FM, FO, FOO);
        FOLD.type = type;
        FOLD.RF = RF;
        FOLD.FR = FR;
        FOLD.EC = EC;
        FOLD.V_sink = V_sink;
        FOLD.V_border = V_border;
    },
    V_P_transform: (V, P, flip) => {
        // rescales V based on normalizing P into a [0, 1] bounding box
        const [p_min, p_max] = M.bounding_box(P);
        const [x_diff, y_diff] = M.sub(p_max, p_min);
        const is_tall = (x_diff < y_diff);
        const diff = is_tall ? y_diff : x_diff;
        const m = [0.5, 0.5];
        const off = M.sub(m, M.div([x_diff, y_diff], 2*diff));
        return V.map(p => M.add(M.div(M.sub(p, p_min), diff), off))
                .map(p => (flip ? M.add(M.refX(M.sub(p, m)), m) : p));
    },
    V_FV_2_FOLD_CELL: (V, FV) => {
        const Ff = FV.map(fV => (M.polygon_area2(fV.map(i => V[i])) < 0));
        const EV_set = new Set();
        for (const fV of FV) {
            let i = fV.length - 1;
            for (let j = 0; j < fV.length; ++j) {
                EV_set.add(M.encode_order_pair([fV[i], fV[j]]));
                i = j;
            }
        }
        const EV = Array.from(EV_set).sort().map(k => M.decode(k));
        const [VV, _] = X.V_EV_2_VV_FV(V, EV);
        const [EF, FE] = X.EV_FV_2_EF_FE(EV, FV);
        const L = EV.map(vs => vs.map(i => V[i]));
        NOTE.time("Constructing points and segments from edges");
        const [P, SP, SE, eps_i] = X.L_2_V_EV_EL(L);
        const eps = M.min_line_length(L)/(2**eps_i);
        NOTE.annotate(P, "points_coords");
        NOTE.annotate(SP, "segments_points");
        NOTE.annotate(SE, "segments_edges");
        NOTE.lap();
        NOTE.time("Constructing cells from segments");
        const [PP,CP] = X.V_EV_2_VV_FV(P, SP);
        NOTE.annotate(CP, "cells_points");
        NOTE.lap();
        NOTE.time("Computing segments_cells");
        const [SC, CS] = X.EV_FV_2_EF_FE(SP, CP);
        NOTE.annotate(SC, "segments_cells");
        NOTE.annotate(CS, "cells_segments");
        NOTE.lap();
        NOTE.time("Making face-cell maps");
        const [CF, FC] = X.EF_FV_P_SP_SE_CP_SC_2_CF_FC(EF, FV, P, SP, SE, CP, SC);
        const BF = X.EF_SP_SE_CP_CF_2_BF(EF, SP, SE, CP, CF);
        const BI = new Map();
        for (const [i, F] of BF.entries()) { BI.set(F, i); }
        NOTE.annotate(BF, "variables_faces");
        NOTE.lap();
        const fixed = false;
        const FOLD = {V, VV, FV, EV, EF, FE, Ff, eps, fixed};
        const CELL = {P, SP, SE, PP, CP, CS, SC, CF, FC, BF, BI};
        return [FOLD, CELL];
    },
    filter_clicked_and_reflect: (V, FV_, FV, F_map_, FG, FO, clicked, line, FR_) => {
        // F_map_ maps uncut faces to their face indices after cutting
        // FV_ is old uncut faces, while FV is new cut faces
        // FO is old orders for uncut faces
        // clicked is set mapping group indices to Bool whether clicked
        const F_map = F_map_.map(() => []);
        // map only clicked regions
        const FVx = []; // new vertices per face
        const FM = [];  // boolean whether face moves
        const FR = [];
        for (let fi = 0; fi < F_map_.length; ++fi) {
            const F = F_map_[fi];
            if (F.length == 2) { // face was cut
                const [f1, f2] = F;
                const g1 = clicked.has(FG[f1]);
                const g2 = clicked.has(FG[f2]);
                if (g1 == g2) {  // same group, face was not cut
                    F_map[fi].push(FVx.length);
                    FVx.push(FV_[fi]);
                    FM.push(g1);
                    FR.push(FR_[f1]);
                } else {         // face was cut, so split them
                    F_map[fi].push(FVx.length);
                    FVx.push(FV[f1]);
                    FM.push(g1);
                    FR.push(FR_[f1]);
                    F_map[fi].push(FVx.length);
                    FVx.push(FV[f2]);
                    FM.push(g2);
                    FR.push(FR_[f2]);
                }
            } else if (F.length == 1) { // face was not cut
                const g1 = clicked.has(FG[F[0]]);
                F_map[fi].push(FVx.length);
                FVx.push(FV_[fi]);
                FM.push(g1);
                FR.push(FR_[F[0]]);
            } else {
                throw new Error("malformed FV");
            }
        }
        const [Vx, FVy] = (() => { // remove unused vertices
            const V_ = [];
            const n = V.length;
            const Vuse = Array(n).fill(false);
            const V_map = Array(n).fill(undefined);
            for (const F of FVx) {
                for (const i of F) {
                    Vuse[i] = true;
                }
            }
            for (let i = 0; i < n; ++i) {
                if (Vuse[i]) {
                    V_map[i] = V_.length;
                    V_.push(V[i]);
                }
            }
            const FV_ = FVx.map(F => F.map(i => V_map[i]));
            return [V_, FV_];
        })();
        const Vy = Vx.map(() => undefined);
        const [u, d] = line;
        for (let fi = 0; fi < FVy.length; ++fi) {
            const F = FVy[fi];
            if (!FM[fi]) { continue; }
            for (const vi of F) {
                if (Vy[vi] != undefined) { continue; }
                const v = Vx[vi];
                const d2 = M.dot(v, u) - d;
                Vy[vi] = M.sub(v, M.mul(u, 2*d2));
            }
        }
        for (let vi = 0; vi < Vy.length; ++vi) {
            if (Vy[vi] == undefined) {
                Vy[vi] = Vx[vi];
            }
        }
        const FOO = []; // old order
        for (const [f, g, o] of FO) {
            for (const f_ of F_map[f]) {
                for (const g_ of F_map[g]) {
                    const pair = M.encode_order_pair([f_, g_]);
                    FOO.push([f_, g_, o]);
                }
            }
        }
        return [Vx, Vy, FVy, FM, FOO, F_map, FR];
    },
    solve: (FOLD, CELL, BA0) => {
        const {EF} = FOLD;
        const {BF, BI, SE, FC, CF, SC} = CELL;
        const BT = X.BF_BI_EF_SE_CF_SC_2_BT(BF, BI, EF, SE, CF, SC);
        const BTn = [0, 0, 0];
        for (const bT of BT) {
            for (let i = 0; i < 3; ++i) { BTn[i] += bT[i].length; }
        }
        for (const [i, d] of [[0, 6], [1, 2], [2, 2]]) { BTn[i] /= d; }
        NOTE.log(`   - Found ${BTn[0]} taco-taco`);
        NOTE.log(`   - Found ${BTn[1]} taco-tortilla`);
        NOTE.log(`   - Found ${BTn[2]} tortilla-tortilla`);
        NOTE.lap();
        NOTE.time("Computing taco-tortilla implied transitivity");
        const CC = X.FC_BF_BI_BT_2_CC(FC, BF, BI, BT);
        NOTE.lap();
        NOTE.time("*** Computing states ***");
        const trans_count = {all: 0, reduced: 0};
        const BA = SOLVER.initial_assignment(BA0, BF, BT, BI,
            FC, CF, CC, trans_count);
        if ((BA.length == 3) && (BA[1].length != undefined)) {
            const [type, F, E] = BA;
            const str = `Unable to resolve ${CON.names[type]} on faces [${F}]`;
            NOTE.log(`   - ${str}`);
            NOTE.log(`   - Faces participating in conflict: [${E}]`);
            NOTE.end();
            return [undefined, undefined];
        }
        NOTE.annotate(BA.map((_, i) => i).filter(i => BA[i] != 0),
            "initially assignable variables");
        NOTE.lap();
        NOTE.time("Finding unassigned components");
        const GB = SOLVER.get_components(BI, BF, BT, BA, FC, CF, CC, trans_count);
        NOTE.count(GB.length - 1, "unassigned components");
        NOTE.log(`   - Found ${trans_count.reduced/3} reduced transitivity`);
        NOTE.log(`   - Found ${trans_count.all/3} total transitivity`);
        NOTE.lap();
        const GA = SOLVER.solve(BI, BF, BT, BA, GB, FC, CF, CC, Infinity);
        return (GA.length == undefined) ? [undefined, undefined] : [GB, GA];
    },
    FO_Ff_EF_2_H_EA: (FO, Ff, EF) => {
        const H = new Map();
        for (const [i, j, o] of FO) {
            const o_ = Ff[j] ? -o : o;
            H.set(M.encode([i, j]),  o_);
            H.set(M.encode([j, i]), -o_);
        }
        const EA = EF.map(F => {
            if (F.length != 2) { return "B"; }
            const [f, g] = F;
            const o = H.get(M.encode(F));
            if (o == undefined) { return "F"; }
            return ((o > 0) == Ff[g]) ? "V" : "M";
        });
        return [H, EA];
    },
    FM_H1_H2_2_HC: (FM, H1, H2) => {
        const HC = new Set();
        for (const [k, o2] of H2) {
            const [f, g] = M.decode(k);
            if ((!FM[f]) || (!FM[g])) { continue; }
            const o1 = H1.get(k);
            if (o1 == ((FM[f] == FM[g]) ? -o2 : o2)) { continue; }
            HC.add(k);
        }
        return HC;
    },
    EF_FM_HC_2_FR_RF: (EF, FM, HC) => {
        const FF = FM.map(() => []);
        for (const F of EF) {
            if (F.length != 2) { continue; }
            const [f, g] = F;
            if (!FM[f] || !FM[g]) { continue; }
            FF[f].push(g);
            FF[g].push(f);
        }
        const FR = FM.map(() => undefined); // faces to initial coloring
        const RF = [];                      // colors to array of faces
        for (let i = 0; i < FM.length; ++i) {
            if (!FM[i]) { continue; } // only moving faces
            if (FR[i] != undefined) { continue; } // already labelled
            const ri = RF.length;
            FR[i] = ri;
            const rF = [i];
            RF.push(rF);
            const Q = [i]; // bfs on component
            let j = 0;
            while (j < Q.length) {
                const f = Q[j]; ++j;
                const v = FR[f];
                for (const g of FF[f]) {
                    if (
                        (FR[g] == undefined) &&
                        !HC.has(M.encode_order_pair([f, g]))
                    ) {
                        Q.push(g);
                        FR[g] = ri;
                        rF.push(g);
                    }
                }
            }
        }
        return [FR, RF];
    },
    FM_FR_H2_HC_2_uniform: (FM, FR, H2, HC) => {
        /*
         *  First, we need to ensure that the layer order within
         *  each color layer is consistent (does not change between states).
         *  We call this property "uniformity": pairs of faces in the same
         *  layer do not change their order from the start to the end state.
         *  (linear in the number of orders among moving faces)
         */
        let uniform = true;
        for (const [k, o] of H2) {
            const [i, j] = M.decode(k);
            if ((!FM[i]) || (!FM[j])) { continue; }
            if (FR[i] != FR[j]) { continue; }
            if (HC.has(k)) {
                uniform = false;
                break;
            }
        }
        return uniform;
    },
    FR_RF_H1_H2_HC_2_separable: (FR, RF, H1, H2, HC) => {
        /*
         *  Next, we check whether pairs of layers interleave,
         *  in either the start or end states. If they do, they must be merged.
         *  Merging two layers is disallowed if doing so would violate
         *  uniformity. If two interleaving layers cannot be merged, the fold
         *  is "non-separable".
         *  Otherwise it is "separable".
         *  We compute in the following way:
         *      - We process one layer at a time.
         *      - We maintian the invariant that the set of layers processed
         *        so far is uniform and separable.
         *      - When processing a layer Ri, we check to see whether it
         *        interleaves an existing layer Rj.
         *        If it does, we reassign the layer of all faces in
         *        Ri to Rj, and reprocess the combination
         */
        let RF_ = RF.map(F => F.map(f => f));
        const RO1 = new Map(); // color layer orders in start state
        const RO2 = new Map(); // color layer orders in end state
        let separable = true;
        for (let Ri = 1; Ri < RF.length; ++Ri) {
            const Fi = RF_[Ri];
            const RiO1 = new Map();
            const RiO2 = new Map();
            // check separability of layer Ri
            // with each previously processed layer
            const Rmerge = new Set();
            for (let Rj = 0; Rj < Ri; ++Rj) {
                const Fj = RF_[Rj];
                let o = undefined;
                for (const fj of Fj) {
                    for (const fi of Fi) {
                        const kji = M.encode([fj, fi]);
                        const oji1 = H1.get(kji);
                        const oji2 = H2.get(kji);
                        if (oji1 == undefined) { continue; }
                        const Rji = M.encode([Rj, Ri]);
                        const oji1_ = RiO1.get(Rji);
                        const oji2_ = RiO2.get(Rji);
                        if (oji1_ == undefined) {
                            const Rij = M.encode([Ri, Rj]);
                            RiO1.set(Rji,  oji1);
                            RiO1.set(Rij, -oji1);
                            RiO2.set(Rji,  oji2);
                            RiO2.set(Rij, -oji2);
                        } else if ((oji1 != oji1_) || (oji2 != oji2_)) {
                            Rmerge.add(Rj);
                        }
                    }
                }
            }
            // merge each interleaved layer into Ri
            // note that after merging, the new layer might newly
            // interleave existing layers, so we reprocess again
            // to check for continued merging
            if (Rmerge.size > 0) {
                for (const Rj of Rmerge) {
                    // remove existing layer orders for Rmerge
                    for (let Rk = 0; Rk < Ri; ++Rk) {
                        if (Rk == Rj) { continue; }
                        RO1.delete(M.encode([Rj, Rk]));
                        RO1.delete(M.encode([Rk, Rj]));
                        RO2.delete(M.encode([Rj, Rk]));
                        RO2.delete(M.encode([Rk, Rj]));
                    }
                    // reassign layer to Ri
                    const Fj = RF_[Rj];
                    while (Fj.length > 0) {
                        const f = Fj.pop();
                        Fi.push(f);
                    }
                }
                // check whether merged layer Ri is uniform
                for (const f of Fi) {
                    for (const g of Fi) {
                        if (f == g) { continue; }
                        if (HC.has(M.encode_order_pair([f, g]))) {
                            separable = false;
                            break;
                        }
                    }
                    if (!separable) { break; }
                }
                --Ri;
            } else {
                for (const [k, v] of RiO1) { RO1.set(k, v); }
                for (const [k, v] of RiO2) { RO2.set(k, v); }
            }
            if (!separable) { break; }
        }
        if (!separable) { return false; }
        RF.length = 0;
        for (const F of RF_) {
            if (F.length == 0) { continue; }
            RF.push(F);
            const r = RF.length - 1;
            for (const f of F) { FR[f] = r; }
        }
        /*
         *  If separable, we check to see whether any pairs of colors can be
         *  merged without losing uniformity and merge them. Here we only merge
         *  if nothing else exists between the layers (in either the start or
         *  end states), so unlike merging interleaved layers merging separated
         *  layers cannot yield additional merging.
         */
        RF_ = RF.map(F => F.map(f => f));
        for (let Ri = 1; Ri < RF_.length; ++Ri) {
            const Fi = RF_[Ri];
            if (Fi.length == 0) { continue; }
            for (let Rj = 0; Rj < Ri; ++Rj) {
                const Fj = RF_[Rj]
                if (Fj.length == 0) { continue; }
                let adjacent = true;
                for (let Rk = 0; Rk < RF_.length; ++Rk) {
                    if ((Rk == Ri) || (Rk == Rj)) { continue; }
                    const Oik1 = RO1.get(M.encode([Ri, Rk]));
                    const Ojk1 = RO1.get(M.encode([Rj, Rk]));
                    const Oik2 = RO2.get(M.encode([Ri, Rk]));
                    const Ojk2 = RO2.get(M.encode([Rj, Rk]));
                    if ((Oik1 == undefined) ||
                        (Ojk1 == undefined) ||
                        ((Oik1 == Ojk1) && (Oik2 == Ojk2))
                    ) { continue; }
                    adjacent = false;
                    break;
                }
                if (!adjacent) { continue; }
                for (let fk = 0; fk < FR.length; ++fk) {
                    if (FR[fk] != undefined) { continue; }
                    for (const fj of Fj) {
                        const ojk = H1.get(M.encode([fj, fk]));
                        if (ojk == undefined) { continue; } // no overlap
                        for (const fi of Fi) {
                            const oik = H1.get(M.encode([fi, fk]));
                            if (oik == undefined) { continue; } // no overlap
                            if (oik == ojk) { continue; } // match
                            adjacent = false;
                            break;
                        }
                        if (!adjacent) { break; }
                    }
                    if (!adjacent) { break; }
                }
                if (!adjacent) { continue; }
                let merge = false;
                for (const fi of Fi) {
                    for (const fj of Fj) {
                        const k = M.encode([fi, fj]);
                        if (H1.has(k)) {
                            merge = true;
                            break;
                        }
                    }
                    if (merge) { break; }
                }
                for (const fi of Fi) {
                    for (const fj of Fj) {
                        const k = M.encode([fi, fj]);
                        if (HC.has(k)) {
                            merge = false;
                            break;
                        }
                    }
                    if (!merge) { break; }
                }
                if (!merge) { continue; }
                for (let Rk = 0; Rk < RF_.length; ++Rk) {
                    if (Rk == Ri) { continue; }
                    if (Rk == Rj) { continue; }
                    const Kjk = M.encode([Rj, Rk]);
                    const Kkj = M.encode([Rk, Rj]);
                    const Ojk = RO1.get(Kjk);
                    const Okj = RO1.get(Kkj);
                    RO1.delete(Kjk);
                    RO1.delete(Kkj);
                    if (Ojk != undefined) {
                        RO1.set(M.encode([Ri, Rk]), Ojk);
                    }
                    if (Okj != undefined) {
                        RO1.set(M.encode([Rk, Ri]), Okj);
                    }
                }
                while (Fj.length > 0) {
                    const f = Fj.pop();
                    Fi.push(f);
                }
                // no need to repeat because old layers already
                // maximally merged
            }
        }
        RF.length = 0;
        for (const F of RF_) {
            if (F.length == 0) { continue; }
            RF.push(F);
            const r = RF.length - 1;
            for (const f of F) { FR[f] = r; }
        }
        return true;
    },
    FM_FE_EF_RF_2_local: (FM, FE, EF, RF) => {
        /*
         *  Require each layer to be local
         *  (adjacent to boundary or non-moving face) O(n)
         */
        const non_local_layers = [];
        for (let i = 0; i < RF.length; ++i) {
            const rF = RF[i];
            let local = false;
            for (const f of rF) {
                for (const e of FE[f]) {
                    const fE = EF[e];
                    if (fE.length == 1) { continue; }
                    console.assert(fE.length == 2);
                    const [a, b] = fE;
                    const g = (a == f) ? b : a;
                    if (!FM[g]) { local = true; }
                    if (local) { break; }
                }
                if (local) { break; }
            }
            if (!local) {
                non_local_layers.push(i);
            }
        }
        return (non_local_layers.length == 0);
    },
    FM_FR_RF_EF_2_connected: (FM, FR, RF, EF) => {
        /*
         *  Require layer adjacencies to be connected
         */
        const adj = FR.map(() => new Set());
        const adj_edges = new Map();
        for (let i = 0; i < EF.length; ++i) {
            const [fi, fj] = EF[i];
            if ((!FM[fi]) || (!FM[fj])) { continue; }
            const ri = FR[fi];
            const rj = FR[fj];
            if (ri == rj) { continue; }
            adj[ri].add(rj);
            adj[rj].add(ri);
            const k = M.encode_order_pair([ri, rj]);
            if (!adj_edges.has(k)) { adj_edges.set(k, []); }
            const Eij = adj_edges.get(k);
            Eij.push(i);
        }
        const Q = [0];
        const seen = new Set([0]);
        for (let qi = 0; qi < Q.length; ++qi) {
            const ri = Q[qi];
            for (const rj of adj[ri]) {
                if (seen.has(rj)) { continue; }
                seen.add(rj);
                Q.push(rj);
            }
        }
        return (seen.size == RF.length);
    },
    classify: (V, EV, EF, FE, Ff1, FM, FO1, FO2) => {
        const Ff2 = Ff1.map((flip, i) => FM[i] ? !flip : flip);
        const [H1, EA1] = COMP.FO_Ff_EF_2_H_EA(FO1, Ff1, EF);
        const [H2, EA2] = COMP.FO_Ff_EF_2_H_EA(FO2, Ff2, EF);
        const EC = EA1.map((a, i) => (a != EA2[i]));
        const HC = COMP.FM_H1_H2_2_HC(FM, H1, H2);
        const [FR, RF] = COMP.EF_FM_HC_2_FR_RF(EF, FM, HC);
        // uniform: orders within color don't change
        const uniform = COMP.FM_FR_H2_HC_2_uniform(FM, FR, H2, HC);
        // separable: orders between pairs of colors are consistent
        const separable = COMP.FR_RF_H1_H2_HC_2_separable(FR, RF, H1, H2, HC);
        // joint: orders between a color and non-moving region are
        // consistent in end state
        const joint = (() => {
            for (let f = 0; f < FM.length; ++f) {
                if (FM[f]) { continue; }
                for (const F of RF) {
                    for (let j = 1; j < F.length; ++j) {
                        const oj = H1.get(M.encode([f, F[j]]));
                        if (oj == undefined) { continue; } // no overlap
                        for (let i = 0; i < j; ++i) {
                            const oi = H1.get(M.encode([f, F[i]]));
                            if (oi == undefined) { continue; } // no overlap
                            if (oi == oj) { continue; } // match
                            return false;
                        }
                    }
                }
            }
            return true;
        })();
        // local: every color borders an unmoving face
        const local = COMP.FM_FE_EF_RF_2_local(FM, FE, EF, RF);
        // connected: colors form a connected component
        const connected = COMP.FM_FR_RF_EF_2_connected(FM, FR, RF, EF);
        const valid = (uniform && separable && joint && local && connected);
        // // OLD: reverse fold if adjacency graph is chain and
        // // no adjacency is closed.
        // equivilently, we can check that no vertex is "popped",
        // i.e., has folded degree >2 and switches majority assignment
        // if (!valid) {
        //     console.log(`fold invalid: ${
        //         uniform   ? "" : "!uniform, "   }${
        //         separable ? "" : "!separable, " }${
        //         joint     ? "" : "!joint, "     }${
        //         local     ? "" : "!local, "     }${
        //         connected ? "" : "!connected"
        //     }`);
        // }
        const VE = V.map(() => []);
        for (let i = 0; i < EV.length; ++i) {
            const [u, v] = EV[i];
            VE[u].push(i);
            VE[v].push(i);
        }
        let is_sink = false;
        let is_open = true;
        let is_closed = true;
        const V_sink = VE.map((E, i) => {
            let n_folded = 0;
            let M_less_V_new = 0;
            let M_less_V_old = 0;
            const R = new Set();
            for (const e of E) {
                const a1 = EA1[e];
                const a2 = EA2[e];
                n_folded += (a2 == "M") || (a2 == "V");
                M_less_V_new += (a1 == "M") ? 1 : ((a1 == "V") ? -1 : 0);
                M_less_V_old += (a2 == "M") ? 1 : ((a2 == "V") ? -1 : 0);
                const [f, g] = EF[e];
                if ((!FM[f]) || (!FM[g])) { return 0; } // not interior
                R.add(FR[f]);
                R.add(FR[g]);
            }
            if ((n_folded < 4) ||
                ((M_less_V_new < 0) == (M_less_V_old < 0))
            ) { return 0; } // not sink
            is_sink = true;
            let _open = true;
            let _closed = true;
            for (const e of E) {
                const a = EA1[e];
                if ((a != "M") && (a != "V")) { continue; }
                const k = M.encode(EF[e]);
                if (!HC.has(k)) { _open = false; }
            }
            if (R.size > 2) { _closed = false; }
            is_open   &&= _open;
            is_closed &&= _closed;
            if (_open) {
                return 1;
            } else if (_closed) {
                return 2;
            } else {
                return 3;
            }
        });
        let is_inside = true;
        let is_outside = true;
        const V_border = VE.map((E, i) => {
            const outer = [];
            const inner = [];
            let MV_diff1 = 0;
            let MV_diff2 = 0;
            for (const e of E) {
                const a_new = EA1[e];
                if (a_new == "M") { ++MV_diff1; }
                if (a_new == "V") { --MV_diff1; }
                const a_old = EA2[e];
                if (a_old == "M") { ++MV_diff2; }
                if (a_old == "V") { --MV_diff2; }
                if (a_new == "B") { return 0; } // on paper boundary
                if (a_new == a_old) { continue; }
                const [f, g] = EF[e];
                if (FM[f] != FM[g]) { outer.push(a_new); }
                else                { inner.push(a_new); }
            }
            if ((inner.length == 0) ||
                (outer.length != 2) || (outer[0] != outer[1])
            ) { return 0; } // not reverse fold
            if (MV_diff1 == MV_diff2) {
                is_outside = false;
                return 1;
            } else {
                is_inside = false;
                return 2;
            }
        });
        let type;
        if (!valid) {
            type = TYPE.COMPLEX;
        } else if (RF.length == 1) {
            type = TYPE.PURELAND;
        } else if (!is_sink) {
            if (is_inside) {
                type = TYPE.INSIDE_REVERSE;
            } else if (is_outside) {
                type = TYPE.OUTSIDE_REVERSE;
            } else {
                type = TYPE.MIXED_REVERSE;
            }
        } else if (is_open) {
            type = TYPE.OPEN_SINK;
        } else if (is_closed) {
            type = TYPE.CLOSED_SINK;
        } else {
            type = TYPE.MIXED_SINK;
        }
        return [type, RF, FR, EC, V_sink, V_border];
    },
    separators_from_edge: (i, j, FOLD, FM, EF_map, E_map, V_boundary) => {
        const {Vf, VV, FV, EV, EF, EA, H} = FOLD;
        const E_sep = EF.map(() => new Set());
        const f_left = EF_map.get(M.encode([i, j]));
        const f_right = EF_map.get(M.encode([j, i]));
        if ((f_left == undefined) || (!FM[f_left])) { return E_sep; }
        if ((f_right == undefined) || (!FM[f_right])) { return E_sep; }
        // found entering edge
        let o = H.get(M.encode([f_left, f_right]));
        if (o == undefined) { return E_sep; }
        const P = [i];
        const processing = Vf.map(() => false);
        processing[i] = true;
        const dfs = (v) => {
            if (processing[v]) { return; }
            P.push(v);
            processing[v] = true;
            if (V_boundary[v] == 3) { // boundary vertex
                const path = [];
                for (let i = 1; i < P.length; ++i) {
                    const p = P[i - 1];
                    const q = P[i];
                    path.push(E_map.get(M.encode_order_pair([p, q])));
                }
                let good = true;
                // const block = new Set(path);
                // const p = P[0];
                // const q = P[1];
                // const f_left = EF_map.get(M.encode([p, q]));
                // const f_right = EF_map.get(M.encode([q, p]));
                // const F_side = new Map();
                // F_side.set(f_left, 0);
                // F_side.set(f_right, 1);
                // const left = [f_left];
                // for (let qi = 0; qi < left.length; ++qi) {
                //     const f = left[qi];
                //     const V_ = FV[f];
                //     let a = V_[V_.length - 1];
                //     for (const b of V_) {
                //         const e = E_map.get(M.encode_order_pair([a, b]));
                //         const g = EF_map.get(M.encode([a, b]));
                //         if (!block.has(e) &&
                //             (g != undefined) &&
                //             FM[g] && (F_side.get(g) != 0)
                //         ) {
                //             left.push(g);
                //             F_side.set(g, 0);
                //         }
                //         a = b;
                //     }
                // }
                // const right = [f_right];
                // for (let qi = 0; qi < right.length; ++qi) {
                //     const f = right[qi];
                //     const V = FV[f];
                //     let a = V[V.length - 1];
                //     for (const b of V) {
                //         const e = E_map.get(M.encode_order_pair([a, b]));
                //         const g = EF_map.get(M.encode([a, b]));
                //         if (!block.has(e) &&
                //             (g != undefined) &&
                //             FM[g] && (F_side.get(g) != 1)
                //         ) {
                //             right.push(g);
                //             F_side.set(g, 1);
                //         }
                //         a = b;
                //     }
                // }
                // const left_set = new Set(left);
                // if (!left_set.has(right[0])) {
                //     for (const a of left) {
                //         for (const b of right) {
                //             const o_ = H.get(M.encode([a, b]));
                //             if ((o_ != undefined) && (o_ != o)) {
                //                 good = false;
                //                 break;
                //             }
                //         }
                //         if (!good) { break; }
                //     }
                // }
                if (good) {
                    for (const e1 of path) {
                        for (const e2 of path) {
                            E_sep[e1].add(e2);
                            E_sep[e2].add(e1);
                        }
                    }
                }
            } else {
                for (const u of VV[v]) {
                    const f_left = EF_map.get(M.encode([v, u]));
                    const f_right = EF_map.get(M.encode([u, v]));
                    if ((f_left == undefined) || (!FM[f_left])) { continue; }
                    const o2 = H.get(M.encode([f_left, f_right]));
                    if (o != o2) { continue; }
                    dfs(u);
                }
            }
            P.pop();
            processing[v] = false;
        };
        dfs(j);
        return E_sep;
    },
    find_separators: (FOLD, FM) => {
        const {Vf, VV, FV, EV, EF, EA, H} = FOLD;
        const EF_map = new Map();
        for (const [i, F] of FV.entries()) {
            for (const [j, v1] of F.entries()) {
                const v2 = F[(j + 1) % F.length];
                EF_map.set(M.encode([v2, v1]), i);
            }
        }
        const E_map = new Map();
        for (let i = 0; i < EV.length; ++i) {
            E_map.set(M.encode_order_pair(EV[i]), i);
        }
        const V_boundary = Vf.map(() => 0);
        for (let i = 0; i < FM.length; ++i) {
            const c = FM[i] ? 1 : 2;
            for (const v of FV[i]) { V_boundary[v] |= c; }
        }
        for (let i = 0; i < Vf.length; ++i) {
            for (const j of VV[i]) {
                if (EA[E_map.get(M.encode_order_pair([i, j]))] == "B") {
                    V_boundary[i] |= 2;
                }
            }
        }
        const E_sep = EF.map(() => new Set());
        for (let i = 0; i < Vf.length; ++i) {
            if (V_boundary[i] != 3) { continue; }
            for (const j of VV[i]) {
                const EE = COMP.separators_from_edge(
                    i, j, FOLD, FM, EF_map, E_map, V_boundary);
                for (let e1 = 0; e1 < EE.length; ++e1) {
                    for (const e2 of EE[e1]) {
                        E_sep[e1].add(e2);
                        E_sep[e2].add(e1);
                    }
                }
            }
        }
        return [V_boundary, E_sep];
    },
    separator_branch: (FOLD, clicked_edges, e) => {
        const {VV, EV} = FOLD;
        const EV_map = new Map();
        for (let i = 0; i < EV.length; ++i) {
            EV_map.set(M.encode_order_pair(EV[i]), i);
        }
        const seen = new Set();
        seen.add(e);
        const dfs = (e) => {
            for (const u of EV[e]) {
                const E = [];
                for (const v of VV[u]) {
                    const k = EV_map.get(M.encode_order_pair([u, v]));
                    if (clicked_edges.has(k)) { E.push(k); }
                }
                if (E.length != 2) { continue; }
                for (const k of E) {
                    if (seen.has(k)) { continue; }
                    seen.add(k);
                    dfs(k);
                }
            }
        }
        dfs(e);
        return seen;
    },
    draw_separators: (FOLD) => {
        const {Vf, VV, FV, EV, EF, FM, EA, H} = FOLD;
        const EF_map = new Map();
        for (const [i, F] of FV.entries()) {
            for (const [j, v1] of F.entries()) {
                const v2 = F[(j + 1) % F.length];
                EF_map.set(M.encode([v2, v1]), i);
            }
        }
        const E_map = new Map();
        for (let i = 0; i < EV.length; ++i) {
            E_map.set(M.encode_order_pair(EV[i]), i);
        }
        const V_boundary = Vf.map(() => 0);
        for (let i = 0; i < FM.length; ++i) {
            const c = FM[i] ? 1 : 2;
            for (const v of FV[i]) { V_boundary[v] |= c; }
        }
        for (let i = 0; i < Vf.length; ++i) {
            for (const j of VV[i]) {
                if (EA[E_map.get(M.encode_order_pair([i, j]))] == "B") {
                    V_boundary[i] |= 2;
                }
            }
        }
        SVG.draw_points(document.getElementById("notes"), Vf, {
            filter: (i) => (V_boundary[i] == 3), r: 5, fill: "green", id: true,
        });
        const E_sep = EF.map(() => false);
        for (let i = 0; i < Vf.length; ++i) {
            if (V_boundary[i] != 3) { continue; }
            for (const j of VV[i]) {
                const f_left = EF_map.get(M.encode([i, j]));
                const f_right = EF_map.get(M.encode([j, i]));
                if ((f_left == undefined) || (!FM[f_left])) { continue; }
                if ((f_right == undefined) || (!FM[f_right])) { continue; }
                // found entering edge
                const o = H.get(M.encode([f_left, f_right]));
                if (o == undefined) { continue; }
                const P = [i];
                const processing = Vf.map(() => false);
                processing[i] = true;
                const dfs = (v) => {
                    if (processing[v]) { return; }
                    P.push(v);
                    processing[v] = true;
                    if (V_boundary[v] == 3) {
                        const path = [];
                        for (let i = 1; i < P.length; ++i) {
                            const p = P[i - 1];
                            const q = P[i];
                            path.push(E_map.get(M.encode_order_pair([p, q])));
                        }
                        const block = new Set(path);
                        const p = P[0];
                        const q = P[1];
                        const f_left = EF_map.get(M.encode([p, q]));
                        const f_right = EF_map.get(M.encode([q, p]));
                        const F_side = new Map();
                        F_side.set(f_left, 0);
                        F_side.set(f_right, 1);
                        const left = [f_left];
                        for (let qi = 0; qi < left.length; ++qi) {
                            const f = left[qi];
                            const V_ = FV[f];
                            let a = V_[V_.length - 1];
                            for (const b of V_) {
                                const e = E_map.get(M.encode_order_pair([a, b]));
                                const g = EF_map.get(M.encode([a, b]));
                                if (!block.has(e) &&
                                    (g != undefined) &&
                                    FM[g] &&
                                    !F_side.has(g)
                                ) {
                                    left.push(g);
                                    F_side.set(g, 0);
                                }
                                a = b;
                            }
                        }
                        const right = [f_right];
                        for (let qi = 0; qi < right.length; ++qi) {
                            const f = right[qi];
                            const V = FV[f];
                            let a = V[V.length - 1];
                            for (const b of V) {
                                const e = E_map.get(M.encode_order_pair([a, b]));
                                const g = EF_map.get(M.encode([a, b]));
                                if (!block.has(e) &&
                                    (g != undefined) &&
                                    FM[g] &&
                                    !F_side.has(g)
                                ) {
                                    right.push(g);
                                    F_side.set(g, 1);
                                }
                                a = b;
                            }
                        }
                        let good = true;
                        for (const a of left) {
                            for (const b of right) {
                                const o_ = H.get(M.encode([a, b]));
                                if ((o_ != undefined) && (o_ != o)) {
                                    good = false;
                                    break;
                                }
                            }
                            if (!good) { break; }
                        }
                        if (good) {
                            for (const e of path) { E_sep[e] = true; }
                        }
                    } else {
                        for (const u of VV[v]) {
                            const f_left = EF_map.get(M.encode([v, u]));
                            const f_right = EF_map.get(M.encode([u, v]));
                            if ((f_left == undefined) || (!FM[f_left])) { continue; }
                            const o2 = H.get(M.encode([f_left, f_right]));
                            if (o != o2) { continue; }
                            dfs(u);
                        }
                    }
                    P.pop();
                    processing[v] = false;
                };
                dfs(j);
                for (let i = 0; i < E_sep.length; ++i) {
                    if (!E_sep[i]) { continue; }
                    const el = document.getElementById(`flat_e${i}`);
                    el.setAttribute("stroke-width", 10);
                }
            }
        }
    },
    slide_LL: (FOLD, CELL, val, flip) => {
        const {Ff, LL} = FOLD;
        const {CD, Ctop, Ccolor} = CELL;
        const F_set = new Set();
        for (let i = 0; i < Ff.length; ++i) { F_set.add(i); }
        if (flip) { LL.reverse(); }
        for (let i = LL.length - 1; i >= val; --i) {
            const layer = LL[i];
            for (const fi of layer) {
                F_set.delete(fi);
            }
        }
        if (flip) { LL.reverse(); }
        for (let i = 0; i < CD.length; ++i) {
            const S = CD[i];
            const D = S.map(a => a);
            if (flip) { D.reverse(); }
            Ctop[i] = undefined;
            while (D.length > 0) {
                const fi = D.pop();
                if (!F_set.has(fi)) { Ctop[i] = fi; break; }
            }
        }
    },
    split_FOLD_on_line: (FOLD, line) => {
        // assumes convex faces (or line divides a face into at most two pieces
        const {FV, V, Vf} = FOLD;
        const eps = FOLD.eps/100;
        const [u, d] = line;
        const [a, b] = COMP.line_2_coords(line);
        const nV = V.length;
        const V2 = V.map(v => v);
        const Vf2 = Vf.map(v => v);
        const VD = V.map(v => { // signed distance from fold line
            const dv = M.dot(u, v) - d;
            return (Math.abs(dv) <= eps) ? 0 : dv;
        });
        const EV_map = new Map();
        const FV2 = FV.map((F) => {
            const pair = [[], []];
            const nF = F.length;
            let [neg, pos] = [false, false];
            for (const v of F) {
                const d = VD[v];
                if (d < 0) { neg = true; }
                if (d > 0) { pos = true; }
            }
            console.assert(neg || pos); // face has zero area?
            if (neg != pos) {
                pair[pos ? 0 : 1] = F.map(i => i);
                return pair;
            }
            let i = 1;
            while ((i < nF) && ((VD[F[i - 1]] < 0) == (VD[F[i]] < 0))) {
                ++i;
            }
            for (let j = 0; j < nF; ++j) {
                const i1 = (i + j) % nF;
                const v1 = F[i1];
                if (Math.abs(VD[v1]) == 0) {
                    pair[0].push(v1);
                    pair[1].push(v1);
                    continue;
                }
                if (VD[v1] > 0) { pair[0].push(v1); }
                if (VD[v1] < 0) { pair[1].push(v1); }
                const i2 = (i1 + 1) % nF;
                const v2 = F[i2];
                if (Math.abs(VD[v2]) == 0) { continue; }
                if ((VD[v1] < 0) != (VD[v2] < 0)) {
                    const s = M.encode_order_pair([v1, v2]);
                    let xi = EV_map.get(s);
                    if (xi == undefined) {
                        const x = M.intersect([V[v1], V[v2]], [a, b], eps);
                        const xf = M.add(Vf[v1],
                            M.mul(
                                M.sub(Vf[v2], Vf[v1]),
                                M.dist(x, V[v1])/M.dist(V[v1], V[v2])
                            )
                        );
                        xi = V2.length;
                        EV_map.set(s, xi);
                        V2.push(x);
                        Vf2.push(xf);
                        VD.push(0);
                    }
                    pair[0].push(xi);
                    pair[1].push(xi);
                }
            }
            return pair;
        });
        const FV_ = [];
        const F_map = FV2.map(() => []);
        for (let fi = 0; fi < FV2.length; ++fi) {
            for (const f of FV2[fi]) {
                if (f.length > 0) {
                    F_map[fi].push(FV_.length);
                    FV_.push(f);
                }
            }
        }
        return [FV_, V2, Vf2, VD, F_map];
    },
    get_groups: (FV, VD) => {
        const EF_map = new Map();
        for (const [i, F] of FV.entries()) {
            for (const [j, v1] of F.entries()) {
                const v2 = F[(j + 1) % F.length];
                EF_map.set(M.encode([v2, v1]), i);
            }
        }
        const FG = FV.map(() => undefined);
        let g = 0;
        const dfs = (i) => {
            if (FG[i] != undefined) { return; }
            FG[i] = g;
            const F = FV[i];
            for (const [j, v1] of F.entries()) {
                const v2 = F[(j + 1) % F.length];
                if ((VD[v1] == 0) && (VD[v2] == 0)) { continue; }
                const fi = EF_map.get(M.encode([v1, v2]));
                if (fi != undefined) { dfs(fi); }
            }
        };
        for (let i = 0; i < FG.length; ++i) {
            if (FG[i] != undefined) { continue; }
            dfs(i);
            ++g;
        }
        return FG;
    },
    map_order: (BI, F_map, FO) => {
        const FO_ = [];
        for (const [f, g, o] of FO) {
            for (const f_ of F_map[f]) {
                for (const g_ of F_map[g]) {
                    const pair = M.encode_order_pair([f_, g_]);
                    if (BI.has(pair)) {
                        FO_.push([f_, g_, o]);
                    }
                }
            }
        }
        return FO_;
    },
    get_lines: (P, eps) => {
        if ((P.length < 2) && (P.length > 4)) { return []; }
        const out = [];
        if (P.length == 2) {
            const [a, b] = P;
            const m = M.div(M.add(a, b), 2);
            const v1 = M.sub(b, a); // line through AB
            const v2 = M.perp(v1);  // perpendicular bisector of AB
            for (let v of [v1, v2]) {
                const u = M.unit(v);
                const d = M.dot(m, u);
                out.push([u, d]);
            }
        } else if (P.length == 3) {
            const [a, b, c] = P;
            if (Math.abs(M.area2(a, b, c)) > eps*eps) {
                {   // angle bisector of ABC
                    const ba = M.unit(M.sub(a, b));
                    const bc = M.unit(M.sub(c, b));
                    const u = M.unit(M.add(ba, bc));
                    const v = M.perp(u);
                    const d = M.dot(b, v);
                    out.push([v, d]);
                }
                {   // through A perpendicular to BC
                    const u = M.unit(M.sub(c, b));
                    const d = M.dot(a, u);
                    out.push([u, d]);
                }
            } else {
                {   // through A perpendicular to BC
                    const u = M.unit(M.sub(c, b));
                    const d = M.dot(a, u);
                    out.push([u, d]);
                }
            }
        } else if (P.length == 4) {
            const [a, b, c, d] = P; // angle bisector of AB and CD
            const ba = M.unit(M.sub(a, b));
            const cd = M.unit(M.sub(d, c));
            const perp = M.perp(cd);
            const parallel = Math.abs(M.dot(perp, ba)) < eps;
            const u = parallel ? ba : M.unit(M.add(ba, cd));
            const U = parallel ? [u] : [u, M.perp(u)];
            for (const v of U) {
                const pc = M.dot(c, v);
                const pa = M.dot(a, v);
                const pb = M.dot(b, v);
                const bb = M.add(a, M.mul(M.sub(b, a), (pc - pa)/(pb - pa)));
                const x = M.div(M.add(bb, c), 2);
                const pv = M.perp(v);
                out.push([pv, M.dot(pv, x)]);
            }
        }
        return out.map(([u, d]) => {
            return (d < 0) ? [M.mul(u, -1), -d] : [u, d];
        });
    },
    line_2_coords: (line) => {
        const [u, d] = line;
        const p = M.mul(u, d);
        const off = M.mul(M.perp(u), 10);
        const p1 = M.add(p, off);
        const p2 = M.sub(p, off);
        return [p1, p2];
    },
};

/*  Axioms:
 *     #pts | abbr  | description
 *  1)  2   | TPTP  | line through two points
 *  2)  2   | P2P   | line folding point to point
 *  3) 3,4  | S2S   | line folding segment to segment
 *  4)  3   | RSTP  | line perpendicular to segment through point
 *  5)  4   | P2STP | line folding point to segment through point
 *  6)  5   | P2SRS | line folding point to segment perpendicular to segment
 *  7)  6   | PP2SS | line folding points to segments
 */
