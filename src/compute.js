import { M } from "./flatfolder/math.js";

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
};
