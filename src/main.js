import { M } from "./math.js";
import { NOTE } from "./note.js";
import { SVG } from "./svg.js";
import { IO } from "./io.js";
import { X } from "./conversion.js";
import { GUI } from "./gui.js";

window.onload = () => { MAIN.startup(); };  // entry point

const MAIN = {
    startup: () => {
        NOTE.clear_log();
        NOTE.start("*** Starting Flat-Folder ***");
        NOTE.time("Initializing interface");
        const [b, s] = [50, SVG.SCALE];
        const main = document.getElementById("main");
        for (const [k, v] of Object.entries({
            xmlns: SVG.NS,
            style: `background: ${GUI.COLORS.background}`,
            viewBox: [0, 0, 2*s, s].join(" "),
        })) {
            main.setAttribute(k, v);
        }
        for (const [i, id] of ["org"].entries()) {
            const svg = document.getElementById(id);
            for (const [k, v] of Object.entries({
                xmlns: SVG.NS,
                height: s,
                width: s,
                x: i*s,
                y: 0,
                viewBox: [-b, -b, s + 2*b, s + 2*b].join(" "),
            })) {
                svg.setAttribute(k, v);
            }
        }
        document.getElementById("import").onchange = (e) => {
            if (e.target.files.length > 0) {
                const file_reader = new FileReader();
                file_reader.onload = MAIN.process_file;
                file_reader.readAsText(e.target.files[0]);
            }
        };
        NOTE.end();
    },
    process_file: (e) => {
        NOTE.clear_log();
        NOTE.start("*** Starting File Import ***");
        const doc = e.target.result;
        const file_name = document.getElementById("import").value;
        const parts = file_name.split(".");
        const type = parts[parts.length - 1].toLowerCase();
        if (type != "fold") {
            console.log(`Found file with extension ${type}, FOLD format required`);
            return;
        }
        NOTE.time(`Importing from file ${file_name}`);
        const [FOLD, CELL] = MAIN.import_state(doc);
        MAIN.update_fold(FOLD, CELL);
        document.getElementById("flip").onchange = (e) => {
            NOTE.start("Flipping model");
            MAIN.update_fold(FOLD, CELL);
            NOTE.end();
        };
    },
    import_state: (doc) => {
        const [V, FV, FO] = MAIN.FOLD_2_V_FV_FO(doc);
        for (const property of [V, FV, FO]) {
            if (property == undefined) { return; }
        }
        const Ff = MAIN.FV_V_2_Ff(FV, V);
        const EV_set = new Set();
        for (const fV of FV) {
            let i = fV.length - 1;
            for (let j = 0; j < fV.length; ++j) {
                EV_set.add(M.encode_order_pair([fV[i], fV[j]]));
                i = j;
            }
        }
        const EV = Array.from(EV_set).sort().map(k => M.decode(k));
        const [EF, FE] = X.EV_FV_2_EF_FE(EV, FV);
        const L = EV.map(vs => vs.map(i => V[i]));
        const eps = M.min_line_length(L) / M.EPS;
        NOTE.time(`Using eps ${eps} from min line length ${
            eps*M.EPS} (factor ${M.EPS})`);
        NOTE.time("Constructing points and segments from edges");
        const [P, SP, SE] = X.L_2_V_EV_EL(L, eps);
        const P_norm = M.normalize_points(P);
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
        const [CF, FC] = X.EF_FV_SP_SE_CP_SC_2_CF_FC(EF, FV, SP, SE, CP, SC);
        const edges = FO.map(([f1, f2, o]) => {
            return M.encode((Ff[f2]*o >= 0) ? [f1, f2] : [f2, f1]);
        });
        const CD = X.CF_edges_flip_2_CD(CF, edges);
        NOTE.count(CF, "face-cell adjacencies");
        NOTE.lap();
        const FOLD = {V, FV, EV, EF, FE, Ff};
        const CELL = {P_norm, SP, SE, PP, CP, CS, SC, CF, FC, CD};
        return [FOLD, CELL];
    },
    update_fold: (FOLD, CELL) => {
        const {EF, Ff} = FOLD;
        const {P_norm, SP, SE, PP, CP, SC, CF, CD} = CELL;
        const svg = SVG.clear("org");
        const flip = document.getElementById("flip").checked;
        const tops = CD.map(S => flip ? S[0] : S[S.length - 1]);
        const SD = X.EF_SE_SC_CF_CD_2_SD(EF, SE, SC, CF, tops);
        const m = [0.5, 0.5];
        const Q = P_norm.map(p => (flip ? M.add(M.refX(M.sub(p, m)), m) : p));
        const cells = CP.map(V => M.expand(V, Q));
        const colors = tops.map(d => {
            if (d == undefined) { return undefined; }
            if (Ff[d] != flip)  { return GUI.COLORS.face.top; }
            else                { return GUI.COLORS.face.bottom; }
        });
        const visible = MAIN.PP_Ctop_CP_SC_2_visible(Q, PP, tops, CP, SC);
        SVG.draw_polygons(svg, cells, {
            id: "fold_c", fill: colors, stroke: colors});
        const lines = SP.map((ps) => M.expand(ps, Q));
        SVG.draw_segments(svg, lines, {
            id: "fold_s_crease", stroke: GUI.COLORS.edge.F,
            filter: (i) => SD[i] == "C"});
        SVG.draw_segments(svg, lines, {
            id: "fold_s_edge", stroke: GUI.COLORS.edge.B,
            filter: (i) => SD[i] == "B"});
        const Lsvg = SVG.append("g", svg, {id: "lines"});
        SVG.draw_points(svg, Q, {
            id: "fold_p", fill: GUI.COLORS.line.N, r: 5,
            filter: (i) => visible[i],
        });
        const clicked = [];
        for (let i = 0; i < Q.length; ++i) {
            const el = document.getElementById(`fold_p${i}`);
            if (el != undefined) {
                el.onclick = () => {
                    NOTE.time(`Clicked point ${i}`);
                    const seen = new Set(clicked.map(([i, el]) => i));
                    if (seen.has(i) || (clicked.length >= 5)) {
                        NOTE.time(`Emptying clicked set`);
                        SVG.clear("lines");
                        while (clicked.length > 0) {
                            const [i, el] = clicked.pop();
                            el.setAttribute("fill", GUI.COLORS.line.N);
                        }
                        return;
                    }
                    clicked.push([i, el]);
                    NOTE.time(`Clicked set is: [${clicked.map(([i, el]) => i)}]`);
                    el.setAttribute("fill", GUI.COLORS.line.A);
                    SVG.clear("lines");
                    const L = MAIN.get_lines(clicked.map(([i, el]) => Q[i]));
                    if (L.length > 0) {
                        SVG.draw_segments(Lsvg, L.map(l => MAIN.line_2_coords(l)), {
                            id: "line", stroke: GUI.COLORS.line.L, stroke_width: 5,
                        });
                    }
                };
            }
        }
        NOTE.lap();
    },
    FOLD_2_V_FV_FO: (doc) => {
        const ex = JSON.parse(doc);
        const properties = ["vertices_coords", "faces_vertices", "faceOrders"];
        return properties.map(property => {
            const val = ex[property];
            if (val == undefined) {
                NOTE.time(`FOLD file must contain ${property}, but not found`);
                return undefined;
            }
            return val;
        });
    },
    FV_V_2_Ff: (FV, V) => FV.map(fV => (M.polygon_area2(fV.map(i => V[i])) < 0)),
    PP_Ctop_CP_SC_2_visible: (P, PP, Ctop, CP, SC) => {
        const SC_map = new Map();
        for (const [i, C] of CP.entries()) {
            for (const [j, p1] of C.entries()) {
                const p2 = C[(j + 1) % C.length];
                SC_map.set(M.encode([p2, p1]), i);
            }
        }
        return PP.map((V, i) => {
            const F = [];
            const A = V.map(j => M.angle(M.sub(P[j], P[i])));
            for (const j of V) {
                const c = SC_map.get(M.encode([i, j]));
                F.push((c == undefined) ? -1 : Ctop[c]);
            }
            const F_set = new Map();
            for (let i = 0; i < F.length; ++i) {
                const f = F[i];
                let ang = F_set.get(f);
                if (ang == undefined) { ang = 0; }
                ang += A[i] - A[(i - 1) % A.length];
                F_set.set(f, ang);
            }
            if (F_set.size > 2) {
                return true;
            }
            for (const [f, ang] of F_set.entries()) {
                if (Math.abs(Math.abs(ang) - Math.PI) > 0.001) {
                    return true;
                }
            }
            return false;
        });
    },
    get_lines: (P) => {
        const out = [];
        if (P.length == 2) {
            const [a, b] = P;
            const m = M.div(M.add(a, b), 2);
            const v1 = M.sub(b, a);
            const v2 = M.perp(v1);
            for (let v of [v1, v2]) {
                const ang = M.angle(v);
                const u = M.unit(v);
                const d = M.dot(m, u);
                out.push([ang, d]);
            }
        } else if (P.length == 3) {
            const [a, b, c] = P;
            if (Math.abs(M.area2(a, b, c)) > M.FLOAT_EPS) {
                const ba = M.unit(M.sub(a, b));
                const bc = M.unit(M.sub(c, b));
                const u = M.unit(M.add(ba, bc));
                const v = M.perp(u);
                const ang = M.angle(v);
                const d = M.dot(b, v);
                out.push([ang, d]);
            }
        } else if (P.length == 4) {
        } else if (P.length == 5) {
        }
        return out;
    },
    line_2_coords: (L) => {
        const [a, d] = L;
        const u = [Math.cos(a), Math.sin(a)];
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
 */
