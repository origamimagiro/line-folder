import { M } from "./math.js";
import { NOTE } from "./note.js";
import { SVG } from "./svg.js";
import { IO } from "./io.js";
import { X } from "./conversion.js";

window.onload = () => { MAIN.startup(); };  // entry point

const MAIN = {
    color: {
        background: "lightgray",
        normal: "black",
        active: "red",
        select: "blue",
        face: {
            top: "gray",
            bottom: "white",
        },
        edge: {
            U: "black",
            F: "lightgray",
            B: "black",
        },
    },
    opacity: {
        normal: 0.01,
        hover: 1,
    },
    radius: {
        normal: 5,
        hover: 10,
    },
    startup: () => {
        NOTE.clear_log();
        NOTE.start("*** Starting Flat-Folder ***");
        NOTE.time("Initializing interface");
        const [b, s] = [50, SVG.SCALE];
        const main = document.getElementById("main");
        for (const [k, v] of Object.entries({
            xmlns: SVG.NS,
            style: `background: ${MAIN.color.background}`,
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
            if (Ff[d] != flip)  { return MAIN.color.face.top; }
            else                { return MAIN.color.face.bottom; }
        });
        const visible = MAIN.PP_Ctop_CP_SC_2_visible(Q, PP, tops, CP, SC);
        SVG.draw_polygons(svg, cells, {
            id: "fold_c", fill: colors, stroke: colors});
        const lines = SP.map((ps) => M.expand(ps, Q));
        SVG.draw_segments(svg, lines, {
            id: "fold_s_crease", stroke: MAIN.color.edge.F,
            filter: (i) => SD[i] == "C"});
        SVG.draw_segments(svg, lines, {
            id: "fold_s_edge", stroke: MAIN.color.edge.B,
            filter: (i) => SD[i] == "B"});
        const Lsvg = SVG.append("g", svg, {id: "lines"});
        SVG.draw_points(svg, Q, {
            id: "fold_p", filter: (i) => visible[i],
            fill: MAIN.color.normal, r: MAIN.radius.normal,
            opacity: MAIN.opacity.normal,
        });
        const clicked = new Map();
        for (let i = 0; i < Q.length; ++i) {
            const el = document.getElementById(`fold_p${i}`);
            if (el != undefined) {
                el.onmouseover = () => MAIN.point_over(el);
                el.onmouseout = () => MAIN.point_out(i, el, clicked);
                el.onclick = () => MAIN.point_click(i, el, clicked, Lsvg, Q);
            }
        }
        NOTE.lap();
    },
    point_over: (el) => {
        el.setAttribute("fill", MAIN.color.select);
        el.setAttribute("r", MAIN.radius.hover);
        el.setAttribute("opacity", MAIN.opacity.hover);
    },
    point_out: (i, el, clicked) => {
        const unclicked = (clicked.get(i) == undefined);
        el.setAttribute("r", MAIN.radius.normal);
        el.setAttribute("opacity", unclicked
            ? MAIN.opacity.normal
            : MAIN.opacity.hover
        );
        el.setAttribute("fill", unclicked
            ? MAIN.color.normal
            : MAIN.color.active
        );
    },
    point_click: (i, el, clicked, svg, Q) => {
        NOTE.time(`Clicked point ${i}`);
        if (clicked.has(i) || (clicked.size > 3)) {
            SVG.clear("lines");
            MAIN.clear_clicked(clicked);
            MAIN.point_over(el);
            return;
        }
        clicked.set(i, el);
        NOTE.time(`Clicked set is: [${Array.from(clicked.keys())}]`);
        el.setAttribute("fill", MAIN.color.select);
        SVG.clear("lines");
        const L = MAIN.get_lines(Array.from(clicked.keys()).map(i => Q[i]));
        if (L.length > 0) {
            const g = SVG.draw_segments(svg, L.map(l => MAIN.line_2_coords(l)), {
                id: "line", stroke: MAIN.color.normal,
                stroke_width: MAIN.radius.normal,
            });
            for (const el of g.children) {
                el.onclick = () => MAIN.line_click(el, clicked);
                el.onmouseover = () => MAIN.line_over(el);
                el.onmouseout = () => MAIN.line_out(el);
            }
        } else if (clicked.size > 1) {
            SVG.clear("lines");
            MAIN.clear_clicked(clicked);
            MAIN.point_over(el);
        }
    },
    clear_clicked: (clicked) => {
        for (const [i, el] of clicked) {
            el.setAttribute("fill", MAIN.color.normal);
            el.setAttribute("opacity", MAIN.opacity.normal);
        }
        clicked.clear();
    },
    line_over: (el) => {
        el.setAttribute("stroke", MAIN.color.select);
        el.setAttribute("stroke-width", MAIN.radius.hover);
    },
    line_out: (el) => {
        el.setAttribute("stroke", MAIN.color.normal);
        el.setAttribute("stroke-width", MAIN.radius.normal);
    },
    line_click: (el, clicked) => {
        SVG.clear("lines");
        MAIN.clear_clicked(clicked);
        document.getElementById("lines").appendChild(el);
        el.onmouseover = undefined;
        el.onmouseout = undefined;
        el.onclick = undefined;
        el.setAttribute("stroke", MAIN.color.active);
        el.setAttribute("stroke-width", MAIN.radius.normal);
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
            const v1 = M.sub(b, a); // line through AB
            const v2 = M.perp(v1);  // perpendicular bisector of AB
            for (let v of [v1, v2]) {
                const u = M.unit(v);
                const d = M.dot(m, u);
                out.push([u, d]);
            }
        } else if (P.length == 3) {
            const [a, b, c] = P;
            if (Math.abs(M.area2(a, b, c)) > 0.000001) {
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
            }
        } else if (P.length == 4) {
            const [a, b, c, d] = P; // angle bisector of AB and CD
            const ba = M.unit(M.sub(a, b));
            const cd = M.unit(M.sub(d, c));
            const perp = M.perp(cd);
            const parallel = Math.abs(M.dot(perp, ba)) < 0.000001;
            const u = parallel ? ba : M.unit(M.add(ba, cd));
            const U = parallel ? [u] : [u, M.perp(u)];
            for (const v of U) {
                const pc = M.dot(c, v);
                const pa = M.dot(a, v);
                const pb = M.dot(b, v);
                const bb = M.add(a, M.mul(M.sub(b, a), (pc - pa)/(pb - pa)));
                const x = M.div(M.add(bb, c), 2);
                const pv = M.perp(v);
                const [aL, bL, cL, dL] = [a, b, c, d].map(p => {
                    return M.dot(M.sub(p, x), pv) < 0;
                });
                if ((aL != bL) || (bL != cL) || (cL != dL)) {
                    out.push([pv, M.dot(pv, x)]);
                }
            }
        }
        return out;
    },
    line_2_coords: (L) => {
        const [u, d] = L;
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
