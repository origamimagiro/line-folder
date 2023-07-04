import { M } from "./flatfolder/math.js";
import { NOTE } from "./flatfolder/note.js";
import { SVG } from "./flatfolder/svg.js";
import { IO } from "./flatfolder/io.js";
import { X } from "./flatfolder/conversion.js";

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
        const [V_org, FV, FO] = MAIN.FOLD_2_V_FV_FO(doc);
        const V = M.normalize_points(V_org);
        for (const property of [V, FV, FO]) {
            if (property == undefined) { return; }
        }
        const [FOLD, CELL] = MAIN.V_FV_2_FOLD_CELL(V, FV);
        FOLD.FO = FO;
        MAIN.update_fold(FOLD, CELL);
        document.getElementById("flip").onchange = (e) => {
            NOTE.start("Flipping model");
            MAIN.update_fold(FOLD, CELL);
            NOTE.end();
        };
    },
    FOLD_CELL_2_STATE: (FOLD, CELL, flip) => {
        const {EF, Ff, FO} = FOLD;
        const {P, SE, PP, CP, SC, CF} = CELL;
        const m = [0.5, 0.5];
        const Q = P.map(p => (flip ? M.add(M.refX(M.sub(p, m)), m) : p));
        const edges = FO.map(([f1, f2, o]) => {
            return M.encode((Ff[f2]*o >= 0) ? [f1, f2] : [f2, f1]);
        });
        const CD = X.CF_edges_flip_2_CD(CF, edges);
        const Ctop = CD.map(S => flip ? S[0] : S[S.length - 1]);
        const Ccolor = Ctop.map(d => {
            if (d == undefined) { return undefined; }
            if (Ff[d] != flip)  { return MAIN.color.face.top; }
            else                { return MAIN.color.face.bottom; }
        });
        const SD = X.EF_SE_SC_CF_CD_2_SD(EF, SE, SC, CF, Ctop);
        const Pvisible = MAIN.PP_Ctop_CP_SC_2_Pvisible(Q, PP, Ctop, CP, SC);
        return {Q, CD, Ctop, Ccolor, SD, Pvisible};
    },
    update_fold: (FOLD, CELL) => {
        SVG.clear("export");
        const fold_button = document.getElementById("fold_button");
        fold_button.style.display = "none";
        const flip = document.getElementById("flip").checked;
        const STATE = MAIN.FOLD_CELL_2_STATE(FOLD, CELL, flip);
        const {CP, SP} = CELL;
        const {Q, SD, Ccolor, Pvisible} = STATE;
        const svg = SVG.clear("org");
        const cells = CP.map(V => M.expand(V, Q));
        SVG.draw_polygons(svg, cells, {
            id: "fold_c", fill: Ccolor, stroke: Ccolor});
        const lines = SP.map((ps) => M.expand(ps, Q));
        SVG.draw_segments(svg, lines, {
            id: "fold_s_crease", stroke: MAIN.color.edge.F,
            filter: (i) => SD[i] == "C"});
        SVG.draw_segments(svg, lines, {
            id: "fold_s_edge", stroke: MAIN.color.edge.B,
            filter: (i) => SD[i] == "B"});
        const Lsvg = SVG.append("g", svg, {id: "lines"});
        SVG.draw_points(svg, Q, {
            id: "fold_p", filter: (i) => Pvisible[i],
            fill: MAIN.color.normal, r: MAIN.radius.normal,
            opacity: MAIN.opacity.normal,
        });
        const clicked = new Map();
        for (let i = 0; i < Q.length; ++i) {
            const el = document.getElementById(`fold_p${i}`);
            if (el != undefined) {
                el.onmouseover = () => MAIN.point_over(el);
                el.onmouseout = () => MAIN.point_out(i, el, clicked);
                el.onclick = () => MAIN.point_click(i, el, clicked,
                    Lsvg, Q, FOLD, CELL);
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
    point_click: (i, el, clicked, svg, Q, FOLD, CELL) => {
        const {P} = CELL;
        NOTE.time(`Clicked point ${i}`);
        if (clicked.has(i) || (clicked.size >= 4)) {
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
            for (let j = 0; j < g.children.length; ++j) {
                const el = g.children[j];
                el.onclick = () => MAIN.line_click(el, clicked, L[j], FOLD, CELL);
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
    FV_VD_2_FG: (FV, VD) => {
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
                if (fi != undefined) {
                    dfs(fi);
                }
            }
        };
        for (let i = 0; i < FG.length; ++i) {
            if (FG[i] != undefined) { continue; }
            dfs(i);
            ++g;
        }
        return FG;
    },
    line_click: (el, clicked, line, FOLD_, CELL_) => {
        const V_ = FOLD_.V;
        const FV_ = FOLD_.FV;
        const FO_ = FOLD_.FO;
        const eps_ = FOLD_.eps;
        const flip = document.getElementById("flip").checked;
        if (flip) {
            const [u, d] = line;
            const o = [0.5, 0.5];
            const p = M.add(M.refX(M.sub(M.mul(u, d), o)), o);
            const v = M.refX(u);
            const d2 = M.dot(p, v);
            line = [v, d2];
        }
        const [FV2, V, VD] = MAIN.FV_V_line_eps_2_FV2_V2_VD2(FV_, V_, line, eps_);
        const FV = [];
        const F_map = FV2.map(() => []);
        for (let fi = 0; fi < FV2.length; ++fi) {
            for (const f of FV2[fi]) {
                if (f.length > 0) {
                    F_map[fi].push(FV.length);
                    FV.push(f);
                }
            }
        }
        const FG = MAIN.FV_VD_2_FG(FV, VD);
        const [FOLD, CELL] = MAIN.V_FV_2_FOLD_CELL(V, FV);
        const {Ff, EF} = FOLD;
        const {P, PP, CP, CF, SE, SC, SP} = CELL;
        const BF_set = new Set(X.CF_2_BF(CF));
        const FO = [];
        for (const [f, g, o] of FO_) {
            for (const f_ of F_map[f]) {
                for (const g_ of F_map[g]) {
                    const pair = M.encode_order_pair([f_, g_]);
                    if (BF_set.has(pair)) {
                        FO.push([f_, g_, o]);
                    }
                }
            }
        }
        FOLD.FO = FO;
        const STATE = MAIN.FOLD_CELL_2_STATE(FOLD, CELL, flip);
        const {Q, SD, Ctop, Ccolor, Pvisible} = STATE;
        const svg = SVG.clear("org");
        const cells = CP.map(V => M.expand(V, Q));
        SVG.draw_polygons(svg, cells, {
            id: "fold_c", fill: Ccolor, stroke: Ccolor});
        const lines = SP.map((ps) => M.expand(ps, Q));
        SVG.draw_segments(svg, lines, {
            id: "fold_s_crease", stroke: MAIN.color.edge.F,
            filter: (i) => SD[i] == "C"});
        SVG.draw_segments(svg, lines, {
            id: "fold_s_edge", stroke: MAIN.color.edge.B,
            filter: (i) => SD[i] == "B"});
        const Lsvg = SVG.append("g", svg, {id: "lines"});
        SVG.draw_points(svg, Q, {
            id: "fold_p", filter: (i) => Pvisible[i],
            fill: MAIN.color.normal, r: MAIN.radius.normal,
            opacity: MAIN.opacity.normal,
        });
        SVG.clear("lines");
        MAIN.clear_clicked(clicked);
        document.getElementById("lines").appendChild(el);
        el.onmouseout = () => {
            el.setAttribute("stroke", MAIN.color.active);
            el.setAttribute("stroke-width", MAIN.radius.normal);
        };
        el.onclick = () => MAIN.update_fold(FOLD_, CELL_);
        el.setAttribute("stroke", MAIN.color.active);
        el.setAttribute("stroke-width", MAIN.radius.normal);
        const clicked_groups = new Set();
        for (let i = 0; i < CF.length; ++i) {
            const el = document.getElementById(`fold_c${i}`);
            el.onmouseover = () => MAIN.cell_over(i, Ctop, FG);
            el.onmouseout = () => MAIN.cell_out(i, Ctop, FG, Ccolor, clicked_groups);
            el.onclick = () => MAIN.cell_click(i, Ctop, FG, Ccolor, clicked_groups);
        }
        const fold_button = document.getElementById("fold_button");
        fold_button.style.display = "inline";
        fold_button.onclick = () => {
            if (clicked_groups.size == 0) {
                MAIN.update_fold(FOLD_, CELL_);
            } else {
                MAIN.write(MAIN.make_fold(V, FV_, FV, F_map, FG, FO_, clicked_groups));
            }
        };
    },
    make_fold: (V, FV_, FV, F_map_old, FG, FO_, clicked_groups) => {
        const FVx = [];
        const F_map = F_map_old.map(() => []);
        for (let fi = 0; fi < F_map_old.length; ++fi) {
            const F = F_map_old[fi];
            if (F.length == 2) {
                const [f1, f2] = F;
                if (clicked_groups.has(FG[f1]) == clicked_groups.has(FG[f2])) {
                    F_map[fi].push(FVx.length);
                    FVx.push(FV_[fi]);
                } else {
                    F_map[fi].push(FVx.length);
                    FVx.push(FV[f1]);
                    F_map[fi].push(FVx.length);
                    FVx.push(FV[f2]);
                }
            } else if (F.length == 1) {
                F_map[fi].push(FVx.length);
                FVx.push(FV_[fi]);
            } else {
                throw new Error("malformed FV2");
            }
        }
        const [FOLD, CELL] = MAIN.V_FV_2_FOLD_CELL(V, FVx);
        const {Ff, EF} = FOLD;
        const {P, PP, CP, CF, SE, SC, SP} = CELL;
        const BF_set = new Set(X.CF_2_BF(CF));
        const FO = [];
        for (const [f, g, o] of FO_) {
            for (const f_ of F_map[f]) {
                for (const g_ of F_map[g]) {
                    const pair = M.encode_order_pair([f_, g_]);
                    if (BF_set.has(pair)) {
                        FO.push([f_, g_, o]);
                    }
                }
            }
        }
        FOLD.FO = FO;
        return FOLD;
    },
    V_FV_2_FOLD_CELL: (V, FV) => {
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
        const FOLD = {V, FV, EV, EF, FE, Ff, eps};
        const CELL = {P, SP, SE, PP, CP, CS, SC, CF, FC};
        return [FOLD, CELL];
    },
    cell_over: (i, Ctop, FG) => {
        const g = FG[Ctop[i]];
        for (let j = 0; j < Ctop.length; ++j) {
            if (g == FG[Ctop[j]]) {
                const el = document.getElementById(`fold_c${j}`);
                el.setAttribute("fill", "lightpink");
            }
        }
    },
    cell_out: (i, Ctop, FG, Ccolor, clicked) => {
        const g1 = FG[Ctop[i]];
        for (let j = 0; j < Ctop.length; ++j) {
            const g2 = FG[Ctop[j]];
            if (g1 == g2) {
                const el = document.getElementById(`fold_c${j}`);
                el.setAttribute("fill", clicked.has(g2) ? "yellow" : Ccolor[j]);
            }
        }
    },
    cell_click: (i, Ctop, FG, Ccolor, clicked) => {
        const g = FG[Ctop[i]];
        if (clicked.has(g)) {
            clicked.delete(g);
        } else {
            clicked.add(g);
        }
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
    PP_Ctop_CP_SC_2_Pvisible: (P, PP, Ctop, CP, SC) => {
        // computes boolean whether each vertex isvisible from top
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
    FV_V_line_eps_2_FV2_V2_VD2: (FV, V, line, eps) => {
        // assumes convex faces (or line divides a face into at most two pieces
        const [u, d] = line;
        const [a, b] = MAIN.line_2_coords(line);
        const nV = V.length;
        const V2 = V.map(v => v);
        const VD = V.map(v => {
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
            if (!neg && !pos) {
                throw new Exception("face has zero area?");
            }
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
                        xi = V2.length;
                        EV_map.set(s, xi);
                        V2.push(x);
                        VD.push(0);
                    }
                    pair[0].push(xi);
                    pair[1].push(xi);
                }
            }
            return pair;
        });
        return [FV2, V2, VD];
    },
    write: (FOLD) => {
        const {V, Vf, EV, EA, FV, FO} = FOLD;
        const path = document.getElementById("import").value.split("\\");
        const name = path[path.length - 1].split(".")[0];
        FOLD = {
            file_spec: 1.1,
            file_creator: "flat-folder",
            file_title: `${name}_state`,
            file_classes: ["singleModel"],
            vertices_coords:  V,
            edges_vertices:   EV,
            edges_assignment: EA,
            faces_vertices:   FV,
            faceOrders:       FO,
        };
        const data = new Blob([JSON.stringify(FOLD, undefined, 2)], {
            type: "application/json"});
        const ex = SVG.clear("export");
        const link = document.createElement("a");
        const button = document.createElement("input");
        ex.appendChild(link);
        link.appendChild(button);
        link.setAttribute("download", `${name}_state.fold`);
        link.setAttribute("href", window.URL.createObjectURL(data));
        button.setAttribute("type", "button");
        button.setAttribute("value", "Output State");
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
