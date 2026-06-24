import { M } from "./flatfolder/math.js";
import { NOTE } from "./flatfolder/note.js";
import { SVG } from "./flatfolder/svg.js";
import { IO } from "./flatfolder/io.js";
import { X } from "./flatfolder/conversion.js";
import { SOLVER } from "./flatfolder/solver.js";
import { CON } from "./flatfolder/constraints.js";

window.onload = () => { MAIN.startup(); };  // entry point

const TYPES = [
    "INVALID",
    "PURELAND",
    "INSIDE_REVERSE",
    "OUTSIDE_REVERSE",
    "MIXED_REVERSE",
    "OPEN_SINK",
    "CLOSED_SINK",
    "MIXED_SINK",
];
const TYPE = Object.fromEntries(TYPES.map((t, i) => [t, i]));

const COLOR = {
    normal: "black", select: "blue", active: "red",
    face: {select: "pink", active: "yellow", top: "#AAA", bottom: "#FFF"},
    edge: {
        U: "gray", F: "lightgray", Fcp: "gray",
        B: "black", V: "blue", M: "red"
    },
};
const OPACITY = {normal: 1, hidden: 0.01};
const LENGTH = {normal: 10, select: 20};
const STYLE = {
    apply: (el, style) => {
        for (const [k, v] of Object.entries(style)) { el.setAttribute(k, v); }
    },
    point_hidden: {fill: COLOR.normal, opacity: OPACITY.hidden, r: LENGTH.normal},
    point_select: {fill: COLOR.select, opacity: OPACITY.normal, r: LENGTH.select},
    point_active: {fill: COLOR.active, opacity: OPACITY.normal, r: LENGTH.normal},
    line_normal: {stroke: COLOR.normal, "stroke-width": LENGTH.normal},
    line_select: {stroke: COLOR.select, "stroke-width": LENGTH.select},
    line_active: {stroke: COLOR.active, "stroke-width": LENGTH.normal},
};

const MAIN = {
    startup: () => {
        CON.build();
        NOTE.clear_log();
        NOTE.start("*** Starting Flat-Folder ***");
        NOTE.time("Initializing interface");
        const [b, s] = [50, SVG.SCALE];
        const main = document.getElementById("main");
        for (const [k, v] of Object.entries({
            xmlns: SVG.NS, style: "background: lightgray",
            viewBox: [0, 0, 3*s, s].join(" "),
        })) { main.setAttribute(k, v); }

        for (const [i, id] of ["cp", "input", "output"].entries()) {
            const svg = document.getElementById(id);
            for (const [k, v] of Object.entries({
                xmlns: SVG.NS, height: s, width: s, x: i*s, y: 0,
                viewBox: [-b, -b, s + 2*b, s + 2*b].join(" "),
            })) { svg.setAttribute(k, v); }
        }
        const type_select = document.getElementById("type_select");
        for (const option of ["all", "pure"]) {
            const el = document.createElement("option");
            el.setAttribute("value", option);
            el.textContent = option;
            type_select.appendChild(el);
        }
        document.getElementById("import").onchange = (e) => {
            if (e.target.files.length > 0) {
                const file_reader = new FileReader();
                file_reader.onload = MAIN.process_file;
                file_reader.readAsText(e.target.files[0]);
            }
        };
        const V = [[0, 0], [0, 1], [1, 0], [1, 1]];
        const FV = [[0, 2, 3, 1]];
        const [FOLD, CELL] = MAIN.V_FV_2_FOLD_CELL(V, FV);
        FOLD.FO = [];
        MAIN.update_fold([[FOLD, CELL]]);
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
        const FS = (() => {
            const ex = JSON.parse(doc);
            const properties = [
                "vertices_coords", "faces_vertices",
                "faceOrders", "file_frames",
            ];
            const [V, FV, FO, frames] = properties.map(property => {
                const val = ex[property];
                if (val == undefined) {
                    NOTE.time(`FOLD file must contain ${property}, but not found`);
                }
                return val;
            });
            const FS = [];
            if (frames == undefined) {
                const [FOLD, CELL] = MAIN.V_FV_2_FOLD_CELL(V, FV);
                FOLD.FO = FO;
                FS.push([FOLD, CELL]);
            } else {
                for (const frame of frames) {
                    const [FOLD, CELL] = MAIN.V_FV_2_FOLD_CELL(
                        frame.vertices_coords,
                        frame.faces_vertices
                    );
                    FOLD.FR = frame["faces_lf:group"];
                    FOLD.FO = frame.faceOrders;
                    FOLD.lfL = frame["lf:line"];
                    FOLD.lfP = frame["lf:points"];
                    FS.push([FOLD, CELL]);
                }
            }
            return FS;
        })();
        MAIN.update_fold(FS);
    },
    update_fold: (FS) => {
        const [FOLD, CELL] = FS[FS.length - 1];
        document.getElementById("back_button").onclick = () => {
            if (FS.length == 1) { return; }
            FS.pop();
            MAIN.update_fold(FS);
        };
        const slider = document.getElementById("slider");
        slider.style.display = "none";
        slider.value = 0;
        for (const id of [
            "cycle", "replace", "fold_button", "state_controls", "state_config"
        ]) { document.getElementById(id).style.display = "none"; }
        const flip_el = document.getElementById("flip");
        SVG.clear("output");
        NOTE.time("Drawing State");
        MAIN.draw_state(SVG.clear("input"), FS);
        flip_el.onchange = () => {
            NOTE.start("Flipping model");
            MAIN.draw_state(SVG.clear("input"), FS);
            NOTE.end();
        };
        NOTE.time("Writing output");
        MAIN.write(FS);
        NOTE.end();
    },
    update_cp: (FOLD, LINE) => {
        const {V, FV, EV, EF, FE, Ff, FO, FOO, FM} = FOLD;
        const [H, EA] = MAIN.FO_Ff_EF_2_H_EA(FO, Ff, EF);
        FOLD.EA = EA;
        FOLD.Vf = X.V_FV_EV_EA_2_Vf_Ff(V, FV, EV, EA)[0];
        if (M.polygon_area2(M.expand(FOLD.FV[0], FOLD.Vf)) < 0) {
            FOLD.Vf = FOLD.Vf.map(v => M.add(M.refY(v), [0, 1]));
        }
        const v0 = FOLD.Vf[0];
        FOLD.Vf = FOLD.Vf.map(p => M.sub(p, v0));
        const [c1, s1] = FOLD.Vf[1];
        FOLD.Vf = FOLD.Vf.map(p => M.rotate_cos_sin(p, c1, -s1));
        FOLD.Vf = FOLD.Vf.map(p => M.rotate_cos_sin(p, 0, 1));
        FOLD.Vf = MAIN.V_P_transform(FOLD.Vf, FOLD.Vf);
        const Vf = FOLD.Vf;
        // rendering
        const cp = SVG.clear("cp");
        const faces = FV.map(F => M.expand(F, Vf));
        const lines = EV.map(E => M.expand(E, Vf));
        const colors = EA.map(a => {
            return (a == "F") ? COLOR.edge.Fcp : COLOR.edge[a];
        });
        const g1 = SVG.append("g", cp, {id: "flat_f"});
        SVG.draw_polygons(g1, faces, {fill: "white", id: true});
        const g2 = SVG.append("g", cp, {id: "flat_e"});
        SVG.draw_segments(g2, lines, {stroke: colors, id: true});
        SVG.append("g", cp, {id: "notes"});
        if (FM != undefined) {
            const [type, RF, FR, EC, V_sink, V_border] = MAIN.classify(
                V, EV, EF, FE, Ff, FM, FO, FOO);
            console.log(` ** Fold type: ${TYPES[type]}`);
            const g = document.getElementById("notes");
            const CLS = ["", "white", "black", "gray"];
            if (V_sink.length > 0) {
                const points = [];
                const colors = [];
                for (let i = 0; i < V_sink.length; ++i) {
                    const sink_type = V_sink[i];
                    if (sink_type == 0) { continue; }
                    points.push(Vf[i]);
                    colors.push(CLS[sink_type]);
                }
                SVG.draw_points(g, points, {r: 12, fill: "black"});
                SVG.draw_points(g, points, {r: 10, fill: colors});
            }
            if (V_border.length > 0) {
                const points = [];
                const colors = [];
                for (let i = 0; i < V_border.length; ++i) {
                    const reverse_type = V_border[i];
                    if (reverse_type == 0) { continue; }
                    points.push(Vf[i]);
                    colors.push(CLS[reverse_type]);
                }
                SVG.draw_points(g, points, {r: 8, fill: "black"});
                SVG.draw_points(g, points, {r: 6, fill: colors});
            }
            for (let i = 0; i < RF.length; ++i) {
                const hue = (i*137) % 360; // Approx Golden Angle Method
                const color = `hsl(${hue}, ${
                    (type == TYPE.INVALID) ? 30 : 100
                }%, 85%)`;
                for (const f of RF[i]) {
                    const el = document.getElementById(`flat_f${f}`);
                    el.setAttribute("fill", color);
                }
            }
            for (let i = 0; i < EC.length; ++i) {
                if (!EC[i]) { continue; }
                const el = document.getElementById(`flat_e${i}`);
                el.setAttribute("stroke-width", LENGTH.normal);
            }
            FOLD.FR = FR.map(l => l ?? -1);
            // if (FOLD.L != undefined) { // draw linearized order
            //     const L = FOLD.L.flat();
            //     const lines = [];
            //     for (let i = 1; i < L.length; ++i) {
            //         const pi = M.centroid(M.expand(FV[L[i - 1]], Vf));
            //         const pj = M.centroid(M.expand(FV[L[i]], Vf));
            //         lines.push([pi, pj]);
            //     }
            //     SVG.draw_segments(g, lines, {});
            // }
        }
        if (LINE == undefined) { return; }
        const clicked = LINE.clicked_groups;
        const FG = LINE.FG;
        for (let i = 0; i < FG.length; ++i) {
            const el = document.getElementById(`flat_f${i}`);
            el.setAttribute("fill", clicked.has(FG[i]) ? "yellow" : "white");
        }
        for (let i = 0; i < EV.length; ++i) {
            if (EF[i].length < 2) { continue; }
            const [f, g] = EF[i];
            if (FG[f] == FG[g]) { continue; }
            const el = document.getElementById(`flat_e${i}`);
            el.setAttribute("stroke-width", LENGTH.normal);
        }
    },
    line_click: (el, lfP, lfL, FS) => {
        const [FOLD_, CELL_] = FS[FS.length - 1];
        const [FV2, V, Vf, VD] = MAIN.FV_V_Vf_lfL_eps_2_FV2_V2_Vf2_VD2(
            FOLD_.FV, FOLD_.V, FOLD_.Vf, lfL, FOLD_.eps);
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
        const FG = (() => {
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
        })();
        const [FOLD, CELL] = MAIN.V_FV_2_FOLD_CELL(V, FV); // fully divided
        FOLD.Vf = Vf;
        const {BF, BI, CF} = CELL;
        const FO = [];
        for (const [f, g, o] of FOLD_.FO) {
            for (const f_ of F_map[f]) {
                for (const g_ of F_map[g]) {
                    const pair = M.encode_order_pair([f_, g_]);
                    if (BI.has(pair)) {
                        FO.push([f_, g_, o]);
                    }
                }
            }
        }
        FOLD.FO = FO;
        const svg = SVG.clear("input");
        const slider = document.getElementById("slider");
        slider.value = 0;
        slider.setAttribute("max", FV.length);
        const clicked_groups = new Set();
        const LINE = {el, FG, clicked_groups};
        FS.push([FOLD, CELL]);
        MAIN.draw_state(svg, FS, LINE);
        const fold_button = document.getElementById("fold_button");
        fold_button.style.display = "inline";
        fold_button.onclick = async () => {
            if (clicked_groups.size == 0) {
                FS.pop();
                MAIN.update_fold(FS);
            } else {
                SVG.clear("input");
                document.getElementById("slider").value = 0;
                MAIN.draw_state(svg, FS);
                svg.appendChild(el);
                MAIN.make_fold(V, FOLD_.FV, FV, F_map,
                    FG, FOLD_.FO, clicked_groups, lfP, lfL, FS);
            }
        };
    },
    draw_state: (svg, FS, LINE) => {
        const [FOLD, CELL] = FS[FS.length - 1];
        const {Ff, EF, FO} = FOLD;
        const {P, PP, CP, CF, SP, SC, SE} = CELL;
        const Q = MAIN.V_P_transform(P, P);
        const edges = FO.map(([f1, f2, o]) => {
            return M.encode(((Ff[f2] ? 1 : -1)*o >= 0) ? [f1, f2] : [f2, f1]);
        });
        const [L, LL] = MAIN.linearize(edges, Ff.length);
        FOLD.L = L;
        const slider = document.getElementById("slider");
        if (LL != undefined) {
            slider.setAttribute("max", LL.length);
        }
        const CD = X.CF_edges_2_CD(CF, edges); // CF ordered in state
        const flip = document.getElementById("flip").checked;
        const Ctop = CD.map(S => flip ? S[0] : S[S.length - 1]);
        const Ccolor = Ctop.map(d => {
            if (d == undefined) { return undefined; }
            if (Ff[d] != flip)  { return COLOR.face.top; }
            else                { return COLOR.face.bottom; }
        });
        MAIN.update_cp(FOLD, LINE);
        {   // update linearized state if exists
            const slider = document.getElementById("slider");
            const cycle = document.getElementById("cycle");
            if (LL == undefined) {
                slider.style.display = "none";
                cycle.style.display = "inline";
            } else {
                slider.style.display = "inline";
                cycle.style.display = "none";
                slider.oninput = () => {
                    SVG.clear(svg.id);
                    MAIN.draw_state(svg, FS, LINE);
                };
                const val = +slider.value;
                const F_set = new Set();
                for (let i = 0; i < Ff.length; ++i) { F_set.add(i); }
                const flip = document.getElementById("flip").checked;
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
                for (let i = 0; i < Ccolor.length; ++i) {
                    const d = Ctop[i];
                    Ccolor[i] = (d == undefined) ? undefined : (
                        (Ff[d] != flip) ? COLOR.face.top : COLOR.face.bottom);
                }
            }
        }
        const SD = X.Ctop_SC_SE_EF_Ff_2_SD(Ctop, SC, SE, EF, Ff);
        const cells = CP.map(V => M.expand(V, Q));
        const fold_c = SVG.append("g", svg, {id: "fold_c"});
        const fold_s_crease = SVG.append("g", svg, {id: "fold_s_crease"});
        const fold_s_edge = SVG.append("g", svg, {id: "fold_s_edge"});
        const Lsvg = SVG.append("g", svg, {id: "lines"});
        const fold_p = SVG.append("g", svg, {id: "fold_p"});
        SVG.draw_polygons(fold_c, cells, {
            id: true, fill: Ccolor, stroke: Ccolor});
        const lines = SP.map((ps) => M.expand(ps, Q));
        SVG.draw_segments(fold_s_crease, lines, {
            id: true, stroke: COLOR.edge.F,
            filter: (i) => SD[i][0] == "C"});
        SVG.draw_segments(fold_s_edge, lines, {
            id: true, stroke: COLOR.edge.B,
            filter: (i) => SD[i][0] == "B"});
        if (svg.id != "input") { return; } // selection interface only on input
        if (LINE == undefined) {
            const Pvisible = MAIN.P_PP_Ctop_CP_SC_2_Pvisible(P, PP, Ctop, CP, SC);
            SVG.draw_points(fold_p, Q, {id: true, filter: (i) => Pvisible[i]});
            const lfP = new Set(); // "lf:points"
            for (let i = 0; i < P.length; ++i) {
                const el = document.getElementById(`fold_p${i}`);
                if (el == undefined) { continue; }
                el.onmouseout  = () => STYLE.apply(el, STYLE["point_" +
                    (lfP.has(i) ? "active" : "hidden")
                ]);
                el.onmouseout();
                el.onmouseover = () => STYLE.apply(el, STYLE.point_select);
                el.onclick     = () => { MAIN.point_click(i, lfP, FS); }
            }
        } else {
            const {el, FG, clicked_groups} = LINE;
            Lsvg.appendChild(el);
            el.onmouseout = () => STYLE.apply(el, STYLE.line_active);
            el.onmouseout();
            el.onclick = () => { MAIN.update_fold(FS); }
            const group_over = (g, Ctop, FG) => {
                for (let j = 0; j < Ctop.length; ++j) {
                    if (g != FG[Ctop[j]]) { continue; }
                    document.getElementById(`fold_c${j}`)
                        .setAttribute("fill", COLOR.face.select);
                }
                for (let j = 0; j < FG.length; ++j) {
                    if (g != FG[j]) { continue; }
                    document.getElementById(`flat_f${j}`)
                        .setAttribute("fill", COLOR.face.select);
                }
            };
            const group_out = (g1, Ctop, FG, Ccolor, clicked) => {
                for (let j = 0; j < Ctop.length; ++j) {
                    const g2 = FG[Ctop[j]];
                    if (g1 != g2) { continue; }
                    const el = document.getElementById(`fold_c${j}`);
                    el.setAttribute("fill",
                        clicked.has(g2) ? COLOR.face.active : Ccolor[j]);
                }
                for (let j = 0; j < FG.length; ++j) {
                    const g2 = FG[j];
                    if (g1 != g2) { continue; }
                    const el = document.getElementById(`flat_f${j}`);
                    el.setAttribute("fill", clicked.has(g2)
                        ? COLOR.face.active : COLOR.face.bottom);
                    // cp always white-side up
                }
            };
            const group_click = (g, clicked) => {
                clicked.has(g) ? clicked.delete(g) : clicked.add(g);
            };
            for (let i = 0; i < CF.length; ++i) {
                const el = document.getElementById(`fold_c${i}`);
                if (el == undefined) { continue; }
                const g = FG[Ctop[i]];
                el.onmouseover = () => group_over(g, Ctop, FG);
                el.onmouseout = () => group_out(
                    g, Ctop, FG, Ccolor, clicked_groups);
                el.onclick = () => group_click(g, clicked_groups);
                el.setAttribute("fill",
                    clicked_groups.has(g) ? "yellow" : Ccolor[i]);
            }
            for (let i = 0; i < FG.length; ++i) {
                const el = document.getElementById(`flat_f${i}`);
                const g = FG[i];
                el.onmouseover = () => group_over(g, Ctop, FG);
                el.onmouseout = () => group_out(
                    g, Ctop, FG, Ccolor, clicked_groups);
                el.onclick = () => group_click(FG[i], clicked_groups);
            }
        }
    },
    point_click: (i, lfP, FS) => {
        const [FOLD, CELL] = FS[FS.length - 1];
        const {eps} = FOLD;
        const {P} = CELL;
        if (lfP.has(i) || (lfP.size > 3)) { lfP.clear(); }
        else                              { lfP.add(i); }
        const L = MAIN.get_lines(Array.from(lfP).map(i => P[i]), eps);
        for (let i = 0; i < P.length; ++i) {    // rendering
            const el = document.getElementById(`fold_p${i}`);
            if (el == undefined) { continue; }
            STYLE.apply(el, lfP.has(i)
                ? STYLE.point_active : STYLE.point_hidden);
        }
        const svg = SVG.clear("lines");
        if (L.length == 0) { return; }
        SVG.draw_segments(svg,
            L.map(l => MAIN.V_P_transform(MAIN.line_2_coords(l), P)),
            {id: true}
        );
        for (let j = 0; j < L.length; ++j) {    // interface
            const el = svg.children[j];
            el.onmouseout = () => STYLE.apply(el, STYLE.line_normal);
            el.onmouseout();
            el.onmouseover = () => STYLE.apply(el, STYLE.line_select);
            el.onclick = () => MAIN.line_click(el, Array.from(lfP), L[j], FS);
        }
    },
    V_P_transform: (V, P) => {
        // rescales V based on normalizing P into a [0, 1] bounding box
        const [p_min, p_max] = M.bounding_box(P);
        const [x_diff, y_diff] = M.sub(p_max, p_min);
        const is_tall = (x_diff < y_diff);
        const diff = is_tall ? y_diff : x_diff;
        const m = [0.5, 0.5];
        const off = M.sub(m, M.div([x_diff, y_diff], 2*diff));
        const flip = document.getElementById("flip").checked;
        return V.map(p => M.add(M.div(M.sub(p, p_min), diff), off))
                .map(p => (flip ? M.add(M.refX(M.sub(p, m)), m) : p));
    },
    draw_separators: (V, FV, FO, FM) => {
        const [FOLD, _] = MAIN.V_FV_2_FOLD_CELL(V, FV);
        const {Ff, EF} = FOLD;
        FOLD.FO = FO;
        MAIN.update_cp(FOLD);
        const [H, EA] = MAIN.FO_Ff_EF_2_H_EA(FO, Ff, EF);
        const EF_map = new Map();
        for (const [i, F] of FV.entries()) {
            for (const [j, v1] of F.entries()) {
                const v2 = F[(j + 1) % F.length];
                EF_map.set(M.encode([v2, v1]), i);
            }
        }
        const E_map = new Map();
        for (let i = 0; i < FOLD.EV.length; ++i) {
            E_map.set(M.encode_order_pair(FOLD.EV[i]), i);
        }
        const V_boundary = V.map(() => 0);
        for (let i = 0; i < FM.length; ++i) {
            const c = FM[i] ? 1 : 2;
            for (const v of FV[i]) { V_boundary[v] |= c; }
        }
        for (let i = 0; i < V.length; ++i) {
            for (const j of FOLD.VV[i]) {
                if (EA[E_map.get(M.encode_order_pair([i, j]))] == "B") {
                    V_boundary[i] |= 2;
                }
            }
        }
        SVG.draw_points(document.getElementById("notes"), FOLD.Vf, {
            filter: (i) => (V_boundary[i] == 3), r: 5, fill: "green", id: true,
        });
        const E_sep = EF.map(() => false);
        for (let i = 0; i < V.length; ++i) {
            if (V_boundary[i] != 3) { continue; }
            for (const j of FOLD.VV[i]) {
                const f_left = EF_map.get(M.encode([i, j]));
                const f_right = EF_map.get(M.encode([j, i]));
                if ((f_left == undefined) || (!FM[f_left])) { continue; }
                if ((f_right == undefined) || (!FM[f_right])) { continue; }
                // found entering edge
                const o = H.get(M.encode([f_left, f_right]));
                if (o == undefined) { continue; }
                const P = [i];
                const processing = V.map(() => false);
                processing[i] = true;
                const dfs = (v) => {
                    if (processing[v]) { return; }
                    console.log("processing: ", v);
                    P.push(v);
                    processing[v] = true;
                    if (V_boundary[v] == 3) {
                        console.log("checking");
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
                            // console.log(P);
                            // console.log("FOUND!");
                        }
                    } else {
                        for (const u of FOLD.VV[v]) {
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
                    el.setAttribute("stroke-width", LENGTH.normal);
                }
            }
        }
    },
    make_fold: (V, FV_, FV, F_map_old, FG, FO_, clicked_groups, lfP, lfL, FS) => {
        // F_map_old maps uncut faces to their face indices after cutting
        // FV_ is old uncut faces, while FV is new cut faces
        // FO_ is old orders for uncut faces
        // clicked_groups is set mapping group indices to Bool whether clicked
        const F_map = F_map_old.map(() => []);
        // map only clicked regions
        const FVx = []; // new vertices per face
        const FM = [];  // boolean whether face moves
        for (let fi = 0; fi < F_map_old.length; ++fi) {
            const F = F_map_old[fi];
            if (F.length == 2) { // face was cut
                const [f1, f2] = F;
                const g1 = clicked_groups.has(FG[f1]);
                const g2 = clicked_groups.has(FG[f2]);
                if (g1 == g2) {  // same group, face was not cut
                    F_map[fi].push(FVx.length);
                    FVx.push(FV_[fi]);
                    FM.push(g1);
                } else {         // face was cut, so split them
                    F_map[fi].push(FVx.length);
                    FVx.push(FV[f1]);
                    FM.push(g1);
                    F_map[fi].push(FVx.length);
                    FVx.push(FV[f2]);
                    FM.push(g2);
                }
            } else if (F.length == 1) { // face was not cut
                const g1 = clicked_groups.has(FG[F[0]]);
                F_map[fi].push(FVx.length);
                FVx.push(FV_[fi]);
                FM.push(g1);
            } else {
                throw new Error("malformed FV2");
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
        const [u, d] = lfL;
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
        for (const [f, g, o] of FO_) {
            for (const f_ of F_map[f]) {
                for (const g_ of F_map[g]) {
                    const pair = M.encode_order_pair([f_, g_]);
                    FOO.push([f_, g_, o]);
                }
            }
        }

        // MAIN.draw_separators(Vx, FVy, FOO, FM);
        // return;

        const [FOLD, CELL] = MAIN.V_FV_2_FOLD_CELL(Vy, FVy);
        const {Ff, EF} = FOLD;
        const {P, PP, CP, FC, CF, SE, SC, SP, BF, BI} = CELL;
        const type = document.getElementById("type_select").value;
        const FO = [];
        for (const [f, g, o] of FO_) {
            for (const f_ of F_map[f]) {
                for (const g_ of F_map[g]) {
                    const pair = M.encode_order_pair([f_, g_]);
                    if (BI.has(pair)) { // new overlaps
                        if (type == "pure") {
                            if (FM[f_] == FM[g_]) {
                                FO.push([f_, g_, o]);
                            }
                        } else {
                            if (!FM[f_] && !FM[g_]) {
                                FO.push([f_, g_, o]);
                            }
                        }
                    }
                }
            }
        }
        FOLD.FM = FM;
        FOLD.FOO = FOO;
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
        const BA0 = (() => {
            const BA_map = new Map();
            for (let i = 0; i < FO.length; ++i) {
                const [f, g, o] = FO[i];
                const a1 = (Ff[g] == (o > 0)) ? 1 : 2;
                const a2 = (f < g) ? a1 : ((a1 == 1) ? 2 : 1);
                const s = M.encode_order_pair([f, g]);
                BA_map.set(s, a2);
            }
            return BF.map(s => {
                const out = BA_map.get(s);
                return (out == undefined) ? 0 : out;
            });
        })();
        const trans_count = {all: 0, reduced: 0};
        const BA = SOLVER.initial_assignment(BA0, BF, BT, BI,
            FC, CF, CC, trans_count);
        if ((BA.length == 3) && (BA[1].length != undefined)) {
            const [type, F, E] = BA;
            const str = `Unable to resolve ${CON.names[type]} on faces [${F}]`;
            NOTE.log(`   - ${str}`);
            NOTE.log(`   - Faces participating in conflict: [${E}]`);
            NOTE.end();
            FS.pop();
            MAIN.update_fold(FS);
            return;
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
        if (GA.length == undefined) {
            const gi = GA;
            const F = new Set();
            for (const bi of GB[gi]) {
                for (const f of M.decode(BF[bi])) { F.add(f); }
            }
            postMessage({
                type: "end",
                arg: ["component_error", [gi, Array.from(F)]]
            });
            return;
        }
        const Gn = GA.map(A => A.length);
        NOTE.time("Solve completed");
        const n = Gn.reduce((s, gn) => s*BigInt(gn), BigInt(1));
        NOTE.count(n, "folded states");
        NOTE.lap();
        if (n == 0) {
            NOTE.end();
            FS.pop();
            MAIN.update_fold(FS);
            return;
        }
        const GI = GB.map(() => 0);
        document.getElementById("state_controls").style.display = "inline";
        const comp_select = SVG.clear("component_select");
        for (const [i, _] of GB.entries()) {
            const el = document.createElement("option");
            el.setAttribute("value", `${i}`);
            el.textContent = `${i}`;
            if (i == 1) {
                el.setAttribute("selected", true);
            }
            comp_select.appendChild(el);
        }
        const SOLUTION = {GB, GA, GI};
        FS.pop();
        FS.push([FOLD, CELL]);
        FOLD.lfP = lfP;
        FOLD.lfL = lfL;
        comp_select.onchange = () => {
            NOTE.start("Changing component");
            MAIN.update_component(FS, SOLUTION);
            NOTE.end();
        };
        MAIN.update_component(FS, SOLUTION);
        NOTE.lap();
        stop = Date.now();
        NOTE.end();
    },
    update_component: (FS, SOLUTION) => {
        const [FOLD, CELL] = FS[FS.length - 1];
        const {BF} = CELL;
        const {GB, GA, GI} = SOLUTION;
        const comp_select = document.getElementById("component_select");
        const c = comp_select.value;
        document.getElementById("state_config").style.display = "none";
        const C = [];
        C.push(c);
        const n = GA[c].length;
        document.getElementById("state_config").style.display = "inline";
        const state_label = document.getElementById("state_label");
        const state_select = document.getElementById("state_select");
        state_label.innerHTML = `${n} State${(n == 1) ? "" : "s"}`;
        state_select.setAttribute("min", 1);
        state_select.setAttribute("max", n);
        state_select.value = GI[c] + 1;
        state_select.onchange = () => {
            NOTE.start("Computing new state");
            let j = +state_select.value;
            if (j < 1) { j = 1; }
            if (j > n) { j = n; }
            state_select.value = j;
            GI[c] = j - 1;
            NOTE.time("Computing state");
            const edges = X.BF_GB_GA_GI_2_edges(BF, GB, GA, GI);
            FOLD.FO = X.edges_Ff_2_FO(edges, FOLD.Ff);
            NOTE.time("Drawing fold");
            const svg = SVG.clear("output");
            MAIN.draw_state(svg, FS);
            const flip_el = document.getElementById("flip");
            flip_el.onchange = () => {
                NOTE.start("Flipping model");
                MAIN.draw_state(SVG.clear("output"), FS);
                NOTE.end();
            };
            const replace = document.getElementById("replace");
            replace.style.display = "inline";
            replace.onclick = () => {
                FOLD.FOO = undefined;
                FOLD.FM = undefined;
                MAIN.update_fold(FS);
            };
        };
        state_select.onchange();
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
        const FOLD = {V, VV, FV, EV, EF, FE, Ff, eps};
        const CELL = {P, SP, SE, PP, CP, CS, SC, CF, FC, BF, BI};
        return [FOLD, CELL];
    },
    linearize: (edges, n) => {
        const Adj = Array(n).fill(0).map(() => []);
        for (const s of edges) {
            const [f1, f2] = M.decode(s);
            Adj[f1].push(f2);
        }
        const L = [];
        const seen = Array(n).fill(false);
        const dfs = (i) => {
            if (seen[i]) { return; }
            seen[i] = true;
            for (const j of Adj[i]) {
                dfs(j);
            }
            L.push(i);
        };
        for (let i = 0; i < n; ++i) {
            dfs(i);
        }
        L.reverse();
        console.assert(L.length == n);
        const idx_map = Array(n).fill(undefined);
        for (let i = 0; i < n; ++i) {
            const fi = L[i];
            idx_map[fi] = i;
        }
        for (const s of edges) {
            const [f1, f2] = M.decode(s);
            if (idx_map[f1] > idx_map[f2]) {
                return [undefined, undefined]; // cycle
            }
        }
        for (let i = 0; i < n; ++i) {
            seen[i] = false;
        }
        const layers = [];
        for (let i = 0; i < n; ++i) {
            const fi = L[i];
            if (seen[fi]) { continue; }
            seen[fi] = true;
            const layer = [fi];
            const Adj_set = new Set();
            for (const fj of Adj[fi]) {
                Adj_set.add(fj);
            }
            for (let j = i + 1; j < L.length; ++j) {
                const fj = L[j];
                if (seen[fj]) { continue; }
                if (!Adj_set.has(fj)) {
                    seen[fj] = true;
                    layer.push(fj);
                }
                for (const fk of Adj[fj]) {
                    Adj_set.add(fk);
                }
            }
            layers.push(layer);
        }
        return [L, layers];
    },
    P_PP_Ctop_CP_SC_2_Pvisible: (P, PP, Ctop, CP, SC) => {
        // computes boolean whether each vertex is visible from top
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
            if (F_set.size > 2) { // degree >= 3
                return true;
            }
            for (const [f, ang] of F_set.entries()) { // bent corner
                if (Math.abs(Math.abs(ang) - Math.PI) > 0.001) {
                    return true;
                }
            }
            return false;
        });
    },
    FV_V_Vf_lfL_eps_2_FV2_V2_Vf2_VD2: (FV, V, Vf, lfL, eps) => {
        // assumes convex faces (or line divides a face into at most two pieces
        const [u, d] = lfL;
        const [a, b] = MAIN.line_2_coords(lfL);
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
        return [FV2, V2, Vf2, VD];
    },
    write: (FS) => {
        const frames = [];
        for (const [FOLD, _] of FS) {
            const {V, Vf, FV, FO, FR, lfP, lfL} = FOLD;
            const frame_FOLD = {
                vertices_coords:  V,
                faces_vertices:   FV,
                faceOrders:       FO,
                "faces_lf:group": FR,
                "lf:points":      lfP,
                "lf:line":        lfL,
            };
            frames.push(frame_FOLD);
        }
        const [FOLD, CELL] = FS[FS.length - 1];
        const {V, Vf, EV, EA, FV, FO, FR} = FOLD;
        const path = document.getElementById("import").value.split("\\");
        const name = path[path.length - 1].split(".")[0];
        const export_seq = {
            file_spec: 1.1,
            file_creator: "line-folder",
            file_title: `${name}_state`,
            file_classes: ["singleModel"],
            vertices_coords:  V,
            faces_vertices:   FV,
            faceOrders:       FO,
            file_frames:  frames,
        };
        const seq_data = new Blob([JSON.stringify(export_seq, undefined, 2)], {
            type: "application/json"});
        const seq_link = document.getElementById("seq_anchor");
        seq_link.setAttribute("download", `state.fold`);
        seq_link.setAttribute("href", window.URL.createObjectURL(seq_data));
        seq_link.style.textDecoration = "none";
        const export_cp = {
            file_spec: 1.1,
            file_creator: "line-folder",
            file_title: `${name}_cp`,
            file_classes: ["singleModel"],
            vertices_coords:  Vf,
            faces_vertices:   FV,
            edges_vertices:   EV,
            edges_assignment: EA,
        };
        const cp_data = new Blob([JSON.stringify(export_cp, undefined, 2)], {
            type: "application/json"});
        const cp_link = document.getElementById("cp_anchor");
        cp_link.setAttribute("download", `cp.fold`);
        cp_link.setAttribute("href", window.URL.createObjectURL(cp_data));
        cp_link.style.textDecoration = "none";
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
        for (const [k, o] of H2) {
            const [f, g] = M.decode(k);
            if ((!FM[f]) || (!FM[g])) { continue; }
            if (H1.get(k) == o) { continue; }
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
    classify: (V, EV, EF, FE, Ff, FM, FO, FOO) => {
        const [H1, EA1] = MAIN.FO_Ff_EF_2_H_EA(FO, Ff, EF);
        const [H2, EA2] = MAIN.FO_Ff_EF_2_H_EA(FOO, Ff, EF);
        const EC = EA1.map((a, i) => (a != EA2[i]));
        const HC = MAIN.FM_H1_H2_2_HC(FM, H1, H2);
        const [FR, RF] = MAIN.EF_FM_HC_2_FR_RF(EF, FM, HC);
        // uniform: orders within color don't change
        const uniform = MAIN.FM_FR_H2_HC_2_uniform(FM, FR, H2, HC);
        // separable: orders between pairs of colors are consistent
        const separable = MAIN.FR_RF_H1_H2_HC_2_separable(FR, RF, H1, H2, HC);
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
        const local = MAIN.FM_FE_EF_RF_2_local(FM, FE, EF, RF);
        // connected: colors form a connected component
        const connected = MAIN.FM_FR_RF_EF_2_connected(FM, FR, RF, EF);
        const valid = (uniform && separable && joint && local && connected);
        // // OLD: reverse fold if adjacency graph is chain and
        // // no adjacency is closed.
        // equivilently, we can check that no vertex is "popped",
        // i.e., has folded degree >2 and switches majority assignment
        if (!valid) {
            console.log(`fold invalid: ${
                uniform   ? "" : "!uniform, "   }${
                separable ? "" : "!separable, " }${
                joint     ? "" : "!joint, "     }${
                local     ? "" : "!local, "     }${
                connected ? "" : "!connected"
            }`);
        }
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
            for (const e of E) {
                const a_new = EA1[e];
                if (a_new == "B") { return 0; } // on paper boundary
                const a_old = EA2[e];
                if (a_new == a_old) { continue; }
                const [f, g] = EF[e];
                if (FM[f] == FM[g]) {
                    inner.push(a_new);
                } else {
                    outer.push(a_new);
                }
            }
            if ((inner.length != 1) ||
                (outer.length != 2) ||
                (outer[0] != outer[1])
            ) { return 0; }
            if (inner[0] != outer[0]) {
                is_outside = false;
                return 1;
            } else {
                is_inside = false;
                return 2;
            }
        });
        let type;
        if (!valid) {
            type = TYPE.INVALID;
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
