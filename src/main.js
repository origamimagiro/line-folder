import { M } from "./flatfolder/math.js";
import { NOTE } from "./flatfolder/note.js";
import { SVG } from "./flatfolder/svg.js";
import { X } from "./flatfolder/conversion.js";
import { CON } from "./flatfolder/constraints.js";

import { TYPE_LABEL, TYPE, COMP } from "./compute.js";
import { IO } from "./io.js";

window.onload = () => { MAIN.startup(); };  // entry point

const COLOR = {
    normal: "black", select: "blue", active: "red",
    face: {
        select: "yellow", active: "hsl(0, 100%, 85%)",
        top: "#AAA", bottom: "#FFF"
    },
    edge: {
        U: "gray", F: "lightgray", Fcp: "gray",
        B: "black", V: "blue", M: "red"
    },
};
const OPACITY = {normal: 1, hidden: 0.01};
const LENGTH = {normal: 1, bold: 3, half: 6, active: 12, select: 18};
const STYLE = {
    apply: (el, style) => {
        for (const [k, v] of Object.entries(style)) { el.setAttribute(k, v); }
    },
    point_hidden: {fill: COLOR.normal, opacity: OPACITY.hidden, r: LENGTH.active},
    point_select: {fill: COLOR.select, opacity: OPACITY.normal, r: LENGTH.select},
    point_active: {fill: COLOR.active, opacity: OPACITY.normal, r: LENGTH.active},
    line_normal: {stroke: COLOR.normal, "stroke-width": LENGTH.active},
    line_select: {stroke: COLOR.select, "stroke-width": LENGTH.select},
    line_active: {stroke: COLOR.active, "stroke-width": LENGTH.active},
};

export const MAIN = {
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
        const mode_select = document.getElementById("mode_select");
        for (const option of ["all", "select"]) {
            const el = document.createElement("option");
            el.setAttribute("value", option);
            el.textContent = option;
            mode_select.appendChild(el);
        }
        const replace_file = (FS) => {
            document.getElementById("back_button").onclick = () => {
                if (FS.length == 1) { return; }
                FS.pop();
                MAIN.update_interface(FS);
            };
            document.getElementById("replace").onclick = () => {
                const [FOLD, CELL] = FS[FS.length - 1];
                FOLD.FOO = undefined;
                FOLD.FM = undefined;
                MAIN.update_interface(FS);
            };
            document.getElementById("flip").onchange = () => {
                NOTE.start("Flipping model");
                if (!FS[FS.length - 1][0].fixed) { FS.pop(); }
                const [FOLD, CELL] = FS[FS.length - 1];
                MAIN.draw_cp(FOLD);
                MAIN.draw_state((FOLD.FM == undefined) ? "input" : "output", FS);
                NOTE.end();
            };
            MAIN.update_interface(FS);
        };
        document.getElementById("import").onchange = (e) => {
            if (e.target.files.length > 0) {
                const file_reader = new FileReader();
                file_reader.onload = (e) => {
                    const FS = IO.process_file(e);
                    replace_file(FS);
                }
                file_reader.readAsText(e.target.files[0]);
            }
        };
        const V = [[0, 0], [0, 1], [1, 0], [1, 1]];
        const FV = [[0, 2, 3, 1]];
        const [FOLD, CELL] = COMP.V_FV_2_FOLD_CELL(V, FV);
        FOLD.FO = [];
        COMP.augment_FOLD_FO(FOLD);
        const FS = [[FOLD, CELL]];
        replace_file(FS);
    },
    update_interface: (FS) => {
        const [FOLD, CELL] = FS[FS.length - 1];
        FOLD.type = undefined;
        FOLD.RF = undefined;
        FOLD.FR = undefined;
        FOLD.V_border = undefined;
        FOLD.V_sink = undefined;
        FOLD.EC = undefined;
        const slider = document.getElementById("slider");
        slider.style.display = "none";
        slider.value = 0;
        for (const id of [
            "cycle", "replace", "fold_button", "state_controls", "state_config"
        ]) { document.getElementById(id).style.display = "none"; }
        NOTE.time("Writing output");
        IO.write(FS);
        NOTE.time("Drawing State");
        SVG.clear("output");
        MAIN.draw_cp(FOLD);
        MAIN.draw_state("input", FS);
        NOTE.end();
    },
    draw_cp: (FOLD, LINE) => {
        const {Vf, FV, EV, EA, EF} = FOLD;
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
        MAIN.update_cp(FOLD, LINE);
    },
    line_click: (el, lfP, lfL, FS) => {
        const [FOLD_, CELL_] = FS[FS.length - 1];
        const [FV, V, Vf, VD, F_map] = COMP.split_FOLD_on_line(FOLD_, lfL);
        const [FOLD, CELL] = COMP.V_FV_2_FOLD_CELL(V, FV); // fully divided
        const FG = COMP.get_groups(FV, VD);
        FOLD.Vf = Vf;
        const {BF, BI, CF} = CELL;
        FOLD.FO = COMP.map_order(BI, F_map, FOLD_.FO);
        COMP.augment_FOLD_FO(FOLD);
        const slider = document.getElementById("slider");
        slider.value = 0;
        slider.setAttribute("max", FV.length);
        const clicked_groups = new Set();
        const clicked_edges = new Map();
        const LINE = {el, FG, clicked_groups, clicked_edges};
        FS.push([FOLD, CELL]);
        MAIN.draw_cp(FOLD, LINE);
        MAIN.draw_state("input", FS, LINE);
        const fold_button = document.getElementById("fold_button");
        fold_button.style.display = "inline";
        fold_button.onclick = async () => {
            if (clicked_groups.size == 0) {
                FS.pop();
                MAIN.update_interface(FS);
                return;
            }
            const mode = document.getElementById("mode_select").value;
            let nFM = 0;
            for (const g of FG) {
                if (clicked_groups.has(g)) { ++nFM; }
            }
            // if ((mode == "all") && (nFM > 20)) {
            //     console.log("Too many faces?");
            //     return;
            // }
            SVG.clear("input");
            document.getElementById("slider").value = 0;
            MAIN.draw_state("input", FS);
            document.getElementById("input").appendChild(el);
            el.onmouseover = undefined;
            el.onmouseout = undefined;
            el.onclick = undefined;
            MAIN.make_fold(V, FOLD_.FV, FV, F_map,
                FG, FOLD_.FO, clicked_groups, lfP, lfL, FS);
        };
    },
    update_cp: (FOLD, LINE) => {
        const {EV, EA, EF} = FOLD;
        const {FG, clicked_groups, clicked_edges, over_group} = (
            LINE ?? {
                FG: Array(FOLD.Ff.length).fill(0),
                clicked_groups: new Set(),
                clicked_edges: new Map(),
                over_group: undefined,
            }
        );
        const cp = document.getElementById("cp");
        for (let i = 0; i < FG.length; ++i) {
            const el = cp.getElementById(`flat_f${i}`);
            const g = FG[i];
            el.setAttribute("fill",
                (over_group == g)     ? COLOR.face.select : (
                clicked_groups.has(g) ? COLOR.face.active
                                      : COLOR.face.bottom)
            );
        }
        for (let i = 0; i < EV.length; ++i) {
            if (EF[i].length < 2) { continue; }
            const [f, g] = EF[i];
            const el = document.getElementById(`flat_e${i}`);
            el.setAttribute("stroke-width",
                (FG[f] == FG[g]) ? LENGTH.normal : LENGTH.half);
        }
        let FR = FOLD.FR;
        let RF = FOLD.RF;
        const mode = document.getElementById("mode_select").value;
        if (LINE != undefined) {
            const FM = FG.map(g => clicked_groups.has(g));
            const HC = new Set();
            for (const [e, _] of clicked_edges) {
                HC.add(M.encode_order_pair(EF[e]));
            }
            [FR, RF] = COMP.EF_FM_HC_2_FR_RF(EF, FM, HC);
            FOLD.FR = FR;
            FOLD.RF = RF;
            if (mode == "select") {
                const [V_boundary, E_sep] = COMP.find_separators(FOLD, FM);
                const reset_edges = () => {
                    for (let e = 0; e < EF.length; ++e) {
                        if (!((EA[e] == "M") || (EA[e] == "V"))) { continue; }
                        const [f, g] = EF[e];
                        if ((!clicked_groups.has(FG[f])) ||
                            (!clicked_groups.has(FG[g]))) { continue; }
                        document.getElementById(`flat_e${e}`)
                            .setAttribute("stroke-width",
                            clicked_edges.has(e) ? LENGTH.half : (
                                (E_sep[e].size > 0) ? LENGTH.bold : LENGTH.normal));
                    }
                };
                for (let i = 0; i < EF.length; ++i) {
                    if (!((EA[i] == "M") || (EA[i] == "V"))) { continue; }
                    const [f, g] = EF[i];
                    if ((!clicked_groups.has(FG[f])) ||
                        (!clicked_groups.has(FG[g]))) { continue; }
                    const el = document.getElementById(`flat_e${i}`);
                    el.onmouseout = () => reset_edges();
                    el.onmouseover = () => {
                        if (clicked_edges.has(i)) {
                            const E = COMP.separator_branch(FOLD, clicked_edges, i);
                            for (const e of E) {
                                document.getElementById(`flat_e${e}`)
                                    .setAttribute("stroke-width", LENGTH.active);
                            }
                        } else {
                            const E = E_sep[i];
                            for (const e of E) {
                                document.getElementById(`flat_e${e}`)
                                    .setAttribute("stroke-width", LENGTH.active);
                            }
                        }
                    };
                    el.onclick = () => {
                        if (clicked_edges.has(i)) {
                            const E = COMP.separator_branch(FOLD, clicked_edges, i);
                            for (const e of E) { clicked_edges.delete(e); }
                        } else {
                            const E = E_sep[i];
                            for (const e of E) { clicked_edges.set(e, 1); }
                        }
                        MAIN.update_cp(FOLD, LINE);
                    };
                }
                reset_edges();
            }
        }
        if (RF != undefined) {
            for (let i = 0; i < RF.length; ++i) {
                const hue = (i*139) % 360; // Approx Golden Angle Method
                const g = FG[RF[i][0]];
                const color = ((FOLD.FOO == undefined) && (mode == "all")) ? (
                        (over_group == g) ? COLOR.face.select : COLOR.face.active
                    ) : `hsl(${hue}, ${(FOLD.type == TYPE.COMPLEX) ? 30 : 100
                }%, 85%)`;
                for (const f of RF[i]) {
                    const el = document.getElementById(`flat_f${f}`);
                    el.setAttribute("fill", color);
                }
            }
        }
        if (FOLD.FM != undefined) {
            const {Vf, type, EC, FR, RF, V_border, V_sink} = FOLD;
            const g = SVG.clear("notes");
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
            const colors = EA.map(a => {
                return (a == "F") ? COLOR.edge.Fcp : COLOR.edge[a];
            });
            for (let i = 0; i < EC.length; ++i) {
                const el = document.getElementById(`flat_e${i}`);
                el.setAttribute("stroke", colors[i]);
                if (!EC[i]) { continue; }
                el.setAttribute("stroke-width", LENGTH.half);
            }
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
    },
    update_state: (svg, FOLD, CELL, LINE) => {
        const {Ccolor, Ctop} = CELL;
        const {FG, clicked_groups, over_group} = LINE ?? {
            FG: Array(FOLD.Ff.length).fill(0),
            clicked_groups: new Set(),
            over_group: undefined,
        };
        for (let i = 0; i < Ctop.length; ++i) {
            const el = svg.getElementById(`fold_c${i}`);
            if (el == undefined) { continue; }
            const f = Ctop[i];
            const g = FG[f];
            el.setAttribute("fill",
                (over_group == g)     ? COLOR.face.select : (
                clicked_groups.has(g) ? COLOR.face.active
                                      : Ccolor[i])
            );
        }
    },
    draw_state: (id, FS, LINE) => {
        const svg = SVG.clear(id);
        const [FOLD, CELL] = FS[FS.length - 1];
        const {Ff, EF, FO, EA, H, LL} = FOLD;
        const {P, PP, CP, CF, SP, SC, SE} = CELL;
        const flip = document.getElementById("flip").checked;
        const Q = COMP.V_P_transform(P, P, flip);
        const slider = document.getElementById("slider");
        const CD = X.CF_edges_2_CD(CF, FO.map(([f1, f2, o]) => {
            return M.encode(((Ff[f2] ? 1 : -1)*o >= 0) ? [f1, f2] : [f2, f1]);
        })); // CF ordered in state
        const Ctop = CD.map(S => flip ? S[0] : S[S.length - 1]);
        const Ccolor = Ctop.map(d => {
            if (d == undefined) { return undefined; }
            if (Ff[d] != flip)  { return COLOR.face.top; }
            else                { return COLOR.face.bottom; }
        });
        CELL.Q = Q;
        CELL.CD = CD;
        CELL.Ctop = Ctop;
        CELL.Ccolor = Ccolor;
        if (LL != undefined) { slider.setAttribute("max", LL.length); }
        const cycle = document.getElementById("cycle");
        if (LL == undefined) {
            slider.style.display = "none";
            cycle.style.display = "inline";
        } else {    // update linearized state if exists
            slider.style.display = "inline";
            cycle.style.display = "none";
            slider.oninput = () => MAIN.draw_state(id, FS, LINE);
            const val = +slider.value;
            const flip = document.getElementById("flip").checked;
            COMP.slide_LL(FOLD, CELL, val, flip);
            for (let i = 0; i < Ccolor.length; ++i) {
                const d = Ctop[i];
                Ccolor[i] = (d == undefined) ? undefined : (
                    (Ff[d] != flip) ? COLOR.face.top : COLOR.face.bottom);
            }
        }
        const SD = X.Ctop_SC_SE_EF_Ff_2_SD(Ctop, SC, SE, EF, Ff);
        const cells = CP.map(V => M.expand(V, Q));
        const fold_c = SVG.append("g", svg, {id: "fold_c"});
        const fold_s_crease = SVG.append("g", svg, {id: "fold_s_crease"});
        const fold_s_edge = SVG.append("g", svg, {id: "fold_s_edge"});
        SVG.append("g", svg, {id: "lines"});
        SVG.append("g", svg, {id: "fold_p"});
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
            MAIN.draw_point_interface(id, FS);
        } else {
            MAIN.draw_line_interface(id, FS, LINE);
        }
    },
    draw_line_interface: (id, FS, LINE) => {
        const [FOLD, CELL] = FS[FS.length - 1];
        const {CF, Ctop} = CELL;
        const {el, FG, clicked_groups, clicked_edges} = LINE;
        const Lsvg = document.getElementById("lines");
        Lsvg.appendChild(el);
        (el.onmouseout = () => STYLE.apply(el, STYLE.line_active))();
        el.onmouseover = () => STYLE.apply(el, STYLE.line_select);
        el.onclick = () => { MAIN.update_interface(FS); }
        const input = document.getElementById("input");
        const group_interface = (el, g) => {
            el.onmouseover = () => {
                LINE.over_group = g;
                MAIN.update_cp(FOLD, LINE);
                MAIN.update_state(input, FOLD, CELL, LINE);
            };
            el.onmouseout = () => {
                LINE.over_group = undefined;
                MAIN.update_cp(FOLD, LINE);
                MAIN.update_state(input, FOLD, CELL, LINE);
            };
            el.onclick = () => {
                if (clicked_groups.has(g)) {
                    clicked_groups.delete(g);
                    clicked_edges.clear();
                } else {
                    clicked_groups.add(g);
                }
                MAIN.draw_cp(FOLD, LINE);
                MAIN.update_state(input, FOLD, CELL, LINE);
                MAIN.draw_line_interface("input", FS, LINE);
            };
        }
        for (let i = 0; i < CF.length; ++i) {
            const el = input.getElementById(`fold_c${i}`);
            const g = FG[Ctop[i]];
            group_interface(el, g);
        }
        const cp = document.getElementById("cp");
        for (let i = 0; i < FG.length; ++i) {
            const el = cp.getElementById(`flat_f${i}`);
            const g = FG[i];
            group_interface(el, g);
        }
    },
    draw_point_interface: (id, FS) => {
        const svg = document.getElementById(id);
        const [FOLD, CELL] = FS[FS.length - 1];
        const {Ff, EF, FO, EA, H, LL} = FOLD;
        const {Q, P, PP, CP, CF, SP, SC, SE, Ctop} = CELL;
        const Pvisible = (() => {
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
        })();
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
    },
    point_click: (i, lfP, FS) => {
        const [FOLD, CELL] = FS[FS.length - 1];
        const {eps} = FOLD;
        const {P} = CELL;
        if (lfP.has(i) || (lfP.size > 3)) { lfP.clear(); }
        else                              { lfP.add(i); }
        const L = COMP.get_lines(Array.from(lfP).map(i => P[i]), eps);
        for (let i = 0; i < P.length; ++i) {    // rendering
            const el = document.getElementById(`fold_p${i}`);
            if (el == undefined) { continue; }
            STYLE.apply(el, lfP.has(i)
                ? STYLE.point_active : STYLE.point_hidden);
        }
        const svg = SVG.clear("lines");
        if (L.length == 0) { return; }
        const flip = document.getElementById("flip").checked;
        SVG.draw_segments(svg,
            L.map(l => COMP.V_P_transform(COMP.line_2_coords(l), P, flip)),
            {id: true}
        );
        for (let j = 0; j < L.length; ++j) {    // interface
            const el = svg.children[j];
            (el.onmouseout = () => STYLE.apply(el, STYLE.line_normal))();
            el.onmouseover = () => STYLE.apply(el, STYLE.line_select);
            el.onclick = () => MAIN.line_click(el, Array.from(lfP), L[j], FS);
        }
    },
    make_fold: (V, FV_, FV, F_map_old, FG, FO_, clicked_groups, lfP, lfL, FS) => {
        const [FOLD_old, _] = FS.pop();
        const {FR} = FOLD_old;
        const [Vx, Vy, FVy, FM, FOO, F_map, FR_] = COMP.filter_clicked_and_reflect(
            V, FV_, FV, F_map_old, FG, FO_, clicked_groups, lfL, FR);

        // const [FOLD_, _] = COMP.V_FV_2_FOLD_CELL(Vx, FVy);
        // FOLD_.FO = FOO;
        // FOLD_.FM = FM;
        // COMP.augment_FOLD_FO(FOLD_);
        // COMP.draw_separators(FOLD_);
        // return;

        const [FOLD, CELL] = COMP.V_FV_2_FOLD_CELL(Vy, FVy);
        FOLD.FM = FM;
        FOLD.FOO = FOO;
        FOLD.lfP = lfP;
        FOLD.lfL = lfL;
        FOLD.fixed = true;
        const {EV, EF, FE, Ff} = FOLD;
        const {P, PP, CP, FC, CF, SE, SC, SP, BF, BI} = CELL;
        const type = document.getElementById("type_select").value;
        const EF_set = new Set();
        for (const F of EF) {
            if (F.length != 2) { continue; }
            EF_set.add(M.encode_order_pair(F));
        }
        const FO = [];
        const mode = document.getElementById("mode_select").value;
        for (const [f, g, o] of FO_) {
            for (const f_ of F_map[f]) {
                for (const g_ of F_map[g]) {
                    const e = M.encode_order_pair([f_, g_]);
                    if (!BI.has(e)) { continue; } // no longer exist
                    const mode = document.getElementById("mode_select").value;
                    if (mode == "all") {
                        if (!FM[f_] && !FM[g_]) { FO.push([f_, g_, o]); }
                    } else {
                        if (FR_[f_] == FR_[g_]) { FO.push([f_, g_, o]); continue; }
                        if (FR_[f_] == undefined) { continue; }
                        if (FR_[g_] == undefined) { continue; }
                        if (!EF_set.has(e)) { continue; }
                        FO.push([f_, g_, (o < 1) ? 1 : -1]);
                    }
                }
            }
        }
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
        const [GB, GA] = COMP.solve(FOLD, CELL, BA0);
        if (GB == undefined) { // failed to solve
            MAIN.update_interface(FS);
            NOTE.end();
            return;
        }
        const Gn = GA.map(A => A.length);
        NOTE.time("Solve completed");
        const n = Gn.reduce((s, gn) => s*BigInt(gn), BigInt(1));
        NOTE.count(n, "folded states");
        NOTE.lap();
        if (n == 0) {
            NOTE.end();
            MAIN.update_interface(FS);
            return;
        }
        const GI = GB.map(() => 0);
        const type_states = TYPE_LABEL.map(() => []);
        NOTE.time("Classifying states");
        NOTE.start_check("state", Number(n));
        for (let i = 0; i < n; ++i) {
            NOTE.check(i);
            const edges = X.BF_GB_GA_GI_2_edges(BF, GB, GA, GI);
            const FO = X.edges_Ff_2_FO(edges, Ff);
            const [type, RF, FR, EC, V_sink, V_border] = COMP.classify(
                V, EV, EF, FE, Ff, FM, FO, FOO);
            type_states[type].push({gi: GI.map(i => i), n: RF.length});
            for (let i = 0; i < GI.length; ++i) {
                if (GI[i] != (Gn[i] - 1)) { GI[i] += 1; break; }
                GI[i] = 0;
            }
        }
        for (const states of type_states) {
            states.sort((a, b) => (a.n - b.n));
            for (let i = 0; i < states.length; ++i) {
                states[i] = states[i].gi;
            }
        }
        NOTE.time("Classified the following states:");
        for (let t = 0; t < type_states.length; ++t) {
            const tn = type_states[t].length;
            if (tn == 0) { continue; }
            NOTE.log(`   - ${tn} ${TYPE_LABEL[t]}`);
        }
        document.getElementById("state_controls").style.display = "inline";
        const type_select = SVG.clear("type_select");
        for (let t = 0, first = true; t < TYPE_LABEL.length; ++t) {
            const tn = type_states[t].length;
            if (tn == 0) { continue; }
            const el = document.createElement("option");
            el.setAttribute("value", t);
            el.textContent = `${TYPE_LABEL[t]} (${tn})`;
            if (first) {
                el.setAttribute("selected", true);
                first = false;
            }
            type_select.appendChild(el);
        }
        const SOLUTION = {GB, GA, GI, type_states};
        FS.push([FOLD, CELL]);
        const replace = document.getElementById("replace");
        replace.style.display = "inline";
        type_select.onchange = () => {
            NOTE.start("Changing type");
            MAIN.update_type(FS, SOLUTION);
            NOTE.end();
        };
        MAIN.update_type(FS, SOLUTION);
        NOTE.lap();
        stop = Date.now();
        NOTE.end();
    },
    update_type: (FS, SOLUTION) => {
        const [FOLD, CELL] = FS[FS.length - 1];
        const {BF} = CELL;
        const {GB, GA, GI, type_states} = SOLUTION;
        const type_select = document.getElementById("type_select");
        const t = +type_select.value;
        const states = type_states[t];
        let state_idx = 0;
        const n = states.length;
        document.getElementById("state_config").style.display = "inline";
        const state_label = document.getElementById("state_label");
        const state_select = document.getElementById("state_select");
        state_label.innerHTML = `State`;
        state_select.setAttribute("min", 1);
        state_select.setAttribute("max", n);
        state_select.value = 1;
        const compute_state = () => {
            let j = +state_select.value;
            if (j < 1) { j = 1; }
            if (j > n) { j = n; }
            state_select.value = j;
            state_idx = j - 1;
            NOTE.time("Computing state");
            const edges = X.BF_GB_GA_GI_2_edges(BF, GB, GA, states[state_idx]);
            FOLD.FO = X.edges_Ff_2_FO(edges, FOLD.Ff);
            COMP.augment_FOLD_FO(FOLD);
            console.log(` ** Fold type: ${TYPE_LABEL[FOLD.type]}`);
        };
        compute_state();
        MAIN.draw_state("output", FS);
        MAIN.draw_cp(FOLD);
        state_select.onchange = () => {
            compute_state();
            MAIN.draw_state("output", FS);
            MAIN.update_cp(FOLD);
        };
    },
};
