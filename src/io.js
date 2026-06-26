import { MAIN } from "./main.js";
import { NOTE } from "./flatfolder/note.js";

import { COMP } from "./compute.js";

export const IO = {
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
        return FS;
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
        if (FOLD.EA == undefined) {
            [FOLD.H, FOLD.EA] = COMP.FO_Ff_EF_2_H_EA(FOLD.FO, FOLD.Ff, FOLD.EF);
        }
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
            edges_assignment: EA.map(a => ( // as seen from colored side
                (a == "M") ? "V" : (
                (a == "V") ? "M" : a
            ))),
        };
        const cp_data = new Blob([JSON.stringify(export_cp, undefined, 2)], {
            type: "application/json"});
        const cp_link = document.getElementById("cp_anchor");
        cp_link.setAttribute("download", `cp.fold`);
        cp_link.setAttribute("href", window.URL.createObjectURL(cp_data));
        cp_link.style.textDecoration = "none";
    },
};
