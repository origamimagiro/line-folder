import { CON } from "../../flatfolder/constraints.js";

import { SMPL } from "../sample.js";
import { PRJ } from "../project.js";
import { GUI_STATE } from "./state.js";
import { GUI_PAGE } from "./page.js";
import { GUI_IO } from "./io.js";
import { LOAD } from "./load.js";

export const GUI = {
    samples: [
        SMPL.hf,
        SMPL.sq,
        SMPL.windmil,
        SMPL.hf,
        SMPL.nonlin,
        SMPL.sq,
        SMPL.hanikamu],

    startup: async () => {
        const defs = document.getElementById("defs");
        await GUI.fetch(defs, './resources/defs.xml');
        CON.build();
        GUI_IO.startup();
        GUI_STATE.startup();
        await LOAD.set(
            GUI.samples.length,
            async () => {
                await GUI.build();
            });
        GUI_PAGE.startup();
    },

    build: async () => {
        for (const [idx, sample] of GUI.samples.entries()) {
            GUI_IO.import_cp("sample", sample, idx == 0);
            await LOAD.report();
        }
    },

    open_close: (id, display_style) => {
        var el = document.getElementById(id);
        document.getElementById(id + "_b").onclick = () => {
            if (el.style.display == display_style) {
                el.style.display = "none";
            }
            else {
                el.style.display = display_style;
            }
        }
    },

    fetch: async (par, url) => {
        try {
            const response = await fetch(url);
            const xml = await response.text();
            par.innerHTML = xml;
        }
        catch (error) {
            console.error("fetch error:", error)
        }
    },


}