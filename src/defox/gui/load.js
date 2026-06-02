
export const LOAD = {
    EL: undefined,

    tasks: 0,
    checked: 0,
    time: 0,


    startup: async () => {
        LOAD.EL = document.getElementById("loading");
        await fetch('./resources/loading.xml')
            .then(response => response.text())
            .then(xml => {
                LOAD.EL.innerHTML = xml;
                LOAD.EL.style.visibility = "visible";
            });
    },


    start: (func, tasks = 1) => {
        LOAD.set(tasks);

        const el = LOAD.EL;
        el.style.visibility = "visible";
        func();
        LOAD.check();
        LOAD.end();
        el.style.visibility = "hidden";
    },

    set: (tasks) => {
        LOAD.tasks = tasks;
        LOAD.checked = -1;
        LOAD.time = Date.now();

        LOAD.check();
    },

    end: () => {
        LOAD.tasks = 0;
        LOAD.checked = -1;
        LOAD.time = 0;

        document.getElementById("progress").innerHTML = "";
        document.getElementById("remains").innerHTML = "";
    },
    check: () => {
        LOAD.checked = LOAD.checked + 1;
        if (LOAD.checked > LOAD.tasks) {
            LOAD.end();
            return
        }
        const el = document.getElementById("progress");
        const done = "|";
        const tbd = "-";
        const len = 30;
        const d = len * (LOAD.checked) / LOAD.tasks;
        el.innerHTML = done.repeat(d) + tbd.repeat(len - d);

        const rem = LOAD.time_str((Date.now() - LOAD.time) * (LOAD.tasks - LOAD.checked) / LOAD.tasks);
        document.getElementById("remains").innerHTML = rem;
    },

    time_str: (msec) => {
        if (msec < 1000) {
            const milli = Math.ceil(msec);
            return `${milli} millisecs`;
        } else if (msec < 60000) {
            const secs = Math.ceil(msec / 1000);
            return `${secs} secs`;
        } else {
            const mins = Math.floor(msec / 60000);
            const secs = Math.ceil((msec - mins * 60000) / 1000);
            return `${mins} mins ${secs} secs`;
        }
    },
}