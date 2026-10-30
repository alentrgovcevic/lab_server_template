// Shared page settings (2026-10-29): light/dark/auto theme and a basic/advanced view, remembered per
// browser. Load it with a plain <script src="/ui.js"> AFTER the page's <style>, so there is no flash.
//
// Theme: each page keeps its own palette in `@media (prefers-color-scheme: dark) { :root {…} }`.
// This script narrows that rule to `:root:not([data-theme=light])` and copies it to
// `:root[data-theme=dark]`, so one palette per page serves auto, light and dark.
// View: elements with class "adv" show only in Advanced (`[data-view=basic] .adv` is hidden).
// Controls: put <span id="uictl"></span> in the header; it gets Basic|Advanced and Auto|Light|Dark.
// pmUI.setView("advanced") switches from page code.
(function () {
  const root = document.documentElement;
  const get = k => { try { return localStorage.getItem(k); } catch (e) { return null; } };
  const set = (k, v) => { try { localStorage.setItem(k, v); } catch (e) { /* private window */ } };

  let darkVars = "";
  for (const sheet of document.styleSheets) {
    let rules;
    try { rules = sheet.cssRules; } catch (e) { continue; }
    for (const r of rules) {
      if (r.conditionText && r.conditionText.includes("prefers-color-scheme: dark")) {
        for (const inner of r.cssRules) {
          if (inner.selectorText === ":root") {
            darkVars += inner.style.cssText;
            inner.selectorText = ":root:not([data-theme=light])";
          }
        }
      }
    }
  }
  const st = document.createElement("style");
  st.textContent = `:root[data-theme=dark]{${darkVars}color-scheme:dark}
:root[data-theme=light]{color-scheme:light}
[data-view=basic] .adv{display:none!important}
.uictl{display:inline-flex;gap:8px;align-items:center}
.uigrp{display:inline-flex}
.uictl button{font:inherit;font-size:12px;min-height:30px;padding:3px 10px;border:1px solid var(--ctl-line,var(--line));margin-left:-1px;background:var(--panel);color:var(--muted);cursor:pointer}
.uigrp button:first-child{border-radius:14px 0 0 14px}.uigrp button:last-child{border-radius:0 14px 14px 0}
@media (pointer:coarse){.uictl button,header button,.tabs button{min-height:44px;padding:6px 12px}header a{display:inline-block;padding:10px 4px}}
.uictl button:hover{color:var(--ink);border-color:var(--accent)}
.uictl button[aria-pressed=true]{background:var(--accent-soft);color:var(--accent);border-color:var(--accent)}
:focus-visible{outline:2px solid var(--accent);outline-offset:2px}
@media print{.uictl{display:none!important}}`;
  document.head.appendChild(st);

  const THEMES = ["auto", "light", "dark"], ICON = {auto: "◐ Auto", light: "☀ Light", dark: "☾ Dark"};
  function applyTheme(t) { if (t === "auto") delete root.dataset.theme; else root.dataset.theme = t; }
  function applyView(v) { root.dataset.view = v; document.dispatchEvent(new CustomEvent("pm-view", {detail: v})); }
  let theme = THEMES.includes(get("pm-theme")) ? get("pm-theme") : "auto";
  // No choice saved yet: a browser that has used the site before (it has seen a guide or What's new)
  // starts in Advanced, so nobody loses a tab they used; a new browser starts in Basic (2026-10-29).
  const returning = ["pmt-home", "pmt-bank", "pm-whatsnew"].some(k => get(k) !== null);
  let view = get("pm-view") ? (get("pm-view") === "advanced" ? "advanced" : "basic") : (returning ? "advanced" : "basic");
  applyTheme(theme);
  applyView(view);
  // Save the first default at once: whatsnew.js writes its marker on the first page, which would make
  // the second page count this browser as "returning" (UX review 2026-10-29).
  if (!get("pm-view")) set("pm-view", view);

  const setView = v => { view = v; set("pm-view", v); applyView(v); paint(); };
  const setTheme = t => { theme = t; set("pm-theme", t); applyTheme(t); paint(); };
  let mounted = null;
  function paint() {
    if (!mounted) return;
    mounted.querySelectorAll("[data-view]").forEach(b => b.setAttribute("aria-pressed", b.dataset.view === view));
    mounted.querySelectorAll("[data-theme-set]").forEach(b => b.setAttribute("aria-pressed", b.dataset.themeSet === theme));
  }
  window.pmUI = {
    get view() { return view; },
    setView,
    // show Advanced for this page only (a deep link to a tool Basic hides), without changing the choice
    reveal(v) { view = v; applyView(v); paint(); },
    mount(el) {
      if (!el) return;
      mounted = el;
      el.classList.add("uictl");
      el.innerHTML = `<span class="uigrp" role="group" aria-label="View" title="Basic shows the essentials; Advanced shows every column and tool">`
        + `<button type="button" data-view="basic">Basic</button><button type="button" data-view="advanced">Advanced</button></span>`
        + `<span class="uigrp" role="group" aria-label="Colour theme">`
        + THEMES.map(t => `<button type="button" data-theme-set="${t}" aria-label="${t} theme" title="${ICON[t]}">${ICON[t].split(" ")[0]}</button>`).join("") + `</span>`;
      el.onclick = e => {
        const b = e.target.closest("button");
        if (b?.dataset.view) setView(b.dataset.view);
        if (b?.dataset.themeSet) setTheme(b.dataset.themeSet);
      };
      if (!document.querySelector(".adv") && !el.dataset.always) el.querySelector(".uigrp").hidden = true;
      paint();
    },
  };
  document.addEventListener("DOMContentLoaded", () => {
    window.pmUI.mount(document.getElementById("uictl"));
  });
})();
