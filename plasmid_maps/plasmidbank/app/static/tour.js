/* Quick guide: a few pointer bubbles over the page, with Back / Next / Skip.
   pmTour(key, steps) shows it once per browser (localStorage "<key>" = done); pmTour(key, steps, true)
   always shows it (the "Guide" link). A step is {sel, title, text, before?}; a missing target centres the bubble.
   `before` (What's new) may be async: it runs first, e.g. to open a map or a tab the step points at. */
(function () {
  const css = `
.pmt-dim { position: fixed; inset: 0; z-index: 9000; }
.pmt-hole { position: fixed; z-index: 9001; border-radius: 8px; pointer-events: none;
  box-shadow: 0 0 0 3px var(--accent, #2b7a78), 0 0 0 9999px rgba(0,0,0,.45); transition: all .2s; }
.pmt-tip { position: fixed; z-index: 9002; max-width: 340px; background: var(--panel, #fff); color: var(--ink, #222);
  border: 1px solid var(--line, #ccc); border-radius: 10px; padding: 12px 14px; font-size: 13px; line-height: 1.45;
  box-shadow: 0 6px 24px rgba(0,0,0,.25); }
.pmt-tip .pmt-title { display: block; font-size: 14px; margin-bottom: 4px; }
.pmt-row { display: flex; gap: 6px; align-items: center; margin-top: 10px; }
.pmt-row span { color: var(--muted, #777); font-size: 12px; margin-right: auto; }
.pmt-row button { padding: 3px 10px; border: 1px solid var(--line, #ccc); border-radius: 6px; background: var(--panel, #fff);
  color: inherit; cursor: pointer; font: inherit; }
.pmt-row .pmt-next { background: var(--accent, #2b7a78); border-color: var(--accent, #2b7a78); color: #fff; }
@media print { .pmt-dim, .pmt-hole, .pmt-tip { display: none !important; } }`;
  const store = {
    get: k => { try { return localStorage.getItem(k); } catch (e) { return null; } },
    set: k => { try { localStorage.setItem(k, "done"); } catch (e) { /* private window: shows again next time */ } },
  };
  const esc = s => String(s).replace(/[&<>"]/g, c => ({"&": "&amp;", "<": "&lt;", ">": "&gt;", '"': "&quot;"}[c]));

  window.pmTour = function (key, steps, force = false) {
    if (!force && store.get(key)) return;
    if (!document.getElementById("pmt-css")) {
      const st = document.createElement("style"); st.id = "pmt-css"; st.textContent = css; document.head.appendChild(st);
    }
    document.querySelectorAll(".pmt-dim, .pmt-hole, .pmt-tip").forEach(e => e.remove());
    const dim = Object.assign(document.createElement("div"), {className: "pmt-dim"});
    const hole = Object.assign(document.createElement("div"), {className: "pmt-hole"});
    const tip = Object.assign(document.createElement("div"), {className: "pmt-tip"});
    tip.setAttribute("role", "dialog");
    document.body.append(dim, hole, tip);
    let i = 0;
    const end = () => { store.set(key); dim.remove(); hole.remove(); tip.remove(); removeEventListener("keydown", onKey); removeEventListener("resize", show); };
    const onKey = e => { if (e.key === "Escape") end(); else if (e.key === "ArrowRight") go(1); else if (e.key === "ArrowLeft") go(-1); };
    let busy = false;
    const go = d => { if (busy) return; i += d; if (i >= steps.length) end(); else if (i >= 0) show(); else i = 0; };
    async function show() {
      const s = steps[i];
      if (s.before) { busy = true; try { await s.before(); } catch (e) { /* the bubble still shows, centred */ } busy = false; }
      const el = s.sel && document.querySelector(s.sel);
      tip.innerHTML = `<b class="pmt-title">${esc(s.title)}</b>${s.text}<div class="pmt-row"><span>${i + 1} / ${steps.length}</span>
        <button class="pmt-skip">Skip</button>${i ? `<button class="pmt-back">Back</button>` : ""}
        <button class="pmt-next">${i === steps.length - 1 ? "Done" : "Next"}</button></div>`;
      tip.querySelector(".pmt-skip").onclick = end;
      tip.querySelector(".pmt-next").onclick = () => go(1);
      tip.querySelector(".pmt-back") && (tip.querySelector(".pmt-back").onclick = () => go(-1));
      const r = el && el.getClientRects().length ? (el.scrollIntoView({block: "nearest"}), el.getBoundingClientRect()) : null;
      const W = innerWidth, H = innerHeight, tw = Math.min(340, W - 32);
      tip.style.width = tw + "px";
      if (!r) {
        hole.style.cssText = "left:50%;top:50%;width:0;height:0";
        tip.style.left = (W - tw) / 2 + "px"; tip.style.top = Math.max(16, H / 2 - tip.offsetHeight / 2) + "px";
      } else {
        const pad = 6;
        Object.assign(hole.style, {left: r.left - pad + "px", top: r.top - pad + "px", width: r.width + 2 * pad + "px", height: r.height + 2 * pad + "px"});
        const below = r.bottom + 14 + tip.offsetHeight < H;
        tip.style.top = (below ? r.bottom + 14 : Math.max(16, r.top - 14 - tip.offsetHeight)) + "px";
        tip.style.left = Math.min(Math.max(16, r.left), W - tw - 16) + "px";
      }
      tip.querySelector(".pmt-next").focus();
    }
    dim.onclick = () => go(1);
    addEventListener("keydown", onKey); addEventListener("resize", show);
    show();
  };
})();
