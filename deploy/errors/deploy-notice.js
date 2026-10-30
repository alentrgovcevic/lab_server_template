// Restart warning for every page behind Caddy (bank, dashboard, admin, genome browser).
// Before a planned restart, write {"deploy_at": <unix seconds>} to deploy-notice.json beside this file
// (a few minutes ahead) and remove it afterwards; Caddy answers 204 while there is none. The countdown
// uses the server's Date header, so a wrong clock on the user's computer does not matter.
// Run bare on a laptop (no Caddy) both URLs are 404s and nothing shows.
(() => {
  const POLL_MS = 30000, STALE_S = 900;
  const baseTitle = document.title;
  let basePad = null;
  let deployAt = 0, skew = 0, dismissedFor = null, commit = null, bar = null, tick = null;

  function el() {
    if (bar) return bar;
    bar = document.createElement("div");
    bar.setAttribute("role", "alert");
    bar.style.cssText = "position:fixed;top:0;left:0;right:0;z-index:2147483647;display:flex;gap:12px;" +
      "align-items:center;justify-content:center;padding:10px 44px 10px 16px;background:#f5b700;color:#1a1300;" +
      "font:600 15px/1.35 system-ui,-apple-system,sans-serif;box-shadow:0 2px 8px rgba(0,0,0,.3);text-align:center";
    bar.innerHTML = '<span></span><button type="button" aria-label="Hide" title="Hide (the countdown stays in the tab title)" ' +
      'style="position:absolute;right:10px;top:50%;transform:translateY(-50%);border:0;background:none;' +
      'color:inherit;font:700 20px/1 system-ui;cursor:pointer;padding:4px 8px">&times;</button>';
    bar.querySelector("button").onclick = () => { dismissedFor = commit; render(); };
    document.body.appendChild(bar);
    return bar;
  }

  function render() {
    if (!deployAt) {
      if (bar) bar.hidden = true;
      if (basePad !== null) document.body.style.paddingTop = basePad;
      document.title = baseTitle;
      clearInterval(tick); tick = null;
      return;
    }
    const left = Math.max(0, Math.round(deployAt - (Date.now() / 1000 + skew)));
    const mmss = `${Math.floor(left / 60)}:${String(left % 60).padStart(2, "0")}`;
    const b = el();
    b.querySelector("span").textContent = left > 0
      ? `\u26a0 The lab server restarts in ${mmss} for an update. Save your work now and close the browser, or unsaved changes will be lost. It is back a few minutes later.`
      : "\u26a0 The lab server is restarting for an update. Anything not saved yet may be lost; it is back in a few minutes.";
    // hidden by the user: show it again for the last minute
    b.hidden = dismissedFor === commit && left > 60;
    // push the page down so the bar never covers its header (and its save buttons)
    if (basePad === null) basePad = document.body.style.paddingTop;
    document.body.style.paddingTop = b.hidden ? basePad : b.offsetHeight + "px";
    document.title = (left > 0 ? `\u26a0 restart in ${mmss} \u00b7 ` : "\u26a0 restarting \u00b7 ") + baseTitle;
  }

  async function poll() {
    try {
      const r = await fetch("/deploy-notice.json", { cache: "no-store", credentials: "same-origin" });
      const date = Date.parse(r.headers.get("Date") || "");
      const serverNow = Number.isNaN(date) ? Date.now() / 1000 : date / 1000;
      skew = serverNow - Date.now() / 1000;
      const n = r.status === 200 ? await r.json() : null;
      // a notice long past its time is left over from a crash, not a restart about to happen
      deployAt = n && n.deploy_at > serverNow - STALE_S ? n.deploy_at : 0;
      commit = n ? n.commit : null;
    } catch { return; }   // server unreachable (mid-restart): keep what we know
    render();
    if (deployAt && !tick) tick = setInterval(render, 1000);
  }

  poll();
  setInterval(poll, POLL_MS);
  document.addEventListener("visibilitychange", () => { if (!document.hidden) poll(); });
})();
