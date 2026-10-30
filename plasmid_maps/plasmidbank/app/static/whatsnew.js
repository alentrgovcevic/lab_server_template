/* What's new: a card listing a release's changes, once per browser per release, like an app's
   release notes. Each item has a "Show me" that runs a guided tour of it (tour.js) on its page;
   an item of the other page links there with #whatsnew/<id>. The page passes run(id), which knows its
   own DOM: pmWhatsNew({page: "home"|"bank", run, force}). Someone who has never seen the page's guide
   is new, not returning: they get the guide, and this release is marked seen.
   To publish a release: add it at the top of RELEASES; the newest one is the one shown. A release with
   `from: "YYYY-MM-DD"` stays hidden until that day (browser's local date), so it can ship earlier. */
(function () {
  // A sample release: replace it with your own. Each item's id is a tour defined on its page
  // (whatsNewTour in index.html for "bank", in home.html for "home").
  const RELEASES = [{
    id: "sample-1", title: "What's new",
    intro: "Click <b>Show me</b> on anything to be walked through it. Every release is on the <a href=\"/changelog\">changelog</a>.",
    items: [
      {id: "findseq", page: "bank", icon: "🔎", title: "Find a sequence, copy with Ctrl+C",
       text: "Paste bases into <b>find a sequence</b> above the map: the match is selected, on either strand and across the origin. Ctrl/⌘+C copies any selection."},
      {id: "pointann", page: "bank", icon: "🧬", title: "Annotate a single base",
       text: "Select a stretch, even one base, and use <b>save selection as annotation</b> → <i>this map</i>, with any label and type."},
      {id: "wstree", page: "bank", icon: "📁", title: "My workspace in folders",
       text: "Drag files and folders between folders, rename, duplicate and move them, and pick or create a folder when you save."},
      {id: "multiup", page: "bank", icon: "📤", title: "Upload several maps at once",
       text: "Drop several GenBank or FASTA files on <b>＋ New map</b>, check each lab ID and name in the table, then upload them in one go."},
      {id: "ggzoom", page: "bank", icon: "🔍", title: "Golden Gate: check every overhang",
       text: "Under each product, “Check overhangs” shows every junction base by base on both strands and flags palindromes and near-identical overhangs."},
      {id: "theme", page: "home", icon: "◐", title: "Light or dark, Basic or Advanced",
       text: "Top right on every page: a light or dark theme (or follow your computer), and Basic or Advanced for more tools."},
    ],
  }];
  const KEY = "pm-whatsnew";
  const store = {
    get: k => { try { return localStorage.getItem(k); } catch (e) { return null; } },
    set: (k, v) => { try { localStorage.setItem(k, v); } catch (e) { /* private window */ } },
  };
  const css = `
.pmw-back { position: fixed; inset: 0; z-index: 8900; background: rgba(0,0,0,.45); display: flex; align-items: center; justify-content: center; padding: 16px; }
.pmw { background: var(--panel, #fff); color: var(--ink, #222); border-radius: 12px; max-width: 620px; width: 100%; max-height: 88vh;
  overflow: auto; padding: 18px 20px; box-shadow: 0 10px 40px rgba(0,0,0,.3); font-size: 14px; line-height: 1.45; }
.pmw h2 { margin: 0 0 4px; font-size: 19px; }
.pmw .pmw-intro { color: var(--muted, #666); margin: 0 0 12px; }
.pmw ul { list-style: none; padding: 0; margin: 0; }
.pmw li { display: flex; gap: 10px; align-items: flex-start; padding: 9px 0; border-top: 1px solid var(--line, #ddd); }
.pmw .pmw-ic { font-size: 20px; width: 28px; text-align: center; flex: none; }
.pmw .pmw-t { flex: 1; }
.pmw .pmw-t > b:first-child { display: block; }
.pmw .pmw-t span { color: var(--muted, #666); font-size: 13px; }
.pmw button { padding: 4px 11px; border: 1px solid var(--line, #ccc); border-radius: 6px; background: var(--panel, #fff); color: inherit; cursor: pointer; font: inherit; font-size: 13px; flex: none; }
.pmw .pmw-go { background: var(--accent, #2b7a78); border-color: var(--accent, #2b7a78); color: #fff; }
.pmw .pmw-log { margin-right: auto; align-self: center; font-size: 13px; color: var(--accent, #2b7a78); }
.pmw .pmw-row { display: flex; gap: 8px; justify-content: flex-end; margin-top: 6px; position: sticky; bottom: -18px; background: var(--panel, #fff); padding: 10px 0 12px; border-top: 1px solid var(--line, #ddd); }
@media print { .pmw-back { display: none !important; } }`;
  const today = new Date().toLocaleDateString("sv");  // YYYY-MM-DD, local
  // ?whatsnew=next previews the newest release before its `from` day, to test it
  const preview = new URLSearchParams(location.search).get("whatsnew") === "next";
  const latest = preview ? RELEASES[0] : RELEASES.find(r => !r.from || r.from <= today) || RELEASES[RELEASES.length - 1];

  window.pmWhatsNewItems = latest.items;
  // every release already out, newest first, for the changelog page (/changelog)
  window.pmReleases = RELEASES.filter(r => preview || !r.from || r.from <= today);
  window.pmWhatsNew = function ({page, run, force = false, guideKey}) {
    if (!force) {
      if (store.get(KEY) === latest.id) return false;
      if (guideKey && !store.get(guideKey)) { store.set(KEY, latest.id); return false; }  // new here: the guide covers it
    }
    if (!document.getElementById("pmw-css")) {
      const st = document.createElement("style"); st.id = "pmw-css"; st.textContent = css; document.head.appendChild(st);
    }
    document.querySelectorAll(".pmw-back").forEach(e => e.remove());
    const back = Object.assign(document.createElement("div"), {className: "pmw-back"});
    const here = latest.items.filter(x => x.page === page), there = latest.items.filter(x => x.page !== page);
    const other = page === "home" ? "/bank" : "/";
    back.innerHTML = `<div class="pmw" role="dialog" aria-label="What's new"><h2>${latest.title}</h2><p class="pmw-intro">${latest.intro}</p>
      <ul>${latest.items.map(x => `<li><span class="pmw-ic">${x.icon}</span><span class="pmw-t"><b>${x.title}</b><span>${x.text}</span></span>
        <button data-w="${x.id}">Show me</button></li>`).join("")}</ul>
      <div class="pmw-row"><a class="pmw-log" href="/changelog">All changes →</a><button class="pmw-close">Close</button>${here.length ? `<button class="pmw-go pmw-all">Take the tour (${here.length} on this page)</button>` : ""}</div></div>`;
    document.body.appendChild(back);
    store.set(KEY, latest.id);
    const close = () => { back.remove(); removeEventListener("keydown", onKey); };
    const onKey = e => { if (e.key === "Escape") close(); };
    addEventListener("keydown", onKey);
    back.onclick = e => { if (e.target === back) close(); };
    back.querySelector(".pmw-close").onclick = close;
    back.querySelector(".pmw-all") && (back.querySelector(".pmw-all").onclick = () => { close(); run(here.map(x => x.id), there.length ? other + "#whatsnew" : null); });
    back.querySelectorAll("[data-w]").forEach(b => b.onclick = () => {
      const it = latest.items.find(x => x.id === b.dataset.w);
      close();
      if (it.page === page) run([it.id]);
      else location.href = `${it.page === "home" ? "/" : "/bank"}#whatsnew/${it.id}`;
    });
    return true;
  };
})();
