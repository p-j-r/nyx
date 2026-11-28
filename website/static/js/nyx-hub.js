(() => {
  "use strict";

  const SLIDE_DURATION_MS = 6000;
  const reducedMotionQuery = window.matchMedia
    ? window.matchMedia("(prefers-reduced-motion: reduce)")
    : { matches: false };

  // ---------------------------------------------------------------------------
  // Theme Toggle
  // ---------------------------------------------------------------------------
  function initThemeToggle() {
    const btn = document.querySelector("[data-theme-toggle]");
    if (!btn) return;

    const syncBtn = (theme) => {
      const isDark = theme === "dark";
      btn.setAttribute("aria-pressed", String(isDark));
      const label = btn.querySelector("[data-theme-label]");
      if (label) {
        label.textContent = isDark ? "Dark" : "Light";
      }
      btn.title = isDark ? "Switch to light monograph theme" : "Switch to dark monograph theme";
    };

    const current = document.documentElement.getAttribute("data-theme") || "light";
    syncBtn(current);

    btn.addEventListener("click", () => {
      const prev = document.documentElement.getAttribute("data-theme") || "light";
      const next = prev === "dark" ? "light" : "dark";
      document.documentElement.setAttribute("data-theme", next);
      try {
        localStorage.setItem("nyx-theme", next);
      } catch (_) {
        // Ignore storage restrictions
      }
      syncBtn(next);
    });
  }

  // ---------------------------------------------------------------------------
  // Mobile Header Navigation
  // ---------------------------------------------------------------------------
  function initMobileNav() {
    const toggle = document.querySelector("[data-nav-toggle]");
    const menu = document.querySelector("[data-nav-menu]");
    if (!toggle || !menu) return;

    toggle.addEventListener("click", () => {
      const open = menu.classList.toggle("is-open");
      toggle.setAttribute("aria-expanded", String(open));
    });
  }

  // ---------------------------------------------------------------------------
  // Accessible Mobile Sidenotes
  // ---------------------------------------------------------------------------
  function initSidenotes() {
    document.querySelectorAll(".sn-wrap").forEach((wrap) => {
      const checkbox = wrap.querySelector(".sn-toggle");
      const label = wrap.querySelector(".sn-number");
      if (!checkbox || !label) return;

      const sync = () => {
        label.setAttribute("aria-expanded", String(checkbox.checked));
      };

      checkbox.addEventListener("change", sync);
      label.addEventListener("keydown", (ev) => {
        if (ev.key === "Enter" || ev.key === " ") {
          ev.preventDefault();
          checkbox.checked = !checkbox.checked;
          sync();
        }
      });
      sync();
    });
  }

  // ---------------------------------------------------------------------------
  // Copy-to-Clipboard
  // ---------------------------------------------------------------------------
  async function copyText(text) {
    if (navigator.clipboard && navigator.clipboard.writeText) {
      await navigator.clipboard.writeText(text);
      return;
    }
    const ta = document.createElement("textarea");
    ta.value = text;
    ta.setAttribute("readonly", "");
    ta.style.position = "fixed";
    ta.style.opacity = "0";
    document.body.appendChild(ta);
    ta.select();
    document.execCommand("copy");
    document.body.removeChild(ta);
  }

  function initCopyButtons() {
    document.addEventListener("click", async (ev) => {
      const btn = ev.target.closest("[data-copy-btn]");
      if (!btn) return;

      let text = btn.getAttribute("data-copy-text") || "";
      if (!text) {
        const targetSel = btn.getAttribute("data-copy-target");
        const target = targetSel ? document.querySelector(targetSel) : null;
        text = target ? target.textContent.trim() : "";
      }
      if (!text) return;

      try {
        await copyText(text);
        const orig = btn.getAttribute("data-orig-label") || btn.textContent;
        btn.setAttribute("data-orig-label", orig);
        btn.textContent = "Copied ✓";
        btn.classList.add("is-copied");
        window.clearTimeout(btn._copyTimer);
        btn._copyTimer = window.setTimeout(() => {
          btn.textContent = orig;
          btn.classList.remove("is-copied");
        }, 1800);
      } catch (_) {
        // Ignore clipboard errors in restricted environments
      }
    });
  }

  // ---------------------------------------------------------------------------
  // <nyx-slideshow> Custom Element (Hero PLATE 01)
  // ---------------------------------------------------------------------------
  class NyxSlideshow extends HTMLElement {
    connectedCallback() {
      if (this._initialized) return;
      this._initialized = true;

      this.slides = Array.from(this.querySelectorAll("[data-slide]"));
      this.segments = Array.from(this.querySelectorAll("[data-slide-dot]"));
      this.stage = this.querySelector("[data-slide-stage]");
      this.rootBox = this.querySelector(".plate-slideshow");
      this.prevBtn = this.querySelector("[data-slide-prev]");
      this.nextBtn = this.querySelector("[data-slide-next]");
      this.toggleBtn = this.querySelector("[data-slide-toggle]");
      this.counterEl = this.querySelector("[data-slide-counter]");

      this.hudPlateId = this.querySelector("[data-hud-plate-id]");
      this.hudBadge = this.querySelector("[data-hud-badge]");
      this.captionTitle = this.querySelector("[data-caption-title]");
      this.captionText = this.querySelector("[data-caption-text]");

      this.logPlateId = this.querySelector("[data-log-plate-id]");
      this.logCmd = this.querySelector("[data-log-cmd]");
      this.logCopy = this.querySelector("[data-log-copy]");
      this.logMetrics = this.querySelector("[data-log-metrics]");
      this.logStdout = this.querySelector("[data-log-stdout]");

      if (this.slides.length === 0) return;

      this.index = 0;
      this.manualPaused = Boolean(reducedMotionQuery.matches);
      this.hoverPaused = false;
      this.focusPaused = false;
      this.hiddenPaused = document.visibilityState === "hidden";
      this._timer = null;

      this.bindEvents();
      this.render(0, true);
      this.syncTimer();
    }

    disconnectedCallback() {
      this.clearTimer();
      if (this._onVisChange) {
        document.removeEventListener("visibilitychange", this._onVisChange);
      }
    }

    isPaused() {
      return this.manualPaused || this.hoverPaused || this.focusPaused || this.hiddenPaused;
    }

    clearTimer() {
      if (this._timer !== null) {
        window.clearTimeout(this._timer);
        this._timer = null;
      }
    }

    syncTimer() {
      this.clearTimer();
      const paused = this.isPaused();
      if (this.rootBox) {
        this.rootBox.classList.toggle("is-paused", paused);
      }
      if (this.toggleBtn) {
        this.toggleBtn.setAttribute("aria-pressed", String(this.manualPaused));
        this.toggleBtn.textContent = this.manualPaused ? "Paused" : "Auto · 6s";
      }
      if (!paused && this.slides.length > 1) {
        this._timer = window.setTimeout(() => this.next(), SLIDE_DURATION_MS);
      }
    }

    goTo(nextIdx) {
      const n = this.slides.length;
      const wrapped = ((nextIdx % n) + n) % n;
      this.index = wrapped;
      this.render(wrapped, true);
      this.syncTimer();
    }

    next() {
      this.goTo(this.index + 1);
    }

    prev() {
      this.goTo(this.index - 1);
    }

    render(idx, restartAnim) {
      const total = this.slides.length;
      const active = this.slides[idx];

      this.slides.forEach((slide, i) => {
        if (i === idx) {
          slide.removeAttribute("hidden");
        } else {
          slide.setAttribute("hidden", "");
        }
      });

      this.segments.forEach((seg, i) => {
        seg.classList.remove("is-active", "is-done");
        seg.setAttribute("aria-current", i === idx ? "true" : "false");
        if (i < idx) {
          seg.classList.add("is-done");
        } else if (i === idx) {
          if (restartAnim) {
            void seg.offsetWidth;
          }
          seg.classList.add("is-active");
        }
      });

      if (this.counterEl) {
        const cur = String(idx + 1).padStart(2, "0");
        const max = String(total).padStart(2, "0");
        this.counterEl.textContent = `${cur} / ${max}`;
      }

      const jsonScript = active.querySelector("script[data-slide-json]");
      let payload = {};
      if (jsonScript) {
        try {
          payload = JSON.parse(jsonScript.textContent);
        } catch (_) {
          payload = {};
        }
      }

      const plateId = payload.plate_id || active.getAttribute("data-plate-id") || "";
      const badge = payload.badge || active.getAttribute("data-badge") || "";
      const title = payload.title || active.getAttribute("data-title") || "";
      const caption = payload.caption || active.getAttribute("data-caption") || "";
      const command = payload.command || active.getAttribute("data-command") || "";
      const metrics = payload.metrics || active.getAttribute("data-metrics") || "";
      const lines = Array.isArray(payload.solver_log) ? payload.solver_log : [];

      if (this.hudPlateId) this.hudPlateId.textContent = plateId;
      if (this.hudBadge) this.hudBadge.textContent = badge;
      if (this.captionTitle) this.captionTitle.textContent = title;
      if (this.captionText) this.captionText.textContent = caption;

      if (this.logPlateId) this.logPlateId.textContent = plateId;
      if (this.logCmd) this.logCmd.textContent = command;
      if (this.logCopy) this.logCopy.setAttribute("data-copy-text", command);
      if (this.logMetrics) {
        this.logMetrics.textContent = metrics;
        this.logMetrics.hidden = !metrics;
      }
      if (this.logStdout) {
        this.logStdout.textContent = lines.join("\n");
      }
    }

    bindEvents() {
      if (this.stage) {
        this.stage.addEventListener("click", () => this.next());
        this.stage.addEventListener("keydown", (ev) => {
          if (ev.key === "Enter" || ev.key === " ") {
            ev.preventDefault();
            this.next();
          }
        });
      }

      if (this.prevBtn) {
        this.prevBtn.addEventListener("click", () => this.prev());
      }
      if (this.nextBtn) {
        this.nextBtn.addEventListener("click", () => this.next());
      }
      if (this.toggleBtn) {
        this.toggleBtn.addEventListener("click", () => {
          this.manualPaused = !this.manualPaused;
          this.syncTimer();
        });
      }

      this.segments.forEach((seg, i) => {
        seg.addEventListener("click", () => this.goTo(i));
      });

      this.addEventListener("keydown", (ev) => {
        if (ev.key === "ArrowLeft") {
          ev.preventDefault();
          this.prev();
        } else if (ev.key === "ArrowRight") {
          ev.preventDefault();
          this.next();
        }
      });

      this.addEventListener("mouseenter", () => {
        this.hoverPaused = true;
        this.syncTimer();
      });
      this.addEventListener("mouseleave", () => {
        this.hoverPaused = false;
        this.syncTimer();
      });
      this.addEventListener("focusin", () => {
        this.focusPaused = true;
        this.syncTimer();
      });
      this.addEventListener("focusout", (ev) => {
        if (!this.contains(ev.relatedTarget)) {
          this.focusPaused = false;
          this.syncTimer();
        }
      });

      this._onVisChange = () => {
        this.hiddenPaused = document.visibilityState === "hidden";
        this.syncTimer();
      };
      document.addEventListener("visibilitychange", this._onVisChange);
    }
  }

  if (window.customElements && !window.customElements.get("nyx-slideshow")) {
    window.customElements.define("nyx-slideshow", NyxSlideshow);
  }

  // ---------------------------------------------------------------------------
  // Live Filters & Search (Specimens, Materials, Benchmarks, API)
  // ---------------------------------------------------------------------------
  function initFilters() {
    const roots = document.querySelectorAll("[data-filter-root]");
    roots.forEach((root) => {
      const items = Array.from(root.querySelectorAll("[data-filter-item]"));
      const buttons = Array.from(root.querySelectorAll("[data-filter-btn]"));
      const searchInput = root.querySelector("[data-filter-search]");
      const countEl = root.querySelector("[data-filter-count]");
      const emptyEl = root.querySelector("[data-filter-empty]");

      const activeGroups = {};
      buttons.forEach((btn) => {
        const group = btn.getAttribute("data-filter-group") || "category";
        if (btn.getAttribute("aria-pressed") === "true") {
          activeGroups[group] = btn.getAttribute("data-filter-value") || "all";
        } else if (!(group in activeGroups)) {
          activeGroups[group] = "all";
        }
      });

      const apply = () => {
        const query = searchInput ? searchInput.value.trim().toLowerCase() : "";
        let visible = 0;

        items.forEach((item) => {
          let match = true;

          for (const [group, selected] of Object.entries(activeGroups)) {
            if (!selected || selected === "all") continue;
            const rawAttr = (item.getAttribute(`data-${group}`) || "").toLowerCase();
            const tokens = rawAttr.split(/[\s,|]+/).filter(Boolean);
            if (!tokens.includes(selected.toLowerCase()) && rawAttr !== selected.toLowerCase()) {
              match = false;
              break;
            }
          }

          if (match && query) {
            const haystack = (item.getAttribute("data-search") || item.textContent || "").toLowerCase();
            if (!haystack.includes(query)) {
              match = false;
            }
          }

          item.hidden = !match;
          if (match) visible += 1;
        });

        if (countEl) {
          countEl.textContent = `${visible} / ${items.length}`;
        }
        if (emptyEl) {
          emptyEl.hidden = visible > 0;
        }
      };

      buttons.forEach((btn) => {
        btn.addEventListener("click", () => {
          const group = btn.getAttribute("data-filter-group") || "category";
          const val = btn.getAttribute("data-filter-value") || "all";
          activeGroups[group] = val;

          buttons
            .filter((b) => (b.getAttribute("data-filter-group") || "category") === group)
            .forEach((b) => {
              b.setAttribute("aria-pressed", String(b === btn));
            });

          apply();
        });
      });

      if (searchInput) {
        searchInput.addEventListener("input", apply);
      }

      apply();
    });

    // Global "/" search shortcut
    document.addEventListener("keydown", (ev) => {
      if (ev.key !== "/" || ev.metaKey || ev.ctrlKey || ev.altKey) return;
      const tag = (ev.target && ev.target.tagName || "").toLowerCase();
      if (tag === "input" || tag === "textarea" || tag === "select") return;
      const searchInput = document.querySelector("[data-filter-search]");
      if (searchInput) {
        ev.preventDefault();
        searchInput.focus();
      }
    });
  }

  // ---------------------------------------------------------------------------
  // Expandable Details / Accordion Rows
  // ---------------------------------------------------------------------------
  function initExpanders() {
    document.addEventListener("click", (ev) => {
      const btn = ev.target.closest("[data-expand-btn]");
      if (!btn) return;
      const targetId = btn.getAttribute("aria-controls") || btn.getAttribute("data-expand-target");
      if (!targetId) return;
      const panel = document.getElementById(targetId);
      if (!panel) return;
      const expanded = btn.getAttribute("aria-expanded") === "true";
      btn.setAttribute("aria-expanded", String(!expanded));
      panel.hidden = expanded;
    });
  }

  // ---------------------------------------------------------------------------
  // Platform Tabs (Install Page)
  // ---------------------------------------------------------------------------
  function initTabs() {
    document.querySelectorAll("[data-tabs]").forEach((root) => {
      const tabs = Array.from(root.querySelectorAll("[role='tab']"));
      const panels = Array.from(root.querySelectorAll("[role='tabpanel']"));
      if (tabs.length === 0) return;

      const activate = (idx) => {
        tabs.forEach((tab, i) => {
          const active = i === idx;
          tab.setAttribute("aria-selected", String(active));
          tab.setAttribute("tabindex", active ? "0" : "-1");
          if (panels[i]) {
            panels[i].hidden = !active;
          }
        });
      };

      tabs.forEach((tab, i) => {
        tab.addEventListener("click", () => activate(i));
        tab.addEventListener("keydown", (ev) => {
          if (ev.key === "ArrowRight") {
            ev.preventDefault();
            const next = (i + 1) % tabs.length;
            activate(next);
            tabs[next].focus();
          } else if (ev.key === "ArrowLeft") {
            ev.preventDefault();
            const prev = (i - 1 + tabs.length) % tabs.length;
            activate(prev);
            tabs[prev].focus();
          }
        });
      });
    });
  }

  // ---------------------------------------------------------------------------
  // Boot
  // ---------------------------------------------------------------------------
  document.addEventListener("DOMContentLoaded", () => {
    initThemeToggle();
    initMobileNav();
    initSidenotes();
    initCopyButtons();
    initFilters();
    initExpanders();
    initTabs();
  });
})();
