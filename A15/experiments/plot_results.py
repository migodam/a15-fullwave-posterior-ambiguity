#!/usr/bin/env python
"""A15 figure generation (bounded, read-only wrt source, configs, protocol, results).

This script ONLY writes into:
    A15/figures/*
    A15/research/delegated/figures/*
It never edits src/, configs/, EXPERIMENT_PROTOCOL.md or results/.

Figure status vocabulary:
    GENERATED    - drawn from a frozen/model input that does not depend on the
                   pending confirmatory posterior reference.
    PROVISIONAL  - drawn, but its inputs are discovery screening data, or the
                   confirmatory reference gate (results/reference_validation_gate.json)
                   is absent / not passed. NOT a confirmatory result.
    NOT_RUN      - a required input file is missing; a clearly labelled
                   placeholder panel is written instead. No data is invented.
    ERROR        - inputs existed but the figure could not be drawn; the error
                   text is recorded and a labelled placeholder is written.

No scientific verdict is emitted here. Every figure carries a status banner.

Run (from A15/):
    ../../Gaussian/.venv_nn/bin/python experiments/plot_results.py
"""
from __future__ import annotations

import json
import sys
import traceback
from datetime import datetime, timezone
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

plt.rcParams["axes.unicode_minus"] = False
plt.rcParams["figure.dpi"] = 120
plt.rcParams["savefig.dpi"] = 140
plt.rcParams["font.size"] = 9

# --------------------------------------------------------------------------
# Paths
# --------------------------------------------------------------------------
HERE = Path(__file__).resolve().parent
A15 = HERE.parent
SRC = A15 / "src"
FIG_DIR = A15 / "figures"
DELEG_DIR = A15 / "research" / "delegated" / "figures"
RESULTS = A15 / "results"

PROTOCOL_JSON = A15 / "configs" / "protocol.json"
DISCOVERY_JSON = A15 / "configs" / "discovery.json"

# Discovery-stage inputs available now.
DISC_SWEEP = RESULTS / "discovery" / "r0.10_e12.npz"          # k-sweep on 24 thetas
DISC_SWEEP_ALL = sorted((RESULTS / "discovery").glob("r0.1*_e*.npz"))
DISC_COARSE = RESULTS / "discovery" / "coarse_r0.10_e12.npz"   # 9^3 screening table
DISC_COARSE_ALL = sorted((RESULTS / "discovery").glob("coarse_r0.1*_e*.npz"))
DISC_SUMMARY = RESULTS / "discovery" / "summary.json"

# Confirmatory inputs (expected to arrive later).
GATE_JSON = RESULTS / "reference_validation_gate.json"
TABLE_DIR = RESULTS / "table_41_coupled"
TABLE_PRIOR = TABLE_DIR / "prior.npz"
TABLE_Q = TABLE_DIR / "Q.npy"
TABLE_Y = TABLE_DIR / "y.npy"
TABLE_STATUS = TABLE_DIR / "status.json"
VAL_DIR = RESULTS / "validation_41_coupled"
VAL_SUMMARY = VAL_DIR / "summary.json"
QUERIES_NPZ = RESULTS / "queries_validation_11_coupled.npz"
UNCOUPLED_SUMMARY = RESULTS / "validation_33_uncoupled" / "summary.json"
MODAL_NPZ = RESULTS / "modal_audit.npz"
MODAL_META = RESULTS / "modal_audit.json"

for d in (FIG_DIR, DELEG_DIR):
    d.mkdir(parents=True, exist_ok=True)

MANIFEST: list[dict] = []


# --------------------------------------------------------------------------
# small helpers
# --------------------------------------------------------------------------
def rel(p: Path) -> str:
    try:
        return str(Path(p).resolve().relative_to(A15.resolve()))
    except Exception:
        return str(p)


def load_json(p: Path):
    return json.loads(Path(p).read_text())


def gate_passed() -> bool:
    """True only if the confirmatory reference gate file exists and passed."""
    if not GATE_JSON.exists():
        return False
    try:
        g = load_json(GATE_JSON)
    except Exception:
        return False
    if isinstance(g, dict):
        if g.get("passed") is True:
            return True
        st = str(g.get("status", "")).upper()
        return st in {"PASS", "PASSED", "OK"}
    return False


def banner(ax, figid: str, title: str, status: str, note: str = ""):
    # Figure-level footer records evidence status without colliding with titles.
    return
    """Status banner drawn in the figure corner (small, non-intrusive)."""
    color = {"GENERATED": "#1a7f37",
             "PROVISIONAL": "#b8860b",
             "NOT_RUN": "#8b1a1a",
             "ERROR": "#8b1a1a"}.get(status, "#444444")
    ax.text(0.0, 1.0, f"{figid}  [{status}]", transform=ax.transAxes,
            va="bottom", ha="left", fontsize=9, color=color, weight="bold")


def placeholder_fig(figid: str, title: str, status: str, missing, note=""):
    """Full-page placeholder with an unmistakable status label."""
    fig = plt.figure(figsize=(8.0, 4.5))
    ax = fig.add_axes([0.05, 0.05, 0.9, 0.9])
    ax.axis("off")
    color = "#8b1a1a" if status in ("NOT_RUN", "ERROR") else "#b8860b"
    ax.add_patch(plt.Rectangle((0.02, 0.02), 0.96, 0.96, fill=False,
                               ec=color, lw=3, transform=ax.transAxes))
    ax.text(0.5, 0.80, status, transform=ax.transAxes, ha="center",
            va="center", fontsize=30, color=color, weight="bold")
    ax.text(0.5, 0.63, f"{figid}: {title}", transform=ax.transAxes,
            ha="center", va="center", fontsize=13)
    if missing:
        ax.text(0.5, 0.44, "missing input(s):", transform=ax.transAxes,
                ha="center", va="center", fontsize=10, weight="bold")
        ax.text(0.5, 0.30, "\n".join(rel(m) for m in missing),
                transform=ax.transAxes, ha="center", va="center",
                fontsize=9, family="monospace")
    if note:
        ax.text(0.5, 0.10, note, transform=ax.transAxes, ha="center",
                va="center", fontsize=9, color="#333333", wrap=True)
    return fig


def save(fig, name: str, status: str, title: str, sources, missing, note=""):
    """Save into figures/ and copy into research/delegated/figures/."""
    out = FIG_DIR / name
    if status == "PROVISIONAL":
        fig.text(.5, -.015, "PROVISIONAL — numerical/posterior evidence not certified", ha="center", color="#9b2c2c", fontsize=10, weight="bold")
    fig.savefig(out, bbox_inches="tight")
    copy = DELEG_DIR / name
    fig.savefig(copy, bbox_inches="tight")
    plt.close(fig)
    MANIFEST.append({
        "id": name.split("_")[0],
        "title": title,
        "status": status,
        "figure": rel(out),
        "copy": rel(copy),
        "sources": [rel(s) for s in sources],
        "missing": [rel(m) for m in missing],
        "note": note,
    })
    print(f"[{status:11s}] {rel(out)}")


def make_figure(figid, name, title, needs, draw, confirmatory=False,
                figsize=(9.5, 6.0)):
    """Wrap drawing with input checks, status rules and placeholder fallbacks."""
    missing = [p for p in needs if not Path(p).exists()]
    if missing:
        fig = placeholder_fig(figid, title, "NOT_RUN", missing,
                              note="figure will be produced when inputs exist")
        save(fig, name, "NOT_RUN", title, [], missing,
             note="missing required input file(s)")
        return
    fig = plt.figure(figsize=figsize)
    try:
        status, note, sources = draw(fig)
    except Exception as exc:  # noqa: BLE001 - we want every failure recorded
        plt.close(fig)
        tb = traceback.format_exc().strip().splitlines()[-1]
        fig = placeholder_fig(figid, title, "ERROR", needs,
                              note=f"{exc.__class__.__name__}: {exc}")
        save(fig, name, "ERROR", title, needs, [], note=tb)
        return
    if confirmatory and status == "GENERATED" and not gate_passed():
        status = "PROVISIONAL"
        note = (note + " | reference_validation_gate.json absent or not passed")
    save(fig, name, status, title, sources, [], note=note)


# --------------------------------------------------------------------------
# shared data helpers
# --------------------------------------------------------------------------
def protocol() -> dict:
    return load_json(PROTOCOL_JSON)


def _unique_axis(values: np.ndarray):
    u = np.unique(np.round(values, 12))
    idx = np.searchsorted(u, np.round(values, 12))
    return u, idx


def grid_from_theta(theta: np.ndarray, q: np.ndarray):
    """Reshape a structured (s,Delta,gamma) table onto a 3-D array."""
    us, is_ = _unique_axis(theta[:, 0])
    ud, id_ = _unique_axis(theta[:, 1])
    ug, ig = _unique_axis(theta[:, 2])
    A = np.full((len(us), len(ud), len(ug)), np.nan)
    A[is_, id_, ig] = q
    return us, ud, ug, A


def stable_logsumexp(x):
    m = np.max(x)
    if not np.isfinite(m):
        return -np.inf
    return m + np.log(np.sum(np.exp(x - m)))


def weighted_quantiles(values, weights, qs):
    order = np.argsort(values)
    v = values[order]
    w = weights[order]
    cw = np.cumsum(w) - 0.5 * w
    cw /= np.sum(w)
    return np.interp(qs, cw, v)


# --------------------------------------------------------------------------
# Figure 1 - geometry / ROI / acquisition
# --------------------------------------------------------------------------
def fig01(fig):
    p = protocol()
    fam = p["family"]
    r = fam["radius"]
    tx = fam["target_x"]
    br_gamma = fam["bright_gamma"]
    ext = np.asarray(p["external_receivers"], float)
    near = np.asarray(p["near_receivers"], float)

    ax = fig.add_subplot(111)
    ax.set_aspect("equal")

    def draw_scene(a, scale):
        """Draw geometry on axes `a`; `scale` sets marker/text sizes."""
        for gap, style in ((fam["gap_min"], ":"), (fam["gap_max"], "--")):
            cx = tx - 2 * r - gap
            a.add_patch(plt.Circle((cx, 0), r, fill=False, ec="#1f4e79",
                                   ls=style, lw=1.6))
            a.plot([cx], [0], marker="+", color="#1f4e79", ms=8 * scale)
        a.add_patch(plt.Circle((tx, 0), r, fc="#f2c14e", ec="#7a5b00", alpha=0.45))
        a.add_patch(plt.Circle((tx, 0), r, fill=False, ec="#7a5b00", lw=1.8))
        a.plot([tx], [0], marker="+", color="#7a5b00", ms=8 * scale)
        a.plot(ext[:, 0], ext[:, 1], ls="none", marker="o", ms=6 * scale,
               mfc="white", mec="#444444")
        a.plot(near[:, 0], near[:, 1], ls="none", marker="s", ms=7 * scale,
               mfc="#c0504d", mec="#7a2b28")
        roles = [(p["primary_channels"], "#2e7d32", "o"),
                 (p["extra_fixed"], "#1565c0", "D"),
                 (p["extra_optimized"], "#8e24aa", "s"),
                 (p["extra_random"], "#ef6c00", "^")]
        for chans, col, mk in roles:
            for c in chans:
                c = int(c)
                if c < 16:
                    a.plot(ext[c, 0], ext[c, 1], marker=mk, ms=9 * scale,
                           mfc="none", mec=col, mew=2.0)
        for c in p["near_channels"]:
            if int(c) >= 16:
                a.plot(near[int(c) - 16, 0], near[int(c) - 16, 1], marker="s",
                       ms=10 * scale, mfc="none", mec="#8e24aa", mew=2.0)

    draw_scene(ax, 1.0)
    # receiver labels
    ax.plot(ext[:, 0], ext[:, 1], ls="none", marker="o", ms=6,
            mfc="white", mec="#444444", label="external channels (0..15)")
    for i, (x, y) in enumerate(ext):
        ax.annotate(str(i), (x, y), fontsize=6, ha="center", va="center",
                    color="#444444")
    ax.plot(near[:, 0], near[:, 1], ls="none", marker="s", ms=7, mfc="#c0504d",
            mec="#7a2b28", label="near-field channels (16,17)")
    for i, (x, y) in enumerate(near):
        ax.annotate(str(16 + i), (x, y), fontsize=6, ha="center", va="center",
                    color="white")

    # incident plane wave +x
    ax.annotate("", xy=(1.9, 0.9), xytext=(-2.0, 0.9),
                arrowprops=dict(arrowstyle="-|>", color="#666666"))
    ax.text(-1.9, 0.97, "incident E_z plane wave, +x", fontsize=8, color="#666666")
    ax.set_xlim(-2.35, 2.45)
    ax.set_ylim(-0.75, 2.15)
    ax.set_xlabel("x")
    ax.set_ylabel("y")

    # zoom inset on the resonators
    axi = ax.inset_axes([0.34, 0.48, 0.30, 0.37])
    draw_scene(axi, 0.6)
    axi.set_xlim(-0.30, 0.36)
    axi.set_ylim(-0.22, 0.22)
    axi.set_title("zoom: bright body + target ROI", fontsize=7)
    axi.tick_params(labelsize=6)
    axi.text(tx, -r - 0.035, "target = ROI (Q region)", ha="center", fontsize=6,
             color="#7a5b00")
    axi.text(tx - 2 * r - fam["gap_max"], r + 0.02,
             "bright body, gap s in [%.2f, %.2f]" % (fam["gap_min"], fam["gap_max"]),
             ha="center", fontsize=6, color="#1f4e79")

    handles = [plt.Line2D([], [], ls="none", marker="o", mfc="white", mec="#444444",
                          label="external channel"),
               plt.Line2D([], [], ls="none", marker="s", mfc="#c0504d", mec="#7a2b28",
                          label="near-field channel"),
               plt.Line2D([], [], ls="none", marker="o", mfc="none", mec="#2e7d32",
                          mew=2, label="base channels"),
               plt.Line2D([], [], ls="none", marker="D", mfc="none", mec="#1565c0",
                          mew=2, label="fixed +2"),
               plt.Line2D([], [], ls="none", marker="s", mfc="none", mec="#8e24aa",
                          mew=2, label="optimized +2"),
               plt.Line2D([], [], ls="none", marker="^", mfc="none", mec="#ef6c00",
                          mew=2, label="random +2")]
    ax.legend(handles=handles, loc="upper left", fontsize=7, framealpha=0.9)
    ax.set_title("Fig 1  Canonical two-resonator geometry, ROI and acquisition\n"
                 "radius=%.3g, epsilon_base=%.3g, bright gamma=%.3g, k0=%.4g"
                 % (r, fam["epsilon_base"], br_gamma, p["k0"]))
    note = ("Model definition only (protocol.json). Q = ln(mean ROI |E|^2 / incident). "
            "Receiver roles frozen before validation.")
    ax.text(0.01, 0.01, note, transform=ax.transAxes, fontsize=7,
            color="#333333")
    banner(ax, "FIG01", "geometry", "GENERATED")
    return "GENERATED", "geometry/ROI/acquisition from protocol.json", [PROTOCOL_JSON]


# --------------------------------------------------------------------------
# Figure 2 - frequency spectrum
# --------------------------------------------------------------------------
def fig02(fig):
    p = protocol()
    k0 = p["k0"]
    freqs = p["frequencies"]
    ax1 = fig.add_subplot(2, 1, 1)
    ax2 = fig.add_subplot(2, 1, 2, sharex=ax1)
    used = []
    for f in DISC_SWEEP_ALL:
        d = np.load(f)
        if "k" not in d or d["y"].ndim != 3:
            continue
        k = d["k"]
        yk = np.abs(d["y"]).mean(axis=(0, 2))
        Qk = d["Q"].mean(axis=0)
        lbl = f.stem
        ax1.plot(k, yk, lw=1.0, alpha=0.8, label=lbl)
        ax2.plot(k, Qk, lw=1.0, alpha=0.8, label=lbl)
        used.append(f)
    for a in (ax1, ax2):
        for kk in freqs:
            a.axvline(kk, color="#b8860b", ls="--", lw=0.8, alpha=0.7)
        a.axvline(k0, color="#1a7f37", lw=1.2, alpha=0.8)
    ax1.set_ylabel("mean |y| over receivers")
    ax1.set_title("Fig 2  Discovery frequency sweep (24 thetas per template)")
    ax2.set_ylabel("mean Q over thetas")
    ax2.set_xlabel("k (wavenumber)")
    ax1.legend(fontsize=7, ncol=3, loc="upper right")
    ax2.text(0.01, 0.05, "vertical dotted = frozen protocol frequencies %s; "
             "green = k0" % [round(x, 4) for x in freqs],
             transform=ax2.transAxes, fontsize=7, color="#333333")
    fig.text(.5, -.045, "Bright/dark identification is not established. Frozen k0 lies near an internal-energy minimum.", ha="center", fontsize=9, color="#8b1a1a")
    banner(ax1, "FIG02", "spectrum", "PROVISIONAL")
    return ("PROVISIONAL",
            "discovery screening sweep; bright/dark labels withheld pending modal audit",
            used + [DISC_SUMMARY])


# --------------------------------------------------------------------------
# Figure 3 - mode diagnostic (needs modal_audit.npz)
# --------------------------------------------------------------------------
def fig03(fig):
    d = np.load(MODAL_NPZ)
    cm = d["current_modes"]
    for i in range(2):
        ax=fig.add_subplot(1,3,i+1)
        im=ax.imshow(np.abs(cm[:,:,i]).T,origin="lower",extent=(-.5,.5,-.5,.5),cmap="magma")
        ax.set_title(f"Small singular direction {i+1}")
        fig.colorbar(im,ax=ax,shrink=.7)
    ax=fig.add_subplot(1,3,3)
    ix=np.argsort(d['singular_values'])[:20]
    ax.semilogy(np.arange(20),d['singular_values'][ix],label='sigma(M)')
    ax.set_xlabel('Smallest-mode rank');ax.set_title('VIE n64 diagnostic');ax.legend()
    fig.suptitle('Fig 3  Current directions; two-mode dominance NOT established')
    fig.tight_layout(rect=(0,0,1,.9))
    return "PROVISIONAL", "Fixed-grid diagnostic; known spatial discretization error", [MODAL_NPZ,MODAL_META]

# --------------------------------------------------------------------------
# Figure 4 - parameter response / Q map (discovery coarse table)
# --------------------------------------------------------------------------
def fig04(fig):
    p = protocol()
    d = np.load(DISC_COARSE)
    theta = d["theta"]
    Q = d["Q"]
    y = d["y"]
    k = d["k"]
    k0_idx = int(np.argmin(np.abs(k - p["k0"])))
    uk, vk, wk, A = grid_from_theta(theta, Q[:, k0_idx])
    jg = int(np.argmin(np.abs(wk - np.median(wk))))
    jd = int(np.argmin(np.abs(vk - 0.0)))

    def plane(ax, xv, yv, Z, xl, yl, cbl):
        pc = ax.pcolormesh(xv, yv, Z.T, shading="nearest", cmap="coolwarm")
        ax.set_xlabel(xl)
        ax.set_ylabel(yl)
        fig.colorbar(pc, ax=ax, label="Q (ln energy ratio)")

    ax1 = fig.add_subplot(2, 2, 1)
    plane(ax1, uk, vk, A[:, :, jg], "s (gap)", "Delta", 0)
    ax1.set_title(f"Q at k0, gamma at center ({wk[jg]:.3g})")
    ax2 = fig.add_subplot(2, 2, 2)
    Z2 = A[:, jd, :]
    pc2 = ax2.pcolormesh(uk, wk, Z2.T, shading="nearest", cmap="coolwarm")
    ax2.set_xlabel("s (gap)")
    ax2.set_ylabel("gamma (loss)")
    ax2.set_title("Q at k0, Delta = 0")
    fig.colorbar(pc2, ax=ax2, label="Q")
    ax3 = fig.add_subplot(2, 2, 3)
    Z3 = A[jd, :, :]
    pc3 = ax3.pcolormesh(vk, wk, Z3.T, shading="nearest", cmap="coolwarm")
    ax3.set_xlabel("Delta (detuning)")
    ax3.set_ylabel("gamma (loss)")
    ax3.set_title("Q at k0, s at center")
    fig.colorbar(pc3, ax=ax3, label="Q")
    ax4 = fig.add_subplot(2, 2, 4)
    c0 = int(p["primary_channels"][0])
    resp = np.abs(y[:, k0_idx, c0])
    sc = ax4.scatter(resp, Q[:, k0_idx], c=theta[:, 2], cmap="plasma", s=10)
    ax4.set_xlabel("external response |y| (base channel)")
    ax4.set_ylabel("internal ROI energy Q at k0")
    ax4.set_title("external response vs internal task")
    fig.colorbar(sc, ax=ax4, label="gamma")
    fig.suptitle("Fig 4  Discovery 9x9x9 screening: parameter -> Q / external response map\n"
                 "(NOT posterior prevalence evidence)")
    fig.tight_layout(rect=(0, 0, 1, 0.94))
    banner(ax1, "FIG04", "param map", "PROVISIONAL")
    return ("PROVISIONAL",
            "discovery coarse screening table; explicitly not confirmatory prevalence",
            [DISC_COARSE, DISC_SUMMARY])


# --------------------------------------------------------------------------
# Confirmatory posterior helpers (used by Fig 5-7)
# --------------------------------------------------------------------------
def _load_table():
    prior = np.load(TABLE_PRIOR)
    theta = np.asarray(prior["theta"], float)
    weights = np.asarray(prior["weights"], float)
    Q = np.load(TABLE_Q)
    y = np.load(TABLE_Y)
    status = load_json(TABLE_STATUS) if TABLE_STATUS.exists() else {}
    return theta, weights, Q, y, status


def _load_queries():
    d = np.load(QUERIES_NPZ, allow_pickle=True)
    return {k: d[k] for k in d.files}


def _representative_query():
    """Query with maximum base-method posterior std in validation (never test)."""
    files = sorted(VAL_DIR.glob("seed11_noise0.02_base.json"))
    if not files:
        files = sorted(VAL_DIR.glob("*.json"))
    rows = []
    for f in files:
        obj = load_json(f)
        rs = (obj if isinstance(obj, list) else obj.get("rows", []))
        for rr in rs:
            if isinstance(rr, dict) and "std" in rr:
                rows.append((float(rr["std"]), f, rr))
    if not rows:
        raise RuntimeError("no validation rows with 'std' found")
    rows.sort(key=lambda t: t[0])
    return rows[-1][1], rows[-1][2]


def _obs_for_query(qdata, qi, chan_idx, freq_idx):
    y = qdata["y"]
    arr = np.asarray(y)[qi]
    if arr.ndim == 2:                       # (nfreq, nchan)
        return arr[freq_idx, chan_idx]
    if arr.ndim == 1:                       # already the channel vector
        return arr
    if arr.ndim == 3:                       # (nfreq, nchan, 2) re/im packed
        real = arr[freq_idx, chan_idx, 0]
        imag = arr[freq_idx, chan_idx, 1]
        return real + 1j * imag
    raise RuntimeError(f"unsupported query y shape {np.asarray(y).shape}")


def _posterior_for_query(table, qdata, qi):
    """Same observation + likelihood as the table; base channels at k0."""
    p = protocol()
    theta_t, w_t, Q_t, y_t, _ = table
    k0_idx = int(np.argmin(np.abs(np.asarray(p["frequencies"]) - p["k0"])))
    chans = [int(c) for c in p["primary_channels"]]
    sigma = float(p["sigma_primary"])
    y_obs = _obs_for_query(qdata, qi, chans, k0_idx) + sigma * qdata["noise"][qi, k0_idx, chans]
    y_pred = np.asarray(y_t)[:, k0_idx, :][:, chans]
    resid = y_obs[None, :] - y_pred
    ll = -0.5 * (np.sum(resid.real ** 2, axis=1) + np.sum(resid.imag ** 2, axis=1)) / sigma ** 2
    logw = np.log(np.maximum(w_t, 1e-300)) + ll
    logw -= stable_logsumexp(logw)
    w = np.exp(logw)
    Q_post = Q_t[:, k0_idx].copy()
    ess = 1.0 / np.sum(w ** 2)
    return theta_t, w, Q_post, y_obs, k0_idx, ess


def _param_panels(fig, theta, w, title):
    names = ["s (gap)", "Delta (detuning)", "gamma (loss)"]
    axes = []
    for i in range(3):
        ax = fig.add_subplot(2, 2, i + 1)
        ax.hist(theta[:, i], bins=30, weights=w, color="#1565c0")
        ax.set_title(names[i], fontsize=9)
        axes.append(ax)
    ax22 = fig.add_subplot(2, 2, 4)
    ax22.hist2d(theta[:,0], theta[:,1], bins=35, weights=w, cmap="Blues")
    ax22.set_xlabel(names[0]);ax22.set_ylabel(names[1])
    ax22.set_title("Joint posterior mass", fontsize=9)
    fig.suptitle(title)
    fig.tight_layout(rect=(0, 0, 1, 0.93))
    return axes[0]


def fig05(fig):
    table = _load_table()
    qdata = _load_queries()
    src_file, row = _representative_query()
    qi = int(row["query"])
    theta_t, w, Q_post, y_obs, k0_idx, ess = _posterior_for_query(table, qdata, qi)
    ax = _param_panels(fig, theta_t, w,
                       f"Fig 5  Weighted parameter marginals p(theta|y), seed={row.get('seed')} "
                       f"query={qi} (max-variance validation query)\neffective quadrature nodes={ess:.1f} of {len(w)}")
    banner(ax, "FIG05", "posterior", "PROVISIONAL")
    return ("GENERATED",
            f"table 41^3, rep query {qi} (max variance validation), sigma from protocol",
            [TABLE_PRIOR, TABLE_Q, TABLE_Y, QUERIES_NPZ, src_file, GATE_JSON])


def fig06(fig):
    table = _load_table()
    qdata = _load_queries()
    src_file, row = _representative_query()
    qi = int(row["query"])
    theta_t, w, Q_post, y_obs, k0_idx, ess = _posterior_for_query(table, qdata, qi)
    ax = fig.add_subplot(111)
    _, w_t, _, _, _ = table
    ax.hist(Q_post, bins=40, weights=np.full_like(w_t, 1.0 / len(w_t)),
            density=True, color="#999999", alpha=0.45, label="prior predictive Q")
    ax.hist(Q_post, bins=40, weights=w, density=True, color="#1565c0",
            alpha=0.75, label="posterior p(Q|y)")
    qs = weighted_quantiles(Q_post, w, [0.025, 0.5, 0.975])
    for q in qs:
        ax.axvline(q, color="#8b1a1a", ls="--", lw=1.0)
    if "Q_true" in row:
        ax.axvline(float(row["Q_true"]), color="#1a7f37", lw=1.4, label="Q_true")
    ax.set_xlabel("Q (internal ROI log-energy)")
    ax.set_ylabel("density")
    ax.legend(fontsize=8)
    ax.set_title(f"Fig 6  Posterior of the internal task Q, rep query {qi}\n"
                 f"mean={np.sum(w*Q_post):.3f}, sd={np.sqrt(np.sum(w*(Q_post-np.sum(w*Q_post))**2)):.3f}, "
                 f"effective quadrature nodes={ess:.1f}")
    banner(ax, "FIG06", "p(Q|y)", "GENERATED")
    return ("GENERATED", f"p(Q|y) for rep query {qi}", [TABLE_PRIOR, TABLE_Q, TABLE_Y,
            QUERIES_NPZ, src_file, GATE_JSON])


def fig07(fig):
    sys.path.insert(0, str(SRC))
    from cylinders import CylinderModel  # noqa: WPS433
    from physics import Family
    proto = protocol()
    table = _load_table()
    qdata = _load_queries()
    src_file, row = _representative_query()
    qi = int(row["query"])
    theta_t, w, Q_post, y_obs, k0_idx, ess = _posterior_for_query(table, qdata, qi)
    qlo, qhi = weighted_quantiles(Q_post, w, [0.10, 0.90])

    idxsort=np.argsort(Q_post);cw=np.cumsum(w[idxsort])
    ia,ib=idxsort[np.searchsorted(cw,[.10,.90])]
    model = CylinderModel(family=Family(**proto["family"]), order=20, radial_order=48)
    k = float(proto["k0"])
    to_channels = sorted({int(c) for c in
                          (proto["primary_channels"] + proto["extra_fixed"] +
                           proto["extra_optimized"] + proto["near_channels"])})
    x = np.linspace(-0.22, 0.42, 240)
    yy = np.linspace(-0.34, 0.34, 220)
    XX, YY = np.meshgrid(x, yy)
    pts = np.c_[XX.ravel(), YY.ravel()]
    panels = [("low-Q node (10%)", ia, "#1565c0"), ("high-Q node (90%)", ib, "#c0504d")]
    for j, (lab, idx, col) in enumerate(panels):
        theta = theta_t[idx]
        sol = model.solve(theta, k, fields=True)
        assert abs(sol["Q"] - Q_post[idx]) < 1e-9, "Figure/table physical model mismatch"
        E = model.field(sol, pts).reshape(XX.shape)
        ax = fig.add_subplot(2, 3, 1 + j)
        im = ax.pcolormesh(XX, YY, np.abs(E), shading="auto", cmap="inferno")
        fig.colorbar(im, ax=ax, fraction=0.046)
        ax.set_aspect("equal")
        ax.set_title(f"{lab}\nQ={Q_post[idx]:.3f}", fontsize=9)
        for b, ctr in enumerate(sol["state"]["centers"]):
            ax.add_patch(plt.Circle(ctr, model.family.radius, fill=False,
                                    ec="white", lw=1.0))
        ax = fig.add_subplot(2, 3, 4 + j)
        eps = np.real(sol["state"]["epsilon"])
        gam = np.imag(sol["state"]["epsilon"])
        ax.bar([0, 1], eps, color=[col, "#f2c14e"], edgecolor="k")
        ax.set_xticks([0, 1])
        ax.set_xticklabels(["bright", "target"])
        ax.set_ylabel("Re epsilon")
        ax2 = ax.twinx()
        ax2.plot([0, 1], gam, "ko--")
        ax2.set_ylabel("Im epsilon")
        ax.set_title("material", fontsize=9)
    # external measurement comparison
    ax = fig.add_subplot(2, 3, 6)
    ya = model.solve(theta_t[ia], k)["y"][to_channels]
    yb = model.solve(theta_t[ib], k)["y"][to_channels]
    ax.plot(np.abs(ya), "o-", color="#1565c0", label="low-Q node")
    ax.plot(np.abs(yb), "s-", color="#c0504d", label="high-Q node")
    ax.set_xticks(np.arange(len(to_channels)));ax.set_xticklabels(to_channels)
    ax.set_xlabel("physical channel ID")
    ax.set_ylabel("|y|")
    ax.legend(fontsize=8)
    ax.set_title("external measurement", fontsize=9)
    ax=fig.add_subplot(2,3,3)
    yall=table[3];freqs=np.asarray(proto['frequencies'])
    ax.plot(freqs,table[2][ia], 'o-',label='low-Q node')
    ax.plot(freqs,table[2][ib], 's-',label='high-Q node')
    kl=.5*np.sum(abs(yall[ia,k0_idx,:2]-yall[ib,k0_idx,:2])**2)/proto['sigma_primary']**2
    ax.set_title(f'Base measurement KL={kl:.3g}\nROI energy ratio={np.exp(Q_post[ib]-Q_post[ia]):.3g}',fontsize=9)
    ax.set_xlabel('k');ax.set_ylabel('Q');ax.legend(fontsize=7)
    fig.suptitle("Fig 7  Parameter/field pair at the 10% and 90% Q quantiles of the "
                 "representative posterior\n(selected nodes, NOT independent branches)")
    fig.tight_layout(rect=(0, 0, 1, 0.93))
    banner(ax, "FIG07", "field pair", "PROVISIONAL")
    return ("GENERATED", f"10/90% Q nodes for rep query {qi}; not a branch claim",
            [TABLE_PRIOR, TABLE_Q, TABLE_Y, QUERIES_NPZ, src_file, GATE_JSON])


# --------------------------------------------------------------------------
# Figure 8 - distribution of per-query posterior SD
# --------------------------------------------------------------------------
def fig08(fig):
    p = protocol()
    rows = []
    for f in sorted(VAL_DIR.glob("*.json")):
        if f.name == "summary.json":
            continue
        obj = load_json(f)
        for rr in (obj if isinstance(obj, list) else obj.get("rows", [])):
            if isinstance(rr, dict) and "std" in rr and rr.get("noise") == .02:
                rows.append((f.name, rr))
    if not rows:
        raise RuntimeError("no per-query rows found in validation_41_coupled/")
    ax = fig.add_subplot(111)
    by_method = {}
    for name, rr in rows:
        m = rr.get("method", name)
        by_method.setdefault(m, []).append(float(rr["std"]))
    methods = sorted(by_method)
    ax.hist([by_method[m] for m in methods], bins=30, label=methods, histtype="step")
    ax.axvline(p["tau_q"], color="#8b1a1a", ls="--", lw=1.4,
               label=f"tau_q={p['tau_q']:.4g}")
    ax.set_xlabel("per-query posterior SD(Q|y)")
    ax.set_ylabel("count")
    ax.legend(fontsize=7, ncol=2)
    ax.set_title("Fig 8  Prior-predictive distribution of posterior task risk "
                 "(2% noise, same 384 queries per method)")
    banner(ax, "FIG08", "query SD", "PROVISIONAL")
    return "GENERATED", "per-query SD from validation rows", [VAL_SUMMARY] + \
        [f for f in sorted(VAL_DIR.glob("*.json"))]


# --------------------------------------------------------------------------
# Figure 9 - risk comparisons across methods / noise
# --------------------------------------------------------------------------
def fig09(fig):
    p = protocol()
    summ = load_json(VAL_SUMMARY)
    rows = summ.get("summary", summ if isinstance(summ, list) else [])
    if not rows:
        raise RuntimeError("validation summary has no 'summary' list")
    noises = sorted({float(r["noise"]) for r in rows})
    methods = sorted({r["method"] for r in rows})
    ax = fig.add_subplot(1, 2, 1)
    width = 0.8 / max(len(methods), 1)
    for i, m in enumerate(methods):
        xs = np.arange(len(noises))
        ys = [next((r["risk"] for r in rows
                    if float(r["noise"]) == n and r["method"] == m), np.nan)
              for n in noises]
        ax.bar(xs + i * width, ys, width=width, label=m)
    ax.set_xticks(np.arange(len(noises)) + 0.4 - width / 2)
    ax.set_xticklabels([f"{n:.2f}" for n in noises])
    ax.set_xlabel("noise level")
    ax.set_ylabel("task risk R_Q = E Var(Q|y)")
    ax.set_title("posterior task risk by method")
    ax.legend(fontsize=7, ncol=2)
    ax2 = fig.add_subplot(1, 2, 2)
    for m in methods:
        xs = [n for n in noises]
        ys = [next((r["prevalence"] for r in rows
                    if float(r["noise"]) == n and r["method"] == m), np.nan)
              for n in noises]
        ax2.plot(xs, ys, marker="o", label=m)
    ax2.axhline(p["prevalence_min"], color="#8b1a1a", ls="--", lw=1.0,
                label="prevalence_min")
    ax2.set_xlabel("noise level")
    ax2.set_ylabel("P[SD(Q|y) > tau_q]")
    ax2.legend(fontsize=7)
    ax2.set_title("ambiguous-query prevalence")
    method_notes = ", ".join(methods)
    fig.suptitle("Fig 9  Same-frequency vs neighboring-frequency vs near-field "
                 "risk reduction\nmethods: " + method_notes)
    fig.tight_layout(rect=(0, 0, 1, 0.9))
    banner(ax, "FIG09", "risk", "PROVISIONAL")
    return "GENERATED", "validation summary risk/prevalence comparison", [VAL_SUMMARY]


# --------------------------------------------------------------------------
# Figure 10 - coupling ablation
# --------------------------------------------------------------------------
def fig10(fig):
    p = protocol()
    coup = load_json(RESULTS / "validation_33_coupled" / "summary.json")
    unc = load_json(UNCOUPLED_SUMMARY)

    def rows(obj):
        return obj.get("summary", obj if isinstance(obj, list) else [])

    def pick(rs, noise, method, key):
        for r in rs:
            if abs(float(r["noise"]) - noise) < 1e-9 and r["method"] == method:
                return float(r[key])
        return np.nan

    noises = sorted({float(r["noise"]) for r in rows(coup)})
    methods = sorted({r["method"] for r in rows(coup)})
    ax = fig.add_subplot(1, 2, 1)
    width = 0.8 / max(len(methods), 1)
    for i, m in enumerate(methods):
        xs = np.arange(len(noises))
        yc = [pick(rows(coup), n, m, "risk") for n in noises]
        yu = [pick(rows(unc), n, m, "risk") for n in noises]
        ax.bar(xs + i * width, yc, width=width * 0.45, color="#1565c0",
               label=f"{m} coupled" if i == 0 else None)
        ax.bar(xs + i * width + width * 0.45, yu, width=width * 0.45,
               color="#c0504d", label=f"{m} uncoupled" if i == 0 else None)
    ax.set_xticks(np.arange(len(noises)) + 0.4 - width / 2)
    ax.set_xticklabels([f"{n:.2f}" for n in noises])
    ax.set_xlabel("noise level")
    ax.set_ylabel("task risk R_Q")
    ax.set_title("coupled (blue) vs cross-block removed (red)")
    ax2 = fig.add_subplot(1, 2, 2)
    for m in methods:
        xs_c = [pick(rows(coup), n, m, "risk") for n in noises]
        xs_u = [pick(rows(unc), n, m, "risk") for n in noises]
        if np.any(np.isfinite(xs_c)):
            ax2.plot(xs_c, xs_u, marker="o", label=m)
    lim = [0, max(1e-6, *[pick(rows(coup), n, m, "risk") for n in noises for m in methods])]
    ax2.plot(lim, lim, "k:", lw=0.8)
    ax2.set_xlabel("coupled risk")
    ax2.set_ylabel("uncoupled risk")
    ax2.set_title("ablation of cross-body coupling")
    ax2.legend(fontsize=7)
    fig.suptitle("Fig 10  Coupling ablation (cross-body feedback removed)")
    fig.tight_layout(rect=(0, 0, 1, 0.9))
    banner(ax, "FIG10", "ablation", "PROVISIONAL")
    return "GENERATED", "equal 33-cubed rules, unvalidated ablation risk", [RESULTS / "validation_33_coupled" / "summary.json", UNCOUPLED_SUMMARY]


# --------------------------------------------------------------------------
# Figure 11 - modal sufficient-condition decomposition
# --------------------------------------------------------------------------
def fig11(fig):
    d=np.load(MODAL_NPZ);ix=np.argsort(d['singular_values'])[:30]
    terms=[d['injection'][ix],1/d['singular_values'][ix],d['measurement_visibility'][ix],d['task_visibility'][ix]]
    labels=['Material injection ||u* D_E J_theta C||','Amplification 1/sigma(M)','Receiver visibility ||S v||','Task visibility |L_Q v|']
    for j,(term,label) in enumerate(zip(terms,labels)):
        ax=fig.add_subplot(2,2,j+1);ax.semilogy(np.arange(len(ix)),np.maximum(term,1e-20),'.-');ax.set_title(label);ax.set_xlabel('Smallest-singular-direction rank')
    fig.suptitle('Fig 11  All four factors; small sigma alone is insufficient')
    fig.tight_layout(rect=(0,0,1,.9))
    return "PROVISIONAL", "VIE n64 physical injection diagnostic; residual modes in JSON, no sufficient-condition certificate", [MODAL_NPZ,MODAL_META]

# --------------------------------------------------------------------------
# registry & main
# --------------------------------------------------------------------------
FIGURES = [
    ("FIG01", "fig01_geometry_roi_acquisition.png",
     "Canonical two-resonator geometry, measurement and ROI",
     [PROTOCOL_JSON], fig01, False),
    ("FIG02", "fig02_frequency_spectrum.png",
     "Frequency spectrum (discovery sweep)",
     [PROTOCOL_JSON, DISC_SWEEP], fig02, False),
    ("FIG03", "fig03_mode_diagnostic.png",
     "Two dominant modes / radiation visibility",
     [MODAL_NPZ], fig03, True),
    ("FIG04", "fig04_parameter_map.png",
     "External response vs internal ROI energy parameter map",
     [PROTOCOL_JSON, DISC_COARSE], fig04, False),
    ("FIG05", "fig05_parameter_posterior.png",
     "Posterior in (s, Delta, gamma) under a single-frequency measurement",
     [TABLE_PRIOR, TABLE_Q, TABLE_Y, QUERIES_NPZ, VAL_SUMMARY], fig05, True),
    ("FIG06", "fig06_Q_posterior.png",
     "Posterior p(Q|y) of the internal task",
     [TABLE_PRIOR, TABLE_Q, TABLE_Y, QUERIES_NPZ, VAL_SUMMARY], fig06, True),
    ("FIG07", "fig07_material_field_pair.png",
     "Material / field pair at 10% and 90% Q nodes",
     [TABLE_PRIOR, TABLE_Q, TABLE_Y, QUERIES_NPZ, VAL_SUMMARY, SRC / "cylinders.py"],
     fig07, True),
    ("FIG08", "fig08_query_sd_distribution.png",
     "Prior-predictive distribution of posterior task risk",
     [VAL_SUMMARY], fig08, True),
    ("FIG09", "fig09_risk_comparison.png",
     "Task-risk reduction: same-frequency vs neighbor vs near-field",
     [VAL_SUMMARY], fig09, True),
    ("FIG10", "fig10_coupling_ablation.png",
     "Coupling ablation",
     [VAL_SUMMARY, UNCOUPLED_SUMMARY], fig10, True),
    ("FIG11", "fig11_modal_decomposition.png",
     "Material injection / resolvent / visibility / task decomposition",
     [MODAL_NPZ], fig11, True),
]


def main():
    ts = datetime.now(timezone.utc).isoformat()
    gate = gate_passed()
    print(f"A15 figure build  {ts}\nreference gate passed: {gate}")
    for figid, name, title, needs, fn, conf in FIGURES:
        make_figure(figid, name, title, needs, fn, confirmatory=conf)

    manifest = {
        "generated_utc": ts,
        "reference_gate_passed": gate,
        "gate_file_present": GATE_JSON.exists(),
        "figures": MANIFEST,
    }
    (DELEG_DIR / "FIGURES_MANIFEST.json").write_text(json.dumps(manifest, indent=2))

    lines = ["# A15 figures manifest", "",
             f"Generated (UTC): {ts}",
             f"reference_validation_gate.json passed: **{gate}** "
             f"(present: {GATE_JSON.exists()})", "",
             "Status: GENERATED = drawn | PROVISIONAL = discovery or gate not passed "
             "| NOT_RUN = input missing | ERROR = inputs present, draw failed.", "",
             "| id | status | figure | sources | note |",
             "| --- | --- | --- | --- | --- |"]
    for m in MANIFEST:
        lines.append("| {} | {} | `{}` | {} | {} |".format(
            m["id"], m["status"], m["figure"],
            ", ".join(f"`{s}`" for s in m["sources"]) or "-",
            m["note"] or "-"))
    lines += ["", "## Missing inputs (figures not yet run)", ""]
    for m in MANIFEST:
        if m["missing"]:
            lines.append(f"- {m['id']}: " + ", ".join(f"`{x}`" for x in m["missing"]))
    if not any(m["missing"] for m in MANIFEST):
        lines.append("- none")
    (DELEG_DIR / "FIGURES_MANIFEST.md").write_text("\n".join(lines) + "\n")
    print(f"\nmanifest: {rel(DELEG_DIR / 'FIGURES_MANIFEST.md')}")


if __name__ == "__main__":
    main()
