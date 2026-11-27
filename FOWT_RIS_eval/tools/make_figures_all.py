"""All publication figures for the FOWT-RIS experimental evaluation."""
from __future__ import annotations

import csv
import os

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import scipy.io as sio
from matplotlib import gridspec
from matplotlib.patches import Circle
from scipy import signal

import fx_optics as fo
import replay as rp

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
FIG = os.path.join(ROOT, "figures")
RES = os.path.join(ROOT, "results")
os.makedirs(FIG, exist_ok=True)
plt.rcParams.update({
    "font.family": "DejaVu Sans", "font.size": 8, "axes.labelsize": 8,
    "axes.titlesize": 8, "legend.fontsize": 6.4, "xtick.labelsize": 7,
    "ytick.labelsize": 7, "lines.linewidth": 1.0, "axes.linewidth": 0.6,
    "savefig.dpi": 300, "figure.dpi": 150, "axes.grid": True,
    "grid.linewidth": 0.4, "grid.alpha": 0.35,
})
COND = list(csv.DictReader(open(os.path.join(ROOT, "data", "processed", "conditions.csv"), encoding="utf-8")))
DECAY = list(csv.DictReader(open(os.path.join(ROOT, "data", "processed", "decay_summary.csv"), encoding="utf-8")))
SUMMARY = list(csv.DictReader(open(os.path.join(RES, "replay_summary.csv"), encoding="utf-8")))
MODEL = list(csv.DictReader(open(os.path.join(RES, "model_motion_replay.csv"), encoding="utf-8")))
LABELS = [f"{c['amp_mm']}/{c['freq_hz']}" for c in COND]
TILT = np.array([float(c["tilt_rms_deg"]) for c in COND])


def motion(i):
    z = np.load(os.path.join(ROOT, "data", "processed", "fowt_motion.npz"), allow_pickle=True)
    return z[f"t_{i}"], z[f"pos_{i}"], z[f"rpy_{i}"]


def replay(i):
    return np.load(os.path.join(RES, f"replay_{i:02d}.npz"))


def style(ax):
    ax.tick_params(direction="in", top=True, right=True)


# ---------------------------------------------------------------- Fig 1
def fig1():
    fig = plt.figure(figsize=(7.0, 4.4))
    gs = gridspec.GridSpec(2, 3, hspace=0.45, wspace=0.38)
    ax = [fig.add_subplot(gs[j]) for j in range(6)]
    t30, _, r30 = motion(6)
    t45, _, r45 = motion(8)
    s30 = (t30 > 20) & (t30 < 25)
    s45 = (t45 > 20) & (t45 < 25)
    ax[0].plot(t30[s30] - 20, r30[s30, 0] - r30[s30, 0].mean(), "C0", label="roll 30 mm")
    ax[0].plot(t45[s45] - 20, r45[s45, 0] - r45[s45, 0].mean(), "C3", label="roll 45 mm")
    ax[0].plot(t45[s45] - 20, r45[s45, 1] - r45[s45, 1].mean(), "C3--", label="pitch 45 mm")
    ax[0].set_xlabel("time (s)"); ax[0].set_ylabel("angle (deg)")
    ax[0].legend(loc="lower left", fontsize=5.8, frameon=True, framealpha=0.9, edgecolor="none")
    ax[0].set_title("(a) attitude, 0.8 Hz")
    for idx, lab, col, ls in [(6, "30 mm", "C0", "-"), (8, "45 mm", "C3", "-")]:
        t, _, rpy = motion(idx); s = t > 15
        f, Pw = signal.welch(rpy[s, 0] - rpy[s, 0].mean(), fs=100, nperseg=1024, noverlap=512)
        ax[1].loglog(f, Pw, col + ls, label=f"roll {lab}")
        f, Pw = signal.welch(rpy[s, 1] - rpy[s, 1].mean(), fs=100, nperseg=1024, noverlap=512)
        ax[1].loglog(f, Pw, col + "--", label=f"pitch {lab}")
    ax[1].set_xlim(0.2, 5); ax[1].set_xlabel("frequency (Hz)"); ax[1].set_ylabel("PSD (deg$^2$/Hz)")
    ax[1].legend(frameon=False, ncol=2, loc="lower left", fontsize=5.8); ax[1].set_title("(b) attitude spectra")
    f = np.array([float(c["freq_hz"]) for c in COND])
    rao = np.array([float(c["heave_rao_at_f"]) for c in COND])
    rr = np.array([float(c["roll_rao_at_f"]) for c in COND])
    ax[2].semilogy(f, rao, "o", ms=4, mfc="none", label="heave")
    ax[2].set_xlabel("wave frequency (Hz)"); ax[2].set_ylabel("heave RAO (m/m)")
    ax2 = ax[2].twinx(); ax2.semilogy(f + 0.01, rr, "s", ms=3, mfc="none", color="C3")
    ax2.set_ylabel("roll RAO (deg/m)", color="C3"); ax2.tick_params(axis="y", colors="C3", direction="in")
    ax[2].legend(frameon=False, loc="lower left"); ax[2].set_title("(c) measured RAO")
    import fowt_data as fd
    for j, (lab, col) in enumerate([("heave", "C0"), ("pitch", "C3")]):
        d = fd.load_motion(os.path.join(fd.RAW_DIR, lab.capitalize() + "Decay.mat"))
        ch = d["pos_mm"][:, 2] if lab == "heave" else d["rpy_deg"][:, 1]
        y = ch - ch.mean()
        sos = signal.butter(4, [0.15, 1.5], "bandpass", fs=d["fs"], output="sos")
        yf = signal.sosfiltfilt(sos, y)
        ax[3].plot(d["t"], yf / np.max(np.abs(yf)), col,
                   label=f"{lab}: $T_n$={float(DECAY[j]['period_s']):.2f} s, $\\zeta$={100*float(DECAY[j]['damping_ratio']):.2f}%")
    ax[3].set_xlim(0, 40); ax[3].set_xlabel("time (s)"); ax[3].set_ylabel("normalised decay")
    ax[3].legend(frameon=False, loc="lower right"); ax[3].set_title("(d) decay tests")
    for fq, col, mk in [(0.4, "C0", "o"), (0.8, "C1", "s"), (1.2, "C2", "^")]:
        cc = sorted([c for c in COND if abs(float(c["freq_hz"]) - fq) < 1e-6], key=lambda c: float(c["amp_mm"]))
        ax[4].plot([float(c["amp_mm"]) for c in cc], [float(c["tilt_rms_deg"]) for c in cc], mk + "-", color=col, ms=4, label=f"{fq} Hz")
    ax[4].set_xlabel("wave amplitude (mm)"); ax[4].set_ylabel("tilt rms (deg)")
    ax[4].legend(frameon=False); ax[4].set_title("(e) tilt severity")
    M = sio.loadmat(os.path.join(RES, "sea_model_envelope.mat"), squeeze_me=True)
    ax[5].plot(M["U"], M["tilt"], "x-", color="C7", ms=4, lw=0.8, label="RIS model")
    ax[5].axhspan(TILT.min(), TILT.max(), color="C0", alpha=0.18, label="experiment range")
    ax[5].set_xlabel("wind speed (m/s)"); ax[5].set_ylabel("tilt rms (deg)")
    ax[5].set_ylim(0, 8.2); ax[5].legend(frameon=False, loc="upper left"); ax[5].set_title("(f) model vs experiment")
    for a in ax:
        style(a)
    fig.savefig(os.path.join(FIG, "Fig1_dataset.png"), bbox_inches="tight")
    plt.close(fig)


# ---------------------------------------------------------------- Fig 2
def coverage_curve(d, qin, P):
    B = d["B"]
    gx = np.linspace(-1, 1, 401)
    dens = np.zeros((401, 401))
    lo = np.maximum(-1, qin - B); hi = np.minimum(1, qin + B)
    i0 = np.searchsorted(gx, lo[:, 0]); i1 = np.searchsorted(gx, hi[:, 0])
    j0 = np.searchsorted(gx, lo[:, 1]); j1 = np.searchsorted(gx, hi[:, 1])
    for i in range(len(qin)):
        dens[j0[i]:j1[i], i0[i]:i1[i]] += P[i]
    QX, QY = np.meshgrid(gx, gx)
    ang = np.degrees(np.arcsin(np.minimum(1, np.hypot(QX, QY))))
    out = []
    for a in np.arange(45, 70.01, 0.5):
        m = ang <= a
        out.append((a, dens[m].max()))
    return np.array(out)


def fig2():
    d = fo.load_design()
    ref = np.asarray(d["referenceXY"], float)
    P = np.asarray(d["unitDirectPowerW"], float); P /= P.sum()
    fig = plt.figure(figsize=(7.0, 6.4))
    gs = gridspec.GridSpec(3, 2, hspace=0.45, wspace=0.38)
    ax = [fig.add_subplot(gs[j, k]) for j in range(3) for k in range(2)]
    # (a) command plane
    B = float(d["maxCmdDeg"])
    ax[0].add_patch(plt.Rectangle((-B, -B), 2 * B, 2 * B, fill=False, ec="C7", ls="--", lw=0.8))
    ax[0].plot(-B, B, "o", color="C3", ms=7)
    ax[0].annotate("1600/1600 units\nat $(-15.414^\\circ, +15.414^\\circ)$", (-B, B),
                   xytext=(0.12, 0.55), textcoords="axes fraction", fontsize=6.5,
                   arrowprops=dict(arrowstyle="-", lw=0.6, color="C7"))
    ax[0].set_xlim(-17, 17); ax[0].set_ylim(-17, 17)
    ax[0].set_xlabel("x command (deg)"); ax[0].set_ylabel("y command (deg)")
    ax[0].set_title("(a) production commands")
    # (b) realised direction vs unit y
    qin = fo.q_in(d, np.eye(3), np.asarray(d["pRis"], float), ref)
    qout = qin + np.asarray(d["currentGradient"], float)
    ang = np.degrees(np.arcsin(np.hypot(qout[:, 0], qout[:, 1])))
    sc = ax[1].scatter(ref[:, 1] * 1e3, ang, c=P, s=5, cmap="viridis")
    ax[1].set_xlabel("unit y (mm)"); ax[1].set_ylabel("realised direction (deg)")
    ax[1].set_title("(b) direction spread")
    # (c,d) spot maps
    from matplotlib.colors import LogNorm
    norm = LogNorm(vmin=0.02, vmax=8)
    for j, (key, title) in enumerate([("current", "(c) production spot"), ("improved", "(d) collimated spot")]):
        a = ax[2 + j]
        x = np.asarray(d[key + "X"]); y = np.asarray(d[key + "Y"]); I = np.asarray(d[key + "I"])
        im = a.imshow(I, origin="lower", extent=[x[0], x[-1], y[0], y[-1]],
                      cmap="inferno", norm=norm, aspect="equal")
        c = d[key + "Center"]
        a.add_patch(Circle(c, 0.05, fill=False, ec="cyan", lw=0.8))
        a.set_xlabel("x (m)"); a.set_ylabel("y (m)")
        a.set_title(title)
        a.set_xlim(2.72, 3.90); a.set_ylim(-6.15, -4.42)
    cb = fig.colorbar(im, ax=ax[2:4], fraction=0.035, pad=0.02)
    cb.set_label("irradiance (W/m$^2$)")
    # (e) aperture power vs radius
    for key, col, lab in [("current", "C3", "production"), ("improved", "C0", "collimated")]:
        rr, pp = [], []
        for r in (0.01, 0.02, 0.03, 0.05, 0.075, 0.1, 0.15):
            rr.append(r)
            x = np.asarray(d[key + "X"]); y = np.asarray(d[key + "Y"]); I = np.asarray(d[key + "I"])
            X, Y = np.meshgrid(x, y); c = d[key + "Center"]
            pp.append(I[(X - c[0]) ** 2 + (Y - c[1]) ** 2 <= r * r].sum() * (x[1] - x[0]) * (y[1] - y[0]) * 1e3)
        ax[4].loglog(np.array(rr) * 1e3, pp, "o-", color=col, ms=3, label=lab)
    ax[4].set_xlabel("aperture radius (mm)"); ax[4].set_ylabel("nominal power (mW)")
    ax[4].legend(frameon=False); ax[4].set_title("(e) aperture coupling")
    # (f) coverage vs target angle
    cv = coverage_curve(d, qin, P)
    ax[5].plot(cv[:, 0], 100 * cv[:, 1], "C2-", lw=1.1)
    ax[5].axhline(99.9, color="C7", ls=":", lw=0.8)
    ax[5].axvline(53.9, color="C7", ls="--", lw=0.8)
    ax[5].set_xlabel("common target angle from zenith (deg)")
    ax[5].set_ylabel("covered power (%)")
    ax[5].set_title("(f) feasibility trade-off")
    for a in ax:
        style(a)
    fig.savefig(os.path.join(FIG, "Fig2_design.png"), bbox_inches="tight")
    plt.close(fig)


# ---------------------------------------------------------------- Fig 3
def fig3():
    fig = plt.figure(figsize=(7.0, 4.5))
    gs = gridspec.GridSpec(2, 3, hspace=0.48, wspace=0.38)
    ax = [fig.add_subplot(gs[j]) for j in range(6)]
    for col, idx in enumerate([8, 6]):
        z = replay(idx); t = z["t"]; s = (t > 20) & (t < 40)
        ax[col].plot(t[s], z["current_static_p"][s] * 1e3, color="C7", lw=0.8, label="current static")
        ax[col].plot(t[s], z["collimated_static_p"][s] * 1e3, color="C1", lw=0.8, label="collimated static")
        ax[col].plot(t[s], z["point_adaptive_p"][s] * 1e3, color="C0", lw=1.0, label="point adaptive")
        ax[col].set_xlabel("time (s)"); ax[col].set_ylabel("power (mW)")
        ax[col].set_title(f"({'ab'[col]}) {LABELS[idx]} (mm/Hz)")
    handles, labels = ax[0].get_legend_handles_labels()
    fig.legend(handles, labels, loc="upper center", ncol=3, frameon=False,
               bbox_to_anchor=(0.5, 0.985))
    fig.subplots_adjust(top=0.9)
    vals = {k: [] for k in ["current_static_p", "collimated_static_p", "point_adaptive_p"]}
    for i in range(9):
        z = replay(i)
        for k in vals:
            vals[k].append(z[k] * 1e3)
    for k, col, lab in [("current_static_p", "C7", "current static"),
                        ("collimated_static_p", "C1", "collimated static"),
                        ("point_adaptive_p", "C0", "point adaptive")]:
        v = np.sort(np.concatenate(vals[k]))
        ax[2].semilogx(v, np.arange(1, len(v) + 1) / len(v), col, label=lab)
    ax[2].axvline(5, color="k", ls=":", lw=0.7)
    ax[2].set_xlabel("received power (mW)"); ax[2].set_ylabel("CDF")
    ax[2].set_title("(c) all conditions")
    x = np.arange(len(LABELS))
    for j, (sc, col) in enumerate([("current static", "C7"), ("collimated static", "C1"), ("point adaptive", "C0")]):
        m = np.array([float(r["p_mean_mw"]) for r in SUMMARY if r["scenario"] == sc])
        ax[3].bar(x + (j - 1) * 0.25, m, width=0.23, color=col, label=sc, linewidth=0)
    ax[3].axhline(5, color="k", ls=":", lw=0.7)
    ax[3].set_xticks(x); ax[3].set_xticklabels(LABELS, rotation=60, ha="right", fontsize=5.6)
    ax[3].set_ylabel("mean power (mW)")
    ax[3].set_ylim(0, 17); ax[3].set_title("(d) per condition")
    for sc, col, mk in [("current static", "C7", "o"), ("collimated static", "C1", "s"),
                        ("collimated adaptive", "C2", "^"), ("point adaptive", "C0", "d")]:
        off = np.array([float(r["offset_rms_mm"]) for r in SUMMARY if r["scenario"] == sc])
        ax[4].plot(TILT, off, mk, color=col, ms=4, label=sc)
    ax[4].set_xlabel("tilt rms (deg)"); ax[4].set_ylabel("centroid rms (mm)")
    ax[4].legend(frameon=False, fontsize=5.8); ax[4].set_title("(e) pointing error")
    for sc, col, mk in [("current static", "C7", "o"), ("collimated static", "C1", "s"), ("point adaptive", "C0", "d")]:
        av = np.array([float(r["avail_3db"]) for r in SUMMARY if r["scenario"] == sc])
        ax[5].plot(TILT, av, mk, color=col, ms=4, label=sc)
    ax[5].set_xlabel("tilt rms (deg)"); ax[5].set_ylabel("availability (-3 dB)")
    ax[5].set_ylim(0.4, 1.02); ax[5].set_title("(f) availability")
    for a in ax:
        style(a)
    fig.savefig(os.path.join(FIG, "Fig3_replay.png"), bbox_inches="tight")
    plt.close(fig)


# ---------------------------------------------------------------- Fig 4
def fig4():
    rows = list(csv.DictReader(open(os.path.join(RES, "control_study.csv"), encoding="utf-8")))
    target = "FOWTmotion_11_3_45mmp8Hz.mat"
    fig = plt.figure(figsize=(7.0, 4.0))
    gs = gridspec.GridSpec(2, 2, hspace=0.42, wspace=0.32)
    ax = [fig.add_subplot(gs[j]) for j in range(4)]
    for pred, col, mk in [("ZOH", "C1", "o"), ("linear", "C0", "s")]:
        rr = sorted([r for r in rows if r["record"] == target and r["predictor"] == pred],
                    key=lambda r: float(r["delay_ms"]))
        ax[0].plot([float(r["delay_ms"]) for r in rr], [float(r["offset_rms_mm"]) for r in rr],
                   mk + "-", color=col, ms=4, label=pred)
        ax[1].plot([float(r["delay_ms"]) for r in rr], [float(r["loss_mean_db"]) for r in rr],
                   mk + "-", color=col, ms=4, label=pred)
    ax[0].set_xlabel("processing latency (ms)"); ax[0].set_ylabel("centroid rms (mm)")
    ax[0].legend(frameon=False); ax[0].set_title("(a) latency, 45 mm/0.8 Hz")
    ax[1].set_xlabel("processing latency (ms)"); ax[1].set_ylabel("mean loss (dB)")
    ax[1].legend(frameon=False); ax[1].set_title("(b) mean power loss")
    for rec, col, mk, lab in [(target, "C0", "o", "45 mm/0.8 Hz"), ("FOWTmotion_11_3_30mmp4Hz.mat", "C2", "s", "30 mm/0.4 Hz")]:
        rr = sorted([r for r in rows if r["record"] == rec and r["predictor"] == "hold"],
                    key=lambda r: float(r["delay_ms"]))
        ax[2].plot([float(r["delay_ms"]) for r in rr], [float(r["offset_rms_mm"]) for r in rr],
                   mk + "-", color=col, ms=4, label=lab)
    ax[2].set_xlabel("command update interval (ms)"); ax[2].set_ylabel("centroid rms (mm)")
    ax[2].legend(frameon=False); ax[2].set_title("(c) sample-and-hold")
    rr = sorted([r for r in rows if r["record"] == target and r["predictor"].startswith("noise")],
                key=lambda r: float(r["predictor"].replace("noise", "").replace("deg", "")))
    ax[3].plot([float(r["predictor"].replace("noise", "").replace("deg", "")) for r in rr],
               [float(r["offset_rms_mm"]) for r in rr], "o-", color="C3", ms=4)
    ax[3].set_xlabel("attitude estimate noise rms (deg)"); ax[3].set_ylabel("centroid rms (mm)")
    ax[3].set_title("(d) attitude-noise sensitivity")
    for a in ax:
        style(a)
    fig.savefig(os.path.join(FIG, "Fig4_control.png"), bbox_inches="tight")
    plt.close(fig)


if __name__ == "__main__":
    fig1(); fig2(); fig3(); fig4()
    print("figures 1-4 done")
