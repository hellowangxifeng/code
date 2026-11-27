"""Figure 5: hydrodynamic-model consistency and measured-wave bridge."""
from __future__ import annotations

import os

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import scipy.io as sio
from matplotlib import gridspec
from scipy import signal

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
FIG = os.path.join(ROOT, "figures")
RES = os.path.join(ROOT, "results")
plt.rcParams.update({
    "font.family": "DejaVu Sans", "font.size": 8, "axes.labelsize": 8,
    "axes.titlesize": 8, "legend.fontsize": 6.4, "xtick.labelsize": 7,
    "ytick.labelsize": 7, "lines.linewidth": 1.0, "axes.linewidth": 0.6,
    "savefig.dpi": 300, "figure.dpi": 150, "axes.grid": True,
    "grid.linewidth": 0.4, "grid.alpha": 0.35,
})


def style(ax):
    ax.tick_params(direction="in", top=True, right=True)


def fld(o, name):
    v = o[name]
    while isinstance(v, np.ndarray) and v.dtype == object and v.ndim == 0:
        v = v.item()
    return np.asarray(v, float)


def fig5():
    fig = plt.figure(figsize=(7.0, 4.3))
    gs = gridspec.GridSpec(2, 2, hspace=0.42, wspace=0.32)
    ax = [fig.add_subplot(gs[j]) for j in range(4)]
    z = np.load(os.path.join(ROOT, "data", "processed", "fowt_motion.npz"), allow_pickle=True)
    # (a) measured FOWT pitch vs RIS model pitch under the same measured wave
    t, _, rpy = z["t_8"], z["pos_8"], z["rpy_8"]
    s = (t > 15) & (t < 25)
    ax[0].plot(t[s] - 15, rpy[s, 1] - rpy[s, 1].mean(), "C3", lw=1.0, label="tank platform (measured)")
    M = sio.loadmat(os.path.join(RES, "measured_wave_ris_8.mat"), squeeze_me=True)
    o = M["o"]
    tt = fld(o, "t"); pp = fld(o, "pitchDeg") - np.mean(fld(o, "pitchDeg"))
    ax[0].plot(tt, pp, "C0", lw=1.0)
    ax[0].set_xlabel("time (s)"); ax[0].set_ylabel("pitch (deg)")
    ax[0].set_xlim(0, 10)
    ax[0].text(5.6, -4.6, "1 m RIS model\n3.77° rms", color="C0", fontsize=6.2, ha="center")
    ax[0].text(5.6, 0.55, "tank platform (measured)\n1.10° rms", color="C3", fontsize=6.2, ha="center")
    ax[0].set_title("(a) pitch, 45 mm/0.8 Hz")
    # (b) tilt envelope
    M2 = sio.loadmat(os.path.join(RES, "sea_model_envelope.mat"), squeeze_me=True)
    ax[1].plot(M2["U"], M2["tilt"], "x-", color="C7", ms=4, lw=0.8, label="wind-sea model")
    exp = [0.43, 0.14, 0.46, 0.41, 0.30, 0.76, 0.75, 0.44, 1.10]
    ax[1].axhspan(min(exp), max(exp), color="C0", alpha=0.18, label="tank experiment")
    ax[1].set_xlabel("wind speed (m/s)"); ax[1].set_ylabel("tilt rms (deg)")
    ax[1].set_ylim(0, 8.5); ax[1].legend(frameon=False, loc="upper left")
    ax[1].set_title("(b) motion severity")
    # (c) measured wave spectra vs Neumann
    for idx, lab, col in [(6, "30 mm/0.8 Hz", "C0"), (8, "45 mm/0.8 Hz", "C3")]:
        eta = z[f"eta_{idx}"][:, 5]; t = z[f"t_{idx}"]; s = t > 15
        f, P = signal.welch(eta[s], fs=100, nperseg=2048, noverlap=1024)
        ax[2].semilogy(f, P, col, label=lab)
    ff = np.linspace(0.2, 3, 500); w = 2 * np.pi * ff
    for U, ls in [(3, "--"), (5, ":")]:
        S = 3.05 * np.pi / 4 * w ** -6 * np.exp(-2 * 9.81 ** 2 / (U ** 2 * w ** 2))
        S = S * 2 * np.pi
        ax[2].semilogy(ff, S / S.max() * 1e-3, "k" + ls, lw=0.9, label=f"Neumann U={U}")
    ax[2].set_xlim(0.2, 3); ax[2].set_ylim(1e-7, 1e-1)
    ax[2].set_xlabel("frequency (Hz)"); ax[2].set_ylabel("PSD (m$^2$/Hz)")
    ax[2].legend(frameon=False, fontsize=5.8); ax[2].set_title("(c) wave spectra")
    # (d) model pitch and link under the 45 mm/0.8 Hz measured wave
    t2 = fld(o, "t"); p2 = fld(o, "pitchDeg") - np.mean(fld(o, "pitchDeg"))
    ax[3].plot(t2, p2, "C3", lw=0.7)
    ax[3].set_xlabel("time (s)"); ax[3].set_ylabel("RIS pitch (deg)")
    ax[3].set_xlim(0, 20)
    ax[3].set_title("(d) RIS motion and link, 45 mm/0.8 Hz")
    ax3 = ax[3].twinx()
    d = np.load(os.path.join(RES, "model_replay_8.npz"))
    ax3.plot(d["t"], d["collimated_static_p"] * 1e3, color="C1", lw=0.9, label="collimated static")
    ax3.plot(d["t"], d["point_adaptive_p"] * 1e3, color="C0", lw=0.9, label="point adaptive")
    ax3.set_ylabel("received power (mW)"); ax3.set_ylim(0, 17)
    ax3.legend(frameon=False, loc="lower right", fontsize=5.8)
    for a in ax:
        style(a)
    fig.savefig(os.path.join(FIG, "Fig5_model.png"), bbox_inches="tight")
    plt.close(fig)


if __name__ == "__main__":
    fig5()
    print("ok")
