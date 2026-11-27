"""Build a clean, documented FOWT motion dataset from the raw wave-tank files.

Outputs (all inside FOWT_RIS_eval/data/processed):
  fowt_motion.npz   arrays: t, pos_m, rpy_deg, fs, condition id, wave eta
  conditions.csv    one row per regular-wave run
  decay_summary.csv identified natural periods / damping ratios
"""
from __future__ import annotations

import csv
import glob
import os

import numpy as np
from scipy import signal

import fowt_data as fd

OUT = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "data", "processed")
os.makedirs(OUT, exist_ok=True)


def condition_table():
    rows = []
    for p in sorted(glob.glob(os.path.join(fd.RAW_DIR, "FOWTmotion_*.mat"))):
        b = os.path.basename(p)
        m = fd.re.search(r"(\d+)mm(\d+)p(\d+)Hz", b)
        if m:
            amp, f = int(m.group(1)), float(m.group(2) + "." + m.group(3))
        else:
            m = fd.re.search(r"(\d+)mmp(\d+)Hz", b)
            amp, f = int(m.group(1)), float("0." + m.group(2))
        rows.append(dict(motion=p, name=b, amp_mm=amp, freq_hz=f))
    # attach the closest wave-gauge file by (amp, freq)
    waves = []
    for p in sorted(glob.glob(os.path.join(fd.RAW_DIR, "*.txt"))):
        b = os.path.basename(p)
        fm, am = fd.FREQ_RE.search(b), fd.COND_RE.search(b)
        if fm and am:
            waves.append(dict(wave=p, wname=b, amp_mm=int(round(float(am.group(1)) * 1000)),
                              freq_hz=float(fm.group(1))))
    for r in rows:
        cands = [w for w in waves if w["amp_mm"] == r["amp_mm"] and abs(w["freq_hz"] - r["freq_hz"]) < 1e-6]
        r["wave"] = cands[0]["wave"] if cands else None
        r["wname"] = cands[0]["wname"] if cands else None
    return rows


def steady(t, t0=15.0):
    return t >= t0


def spectrum(y, fs, nperseg=None, detrend="linear"):
    if nperseg is None:
        nperseg = min(len(y), 4096)
    f, P = signal.welch(y, fs=fs, nperseg=nperseg, noverlap=nperseg // 2, detrend=detrend)
    return f, P


def main():
    rows = condition_table()
    store = {}
    cond_rows = []
    for i, r in enumerate(rows):
        mo = fd.load_motion(r["motion"])
        t = mo["t"] - mo["t"][0]
        pos = mo["pos_m"].copy()
        pos -= pos[0]  # relative to first sample; a fixed offset is irrelevant
        rpy = mo["rpy_deg"].copy()
        # remove the first-sample value only for reporting; keep the physical trim
        fs = mo["fs"]
        eta = None
        if r["wave"]:
            w = fd.load_wave(r["wave"])
            eta = w
        # resample motion to a uniform 100 Hz grid for the replay
        fs_out = 100.0
        t_out = np.arange(0.0, t[-1], 1.0 / fs_out)
        pos_out = np.column_stack([np.interp(t_out, t, pos[:, j]) for j in range(3)])
        rpy_out = np.column_stack([np.interp(t_out, t, rpy[:, j]) for j in range(3)])
        store[f"pos_{i}"] = pos_out
        store[f"rpy_{i}"] = rpy_out
        store[f"t_{i}"] = t_out
        if eta is not None:
            tw = eta["t"]
            eta_out = np.column_stack([np.interp(t_out, tw, eta["eta"][:, j]) for j in range(6)])
            store[f"eta_{i}"] = eta_out
        sel = t_out >= 15.0
        stats = {}
        for j, nm in enumerate(["surge", "sway", "heave"]):
            stats[f"{nm}_std_mm"] = float(np.std(pos_out[sel, j] * 1e3))
            stats[f"{nm}_ptp_mm"] = float(np.ptp(pos_out[sel, j] * 1e3))
        for j, nm in enumerate(["roll", "pitch", "yaw"]):
            stats[f"{nm}_std_deg"] = float(np.std(rpy_out[sel, j]))
            stats[f"{nm}_ptp_deg"] = float(np.ptp(rpy_out[sel, j]))
            stats[f"{nm}_mean_deg"] = float(np.mean(rpy_out[sel, j]))
        stats["tilt_rms_deg"] = float(np.sqrt(stats["roll_std_deg"] ** 2 + stats["pitch_std_deg"] ** 2))
        if eta is not None:
            ew = eta["eta"][:, 5]
            # steady wave and RAO estimate against gauge 8
            e = np.interp(t_out, eta["t"], ew)
            selw = t_out >= 15.0
            stats["wave_std_mm"] = float(np.std(e[selw] * 1e3))
            nperseg = 2048
            f, Pxy = signal.csd(e[selw], pos_out[selw, 2], fs=fs_out, nperseg=nperseg,
                                 noverlap=nperseg // 2, detrend="linear")
            f2, Pxx = signal.welch(e[selw], fs=fs_out, nperseg=nperseg,
                                   noverlap=nperseg // 2, detrend="linear")
            f3, Pyy = signal.welch(pos_out[selw, 2], fs=fs_out, nperseg=nperseg,
                                   noverlap=nperseg // 2, detrend="linear")
            coh = np.abs(Pxy) ** 2 / np.maximum(Pxx * Pyy, 1e-30)
            band = np.abs(f - r["freq_hz"]) < 0.08
            stats["heave_rao_at_f"] = float(np.sqrt(Pyy[band].mean() / max(Pxx[band].mean(), 1e-30)))
            stats["heave_wave_coh_at_f"] = float(coh[band].mean())
            for j, nm in enumerate(["roll", "pitch"]):
                _, Pxy2 = signal.csd(e[selw], rpy_out[selw, j], fs=fs_out, nperseg=nperseg,
                                     noverlap=nperseg // 2, detrend="linear")
                _, Pyy2 = signal.welch(rpy_out[selw, j], fs=fs_out, nperseg=nperseg,
                                       noverlap=nperseg // 2, detrend="linear")
                stats[f"{nm}_rao_at_f"] = float(np.sqrt(Pyy2[band].mean() / max(Pxx[band].mean(), 1e-30)))
                stats[f"{nm}_coh_at_f"] = float((np.abs(Pxy2) ** 2 / np.maximum(Pxx * Pyy2, 1e-30))[band].mean())
        cond_rows.append(dict(index=i, name=r["name"], amp_mm=r["amp_mm"], freq_hz=r["freq_hz"],
                              wave_name=r["wname"] or "", **stats))
        print(f"{r['name']:45s} heave std {stats['heave_std_mm']:7.2f} mm  tilt rms {stats['tilt_rms_deg']:5.3f} deg"
              + (f"  RAO {stats['heave_rao_at_f']:.2f} coh {stats['heave_wave_coh_at_f']:.2f}" if 'heave_rao_at_f' in stats else ""))
    np.savez_compressed(os.path.join(OUT, "fowt_motion.npz"),
                        n=len(rows), meta=np.array([r["name"] for r in rows]), **store)
    with open(os.path.join(OUT, "conditions.csv"), "w", newline="", encoding="utf-8") as fh:
        w = csv.DictWriter(fh, fieldnames=list(cond_rows[0].keys()))
        w.writeheader()
        w.writerows(cond_rows)

    # decay tests
    decays = []
    for p in sorted(glob.glob(os.path.join(fd.RAW_DIR, "*Decay.mat"))):
        mo = fd.load_motion(p)
        t = mo["t"]
        mode = os.path.basename(p).replace("Decay.mat", "").lower()
        ch = {"heave": (mo["pos_mm"][:, 2], "mm"), "pitch": (mo["rpy_deg"][:, 1], "deg"),
              "surge": (mo["pos_mm"][:, 0], "mm")}[mode]
        y = ch[0] - ch[0].mean()
        fs = mo["fs"]
        sos = signal.butter(4, [0.15, 1.5], "bandpass", fs=fs, output="sos")
        yf = signal.sosfiltfilt(sos, y)
        zc = np.where((yf[:-1] < 0) & (yf[1:] >= 0))[0]
        if len(zc) < 4:
            continue
        per = np.diff(t[zc])
        pk, _ = signal.find_peaks(yf, distance=int(0.3 * fs))
        amp = np.abs(yf[pk])
        tt = t[pk]
        good = amp > 0.15 * np.median(amp[:5])
        slope = np.polyfit(tt[good], np.log(amp[good]), 1)[0]
        T = float(np.median(per))
        zeta = float(-slope * T / (2 * np.pi))
        decays.append(dict(mode=mode, period_s=T, damping_ratio=zeta,
                           initial_amp=float(np.ptp(ch[0])), unit=ch[1], n=len(t)))
    with open(os.path.join(OUT, "decay_summary.csv"), "w", newline="", encoding="utf-8") as fh:
        w = csv.DictWriter(fh, fieldnames=list(decays[0].keys()))
        w.writeheader()
        w.writerows(decays)
    for d in decays:
        print("decay", d)
    np.savez(os.path.join(OUT, "decay.npz"), **{d["mode"]: np.array([d["period_s"], d["damping_ratio"]]) for d in decays})


if __name__ == "__main__":
    main()
