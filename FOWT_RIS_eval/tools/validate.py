"""End-to-end verification of the FOWT-RIS evaluation.

Run from FOWT_RIS_eval/tools:  python validate.py
Prints PASS/FAIL for every quantitative claim used in the paper and exits
non-zero on the first failure.
"""
from __future__ import annotations

import csv
import os
import sys

import numpy as np
import scipy.io as sio

import fx_optics as fo
import replay as rp

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
RES = os.path.join(ROOT, "results")
DATA = os.path.join(ROOT, "data", "processed")

FAIL = []


def check(name, condition, detail=""):
    ok = bool(condition)
    print(("PASS  " if ok else "FAIL  ") + name + ("  " + detail if detail else ""))
    if not ok:
        FAIL.append(name)


def main():
    # 1. dataset
    z = np.load(os.path.join(DATA, "fowt_motion.npz"), allow_pickle=True)
    n = int(z["n"])
    check("dataset has 9 regular-wave records", n == 9, f"n={n}")
    cond = list(csv.DictReader(open(os.path.join(DATA, "conditions.csv"), encoding="utf-8")))
    check("condition table has 9 rows", len(cond) == 9)
    for i in range(n):
        t = z[f"t_{i}"]
        check(f"record {i} uniform 100 Hz and finite",
              abs(np.median(np.diff(t)) - 0.01) < 1e-6 and np.isfinite(z[f"rpy_{i}"]).all())
    tilts = np.array([float(c["tilt_rms_deg"]) for c in cond])
    check("tilt rms in 0.14-1.10 deg", tilts.min() > 0.13 and tilts.max() < 1.11,
          f"{tilts.min():.3f}-{tilts.max():.3f} deg")
    decay = list(csv.DictReader(open(os.path.join(DATA, "decay_summary.csv"), encoding="utf-8")))
    check("decay table has 3 modes", len(decay) == 3)
    check("heave decay ~2.16 s", abs(float(decay[0]["period_s"]) - 2.158) < 0.01)

    # 2. design
    d = fo.load_design()
    ref = np.asarray(d["referenceXY"], float)
    P = np.asarray(d["unitDirectPowerW"], float); P /= P.sum()
    check("1600 units", len(ref) == 1600)
    cmd = np.asarray(d["currentCommandDeg"], float)
    check("production commands all at the corner",
          np.allclose(cmd, cmd[0]) and abs(cmd[0, 0] + 15.414) < 1e-6 and abs(cmd[0, 1] - 15.414) < 1e-6,
          f"command={cmd[0]}")
    qin0 = fo.q_in(d, np.eye(3), np.asarray(d["pRis"], float), ref)
    qs, cov = fo.power_weighted_target(d, qin0, P, 0.999)
    check("power-weighted target covers >=99.9%", cov >= 0.999, f"q*={qs}, coverage={cov:.4f}")
    check("target angle is 53.9 +/- 0.5 deg",
          abs(np.degrees(np.arcsin(np.hypot(*qs))) - 53.9) < 0.5,
          f"{np.degrees(np.arcsin(np.hypot(*qs))):.2f} deg")

    # 3. reduced geometry matches MATLAB anchors
    R = np.asarray(d["RRis"], float); p = np.asarray(d["pRis"], float)
    qin = fo.q_in(d, R, p, ref)
    qout = qin + np.asarray(d["currentGradient"], float)
    land, _ = fo.unit_landing(d, R, p, ref, qout)
    c_cur = np.sum(land[:, :2] * P[:, None], axis=0)
    check("current centroid matches MATLAB", np.linalg.norm(c_cur - np.asarray(d["currentCenter"])) < 2e-3,
          f"Python {c_cur}, MATLAB {np.asarray(d['currentCenter'])}")
    u = np.asarray(d["improvedU"], float)
    qt = float(d["nAir"]) * (R.T @ u)[:2]
    _, qout2, _ = fo.snap_commands(d, qin, np.tile(qt, (len(qin), 1)))
    land2, _ = fo.unit_landing(d, R, p, ref, qout2)
    c_imp = np.sum(land2[:, :2] * P[:, None], axis=0)
    check("collimated centroid matches MATLAB", np.linalg.norm(c_imp - np.asarray(d["improvedCenter"])) < 2e-3,
          f"Python {c_imp}, MATLAB {np.asarray(d['improvedCenter'])}")

    # 4. aperture
    apc = rp.aperture_map(d, "current", d["currentCenter"], 0.05)
    api = rp.aperture_map(d, "improved", d["improvedCenter"], 0.05)
    pc = rp.aperture_at(apc, 0, 0) * 1e3
    pi = rp.aperture_at(api, 0, 0) * 1e3
    check("production 5 cm nominal 4.5-4.7 mW", 4.5 < pc < 4.7, f"{pc:.3f} mW")
    check("collimated 5 cm nominal 14.5-15.3 mW", 14.5 < pi < 15.3, f"{pi:.3f} mW")
    check("nominal gain >= 5 dB", 10 * np.log10(pi / pc) >= 5.0, f"{10*np.log10(pi/pc):.2f} dB")

    # 5. replay aggregates
    rows = list(csv.DictReader(open(os.path.join(RES, "replay_summary.csv"), encoding="utf-8")))
    check("replay has 45 rows", len(rows) == 45)
    def agg(sc):
        rr = [r for r in rows if r["scenario"] == sc]
        return (np.mean([float(r["p_mean_mw"]) for r in rr]),
                np.mean([float(r["p_p10_mw"]) for r in rr]),
                np.min([float(r["avail_3db"]) for r in rr]),
                np.mean([float(r["offset_rms_mm"]) for r in rr]))
    a_cur = agg("current static")
    a_pt = agg("point adaptive")
    a_cs = agg("collimated static")
    check("production mean 4.4-4.7 mW", 4.4 < a_cur[0] < 4.7, f"{a_cur[0]:.2f}")
    check("point-adaptive mean 14.7-15.0 mW", 14.7 < a_pt[0] < 15.0, f"{a_pt[0]:.2f}")
    check("point-adaptive availability 1.0", a_pt[2] > 0.999)
    check("frozen collimated worst availability < 0.7", a_cs[2] < 0.7, f"{a_cs[2]:.2f}")
    check("point-adaptive centroid < 5 mm", a_pt[3] < 5.0, f"{a_pt[3]:.2f} mm")
    gain = 10 * np.log10(a_pt[0] / a_cur[0])
    check("total mean gain 5.0-5.4 dB", 5.0 < gain < 5.4, f"{gain:.2f} dB")
    gain10 = 10 * np.log10(a_pt[1] / a_cur[1])
    check("total p10 gain >= 5.3 dB", gain10 >= 5.3, f"{gain10:.2f} dB")

    # 6. control study
    ctl = list(csv.DictReader(open(os.path.join(RES, "control_study.csv"), encoding="utf-8")))
    target = "FOWTmotion_11_3_45mmp8Hz.mat"
    def get(pred, delay):
        for r in ctl:
            if r["record"] == target and r["predictor"] == pred and abs(float(r["delay_ms"]) - delay) < 1e-6:
                return float(r["offset_rms_mm"]), float(r["loss_mean_db"])
        raise KeyError((pred, delay))
    zoh = get("ZOH", 100)
    lin = get("linear", 100)
    check("first-order predictor halves 100 ms error", lin[0] < 0.6 * zoh[0],
          f"ZOH {zoh[0]:.1f} mm, linear {lin[0]:.1f} mm")
    check("linear predictor loss at 100 ms > -1 dB", lin[1] > -1.0, f"{lin[1]:.2f} dB")
    noise = get("noise0.05deg", 0)
    check("0.05 deg noise loss < 0.1 dB", noise[1] > -0.1, f"{noise[1]:.3f} dB")

    # 7. model bridge
    M = sio.loadmat(os.path.join(RES, "measured_wave_ris_8.mat"), squeeze_me=True)
    o = M["o"]
    def fld(name):
        v = o[name]
        while isinstance(v, np.ndarray) and v.dtype == object and v.ndim == 0:
            v = v.item()
        return np.asarray(v, float)
    pitch = fld("pitchDeg") - np.mean(fld("pitchDeg"))
    check("RIS model pitch > 3 deg under 45 mm/0.8 Hz", np.std(pitch) > 3.0, f"{np.std(pitch):.2f} deg")
    mrows = list(csv.DictReader(open(os.path.join(RES, "model_motion_replay.csv"), encoding="utf-8")))
    ms = [r for r in mrows if r["record"] == "RISmodel_8" and r["scenario"] == "collimated static"][0]
    mp = [r for r in mrows if r["record"] == "RISmodel_8" and r["scenario"] == "point adaptive"][0]
    check("model replay static availability < 0.5", float(ms["avail_3db"]) < 0.5, ms["avail_3db"])
    check("model replay point > 14 mW", float(mp["p_mean_mw"]) > 14.0, mp["p_mean_mw"])

    print()
    if FAIL:
        print(f"{len(FAIL)} CHECK(S) FAILED:")
        for f in FAIL:
            print("  -", f)
        sys.exit(1)
    print("ALL CHECKS PASSED")


if __name__ == "__main__":
    main()
