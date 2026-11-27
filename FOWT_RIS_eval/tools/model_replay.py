"""Replay the RIS-platform motion produced by sea_wave under the measured wave
record (30 mm and 45 mm at 0.8 Hz) through the optical chain."""
from __future__ import annotations

import csv
import os

import numpy as np
import scipy.io as sio

import fx_optics as fo
import replay as rp

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
RES = os.path.join(ROOT, "results")


def main():
    d = fo.load_design()
    ref = np.asarray(d["referenceXY"], float)
    P = np.asarray(d["unitDirectPowerW"], float); P /= P.sum()
    states = rp.design_setup(d, P, ref)
    ap_cur = rp.aperture_map(d, "current", d["currentCenter"], 0.05)
    ap_imp = rp.aperture_map(d, "improved", d["improvedCenter"], 0.05)
    pnom_cur = rp.aperture_at(ap_cur, 0, 0)
    pnom_imp = rp.aperture_at(ap_imp, 0, 0)
    rows = []
    for idx in (6, 8):
        M = sio.loadmat(os.path.join(RES, f"measured_wave_ris_{idx}.mat"), squeeze_me=True)
        o = M["o"]
        def fld(name):
            v = o[name]
            while isinstance(v, np.ndarray) and v.dtype == object and v.ndim == 0:
                v = v.item()
            return np.asarray(v, float)
        rec = dict(index=idx, name=f"RIS-model, measured wave {idx}",
                   t=fld("t"), pos=np.column_stack([fld("surge"), fld("sway"), fld("heave")]),
                   rpy=np.column_stack([fld("rollDeg"), fld("pitchDeg"), fld("yawDeg")]))
        out = rp.run_record(d, rec, P, ref, states, ap_cur, ap_imp, decimate=1, trim=True)
        np.savez_compressed(os.path.join(RES, f"model_replay_{idx}.npz"),
                            **{f"{k}_{f}": out[k][f] for k in out if isinstance(out[k], dict) for f in out[k]},
                            t=out["t"], rpy=out["rpy"], pos=out["pos"])
        for key, label, pn in [("current_static", "current static", pnom_cur),
                               ("collimated_static", "collimated static", pnom_imp),
                               ("collimated_adaptive", "collimated adaptive", pnom_imp),
                               ("point_adaptive", "point adaptive", pnom_imp)]:
            r = rp.summary_row(rec, label, out[key]["p"], out[key]["off"], out[key]["ang"], pn)
            r["record"] = f"RISmodel_{idx}"
            rows.append(r)
            print(f"{rec['name']:28s} {label:20s} pmean {r['p_mean_mw']:7.2f} p10 {r['p_p10_mw']:7.2f} "
                  f"loss {r['loss_mean_db']:6.2f} avail3 {r['avail_3db']:.2f} off {r['offset_rms_mm']:6.1f} mm")
    with open(os.path.join(RES, "model_motion_replay.csv"), "w", newline="", encoding="utf-8") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0].keys()))
        w.writeheader(); w.writerows(rows)


if __name__ == "__main__":
    main()
