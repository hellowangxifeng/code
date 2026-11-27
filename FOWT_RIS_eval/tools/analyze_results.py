"""Aggregate the replay results: aperture trade-off, per-condition tables, figures data."""
from __future__ import annotations

import csv
import os

import numpy as np

import fx_optics as fo
import replay as rp

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
OUT = os.path.join(ROOT, "results")


def per_radius():
    d = fo.load_design()
    maps = {tag: {r: rp.aperture_map(d, tag, d[tag + "Center"], r) for r in (0.025, 0.05, 0.1)}
            for tag in ("current", "improved")}
    rows = []
    for idx in range(9):
        z = np.load(os.path.join(OUT, f"replay_{idx:02d}.npz"))
        for key, tag in [("current_static", "current"), ("current_point", "current"),
                         ("collimated_static", "improved"), ("collimated_adaptive", "improved"),
                         ("point_adaptive", "improved")]:
            off = z[key + "_offv"]
            for r in (0.025, 0.05, 0.1):
                ap = maps[tag][r]
                p = np.array([rp.aperture_at(ap, *o) for o in off])
                pnom = rp.aperture_at(ap, 0, 0)
                rows.append(dict(record=str(z["t"].shape), index=idx, scenario=key, radius_m=r,
                                 p_nom_mw=pnom * 1e3, p_mean_mw=float(np.mean(p) * 1e3),
                                 p_p10_mw=float(np.percentile(p, 10) * 1e3),
                                 loss_mean_db=float(10 * np.log10(np.mean(p) / pnom)),
                                 loss_p10_db=float(10 * np.log10(np.percentile(p, 10) / pnom)),
                                 avail_3db=float(np.mean(p >= pnom * 10 ** (-3 / 10)))))
    with open(os.path.join(OUT, "aperture_tradeoff.csv"), "w", newline="", encoding="utf-8") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0].keys()))
        w.writeheader(); w.writerows(rows)
    return rows


def aggregate():
    rows = list(csv.DictReader(open(os.path.join(OUT, "replay_summary.csv"), encoding="utf-8")))
    out = []
    for sc in ["current static", "current point-adaptive", "collimated static",
               "collimated adaptive", "point adaptive"]:
        rr = [r for r in rows if r["scenario"] == sc]
        for f in ["p_nom_mw", "p_mean_mw", "p_p10_mw", "p_min_mw", "loss_mean_db", "loss_p10_db",
                  "avail_3db", "offset_rms_mm", "angle_rms_mrad"]:
            v = np.array([float(r[f]) for r in rr])
            if f == "loss_mean_db":
                v = np.array([10 * np.log10(np.mean([float(r["p_mean_mw"]) for r in rr]) /
                                            np.mean([float(r["p_nom_mw"]) for r in rr]))])
            if f == "loss_p10_db":
                v = np.array([10 * np.log10(np.mean([float(r["p_p10_mw"]) for r in rr]) /
                                            np.mean([float(r["p_nom_mw"]) for r in rr]))])
            out.append(dict(scenario=sc, field=f, mean=float(np.mean(v)), min=float(np.min(v)),
                            max=float(np.max(v)), worst=float(np.min(v)) if "avail" not in f and "p_" in f else float(np.max(v))))
    with open(os.path.join(OUT, "aggregate.csv"), "w", newline="", encoding="utf-8") as fh:
        w = csv.DictWriter(fh, fieldnames=list(out[0].keys()))
        w.writeheader(); w.writerows(out)
    for sc in ["current static", "current point-adaptive", "collimated static",
               "collimated adaptive", "point adaptive"]:
        rr = [r for r in rows if r["scenario"] == sc]
        p = np.array([float(r["p_mean_mw"]) for r in rr])
        p10 = np.array([float(r["p_p10_mw"]) for r in rr])
        off = np.array([float(r["offset_rms_mm"]) for r in rr])
        av = np.array([float(r["avail_3db"]) for r in rr])
        print(f"{sc:24s} mean {p.mean():6.2f} mW [{p.min():5.2f},{p.max():5.2f}]  p10 {p10.mean():6.2f}  off {off.mean():6.1f} mm  avail3 min {av.min():.2f}")


if __name__ == "__main__":
    print("== aperture trade-off ==")
    rows = per_radius()
    for sc in ["current_static", "current_point", "collimated_static", "collimated_adaptive", "point_adaptive"]:
        for r in (0.025, 0.05, 0.1):
            rr = [x for x in rows if x["scenario"] == sc and x["radius_m"] == r]
            print(f"{sc:20s} r={r:5.3f} pnom {np.mean([x['p_nom_mw'] for x in rr]):6.2f} "
                  f"pmean {np.mean([x['p_mean_mw'] for x in rr]):6.2f} "
                  f"loss {np.mean([x['loss_mean_db'] for x in rr]):6.2f} "
                  f"avail3 {np.mean([x['avail_3db'] for x in rr]):.2f}")
    print("== aggregate ==")
    aggregate()
