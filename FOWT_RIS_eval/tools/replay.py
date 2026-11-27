"""Experimental 6-DOF motion replay through the reduced RIS link model.

Scenarios
  current_static   : production +z-closest phase frozen at the design pose
  current_point    : production phase re-aimed at its own receiver every frame
  collimated_static: power-weighted common-direction phase frozen
  collimated_adapt : common world direction updated every frame
  point_adapt      : collimated bias + per-unit focus on the fixed receiver
"""
from __future__ import annotations

import csv
import os

import numpy as np
from scipy import ndimage, signal

import fx_optics as fo

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
OUT = os.path.join(ROOT, "results")
os.makedirs(OUT, exist_ok=True)


def load_motion():
    z = np.load(os.path.join(ROOT, "data", "processed", "fowt_motion.npz"), allow_pickle=True)
    n = int(z["n"])
    return [dict(index=i, name=str(z["meta"][i]), t=z[f"t_{i}"],
                 pos=z[f"pos_{i}"], rpy=z[f"rpy_{i}"]) for i in range(n)]


def aperture_map(d, key, center, radius):
    x = np.asarray(d[key + "X"], float)
    y = np.asarray(d[key + "Y"], float)
    I = np.asarray(d[key + "I"], float)
    dx, dy = x[1] - x[0], y[1] - y[0]
    nx, ny = len(x), len(y)
    yy, xx = np.mgrid[0:ny, 0:nx]
    cx = (center[0] - x[0]) / dx
    cy = (center[1] - y[0]) / dy
    mask = ((xx - cx) * dx) ** 2 + ((yy - cy) * dy) ** 2 <= radius ** 2
    corr = signal.fftconvolve(I, mask[::-1, ::-1], mode="same") * dx * dy
    return dict(corr=corr, dx=dx, dy=dy, nx=nx, ny=ny, cx=cx, cy=cy, radius=radius)


def aperture_at(ap, ox, oy):
    ix = ap["cx"] + ox / ap["dx"]
    iy = ap["cy"] + oy / ap["dy"]
    if ix < 0 or iy < 0 or ix > ap["nx"] - 1 or iy > ap["ny"] - 1:
        return 0.0
    return float(ndimage.map_coordinates(ap["corr"], [[iy], [ix]], order=1, mode="nearest")[0])


def point_targets(d, R, p, ref, qin, p_rx, iterations=2):
    q = qin.copy()
    cmd = None
    for _ in range(iterations):
        top = fo.unit_top(d, R, p, ref, q)
        dvec = np.asarray(p_rx, float)[None, :] - top
        dvec /= np.maximum(np.linalg.norm(dvec, axis=1, keepdims=True), 1e-12)
        qt = float(d["nAir"]) * (dvec @ R)[:, :2]
        cmd, q, valid = fo.snap_commands(d, qin, qt)
    return cmd, q, valid


def design_setup(d, power, ref, coverage=0.999):
    R = np.eye(3)
    p = np.asarray(d["pRis"], float)
    qin = fo.q_in(d, R, p, ref)
    qs, cov = fo.power_weighted_target(d, qin, power, coverage)
    u = np.array([qs[0], qs[1], np.sqrt(1 - qs @ qs)])
    B = float(d["B"])
    states = {}
    # current frozen corner
    qout = qin + np.array([-B, B])
    land, _ = fo.unit_landing(d, R, p, ref, qout)
    rc = np.sum(land[:, :2] * power[:, None], axis=0)
    states["current"] = dict(b=np.array([-B, B]), c0=rc, receiver=d["currentCenter"])
    # point-aimed production at its own receiver
    cmd, qout, _ = point_targets(d, R, p, ref, qin, np.array([*d["currentCenter"], d["rxPlaneZ"]]))
    land, _ = fo.unit_landing(d, R, p, ref, qout)
    states["current_point"] = dict(cmd0=cmd, c0=np.sum(land[:, :2] * power[:, None], axis=0),
                                   receiver=d["currentCenter"])
    # collimated
    qt = float(d["nAir"]) * (R.T @ u)[:2]
    cmd, qout, _ = fo.snap_commands(d, qin, np.tile(qt, (len(qin), 1)))
    land, _ = fo.unit_landing(d, R, p, ref, qout)
    states["collimated"] = dict(b=qt - qin, c0=np.sum(land[:, :2] * power[:, None], axis=0),
                                receiver=d["improvedCenter"], u=u, qs=qs)
    # point-aimed collimated
    cmd, qout, _ = point_targets(d, R, p, ref, qin, np.array([*d["improvedCenter"], d["rxPlaneZ"]]))
    land, _ = fo.unit_landing(d, R, p, ref, qout)
    states["point"] = dict(cmd0=cmd, c0=np.sum(land[:, :2] * power[:, None], axis=0),
                           receiver=d["improvedCenter"])
    return states


def run_record(d, rec, power, ref, states, ap_cur, ap_imp, decimate=2, trim=True):
    t = rec["t"][::decimate]
    pos = rec["pos"][::decimate].copy()
    rpy = rec["rpy"][::decimate].copy()
    if trim:
        rpy -= rpy.mean(axis=0)
    pos -= pos.mean(axis=0)
    R0, p0 = np.eye(3), np.asarray(d["pRis"], float)
    B = float(d["B"])
    u = states["collimated"]["u"]
    qs = np.asarray(states["collimated"]["qs"], float)
    dstar = np.array([qs[0], qs[1], np.sqrt(1 - qs @ qs)])
    qin0 = fo.q_in(d, R0, p0, ref)
    _, dw_c = fo.unit_landing(d, R0, p0, ref, qin0 + np.array([-B, B]))
    dstar_current = np.sum(dw_c * power[:, None], axis=0)
    dstar_current /= np.linalg.norm(dstar_current)
    keys = ["current_static", "current_point", "collimated_static", "collimated_adaptive", "point_adaptive"]
    out = {k: dict(p=np.zeros(len(t)), off=np.zeros(len(t)), ang=np.zeros(len(t)),
                   offv=np.zeros((len(t), 2))) for k in keys}
    # nominal directions
    _, dw0 = fo.unit_landing(d, R0, p0, ref, fo.q_in(d, R0, p0, ref) + np.array([-B, B]))
    d0_current = np.sum(dw0 * power[:, None], axis=0)
    d0_current /= np.linalg.norm(d0_current)
    for k in range(len(t)):
        R = fo.euler_matrix(rpy[k])
        p = p0 + pos[k]
        qin = fo.q_in(d, R, p, ref)
        def metrics(qout, key, c_ref, ap, dref):
            land, dw = fo.unit_landing(d, R, p, ref, qout)
            c = np.sum(land[:, :2] * power[:, None], axis=0)
            mdir = np.sum(dw * power[:, None], axis=0)
            mdir /= max(np.linalg.norm(mdir), 1e-12)
            out[key]["off"][k] = np.linalg.norm(c - c_ref)
            out[key]["offv"][k] = c - c_ref
            out[key]["ang"][k] = np.degrees(np.arccos(np.clip(mdir @ dref, 0, 1))) * 1e3
            out[key]["p"][k] = aperture_at(ap, *(c - c_ref))
        # static current
        metrics(qin + np.array([-B, B]), "current_static",
                states["current"]["c0"], ap_cur, dstar_current)
        # current point-aimed
        _, qo, _ = point_targets(d, R, p, ref, qin,
                                 np.array([*states["current_point"]["receiver"], d["rxPlaneZ"]]))
        metrics(qo, "current_point", states["current_point"]["c0"], ap_cur, dstar_current)
        # collimated static
        metrics(qin + states["collimated"]["b"], "collimated_static",
                states["collimated"]["c0"], ap_imp, dstar)
        # adaptive common world direction
        qt = float(d["nAir"]) * (R.T @ u)[:2]
        _, qo, _ = fo.snap_commands(d, qin, np.tile(qt, (len(qin), 1)))
        metrics(qo, "collimated_adaptive", states["collimated"]["c0"], ap_imp, dstar)
        # point-aimed collimated
        _, qo, _ = point_targets(d, R, p, ref, qin,
                                 np.array([*states["point"]["receiver"], d["rxPlaneZ"]]))
        metrics(qo, "point_adaptive", states["point"]["c0"], ap_imp, dstar)
    out["t"] = t
    out["rpy"] = rpy
    out["pos"] = pos
    return out


def summary_row(rec, scenario, p, off, ang, pnom):
    p = p * 1e3
    r = dict(record=rec["name"], scenario=scenario, p_nom_mw=pnom * 1e3,
             p_mean_mw=float(np.mean(p)), p_min_mw=float(np.min(p)),
             p_p10_mw=float(np.percentile(p, 10)), p_p90_mw=float(np.percentile(p, 90)),
             loss_mean_db=float(10 * np.log10(max(np.mean(p), 1e-12) / max(pnom * 1e3, 1e-12))),
             loss_p10_db=float(10 * np.log10(max(np.percentile(p, 10), 1e-12) / max(pnom * 1e3, 1e-12))),
             avail_1db=float(np.mean(p >= pnom * 1e3 * 10 ** (-1 / 10))),
             avail_3db=float(np.mean(p >= pnom * 1e3 * 10 ** (-3 / 10))),
             avail_10db=float(np.mean(p >= pnom * 1e3 * 10 ** (-10 / 10))),
             offset_rms_mm=float(np.sqrt(np.mean(off ** 2)) * 1e3),
             offset_p95_mm=float(np.percentile(off, 95) * 1e3),
             angle_rms_mrad=float(np.sqrt(np.mean(ang ** 2))))
    bad = p < pnom * 1e3 * 10 ** (-3 / 10)
    edges = np.diff(np.concatenate([[0], bad.astype(int), [0]]))
    starts, stops = np.where(edges == 1)[0], np.where(edges == -1)[0]
    r["outage_3db_fraction"] = float(np.mean(bad))
    r["outage_3db_mean_s"] = float(np.mean((stops - starts) * 0.02)) if len(starts) else 0.0
    r["outage_3db_max_s"] = float(np.max((stops - starts) * 0.02)) if len(starts) else 0.0
    return r


def main():
    d = fo.load_design()
    ref = np.asarray(d["referenceXY"], float)
    P = np.asarray(d["unitDirectPowerW"], float)
    P /= P.sum()
    states = design_setup(d, P, ref)
    print("collimated q*", states["collimated"]["qs"], "u", states["collimated"]["u"])
    ap_cur = aperture_map(d, "current", d["currentCenter"], 0.05)
    ap_imp = aperture_map(d, "improved", d["improvedCenter"], 0.05)
    pnom_cur = aperture_at(ap_cur, 0, 0)
    pnom_imp = aperture_at(ap_imp, 0, 0)
    print("nominal mW: current %.3f improved %.3f" % (pnom_cur * 1e3, pnom_imp * 1e3))
    recs = load_motion()
    rows = []
    for rec in recs:
        out = run_record(d, rec, P, ref, states, ap_cur, ap_imp, decimate=2)
        mapping = [("current_static", "current static", pnom_cur),
                   ("current_point", "current point-adaptive", pnom_cur),
                   ("collimated_static", "collimated static", pnom_imp),
                   ("collimated_adaptive", "collimated adaptive", pnom_imp),
                   ("point_adaptive", "point adaptive", pnom_imp)]
        for key, label, pn in mapping:
            r = summary_row(rec, label, out[key]["p"], out[key]["off"], out[key]["ang"], pn)
            rows.append(r)
            print(f"{rec['name'][:34]:36s} {label:24s} pmean {r['p_mean_mw']:7.3f} p10 {r['p_p10_mw']:7.3f} "
                  f"loss {r['loss_mean_db']:6.2f} avail3 {r['avail_3db']:5.2f} off {r['offset_rms_mm']:6.1f}mm "
                  f"ang {r['angle_rms_mrad']:6.2f}mrad")
        np.savez_compressed(os.path.join(OUT, f"replay_{rec['index']:02d}.npz"),
                            **{f"{k}_{f}": out[k][f] for k in out if isinstance(out[k], dict) for f in out[k]},
                            t=out["t"], rpy=out["rpy"], pos=out["pos"])
    with open(os.path.join(OUT, "replay_summary.csv"), "w", newline="", encoding="utf-8") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0].keys()))
        w.writeheader()
        w.writerows(rows)


if __name__ == "__main__":
    main()
