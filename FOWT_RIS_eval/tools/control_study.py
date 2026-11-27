"""Control-bandwidth, latency, prediction and noise study using the measured
6-DOF records and the linearised beam-centroid sensitivity matrix.

The point-adaptive controller places the beam on the receiver for the *estimated*
pose.  If the estimate differs from the true pose by e, the beam centroid moves
by J e to first order, where J is the static pose sensitivity of the design.
J is verified against exact nonlinear ray traces.
"""
from __future__ import annotations

import csv
import os

import numpy as np

import fx_optics as fo
import replay as rp

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
OUT = os.path.join(ROOT, "results")


def sensitivity(d, power, ref, state, eps_deg=0.05, eps_m=0.005):
    """2x6 Jacobian of the beam centroid wrt [roll,pitch,yaw,x,y,z].

    The phase command is frozen at its design value: this is the open-loop
    sensitivity that maps a pose *estimation* error into a beam displacement
    when the controller has inverted the estimated pose.
    """
    R0, p0 = np.eye(3), np.asarray(d["pRis"], float)
    qin0 = fo.q_in(d, R0, p0, ref)
    if "b" in state:
        qout0 = qin0 + state["b"]
    else:
        _, qout0, _ = rp.point_targets(d, R0, p0, ref, qin0,
                                       np.array([*state["receiver"], d["rxPlaneZ"]]))
    def centroid(rpy, dp):
        R = fo.euler_matrix(rpy)
        land, _ = fo.unit_landing(d, R, p0 + dp, ref, qout0)
        return np.sum(land[:, :2] * power[:, None], axis=0)
    c0 = centroid(np.zeros(3), np.zeros(3))
    J = np.zeros((2, 6))
    for j in range(3):
        rpy = np.zeros(3); rpy[j] = eps_deg
        J[:, j] = (centroid(rpy, np.zeros(3)) - c0) / np.deg2rad(eps_deg)
    for j in range(3):
        dp = np.zeros(3); dp[j] = eps_m
        J[:, 3 + j] = (centroid(np.zeros(3), dp) - c0) / eps_m
    return J, c0


def pose_series(rec, decimate=1):
    t = rec["t"][::decimate]
    pos = rec["pos"][::decimate].copy()
    rpy = rec["rpy"][::decimate].copy()
    rpy -= rpy.mean(axis=0)
    pos -= pos.mean(axis=0)
    return t, np.column_stack([np.deg2rad(rpy), pos])


def predict(pose, delay, order):
    """Estimate pose[k] from samples up to k-delay."""
    n, m = pose.shape
    out = np.empty_like(pose)
    out[:delay + 1] = pose[:delay + 1]
    for k in range(delay + 1, n):
        if order == 0:
            out[k] = pose[k - delay]
        elif order == 1:
            out[k] = pose[k - delay] + delay * (pose[k - delay] - pose[k - delay - 1])
        else:
            v1 = pose[k - delay] - pose[k - delay - 1]
            v2 = pose[k - delay - 1] - pose[k - delay - 2]
            out[k] = pose[k - delay] + delay * v1 + 0.5 * delay * (delay + 1) * (v1 - v2)
    return out


def metrics(t, off, p, pnom, dt):
    bad = p < pnom * 10 ** (-3 / 10)
    e = np.diff(np.concatenate([[0], bad.astype(int), [0]]))
    starts, stops = np.where(e == 1)[0], np.where(e == -1)[0]
    return dict(offset_rms_mm=float(np.sqrt(np.mean(off ** 2)) * 1e3),
                offset_p95_mm=float(np.percentile(off, 95) * 1e3),
                p_mean_mw=float(np.mean(p) * 1e3),
                p_p10_mw=float(np.percentile(p, 10) * 1e3),
                loss_mean_db=float(10 * np.log10(np.mean(p) / pnom)),
                avail_3db=float(np.mean(p >= pnom * 10 ** (-3 / 10))),
                outage_frac=float(np.mean(bad)),
                outage_mean_ms=float(np.mean(stops - starts) * dt * 1e3) if len(starts) else 0.0)


def main():
    d = fo.load_design()
    ref = np.asarray(d["referenceXY"], float)
    P = np.asarray(d["unitDirectPowerW"], float); P /= P.sum()
    states = rp.design_setup(d, P, ref)
    ap = rp.aperture_map(d, "improved", d["improvedCenter"], 0.05)
    pnom = rp.aperture_at(ap, 0, 0)
    J, c0 = sensitivity(d, P, ref, states["point"])
    print("J (m per rad / m):\n", J)
    # exact validation of the open-loop sensitivity at +-0.4 deg roll / pitch
    R0, p0 = np.eye(3), np.asarray(d["pRis"], float)
    qin0 = fo.q_in(d, R0, p0, ref)
    _, qout0, _ = rp.point_targets(d, R0, p0, ref, qin0,
                                   np.array([*states["point"]["receiver"], d["rxPlaneZ"]]))
    for j, nm in [(0, "roll"), (1, "pitch")]:
        for a in (-0.4, -0.2, 0.2, 0.4):
            rpy = np.zeros(3); rpy[j] = a
            R = fo.euler_matrix(rpy)
            land, _ = fo.unit_landing(d, R, p0, ref, qout0)
            c = np.sum(land[:, :2] * P[:, None], axis=0)
            print(f"{nm} {a:+.1f} deg exact {np.linalg.norm(c-c0)*1e3:.3f} mm linear {np.linalg.norm(J[:, j]*np.deg2rad(a))*1e3:.3f} mm")
    recs = rp.load_motion()
    sel = [0, 5, 6, 8]  # 10mm 3min, 30mm 1.2, 30mm 0.4, 45mm 0.8
    rows = []
    for idx in sel:
        rec = recs[idx]
        for dec, dt in [(1, 0.01)]:
            t, pose = pose_series(rec, dec)
            for order, pname in [(0, "ZOH"), (1, "linear"), (2, "quadratic")]:
                for delay_ms in [0, 10, 20, 50, 100]:
                    delay = int(round(delay_ms / 1000 / dt))
                    if delay == 0:
                        est = pose
                    else:
                        est = predict(pose, delay, order)
                    e = pose - est
                    off = np.linalg.norm(e @ J.T, axis=1)
                    vec = e @ J.T
                    p = np.array([rp.aperture_at(ap, vec[k, 0], vec[k, 1]) for k in range(len(t))])
                    r = metrics(t, off, p, pnom, dt)
                    r.update(record=rec["name"], predictor=pname, delay_ms=delay_ms)
                    rows.append(r)
                    print(f"{rec['name'][:32]:34s} {pname:9s} {delay_ms:4d} ms  "
                          f"rms {r['offset_rms_mm']:6.2f} mm p95 {r['offset_p95_mm']:6.2f} "
                          f"loss {r['loss_mean_db']:6.2f} dB avail {r['avail_3db']:5.2f}")
        # update-rate: sample and hold every n samples, ZOH
        for hold_ms in [10, 20, 50, 100, 200]:
            hold = int(round(hold_ms / 1000 / dt))
            est = np.repeat(pose[::hold], hold, axis=0)[:len(pose)]
            e = pose - est
            vec = e @ J.T
            off = np.linalg.norm(vec, axis=1)
            p = np.array([rp.aperture_at(ap, vec[k, 0], vec[k, 1]) for k in range(len(t))])
            r = metrics(t, off, p, pnom, dt)
            r.update(record=rec["name"], predictor="hold", delay_ms=hold_ms)
            rows.append(r)
            print(f"{rec['name'][:32]:34s} hold     {hold_ms:4d} ms  rms {r['offset_rms_mm']:6.2f} p95 {r['offset_p95_mm']:6.2f} loss {r['loss_mean_db']:6.2f} avail {r['avail_3db']:5.2f}")
        # measurement noise on the estimated pose
        rng = np.random.default_rng(1234)
        for noise_deg in [0.01, 0.05, 0.1, 0.2]:
            for noise_mm in [1.0]:
                est = pose.copy()
                est[:, :3] += np.deg2rad(noise_deg) * rng.standard_normal((len(pose), 3))
                est[:, 3:] += (noise_mm * 1e-3) * rng.standard_normal((len(pose), 3))
                e = pose - est
                vec = e @ J.T
                off = np.linalg.norm(vec, axis=1)
                p = np.array([rp.aperture_at(ap, vec[k, 0], vec[k, 1]) for k in range(len(t))])
                r = metrics(t, off, p, pnom, dt)
                r.update(record=rec["name"], predictor=f"noise{noise_deg}deg", delay_ms=0)
                rows.append(r)
                print(f"{rec['name'][:32]:34s} noise {noise_deg:5.2f} deg rms {r['offset_rms_mm']:6.2f} p95 {r['offset_p95_mm']:6.2f} loss {r['loss_mean_db']:6.2f} avail {r['avail_3db']:5.2f}")
    with open(os.path.join(OUT, "control_study.csv"), "w", newline="", encoding="utf-8") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0].keys()))
        w.writeheader(); w.writerows(rows)
    np.save(os.path.join(OUT, "sensitivity.npy"), J)


if __name__ == "__main__":
    main()
