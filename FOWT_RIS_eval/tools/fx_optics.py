"""Reduced-order optical/link model for the floating RIS replay.

The geometry reproduces RIS_phase exactly at the unit-reference level:
  layered Fermat path (water -> glass) -> discrete per-axis command ->
  glass refraction -> world direction -> receiver-plane landing.
The full-wave anchor maps are used only for the aperture integral, shifted by
the change of the power-weighted ray centroid.
"""
from __future__ import annotations

import os

import numpy as np
import scipy.io as sio

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)


def load_design():
    d = sio.loadmat(os.path.join(ROOT, "data", "processed", "nominal_design.mat"),
                    squeeze_me=True, struct_as_record=False)
    out = {}
    for k, v in d.items():
        if k.startswith("__"):
            continue
        out[k] = np.asarray(v).squeeze() if isinstance(v, np.ndarray) else v
    # scipy turns some structs into object arrays; keep numeric arrays numeric
    return out


def layer_path(delta_xy, heights, indices, iterations=48):
    """Exact parallel-layer Fermat path; vectorised over rows.

    delta_xy : (N,2) lateral displacement from source to target
    heights  : (2,) layer thicknesses, or (N,2)
    indices  : (2,) refractive indices
    """
    delta_xy = np.atleast_2d(delta_xy)
    heights = np.asarray(heights, float)
    if heights.ndim == 1:
        heights = np.broadcast_to(heights, (delta_xy.shape[0], heights.size))
    rho = np.linalg.norm(delta_xy, axis=1)
    unit = delta_xy / np.maximum(rho[:, None], np.finfo(float).tiny)
    lo = np.zeros_like(rho)
    hi = np.full_like(rho, np.min(indices) * (1 - 1e-13))
    idx = np.asarray(indices, float)
    for _ in range(iterations):
        q = 0.5 * (lo + hi)
        kz = np.sqrt(np.maximum(idx[None, :] ** 2 - q[:, None] ** 2, np.finfo(float).tiny))
        dist = np.sum(heights * q[:, None] / kz, axis=1)
        lower = dist < rho
        lo[lower] = q[lower]
        hi[~lower] = q[~lower]
    q = 0.5 * (lo + hi)
    q[rho == 0] = 0.0
    kz = np.sqrt(np.maximum(idx[None, :] ** 2 - q[:, None] ** 2, np.finfo(float).tiny))
    return dict(q=unit * q[:, None], kz=kz,
                length=heights * idx[None, :] / kz,
                opl=np.sum(heights * idx[None, :] ** 2 / kz, axis=1))


def rot_x(a):
    c, s = np.cos(a), np.sin(a)
    return np.array([[1, 0, 0], [0, c, -s], [0, s, c]])


def rot_y(a):
    c, s = np.cos(a), np.sin(a)
    return np.array([[c, 0, s], [0, 1, 0], [-s, 0, c]])


def rot_z(a):
    c, s = np.cos(a), np.sin(a)
    return np.array([[c, -s, 0], [s, c, 0], [0, 0, 1]])


def euler_matrix(rpy_deg):
    """ZYX: Rz(yaw) Ry(pitch) Rx(roll); input degrees [roll, pitch, yaw]."""
    r, p, y = np.deg2rad(rpy_deg)
    return rot_z(y) @ rot_y(p) @ rot_x(r)


def q_in(d, R, p, ref=None):
    """Tangential wave vector at the unit reference points (body frame)."""
    if ref is None:
        ref = d["referenceXY"]
    src_local = R.T @ (np.asarray(d["pTx"], float) - p)
    delta = ref - src_local[None, :2]
    path = layer_path(delta, [-src_local[2], d["hIn"]],
                      [d["nWater"], d["nGlass"]], 48)
    return path["q"]


def snap_commands(d, qin, qtarget, repair=True):
    nAir = d["nAir"]
    step = d["angleStepDeg"]
    maxcmd = d["maxCmdDeg"]
    b = qtarget - qin
    cmd = np.round(np.degrees(np.arcsin(np.clip(b / nAir, -1, 1))) / step) * step
    cmd = np.clip(cmd, -maxcmd, maxcmd)
    qout = qin + nAir * np.sin(np.deg2rad(cmd))
    bad = np.sum(qout ** 2, axis=1) >= 1 - 1e-12
    if repair and np.any(bad):
        offsets = np.arange(-2, 3) * step
        for j in np.where(bad)[0]:
            best = None
            for dx in offsets:
                for dy in offsets:
                    c2 = np.array([cmd[j, 0] + dx, cmd[j, 1] + dy])
                    if np.any(np.abs(c2) > maxcmd + 1e-12):
                        continue
                    qq = qin[j] + nAir * np.sin(np.deg2rad(c2))
                    if np.sum(qq ** 2) >= 1 - 1e-12:
                        continue
                    tgt = np.array([*qtarget[j], np.sqrt(max(0, 1 - np.sum(qtarget[j] ** 2)))])
                    dd = np.array([*qq, np.sqrt(max(0, 1 - np.sum(qq ** 2)))])
                    score = float(dd @ tgt)
                    if best is None or score > best[0]:
                        best = (score, c2, qq)
            if best is not None:
                cmd[j] = best[1]
                qout[j] = best[2]
    valid = np.sum(qout ** 2, axis=1) < 1 - 1e-12
    return cmd, qout, valid


def unit_top(d, R, p, ref, qout):
    nG = d["nGlass"]
    gk = np.sqrt(np.maximum(nG ** 2 - np.sum(qout ** 2, axis=1), np.finfo(float).tiny))
    top_b = ref + d["hOut"] * qout / gk[:, None]
    return (R @ np.column_stack([top_b, np.full(len(ref), d["zTop"])]).T).T + p


def unit_landing(d, R, p, ref, qout):
    nA = d["nAir"]
    top_w = unit_top(d, R, p, ref, qout)
    d_world = (R @ np.column_stack([qout, np.sqrt(np.maximum(nA ** 2 - np.sum(qout ** 2, axis=1), 0))]).T).T / nA
    dist = (d["rxPlaneZ"] - top_w[:, 2]) / np.maximum(d_world[:, 2], 1e-12)
    land = top_w + dist[:, None] * d_world
    return land, d_world


def direction_from_q(d, R, qout):
    nA = d["nAir"]
    return (R @ np.column_stack([qout, np.sqrt(np.maximum(nA ** 2 - np.sum(qout ** 2, axis=1), 0))]).T).T / nA


def power_weighted_target(d, qin, power, coverage=0.999, ngrid=401):
    """Minimum-angle body-frame common q that covers `coverage` of the power."""
    B = d["B"]
    gx = np.linspace(-1, 1, ngrid)
    dens = np.zeros((ngrid, ngrid))
    lo = np.maximum(-1, qin - B)
    hi = np.minimum(1, qin + B)
    i0 = np.searchsorted(gx, lo[:, 0] - 1e-12)
    i1 = np.searchsorted(gx, hi[:, 0] + 1e-12)
    j0 = np.searchsorted(gx, lo[:, 1] - 1e-12)
    j1 = np.searchsorted(gx, hi[:, 1] + 1e-12)
    for i in range(len(qin)):
        dens[j0[i]:j1[i], i0[i]:i1[i]] += power[i]
    QX, QY = np.meshgrid(gx, gx)
    ANG = np.degrees(np.arcsin(np.minimum(1, np.hypot(QX, QY))))
    mask = dens >= coverage - 1e-9
    if not np.any(mask):
        raise RuntimeError("coverage infeasible")
    aa = ANG.copy()
    aa[~mask] = np.inf
    j = np.unravel_index(np.argmin(aa), aa.shape)
    return np.array([QX[j], QY[j]]), float(dens[j])


def aperture_power_lut(d, key, center, radii=(0.05,), extent=0.6, step=0.01):
    """Precompute aperture power vs centroid offset for one anchor map."""
    x = d[key + "X"]
    y = d[key + "Y"]
    I = d[key + "I"]
    X, Y = np.meshgrid(x, y)
    cell = (x[1] - x[0]) * (y[1] - y[0])
    off = np.arange(-extent, extent + 1e-12, step)
    OX, OY = np.meshgrid(off, off)
    lut = {}
    dx = X - center[0]
    dy = Y - center[1]
    for r in radii:
        out = np.zeros_like(OX)
        for j in range(OX.shape[0]):
            for i in range(OX.shape[1]):
                m = (dx - OX[j, i]) ** 2 + (dy - OY[j, i]) ** 2 <= r * r
                out[j, i] = np.sum(I[m]) * cell
        lut[r] = (off, out)
    return lut


def lut_lookup(lut, r, ox, oy):
    off, tab = lut[r]
    ix = np.interp(ox, off, np.arange(len(off)), left=np.nan, right=np.nan)
    iy = np.interp(oy, off, np.arange(len(off)), left=np.nan, right=np.nan)
    if np.isnan(ix) or np.isnan(iy):
        return 0.0
    x0, y0 = int(np.floor(ix)), int(np.floor(iy))
    x1, y1 = min(x0 + 1, len(off) - 1), min(y0 + 1, len(off) - 1)
    fx, fy = ix - x0, iy - y0
    return float((1 - fx) * (1 - fy) * tab[y0, x0] + fx * (1 - fy) * tab[y0, x1]
                 + (1 - fx) * fy * tab[y1, x0] + fx * fy * tab[y1, x1])
