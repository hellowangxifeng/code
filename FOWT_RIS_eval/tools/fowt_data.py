"""Loader for the FOWT wave-tank 6-DOF dataset (Qualisys .mat + EDL wave gauges).

All raw files live in
  <workspace>/浮式平台六自由度实验数据/FOWT Platform Testing/

The loader returns SI units: positions in metres, Euler angles in radians,
wave elevation in metres, time in seconds from the first valid sample.

Qualisys layout (verified against MainCode.m and the raw structs):
  Measurement.Qualisys_Data.RigidBodies.Positions  (1,3,N)  [mm]
  Measurement.Qualisys_Data.RigidBodies.RPYs       (1,3,N)  [deg]
  Measurement.FrameRate                            scalar   [Hz]
Wave-gauge file layout: 6 gauge channels, each block of 4 columns
  Value, Peak, Offset, Gauge[m]; the gauge elevation is column 4+4*k.
"""
from __future__ import annotations

import glob
import os
import re

import numpy as np
import scipy.io as sio

RAW_DIR = os.path.join(
    os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))),
    "浮式平台六自由度实验数据", "FOWT Platform Testing",
)


def _squeeze_first(x):
    x = np.asarray(x)
    while x.ndim > 1 and x.shape[0] == 1:
        x = x[0]
    return x


def fill_nan_linear(t, y):
    """Linear interpolation over NaN samples; leading/trailing held constant."""
    y = np.array(y, dtype=float, copy=True)
    good = np.isfinite(y)
    if good.all():
        return y
    if not good.any():
        return np.zeros_like(y)
    idx = np.arange(y.size)
    y[~good] = np.interp(idx[~good], idx[good], y[good])
    return y


def load_motion(path):
    d = sio.loadmat(path, squeeze_me=False, struct_as_record=False)
    key = [k for k in d if not k.startswith("__")][0]
    m = d[key][0, 0]
    rb = m.RigidBodies[0, 0]
    pos = np.atleast_2d(_squeeze_first(rb.Positions)).T  # N x 3, mm
    rpy = np.atleast_2d(_squeeze_first(rb.RPYs)).T      # N x 3, deg
    fs = float(np.asarray(m.FrameRate).ravel()[0])
    frames = int(np.asarray(m.Frames).ravel()[0])
    n = min(frames, pos.shape[0], rpy.shape[0])
    pos, rpy = pos[:n], rpy[:n]
    n_nan = int(np.isnan(pos).any(axis=1).sum())
    for j in range(3):
        pos[:, j] = fill_nan_linear(None, pos[:, j])
        rpy[:, j] = fill_nan_linear(None, rpy[:, j])
    t = np.arange(n) / fs
    return dict(
        name=os.path.basename(path),
        key=key,
        rigid_body=str(np.asarray(rb.Name).ravel()[0]),
        fs=fs,
        n=n,
        n_nan=n_nan,
        t=t,
        pos_mm=pos,
        pos_m=pos * 1e-3,
        rpy_deg=rpy,
        rpy_rad=np.deg2rad(rpy),
        residual=np.asarray(_squeeze_first(rb.Residual), dtype=float)[:n],
    )


def _parse_clock(s):
    hh, mm, ss = s.split(":")
    return int(hh) * 3600.0 + int(mm) * 60.0 + float(ss)


def load_wave(path):
    eta, ts = [], []
    with open(path, "r", errors="replace") as fh:
        lines = fh.readlines()
    for line in lines[5:]:
        p = line.rstrip("\n").split("\t")
        if len(p) < 24:
            continue
        try:
            t = _parse_clock(p[0])
            v = [float(p[4 + 4 * k]) for k in range(6)]
        except (ValueError, IndexError):
            continue
        ts.append(t)
        eta.append(v)
    ts = np.asarray(ts)
    eta = np.asarray(eta)
    # unwrap midnight
    dt = np.diff(ts)
    for i in np.where(dt < -1000)[0]:
        ts[i + 1:] += 86400.0
    t = ts - ts[0]
    # EDL sometimes emits duplicated clock strings; keep strictly increasing
    keep = np.concatenate([[True], np.diff(t) > 0])
    t, eta = t[keep], eta[keep]
    return dict(name=os.path.basename(path), t=t, eta=eta, fs=1.0 / np.median(np.diff(t)))


COND_RE = re.compile(r"Amp=([0-9.]+)m")
FREQ_RE = re.compile(r"Freq=([0-9.]+)Hz")


def catalog():
    """Return a list of records describing every raw file and its condition."""
    out = []
    for p in sorted(glob.glob(os.path.join(RAW_DIR, "FOWTmotion_*.mat"))):
        b = os.path.basename(p)
        rec = dict(kind="regular", motion=p, wave=None, amp_mm=None, freq_hz=None)
        m = re.search(r"(\d+)mm", b)
        if m:
            rec["amp_mm"] = int(m.group(1))
        m = re.search(r"p(\d+)Hz", b)
        if m:
            rec["freq_hz"] = int(m.group(1)) / 10.0
        m = re.search(r"(\d+)p(\d+)Hz", b)
        if m and rec["freq_hz"] is None:
            rec["freq_hz"] = float(m.group(1) + "." + m.group(2))
        m = re.search(r"(\d+)_(\d+)_(\d+)mm(\d+)p(\d+)Hz", b)
        if m:
            rec["amp_mm"] = int(m.group(3))
            rec["freq_hz"] = float(m.group(4) + "." + m.group(5))
        rec["name"] = b
        out.append(rec)
    for p in sorted(glob.glob(os.path.join(RAW_DIR, "*.txt"))):
        b = os.path.basename(p)
        fm, am = FREQ_RE.search(b), COND_RE.search(b)
        out.append(dict(
            kind="wave",
            motion=None,
            wave=p,
            name=b,
            freq_hz=float(fm.group(1)) if fm else None,
            amp_mm=int(round(float(am.group(1)) * 1000)) if am else None,
        ))
    for p in sorted(glob.glob(os.path.join(RAW_DIR, "*Decay.mat"))):
        b = os.path.basename(p)
        out.append(dict(kind="decay", motion=p, wave=None, name=b,
                        mode=b.replace("Decay.mat", "").lower()))
    return out
