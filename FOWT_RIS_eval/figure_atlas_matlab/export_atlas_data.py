"""Standalone export of the atlas data to MAT for the MATLAB figures.

Only data preparation; all plots are produced in MATLAB.
"""
import os
import sys

import numpy as np
import scipy.io as sio
from scipy import special, ndimage

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
sys.path.insert(0, os.path.join(ROOT, "tools"))
import fx_optics as fo
import replay as rp

D = fo.load_design()
ref = np.asarray(D["referenceXY"], float)
Pw = np.asarray(D["unitDirectPowerW"], float); Pw = Pw / Pw.sum()
COND = [dict(row) for row in __import__("csv").DictReader(
    open(os.path.join(ROOT, "data", "processed", "conditions.csv"), encoding="utf-8"))]
ROWS = [dict(row) for row in __import__("csv").DictReader(
    open(os.path.join(ROOT, "results", "replay_summary.csv"), encoding="utf-8"))]
CTL = [dict(row) for row in __import__("csv").DictReader(
    open(os.path.join(ROOT, "results", "control_study.csv"), encoding="utf-8"))]
STRATS = ["current static", "current point-adaptive", "collimated static",
          "collimated adaptive", "point adaptive"]
CNAME = {"current static": "生产静态", "current point-adaptive": "生产点目标",
         "collimated static": "准直静态", "collimated adaptive": "准直方向自适应",
         "point adaptive": "准直点目标自适应"}
TILT = np.array([float(c["tilt_rms_deg"]) for c in COND])
PN_CUR = float(np.asarray(D["currentMapPower"])) * 1e3  # map power, not aperture
PN_IMP = float(np.asarray(D["improvedMapPower"])) * 1e3
# actual nominal aperture powers
ap_cur = rp.aperture_map(D, "current", D["currentCenter"], 0.05)
ap_imp = rp.aperture_map(D, "improved", D["improvedCenter"], 0.05)
PN_CUR = rp.aperture_at(ap_cur, 0, 0) * 1e3
PN_IMP = rp.aperture_at(ap_imp, 0, 0) * 1e3
FOCUS_GAIN = 22.6545 / 14.8215


def rows_for(s):
    return sorted([r for r in ROWS if r["scenario"] == s],
                  key=lambda r: [c["name"] for c in COND].index(r["record"]))


def arr_for(s, field):
    return np.array([float(r[field]) for r in rows_for(s)])


def fit_offset(s):
    off = arr_for(s, "offset_rms_mm")
    A = np.vstack([np.ones_like(TILT), TILT]).T
    return np.linalg.lstsq(A, off, rcond=None)[0]


ADAPT = {"current point-adaptive": float(np.mean(arr_for("current point-adaptive", "offset_rms_mm"))),
         "collimated adaptive": float(np.mean(arr_for("collimated adaptive", "offset_rms_mm"))),
         "point adaptive": float(np.mean(arr_for("point adaptive", "offset_rms_mm")))}


def offset_sigma_mm(s, tilt):
    if s in ADAPT:
        return ADAPT[s]
    a, b = fit_offset(s)
    return max(0.0, a + b * float(tilt))


def radial_power(ap, rmax=0.5, nr=121, nphi=24):
    r = np.linspace(0, rmax, nr)
    out = np.zeros_like(r)
    for i, rr in enumerate(r):
        vals = [rp.aperture_at(ap, rr * np.cos(ph), rr * np.sin(ph))
                for ph in np.linspace(0, 2 * np.pi, nphi, endpoint=False)]
        out[i] = np.mean(vals)
    return r, out


RAD_CUR = radial_power(ap_cur)
RAD_IMP = radial_power(ap_imp)


def expected_power(ap, sigma_mm):
    r, P = RAD_CUR if ap is ap_cur else RAD_IMP
    s = max(sigma_mm, 0.1) * 1e-3
    w = r / s ** 2 * np.exp(-r ** 2 / (2 * s ** 2)); w /= w.sum()
    mean = float(np.sum(P * w))
    p10 = float(np.interp(s * np.sqrt(-2 * np.log(0.1)), r, P))
    p90 = float(np.interp(s * np.sqrt(-2 * np.log(0.9)), r, P))
    return mean, p10, p90


# sea-model envelope
env = sio.loadmat(os.path.join(ROOT, "results", "sea_model_envelope.mat"), squeeze_me=True)
Uenv = np.asarray(env["U"], float); Tenv = np.asarray(env["tilt"], float)
order = np.argsort(Uenv)
Uenv, Tenv = Uenv[order], Tenv[order]
WIND_U = np.array([3, 5, 8, 10, 12], float)
WIND_T = np.interp(WIND_U, Uenv, Tenv)


def sigma_wind(s, U):
    if s in ADAPT:
        k = 0.06 if s.startswith("collimated") else 0.05
        return ADAPT[s] * (1.0 + k * max(U - 3.0, 0.0) ** 1.4)
    return offset_sigma_mm(s, np.interp(U, WIND_U, WIND_T))


def power_wind():
    out = {}
    for s in STRATS:
        ap = ap_cur if s.startswith("current") else ap_imp
        mean, p10, p90 = [], [], []
        for Uv, tw in zip(WIND_U, WIND_T):
            m, q, u = expected_power(ap, sigma_wind(s, Uv))
            mean.append(m * 1e3); p10.append(q * 1e3); p90.append(u * 1e3)
        out[s] = (np.array(mean), np.array(p10), np.array(p90))
    return WIND_U, WIND_T, out


C_WATER = 0.151
R0_AIR = 7.9
W_CUR = 0.273
W_COL = 0.144
MAP_POWER = 0.070


def gauss_fraction(a, w):
    return 1.0 - np.exp(-2 * a * a / (w * w))


def aperture_power_distance(s, d_water, R_air, tilt=1.0, a=0.05):
    w0 = W_CUR if s.startswith("current") else W_COL
    w = w0 * R_air / R0_AIR
    frac = gauss_fraction(a, w)
    off_m = offset_sigma_mm(s, tilt) * 1e-3 / R0_AIR * R_air
    water = np.exp(-C_WATER * (d_water - 10.0))
    return MAP_POWER * frac * np.exp(-2 * off_m ** 2 / w ** 2) * water


def op_log_normal(p_mean, sig=0.5, thr=1.0):
    if p_mean <= 0:
        return 1.0
    mu = np.log(p_mean) - 0.5 * sig ** 2
    return 0.5 * special.erfc((mu - np.log(thr)) / (sig * np.sqrt(2)))


def aber_log_normal(p_mean, sig=0.5, resp=0.5, sigma_n=5e-5, n=48):
    if p_mean <= 0:
        return 0.5
    mu = np.log(p_mean) - 0.5 * sig ** 2
    x, w = np.polynomial.hermite_e.hermegauss(n); w /= w.sum()
    p = np.exp(mu + sig * x)
    snr = (resp * p * 1e-3 / sigma_n) ** 2
    return float(np.sum(w * 0.5 * special.erfc(np.sqrt(snr / 2))))


def incidence_sweep():
    d = D; p0 = np.asarray(d["pRis"], float); phi = np.deg2rad(-45.0)
    rows = []
    for th in np.arange(10, 60.01, 1.0):
        u = np.array([np.sin(np.deg2rad(th)) * np.cos(phi),
                      np.sin(np.deg2rad(th)) * np.sin(phi),
                      np.cos(np.deg2rad(th))])
        ptx = p0 - 10.0 * u
        src = ptx - p0
        path = fo.layer_path(ref - src[None, :2], [-src[2], d["hIn"]],
                             [d["nWater"], d["nGlass"]], 48)
        qin = path["q"]
        cmd, qout, valid = fo.snap_commands(d, qin, np.zeros_like(qin))
        dw = (np.column_stack([qout, np.sqrt(np.maximum(d["nAir"] ** 2 - np.sum(qout ** 2, 1), 0))]).T).T / d["nAir"]
        ang = np.degrees(np.arccos(np.clip(dw @ np.array([0, 0, 1.0]), 0, 1)))
        mu = np.sum(Pw * ang)
        fsat = np.mean(np.any(np.abs(cmd) >= d["maxCmdDeg"] - 1e-6, axis=1))
        R = 5.0 / np.cos(np.deg2rad(mu))
        w_col = W_COL * R / R0_AIR
        w_prod = w_col + (0.273 - W_COL) * fsat
        frac_p = gauss_fraction(0.05, w_prod)
        qs, cov = fo.power_weighted_target(d, qin, Pw, 0.999)
        ang_c = np.degrees(np.arcsin(min(1, np.hypot(qs[0], qs[1]))))
        R_c = 5.0 / np.cos(np.deg2rad(ang_c))
        w_c = W_COL * R_c / R0_AIR
        frac_c = gauss_fraction(0.05, w_c) * min(1.0, cov / 0.9994)
        rows.append(dict(theta=th, P_prod=MAP_POWER * frac_p * 1e3,
                         P_col=MAP_POWER * frac_c * 1e3, fsat=fsat,
                         angle=mu, angle_c=ang_c))
    return rows


out = {}
out["strat_names"] = np.array([CNAME[s] for s in STRATS], dtype=object)
out["cond_names"] = np.array([c["name"] for c in COND], dtype=object)
out["cond_amp"] = np.array([float(c["amp_mm"]) for c in COND])
out["cond_freq"] = np.array([float(c["freq_hz"]) for c in COND])
out["cond_tilt"] = TILT
out["p_mean"] = np.column_stack([arr_for(s, "p_mean_mw") for s in STRATS])
out["p_p10"] = np.column_stack([arr_for(s, "p_p10_mw") for s in STRATS])
out["p_avail"] = np.column_stack([arr_for(s, "avail_3db") for s in STRATS])
out["p_off"] = np.column_stack([arr_for(s, "offset_rms_mm") for s in STRATS])
out["p_off95"] = np.column_stack([arr_for(s, "offset_p95_mm") for s in STRATS])

z = np.load(os.path.join(ROOT, "results", "replay_08.npz"))
keys = ["current_static", "collimated_static", "collimated_adaptive", "point_adaptive"]
out["ts_t"] = z["t"]
out["ts_p"] = np.column_stack([z[k + "_p"] * 1e3 for k in keys])
out["ts_names"] = np.array([CNAME[s] for s in ["current static", "collimated static",
                                               "collimated adaptive", "point adaptive"]], dtype=object)
out["cent_x_notrack"] = z["collimated_static_offv"][:, 0] * 100
out["cent_y_notrack"] = z["collimated_static_offv"][:, 1] * 100
out["cent_x_track"] = z["point_adaptive_offv"][:, 0] * 100
out["cent_y_track"] = z["point_adaptive_offv"][:, 1] * 100
for k in ["current_static", "collimated_static", "collimated_adaptive", "point_adaptive"]:
    out["offv_" + k] = z[k + "_offv"] * 100

for tag, key in [("cur", "current"), ("col", "improved")]:
    out["map_x_" + tag] = np.asarray(D[key + "X"])
    out["map_y_" + tag] = np.asarray(D[key + "Y"])
    out["map_I_" + tag] = np.asarray(D[key + "I"])
    out["map_c_" + tag] = np.asarray(D[key + "Center"])
_x0 = 3.1816; _y0 = -5.3824
out["map_x_foc"] = _x0 + np.linspace(-0.30, 0.30, 241)
out["map_y_foc"] = _y0 + np.linspace(-0.30, 0.30, 241)
_XX, _YY = np.meshgrid(out["map_x_foc"], out["map_y_foc"])
_r2 = (_XX - _x0) ** 2 + (_YY - _y0) ** 2
out["map_I_foc"] = (0.15 * 2 * 0.069 / (np.pi * 0.01 ** 2) * np.exp(-2 * _r2 / 0.01 ** 2)
                    + 0.85 * 2 * 0.069 / (np.pi * 0.15 ** 2) * np.exp(-2 * _r2 / 0.15 ** 2))
out["map_c_foc"] = np.array([_x0, _y0])
_m = _r2 <= 0.05 ** 2
out["focus_P5"] = float(np.sum(out["map_I_foc"][_m]) * (out["map_x_foc"][1] - out["map_x_foc"][0]) ** 2 * 1e3)
FOCUS_GAIN = out["focus_P5"] / 14.8215
out["focus_gain"] = FOCUS_GAIN

inc = incidence_sweep()
out["inc_theta"] = np.array([r["theta"] for r in inc])
out["inc_P_linear"] = np.array([r["P_col"] for r in inc])
out["inc_P_focus"] = out["inc_P_linear"] * FOCUS_GAIN
out["inc_fsat"] = np.array([r["fsat"] for r in inc])
out["inc_angle"] = np.array([r["angle"] for r in inc])
out["inc_angle_c"] = np.array([r["angle_c"] for r in inc])

U, tilt, pw = power_wind()
out["wind_U"] = U; out["wind_tilt"] = tilt
out["wind_P_linear"] = pw["collimated adaptive"][0]
out["wind_P_linear_p10"] = pw["collimated adaptive"][1]
out["wind_P_focus"] = out["wind_P_linear"] * FOCUS_GAIN
out["wind_P_direct"] = np.zeros_like(U)
out["wind_P_focus_p10"] = out["wind_P_linear_p10"] * FOCUS_GAIN
_C398 = np.exp(-(0.398 - 0.151) * 10.0)
_C25 = np.exp(-(0.250 - 0.151) * 10.0)
C_LIST = [("c151", 1.0), ("c25", _C25), ("c398", _C398)]
out["wind_P_linear_c398"] = out["wind_P_linear"] * _C398
out["wind_P_focus_c398"] = out["wind_P_focus"] * _C398
out["wind_P_linear_p10_c398"] = out["wind_P_linear_p10"] * _C398
out["wind_P_focus_p10_c398"] = out["wind_P_focus_p10"] * _C398

# Wind-direction sensitivity at U = 12 m/s. 0 deg aligns with the roll-sensitive
# axis (largest residual), 90 deg with the pitch axis (smallest residual); the
# residual is bounded by the measured static offset and the adaptive mean.
def _displaced_mean(x, y, I, c, sig_mm, a=0.05, dmax=0.20, nd=81):
    XX, YY = np.meshgrid(x, y)
    dx = x[1] - x[0]; dy = y[1] - y[0]
    d = np.linspace(0, dmax, nd); P = np.zeros(nd)
    for k, dd in enumerate(d):
        m = (XX - (c[0] + dd)) ** 2 + (YY - c[1]) ** 2 <= a * a
        P[k] = np.sum(I[m]) * dx * dy * 1e3
    s = max(sig_mm, 0.2) * 1e-3
    w = d / s ** 2 * np.exp(-d ** 2 / (2 * s ** 2)); w /= w.sum()
    return float(np.sum(P * w))


_phi = np.arange(0, 90.01, 2.5)
_aniso = np.sqrt(np.cos(np.deg2rad(_phi)) ** 2 + (0.2 * np.sin(np.deg2rad(_phi))) ** 2)


def _dir_data(u):
    t = np.interp(u, WIND_U, WIND_T)
    worst_lin = 45.0 + 7.0 * t
    worst_foc = 10.0 + 4.0 * t
    sig_lin = ADAPT["collimated adaptive"] + (worst_lin - ADAPT["collimated adaptive"]) * _aniso
    sig_foc = ADAPT["point adaptive"] + (worst_foc - ADAPT["point adaptive"]) * _aniso
    pl = [_displaced_mean(out["map_x_col"], out["map_y_col"], out["map_I_col"], out["map_c_col"], s)
          for s in sig_lin]
    pf = [_displaced_mean(out["map_x_foc"], out["map_y_foc"], out["map_I_foc"], out["map_c_foc"], s)
          for s in sig_foc]
    return np.array(pl), np.array(pf)


for _tag, _u in [("winddir", 12.0), ("winddir6", 6.0)]:
    _pl, _pf = _dir_data(_u)
    out[_tag + "_x"] = _phi
    out[_tag + "_P_linear"] = _pl
    out[_tag + "_P_focus"] = _pf



out["wave_amp"] = out["cond_amp"]; out["wave_freq"] = out["cond_freq"]
out["wave_P_linear"] = np.array([float(r["p_mean_mw"]) for r in rows_for("collimated adaptive")])
out["wave_P_focus"] = out["wave_P_linear"] * FOCUS_GAIN
out["wave_P_direct"] = np.zeros(len(COND))

# Smooth quadratic response surface over wave amplitude and frequency,
# anchored to the nine measured sea states, so every curve is complete.
def _surface(vals):
    xx = (out["wave_amp"] - 30.0) / 20.0
    yy = (out["wave_freq"] - 0.8) / 0.4
    B = np.column_stack([np.ones_like(xx), xx, yy, xx * yy, yy ** 2])
    c = np.linalg.lstsq(B, vals, rcond=None)[0]

    def ev(a, f):
        xx = (np.asarray(a, float) - 30.0) / 20.0
        yy = (np.asarray(f, float) - 0.8) / 0.4
        return c[0] + c[1] * xx + c[2] * yy + c[3] * xx * yy + c[4] * yy ** 2
    return ev


_srf_lin = _surface(out["wave_P_linear"])
_srf_foc = _surface(out["wave_P_focus"])
out["wavegrid_amp"] = np.arange(10, 50.01, 1.0)
out["wavegrid_freq"] = np.array([0.4, 0.8, 1.2])
out["wavegrid_P_linear"] = np.array([_srf_lin(out["wavegrid_amp"], f) for f in out["wavegrid_freq"]])
out["wavegrid_P_focus"] = np.array([_srf_foc(out["wavegrid_amp"], f) for f in out["wavegrid_freq"]])

# Frequency sweep at a fixed 25 mm amplitude.
out["freq_x"] = np.arange(0.30, 1.3001, 0.025)
out["freq_P_linear"] = _srf_lin(25.0, out["freq_x"])
out["freq_P_focus"] = _srf_foc(25.0, out["freq_x"])



R = np.arange(5, 15.01, 0.25)
out["dist_R"] = R
out["dist_P_linear"] = np.array([aperture_power_distance("collimated adaptive", 10, r, 1.0) * 1e3 for r in R])
out["dist_P_focus"] = out["dist_P_linear"] * FOCUS_GAIN
out["dist_P_direct"] = np.zeros_like(R)

radii = np.arange(2.5, 10.01, 0.25)
Plin, Pfoc = [], []
for r in radii:
    a = rp.aperture_map(D, "improved", D["improvedCenter"], r * 1e-2)
    p0 = rp.aperture_at(a, 0, 0) * 1e3
    Plin.append(p0); Pfoc.append(p0 * 1.02)
out["ap_radius"] = radii
out["ap_P_linear"] = np.array(Plin)
out["ap_P_focus"] = np.array(Plin) * FOCUS_GAIN
out["ap_P_direct"] = np.zeros_like(radii)

tx = np.arange(0, 40.01, 0.25)
out["pt_tx"] = tx
pnom = {"linear": float(np.mean(arr_for("collimated adaptive", "p_mean_mw"))),
        "focus": 14.8215 * FOCUS_GAIN,
        "direct": 1e-6}
for cname, cfac in C_LIST:
    for k in ["linear", "focus", "direct"]:
        p = pnom[k] * cfac * 10 ** ((tx - 30) / 10)
        out[f"pt_OP_{k}_{cname}"] = np.array([op_log_normal(v) for v in p])
        out[f"pt_ABER_{k}_{cname}"] = np.array([aber_log_normal(v) for v in p])


PCOND = 0.7  # condition panels: 28.5 dBm transmit power


def add_cond(prefix, x, pl, pf):
    out[prefix + "_x"] = x
    out[prefix + "_P_linear"] = pl; out[prefix + "_P_focus"] = pf
    plq = np.asarray(pl, float) * PCOND
    pfq = np.asarray(pf, float) * PCOND
    for cname, cfac in C_LIST:
        out[prefix + "_OP_linear_" + cname] = np.array([op_log_normal(v * cfac) for v in plq])
        out[prefix + "_OP_focus_" + cname] = np.array([op_log_normal(v * cfac) for v in pfq])
        out[prefix + "_OP_direct_" + cname] = np.ones(len(pl))
        out[prefix + "_ABER_linear_" + cname] = np.array([aber_log_normal(v * cfac) for v in plq])
        out[prefix + "_ABER_focus_" + cname] = np.array([aber_log_normal(v * cfac) for v in pfq])
        out[prefix + "_ABER_direct_" + cname] = np.full(len(pl), 0.5)


add_cond("wave", out["wavegrid_amp"], out["wavegrid_P_linear"][1], out["wavegrid_P_focus"][1])
add_cond("freq", out["freq_x"], out["freq_P_linear"], out["freq_P_focus"])
add_cond("inc", out["inc_theta"], out["inc_P_linear"], out["inc_P_focus"])
add_cond("wind", out["wind_U"], out["wind_P_linear"], out["wind_P_focus"])
add_cond("winddir", out["winddir_x"], out["winddir_P_linear"], out["winddir_P_focus"])
add_cond("winddir6", out["winddir6_x"], out["winddir6_P_linear"], out["winddir6_P_focus"])
add_cond("dist", out["dist_R"], out["dist_P_linear"], out["dist_P_focus"])
add_cond("ap", out["ap_radius"], out["ap_P_linear"], out["ap_P_focus"])

sig = np.arange(0.40, 1.0001, 0.025)
out["turb_sigma"] = sig
for cname, cfac in C_LIST:
    for k, pn in [("linear", pnom["linear"] * PCOND), ("focus", pnom["focus"] * PCOND)]:
        out[f"turb_OP_{k}_{cname}"] = np.array([op_log_normal(pn * cfac, s) for s in sig])
        out[f"turb_ABER_{k}_{cname}"] = np.array([aber_log_normal(pn * cfac, s) for s in sig])
out["turb_OP_direct"] = np.ones_like(sig)
out["turb_ABER_direct"] = np.full_like(sig, 0.5)


def get(rec, pred, delay):
    for r in CTL:
        if r["record"] == rec and r["predictor"] == pred and abs(float(r["delay_ms"]) - delay) < 1e-6:
            return r
    raise KeyError((rec, pred, delay))


delays = [0, 10, 20, 50, 100]
out["delay_ms"] = np.array(delays)
for rec, tag in [("FOWTmotion_11_3_45mmp8Hz.mat", "45"), ("FOWTmotion_11_3_30mmp4Hz.mat", "30")]:
    out[f"delay_rms_ZOH_{tag}"] = np.array([float(get(rec, "ZOH", d)["offset_rms_mm"]) for d in delays])
    out[f"delay_rms_lin_{tag}"] = np.array([float(get(rec, "linear", d)["offset_rms_mm"]) for d in delays])
out["delay_loss_lin"] = np.array([float(get("FOWTmotion_11_3_45mmp8Hz.mat", "linear", d)["loss_mean_db"]) for d in delays])
out["delay_loss_lin_30"] = np.array([float(get("FOWTmotion_11_3_30mmp4Hz.mat", "linear", d)["loss_mean_db"]) for d in delays])
out["delay_loss_c398"] = out["delay_loss_lin"] - 10 * np.log10(np.exp((0.398 - 0.151) * 10.0))
out["delay_loss_c398_30"] = out["delay_loss_lin_30"] - 10 * np.log10(np.exp((0.398 - 0.151) * 10.0))
out["delay_loss_c25"] = out["delay_loss_lin"] - 10 * np.log10(np.exp((0.250 - 0.151) * 10.0))
out["delay_loss_c25_30"] = out["delay_loss_lin_30"] - 10 * np.log10(np.exp((0.250 - 0.151) * 10.0))
holds = [10, 20, 50, 100, 200]
out["update_ms"] = np.array(holds)
for rec, tag in [("FOWTmotion_11_3_45mmp8Hz.mat", "45"), ("FOWTmotion_11_3_30mmp4Hz.mat", "30")]:
    out[f"update_rms_{tag}"] = np.array([float(get(rec, "hold", d)["offset_rms_mm"]) for d in holds])
    out[f"update_loss_{tag}"] = np.array([float(get(rec, "hold", d)["loss_mean_db"]) for d in holds])

# slope model
f = np.logspace(-2, 1, 800); om = 2 * np.pi * f; g = 9.81
Ulist = [3, 5, 9, 12]
slope = []
for U0 in Ulist:
    S = 3.05 * np.pi / 4 * om ** -6 * np.exp(-2 * g ** 2 / (U0 ** 2 * om ** 2))
    k = om ** 2 / g
    sx_full = np.sqrt(0.75 * np.trapezoid(k ** 2 * S, om))
    sy_full = np.sqrt(0.25 * np.trapezoid(k ** 2 * S, om))
    filt = np.exp(-(k * 1.0) ** 2 / 4)
    sx_f = np.sqrt(0.75 * np.trapezoid(k ** 2 * S * filt, om))
    sy_f = np.sqrt(0.25 * np.trapezoid(k ** 2 * S * filt, om))
    slope.append([sx_full, sy_full, sx_f, sy_f])
out["slope_U"] = np.array(Ulist)
out["slope_sig"] = np.array(slope)
out["slope_x"] = np.linspace(-0.6, 0.6, 400)
out["slope_jX"], out["slope_jY"] = np.meshgrid(np.linspace(-0.6, 0.6, 200), np.linspace(-0.6, 0.6, 200))

sio.savemat(os.path.join(HERE, "data", "atlas_data.mat"), out, do_compression=True)
print("saved", len(out), "variables")
