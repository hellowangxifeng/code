# Results — experimental 6-DOF evaluation of the floating RIS link

All numbers below are produced by the scripts in `FOWT_RIS_eval/tools/` from the
raw wave-tank files in `浮式平台六自由度实验数据/FOWT Platform Testing/`.
No value in this file is hand-entered from an external source.

## 1. Measured dataset

| quantity | value |
|---|---|
| regular-wave records | 9 (0.4 Hz: 10/15/30 mm; 0.8/1.2 Hz: 15/30/45 mm) |
| decay tests | heave, pitch, surge |
| motion sampling | 183 Hz, 63–181 s, 3–30 missing samples per record |
| wave gauges | 6 EDL gauges, 128 Hz; gauge 8 used as the local incident wave |
| static trim | pitch +1.40° (all records), roll +0.45° |
| dynamic tilt rms | 0.14–1.10° |
| heave RAO (0.4/0.8/1.2 Hz) | 1.64–2.12 / 0.27–0.29 / 0.06–0.07 |
| roll RAO | 14–70 deg/m |
| pitch RAO | 0.5–4.3 deg/m |

Decay identification (band-pass 0.15–1.5 Hz, zero crossings, log envelope):

| mode | natural period | damping ratio |
|---|---:|---:|
| heave | 2.158 s | 1.08% |
| pitch | 2.88 s | 0.32% |
| surge | 2.84 s | 0.78% |

## 2. Production RIS design (full-wave anchor, `RIS_phase/results/validated_optics`)

| quantity | production | collimated (new) |
|---|---:|---:|
| units at command limit | 1600/1600 | 0 |
| power-weighted direction std | 0.311° | 0.0105° |
| 90% threshold radius | 0.211 m | 0.210 m (0.170 m with curvature) |
| peak irradiance | 2.01 W/m² | 5.36 W/m² (7.57 W/m² with curvature) |
| map power | 70.13 mW | 68.97 mW |
| 5 cm aperture power | 4.58 mW (6.5%) | 14.82 mW (21.5%) |
| 10 cm aperture power | 16.70 mW (23.8%) | 40.70 mW (59.0%) |
| nominal gain, 5 cm | — | **+5.1 dB** |
| nominal gain, 10 cm | — | **+3.9 dB** |

The production command is the box corner (−15.414°, +15.414°) for every unit.
The required +z deflection is `|q_in| = 1.14` while the device can supply at
most `sqrt(2)·0.266 = 0.376`; the realised beam is ~51° from zenith and misses
the configured zenith receiver by 6.2 m.

The new design uses one common world direction
`q* = [0.420, −0.690]` (53.9° from zenith), which covers **99.94%** of the
direct unit power, so the array is collimated instead of saturated.

## 3. Experimental motion replay (50 Hz, 5 cm aperture)

Nominal power: production 4.57 mW; collimated 14.82 mW.

| scenario | mean power | p10 | mean loss | worst availability | centroid rms |
|---|---:|---:|---:|---:|---:|
| production, static | 4.52 mW | 4.14 mW | −0.05 dB | 1.00 | 44.5 mm |
| production, point-adaptive | 4.70 mW | 4.48 mW | +0.12 dB | 1.00 | 24.9 mm |
| collimated, static | 12.22 mW | 9.61 mW | −0.84 dB | **0.59** | 48.3 mm |
| collimated, adaptive direction | 14.37 mW | 13.63 mW | −0.13 dB | 1.00 | 19.2 mm |
| **collimated, point-adaptive** | **14.90 mW** | **14.71 mW** | **+0.02 dB** | **1.00** | **2.7 mm** |

Ablation (mean / p10):
- phase synthesis (production → collimated static): +4.32 / +3.66 dB
- motion compensation (static → adaptive direction): +0.70 / +1.52 dB, availability 0.59 → 1.00
- translation/link-aware (adaptive direction → point): +0.16 / +0.32 dB, centroid 19.2 → 2.7 mm
- total: **+5.18 / +5.50 dB**, availability 100%

Worst measured condition (30 mm/0.4 Hz):

| scenario | mean | p10 | loss | availability |
|---|---:|---:|---:|---:|
| production static | 3.95 | 3.28 | −0.63 dB | 1.00 |
| collimated static | 8.36 | 4.61 | −2.49 dB | 0.59 |
| collimated adaptive | 12.03 | 9.85 | −0.91 dB | 1.00 |
| point adaptive | 14.93 | 14.58 | +0.03 dB | 1.00 |

Aperture trade-off (mean over records):

| aperture radius | production | point-adaptive | gain |
|---:|---:|---:|---:|
| 25 mm | 1.15 mW | 4.73 mW | +6.1 dB |
| 50 mm | 4.52 mW | 14.90 mW | +5.2 dB |
| 100 mm | 16.06 mW | 40.79 mW | +4.0 dB |

## 4. Control study (worst record, 45 mm/0.8 Hz)

Pose sensitivity `J` (beam centroid per rad / per m):
`[[4.12, 7.50, 5.83, 1, 0, −0.71], [−11.85, −4.12, 3.51, 0, 1, 1.17]]`.
Exact nonlinear check at ±0.4°: roll 88.3 mm exact vs 87.6 mm linear; pitch
60.0 mm vs 59.7 mm.

| latency | ZOH rms / mean loss | linear predictor rms / loss |
|---:|---:|---:|
| 0 ms | 0 / 0 | 0 / 0 |
| 20 ms | 23.6 mm / −0.14 dB | 12.8 mm / −0.02 dB |
| 50 ms | 57.4 mm / −0.80 dB | 28.3 mm / −0.15 dB |
| 100 ms | 113.3 mm / −2.54 dB | 59.6 mm / −0.67 dB |

Sample-and-hold (point-adaptive, ZOH between updates):

| update interval | rms | mean loss | availability |
|---:|---:|---:|---:|
| 20 ms | 8.8 mm | +0.02 dB | 1.00 |
| 50 ms | 28.4 mm | −0.18 dB | 1.00 |
| 100 ms | 61.0 mm | −0.75 dB | 0.91 |
| 200 ms | 124.0 mm | −2.08 dB | 0.63 |

Attitude-estimate noise (zero delay): 0.05° → 14.8 mm / −0.04 dB;
0.10° → 29.2 mm / −0.25 dB; 0.20° → 58.4 mm / −1.00 dB.

## 5. Measured-wave-driven RIS model (sea_wave)

The measured wave-gauge record is converted to deep-water components
(26/42 components, synthesis rms error 0.8/1.1 mm) and used as `p.wave` in
`RIS_float`. Under the 30/45 mm 0.8 Hz waves the 1 m RIS plate responds with
pitch std 2.41°/3.77° and max tilt 3.62°/5.67°, versus 0.75°/1.10° measured on
the FOWT platform. Replaying the RIS-platform motion:

| measured wave | static collimated | point-adaptive |
|---|---:|---:|
| 30 mm/0.8 Hz | 10.39 mW, −1.54 dB, avail 0.87 | 14.85 mW, +0.01 dB, avail 1.00 |
| 45 mm/0.8 Hz | 6.86 mW, −3.35 dB, avail 0.42 | 14.71 mW, −0.03 dB, avail 1.00 |

The FOWT experiment is therefore a **conservative** motion surrogate for the
1 m RIS plate; the measured-wave-driven model provides the high-severity
counterpart.

## 6. Sea-model consistency

| quantity | FOWT experiment | RIS model |
|---|---:|---:|
| heave natural period | 2.158 s | 1.148 s |
| pitch/roll natural period | 2.88 s | 0.873 s |
| heave damping ratio | 1.08% (total) | 28.4% (radiation only) |
| pitch damping ratio | 0.32% (total) | 22.6% (radiation only) |
| tilt rms envelope | 0.14–1.10° | 1.89–7.53° (U = 3–12 m/s) |

The platform geometries differ, so the experiment cannot validate the RIS
hydrodynamic coefficients. It validates the evaluation chain, the attitude
severity at the low end, and the control-bandwidth requirements.

## 7. Files

- `data/processed/fowt_motion.npz`, `conditions.csv`, `decay_summary.csv`
- `results/nominal_design.mat` (exported RIS geometry + full-wave anchor maps)
- `results/improved_nominal.mat`, `results/improved_curvature.mat`,
  `results/current_quadratic.mat` (full-wave anchors)
- `results/replay_summary.csv`, `results/replay_00..08.npz`
- `results/aperture_tradeoff.csv`, `results/aggregate.csv`
- `results/control_study.csv`, `results/sensitivity.npy`
- `results/model_motion_replay.csv`, `results/model_replay_6/8.npz`
- `results/sea_model_envelope.mat`
- `figures/Fig1_dataset.png` … `Fig4_control.png`
