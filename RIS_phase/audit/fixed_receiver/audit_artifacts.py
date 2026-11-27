"""Verify revised fixed-RX curves, their actual phase schedules and provenance."""
from pathlib import Path
import csv
import hashlib
import json
import numpy as np
import h5py
from PIL import Image

root = Path(__file__).resolve().parents[2]
out = root / 'plots' / 'stage2'
fixed = out / 'fixed_receiver'
sha = lambda p: hashlib.sha256(p.read_bytes()).hexdigest()
v = json.loads((fixed / 'validation.json').read_text(encoding='utf-8'))
assert v['status'] == 'passed'
with h5py.File(fixed / 'comparison.mat', 'r') as f:
    d = f['data']
    sig = ''.join(chr(x) for x in d['fixedReceiverSignature'][...].ravel())
    assert sig == v['signature']
    assert abs(float(d['fovHalfAngleRad'][0, 0]) - .130) < 1e-12
    assert abs(float(d['positionSigmaM'][0, 0]) - .080) < 1e-12
    assert np.max(np.abs(d['receiverCenterXY'][...].ravel() - np.array(v['receiverCenterXY']))) < 1e-12
    selected = set(np.unique(d['choice'][...]).astype(int))
with (fixed / 'curves_and_phase_schedule.csv').open(encoding='utf-8-sig', newline='') as f:
    schedule = list(csv.DictReader(f))
assert len(schedule) == 2412
assert len({(r['receiverXM'], r['receiverYM']) for r in schedule}) == 1
programs = []
for p in sorted((fixed / 'phase_programs').glob('candidate_*.csv')):
    candidate = int(p.stem.split('_')[-1])
    with p.open(encoding='utf-8-sig', newline='') as f:
        rows = list(csv.DictReader(f))
    assert len(rows) == 1600
    command = np.array([[float(r['commandXDeg']), float(r['commandYDeg'])] for r in rows])
    gradient = np.array([[float(r['gradientXOverK0']), float(r['gradientYOverK0'])] for r in rows])
    piston = np.array([float(r['pistonRad']) for r in rows])
    assert np.max(np.abs(command)) <= 15.414 + 1e-12
    assert np.max(np.abs(command / .021 - np.rint(command / .021))) < 1e-9
    assert np.max(np.abs(gradient - np.sin(np.deg2rad(command)))) < 1e-12
    assert np.all((piston >= 0) & (piston < 2 * np.pi))
    if candidate == 0:
        assert np.max(np.abs(command)) == np.max(np.abs(piston)) == 0
    if candidate == 1:
        assert np.max(np.ptp(command, axis=0)) < 1e-12
    programs.append({'candidate': candidate, 'file': p.name, 'sha256': sha(p)})
assert {r['candidate'] for r in programs} == selected
assert len(programs) == v['programCount']
with (fixed / 'graphics_audit.csv').open(encoding='utf-8-sig', newline='') as f:
    graphics = list(csv.DictReader(f))
assert len(graphics) == 2
for r in graphics:
    assert float(r['fontPt']) == 10 and float(r['legendOverlap']) == 0
    assert r['legendCorner'] in {'northwest', 'northeast', 'southwest', 'southeast'}
    assert float(r['xMin']) == 0 and float(r['xMax']) == 100
    assert float(r['widthIn']) == 2.67 and float(r['heightIn']) == 2.37
figures = []
for stem in ['Fig2_aber', 'Fig3_op']:
    png = out / (stem + '.png')
    assert Image.open(png).size == (801, 711)
    figures.append({'stem': stem, 'pngSHA256': sha(png), 'figSHA256': sha(out / (stem + '.fig'))})
with (root / 'audit' / 'source_inventory.csv').open(encoding='utf-8-sig', newline='') as f:
    protected = [r for r in csv.DictReader(f) if r['path'].replace('\\', '/').startswith('W_RIS/')]
assert all(sha(root.parent / r['path']) == r['sha256'] for r in protected)
readme = json.loads((root / 'audit' / 'fixed_receiver' / 'readme_observed.json').read_text(encoding='utf-8-sig'))
assert sha(root / 'README.md') == readme['sha256']
manifest = {'status': 'passed', 'signature': sig, 'phasePolicy': 'per_power',
            'fixedReceiver': v['receiverCenterXY'], 'sourcePhotons': v['scatterPhotons'],
            'figures': figures, 'programs': programs, 'scheduleRows': len(schedule),
            'readmeUnchangedSinceObservation': True, 'readmeCurrentSHA256': readme['sha256'],
            'readmeOriginalStage2SHA256': readme['originalStage2SHA256'],
            'readmeMatchesOriginalStage2': readme['sha256'] == readme['originalStage2SHA256'],
            'readmeWrittenByThisRevision': False, 'wrisUnchanged': True,
            'datasets': {p.name: sha(p) for p in sorted(fixed.glob('*.mat'))},
            'sources': {str(p.relative_to(root)): sha(p) for p in sorted(
                list(root.glob('rp_fixed*.m')) + list(root.glob('rp_nonstationary_rect_bound.m')) +
                list((root / 'plots').glob('rp_fixed*.m')) + list((root / 'plots').glob('rp_stage2_plot_link*.m')))}}
(fixed / 'artifact_validation.json').write_text(json.dumps(manifest, ensure_ascii=False, indent=2), encoding='utf-8')
print(f"PASS: fixed RX, {len(schedule)} schedule rows, {len(programs)} realized programs, two 801x711 figures, current README/W_RIS protected checks passed.")
