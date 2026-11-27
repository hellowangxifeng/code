"""Independent second-stage file, LC-program, protected-input and math audit."""
from pathlib import Path
import csv
import hashlib
import json
import re
import numpy as np
from scipy.io import loadmat

root = Path(__file__).resolve().parents[2]
workspace = root.parent
out = root / 'plots' / 'stage2'
audit = root / 'audit' / 'stage2'
sha = lambda p: hashlib.sha256(p.read_bytes()).hexdigest()
assert sha(root / 'README.md') == (audit / 'readme_sha256.txt').read_text().strip()
with (root / 'audit' / 'source_inventory.csv').open(encoding='utf-8-sig', newline='') as f:
    protected = [r for r in csv.DictReader(f) if r['path'].replace('\\', '/').startswith('W_RIS/')]
assert protected
assert all(sha(workspace / r['path']) == r['sha256'] for r in protected)
stems = ['Fig1_irradiance', 'Fig2_aber', 'Fig3_op', 'Fig4_width_tradeoff',
         'Fig5_sensitivity', 'Fig6_objective_phase', 'Fig8_radial_profiles',
         'Fig9_turbulence', 'Fig10_incidence',
         'fig12_revision/Fig12a_radius', 'fig12_revision/Fig12b_aber', 'fig12_revision/Fig12c_op',
         'fig12_revision/Fig13a_aber', 'fig12_revision/Fig13b_op',
         'fig12_revision/Fig14a_aber', 'fig12_revision/Fig14b_op']
assert not list((root / 'plots').rglob('Fig11*.fig'))
assert not list((root / 'plots').rglob('Fig11*.png'))
assert not (out / 'width_jitter_display.mat').exists()
assert not list((root / 'plots').rglob('Fig7*.fig'))
assert not list((root / 'plots').rglob('Fig7*.png'))
figures = [{'stem': s, 'png_sha256': sha(out / (s + '.png')),
            'fig_sha256': sha(out / (s + '.fig'))} for s in stems]
programs = []
for p in sorted(out.glob('phase_z*.csv')):
    with p.open(encoding='utf-8-sig', newline='') as f:
        rows = list(csv.DictReader(f))
    assert len(rows) == 1600
    commands = np.array([[float(r['commandXDeg']), float(r['commandYDeg'])] for r in rows])
    gradients = np.array([[float(r['gradientXOverK0']), float(r['gradientYOverK0'])] for r in rows])
    pistons = np.array([float(r['pistonRad']) for r in rows])
    assert np.max(np.abs(commands)) <= 15.414 + 1e-12
    assert np.max(np.abs(commands / .021 - np.rint(commands / .021))) < 1e-9
    assert np.max(np.abs(gradients - np.sin(np.deg2rad(commands)))) < 1e-12
    assert np.all((pistons >= 0) & (pistons < 2 * np.pi))
    programs.append({'file': p.name, 'units': len(rows), 'sha256': sha(p)})
assert len(programs) == 9
luts = sorted(out.glob('phase_lut_z*.mat'))
assert len(luts) == 9
for p in luts:
    v = loadmat(p, simplify_cells=True)['phaseLUT']
    assert np.size(v['electrodeCoordinateM']) == 14400
    assert np.size(v['pistonRad']) == 1600
    assert np.max(np.abs(np.diff(v['electrodeCoordinateM']) - 1e-6)) < 1e-14
docs = list(root.glob('*.md')) + [out / 'fig12_revision' / 'SIMULATION_CONDITIONS.md']
for p in docs:
    text = p.read_text(encoding='utf-8-sig')
    assert not re.search(r'\\[\[\]()]', text), p
tables = {}
for name in ['all_candidates', 'objective_comparison', 'sensitivity', 'width_jitter',
             'adaptive_widths', 'turbulence_scan', 'incidence_z7.5', 'incidence_z10',
             'validation_adaptive']:
    with (out / (name + '.csv')).open(encoding='utf-8-sig', newline='') as f:
        tables[name] = sum(1 for _ in csv.DictReader(f))
assert tables['all_candidates'] == 149 and tables['objective_comparison'] == 9
assert tables['sensitivity'] == 924 and tables['adaptive_widths'] == 155
assert tables['turbulence_scan'] == 2424
assert tables['incidence_z7.5'] == tables['incidence_z10'] == 72
assert tables['validation_adaptive'] == 155
manifest = {'status': 'passed', 'readme_unchanged': True, 'wris_unchanged': True,
            'protected_wris_files': len(protected), 'typora_documents': len(docs),
            'figures': figures, 'programs': programs, 'tables': tables,
            'sources': {str(p.relative_to(root)): sha(p) for p in
                        sorted(list(root.glob('*.m')) + list((root / 'plots').glob('*.m')))}}
(audit / 'artifact_validation.json').write_text(json.dumps(manifest, ensure_ascii=False, indent=2), encoding='utf-8')
print(json.dumps({k: v for k, v in manifest.items() if k not in {'sources', 'figures', 'programs'}}, ensure_ascii=False))
