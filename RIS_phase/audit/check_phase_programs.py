"""Independently audit every saved realized LC program against README limits."""
from pathlib import Path
import json
import numpy as np
import h5py

root = Path(__file__).resolve().parents[1]
records = []
bound = np.sin(np.deg2rad(15.414))
for path in sorted((root/'results'/'finite_width').glob('z_*/*.mat')):
    with h5py.File(path, 'r') as f:
        if 'item' not in f:
            continue
        p = f['item/phase']
        command = np.asarray(p['commandDeg'])
        gradient = np.asarray(p['gradientOverK0'])
        piston = np.asarray(p['phaseAtReferenceRad'])
        curvature = np.asarray(p['phaseCurvature'])
        step_error = np.max(np.abs(command/.021-np.rint(command/.021)))
        gradient_error = np.max(np.abs(gradient-np.sin(np.deg2rad(command))))
        max_increment = 2*np.pi/532e-9*1e-6*np.max(np.abs(gradient))
        passed = (command.shape == (2,1600) and np.max(np.abs(command)) <= 15.414+1e-12
                  and step_error < 1e-9 and gradient_error < 1e-12
                  and np.max(np.abs(gradient)) <= bound+1e-12
                  and np.all(piston >= 0) and np.all(piston < 2*np.pi)
                  and np.max(np.abs(curvature)) == 0 and max_increment <= np.pi+1e-12)
        records.append({'file':str(path.relative_to(root)), 'passed':bool(passed),
                        'max_command_deg':float(np.max(np.abs(command))),
                        'command_grid_error':float(step_error),
                        'gradient_identity_error':float(gradient_error),
                        'max_phase_increment_per_1um_rad':float(max_increment)})
result = {'program_count':len(records),'passed':all(r['passed'] for r in records),
          'scope':'affine phase-only envelope, 1 um electrode centers, 0.021 degree command grid; no voltage/fringing calibration',
          'records':records}
(root/'audit'/'phase_program_validation.json').write_text(json.dumps(result,indent=2),encoding='utf-8')
assert len(records) == 93 and result['passed']
print(f"{len(records)} realized phase programs passed independent device-limit audit.")
