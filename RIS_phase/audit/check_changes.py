"""Compare the pre-change inventory; verify protected inputs and list outputs."""
from pathlib import Path
import csv
import hashlib
import json

root = Path(__file__).resolve().parents[1]
workspace = root.parent
audit = root / 'audit'
with (audit/'source_inventory.csv').open(encoding='utf-8-sig', newline='') as f:
    before = {Path(row['path']).as_posix(): row['sha256'] for row in csv.DictReader(f)}
after = {}
for folder in [root, workspace/'W_RIS']:
    for file in folder.rglob('*'):
        if file.is_file() and not any(p in {'cache','finite_width','audit'} for p in file.relative_to(folder).parts):
            after[file.relative_to(workspace).as_posix()] = hashlib.sha256(file.read_bytes()).hexdigest()
removed = sorted(set(before)-set(after))
added = sorted(set(after)-set(before))
modified = sorted(name for name in set(before)&set(after) if before[name] != after[name])
wris_changed = [p for p in removed+modified if p.startswith('W_RIS/')]
result = {'readme_unchanged': before.get('RIS_phase/README.md') == after.get('RIS_phase/README.md'),
          'wris_unchanged': not wris_changed, 'wris_changed': wris_changed,
          'removed': removed, 'modified': modified, 'added': added}
(audit/'change_manifest.json').write_text(json.dumps(result, ensure_ascii=False, indent=2), encoding='utf-8')
print(f"README unchanged={result['readme_unchanged']}; W_RIS unchanged={result['wris_unchanged']}; removed={len(removed)}, modified={len(modified)}, added={len(added)}")
assert result['readme_unchanged'] and result['wris_unchanged']
