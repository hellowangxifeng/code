"""Read original source/artifact inventory and record reproducible audit evidence."""
from pathlib import Path
import csv
import hashlib
import re
import json
import scipy.io
import h5py
import pymupdf
from PIL import Image, ImageOps, ImageDraw

root = Path(__file__).resolve().parents[1]
workspace = root.parent
out = root / 'audit'
records = []
notes = []
thumbs = []
for folder in [root, workspace / 'W_RIS']:
    for path in sorted(folder.rglob('*')):
        if not path.is_file() or any(p in {'cache', 'finite_width', 'audit'} for p in path.relative_to(folder).parts):
            continue
        data = path.read_bytes()
        record = {'path': str(path.relative_to(workspace)), 'bytes': len(data),
                  'sha256': hashlib.sha256(data).hexdigest(), 'type': path.suffix}
        records.append(record)
        if path.suffix in {'.m', '.md', '.py', '.txt'}:
            text = data.decode('utf-8-sig', errors='replace')
            lines = text.splitlines()
            headings = [line for line in lines if re.match(r'^function|^classdef|^#{1,4} |^%[A-Z_]{3,}', line)]
            calls = sorted(set(re.findall(r'\b((?:rp_|wris_|bsf_)[A-Za-z_0-9]+)\s*\(', text)))
            notes.append({'path': record['path'], 'lines': len(lines), 'headings': headings,
                          'project_calls': calls})
        elif path.suffix == '.mat':
            try:
                variables = scipy.io.whosmat(path)
            except NotImplementedError:
                with h5py.File(path, 'r') as f:
                    variables = [(key, str(f[key].attrs.get('MATLAB_class', b''))) for key in f if key != '#refs#']
            except Exception as error:
                variables = str(error)
            notes.append({'path': record['path'], 'mat_variables': variables})
        elif path.suffix == '.png':
            with Image.open(path) as image:
                image.load()
                record['dimensions'] = f'{image.width}x{image.height}'
                if 'comparison_r' in path.name or path.name == 'precision_overview.png':
                    tile = Image.new('RGB', (440, 290), 'white')
                    preview = ImageOps.contain(image.convert('RGB'), (440, 265))
                    tile.paste(preview, ((440-preview.width)//2, 20))
                    ImageDraw.Draw(tile).text((5, 3), path.name, fill='black')
                    thumbs.append(tile)
        elif path.suffix == '.pdf':
            doc = pymupdf.open(path)
            record['pages'] = len(doc)
            (out / f'{path.stem}_fulltext.txt').write_text('\n\f\n'.join(p.get_text() for p in doc), encoding='utf-8')
            notes.append({'path': record['path'], 'pdf_pages': len(doc)})
if thumbs:
    montage = Image.new('RGB', (440*4, 290*((len(thumbs)+3)//4)), 'white')
    for i, tile in enumerate(thumbs):
        montage.paste(tile, ((i%4)*440, (i//4)*290))
    montage.save(out / 'existing_bsf_montage.png')
for source in [workspace/'yimiao'/'4yimiao.pdf']:
    doc = pymupdf.open(source)
    (out / '4yimiao_fulltext.txt').write_text('\n\f\n'.join(p.get_text() for p in doc), encoding='utf-8')
    for page_id in [0, 3, 4, 12]:
        if page_id < len(doc):
            doc[page_id].get_pixmap(matrix=pymupdf.Matrix(1.3, 1.3)).save(out/f'4yimiao_page_{page_id+1}.png')
with (out/'source_inventory.csv').open('w', encoding='utf-8-sig', newline='') as stream:
    writer = csv.DictWriter(stream, fieldnames=['path','bytes','sha256','type','dimensions','pages'])
    writer.writeheader(); writer.writerows(records)
(out/'source_review.json').write_text(json.dumps(notes, ensure_ascii=False, indent=2, default=str), encoding='utf-8')
print(f'Read {len(records)} files; {len(notes)} source, PDF, or MAT records; inspected {len(thumbs)} old figures.')
