"""Keep MATLAB PDF page boxes at the requested exact inch dimensions."""
from pathlib import Path
import hashlib
import json
import os
import sys
import pymupdf
from PIL import Image

folder = Path(sys.argv[1]).resolve()
stems = ['Fig12a_radius', 'Fig12b_aber', 'Fig12c_op', 'Fig13a_aber',
         'Fig13b_op', 'Fig14a_aber', 'Fig14b_op']
bundle = pymupdf.open()
records = []
for stem in stems:
    path = folder / (stem + '.pdf')
    doc = pymupdf.open(path)
    assert len(doc) == 1
    page = doc[0]
    # MATLAB print rounds the page box to whole points. Correct the box only;
    # leave every vector coordinate and the 10 pt text size unchanged.
    page.set_mediabox(pymupdf.Rect(0, 0, 2.67 * 72, 2.37 * 72))
    # MATLAB maps Times New Roman to the unembedded PDF base font Times-Roman.
    # Embed the requested installed TrueType font, preserving text operators.
    original_fonts = page.get_fonts()
    if any(row[3] == 'Times-Roman' for row in original_fonts):
        font_file = Path(os.environ.get('WINDIR', 'C:/Windows')) / 'Fonts' / 'times.ttf'
        assert font_file.is_file(), 'Times New Roman font file is required.'
        embedded = page.insert_font(fontname='EmbeddedTimesNewRoman', fontfile=str(font_file), set_simple=True)
        for row in original_fonts:
            if row[3] == 'Times-Roman':
                doc.update_object(row[0], doc.xref_object(embedded))
    temporary = folder / (stem + '.sized.pdf')
    doc.save(temporary)
    doc.close()
    temporary.replace(path)
    doc = pymupdf.open(path)
    page = doc[0]
    assert abs(page.rect.width / 72 - 2.67) < 1e-6
    assert abs(page.rect.height / 72 - 2.37) < 1e-6
    assert page.get_text().strip()
    fonts = [row[3] for row in page.get_fonts()]
    assert any('TimesNewRoman' in name for name in fonts)
    bundle.insert_pdf(doc)
    with Image.open(folder / (stem + '.png')) as image:
        assert image.size == (1602, 1422)
    records.append({'figure': stem, 'widthIn': page.rect.width / 72,
                    'heightIn': page.rect.height / 72, 'fonts': fonts,
                    'pngSHA256': hashlib.sha256((folder / (stem + '.png')).read_bytes()).hexdigest(),
                    'pdfSHA256': hashlib.sha256(path.read_bytes()).hexdigest()})
    doc.close()
merged = folder / 'Fig12_14_compact.pdf'
temporary = folder / 'Fig12_14_compact.tmp.pdf'
bundle.save(temporary)
temporary.replace(merged)
(folder / 'figure_exports.json').write_text(json.dumps(
    {'status': 'passed', 'frames': records, 'mergedPages': len(bundle)}, indent=2), encoding='utf-8')
print('Seven exact-size vector PDF pages, PNG exports and review bundle verified.')
