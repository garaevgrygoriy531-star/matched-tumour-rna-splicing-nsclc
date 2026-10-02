#!/usr/bin/env python3
"""Filter the copied predictor variant universe to the frozen Strict68 genes."""
from __future__ import annotations
import csv, hashlib, json, sys
from pathlib import Path

def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open('rb') as fh:
        for block in iter(lambda: fh.read(8 * 1024 * 1024), b''):
            h.update(block)
    return h.hexdigest()

def main() -> int:
    if len(sys.argv) != 2:
        raise SystemExit('usage: strict68_filter_m4_m5_universe.py STRICT_ROOT')
    root = Path(sys.argv[1])
    panel = root / '00_project_freeze/STRICT68_GENE_SET.tsv'
    path = root / 'v2.3/results/inputs/m4_m5/V23_M4_M5_VARIANT_UNIVERSE.tsv'
    out = path.with_name('V23_M4_M5_VARIANT_UNIVERSE_STRICT68.tsv')
    genes = {r['gene'] for r in csv.DictReader(panel.open(encoding='utf-8'), delimiter='\t')}
    with path.open(encoding='utf-8', newline='') as inf:
        reader = csv.DictReader(inf, delimiter='\t')
        if not reader.fieldnames or 'gene' not in reader.fieldnames:
            raise RuntimeError('predictor universe lacks gene column')
        rows = [r for r in reader if r.get('gene', '') in genes]
        fields = list(reader.fieldnames)
    with out.open('w', encoding='utf-8', newline='') as outf:
        writer = csv.DictWriter(outf, fieldnames=fields, delimiter='\t', lineterminator='\n')
        writer.writeheader(); writer.writerows(rows)
    # Replace the working path only after the filtered file is complete.
    path.replace(path.with_name('V23_M4_M5_VARIANT_UNIVERSE_BROAD74_COPY.tsv'))
    out.replace(path)
    audit = {
        'status': 'PASS', 'panel': str(panel), 'panel_sha256': sha256(panel),
        'source_broad_copy': str(path.with_name('V23_M4_M5_VARIANT_UNIVERSE_BROAD74_COPY.tsv')),
        'strict_output': str(path), 'strict_gene_count': len(genes), 'strict_variant_count': len(rows),
        'selection_rule': 'gene in frozen Strict68 panel; no outcome fields used',
        'output_sha256': sha256(path),
    }
    (root / 'v2.3/results/audit/STRICT68_M4_M5_VARIANT_UNIVERSE_AUDIT.json').write_text(json.dumps(audit, indent=2) + '\n', encoding='utf-8')
    print(json.dumps(audit, indent=2))
    return 0

if __name__ == '__main__':
    raise SystemExit(main())
