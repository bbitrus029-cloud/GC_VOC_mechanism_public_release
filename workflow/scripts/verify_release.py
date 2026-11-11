#!/usr/bin/env python3
"""Offline release integrity and complete frozen result coverage; never fits models."""
from pathlib import Path
import csv, hashlib, json, os, stat, sys
sys.dont_write_bytecode = True
ROOT = Path.cwd()

def require(ok, message):
    if not ok:
        raise ValueError(message)

def table(path):
    with (ROOT / path).open(newline='', encoding='utf-8') as f:
        return list(csv.DictReader(f, delimiter='\t'))

def sha(path):
    h = hashlib.sha256()
    with path.open('rb') as f:
        for chunk in iter(lambda: f.read(1024 * 1024), b''):
            h.update(chunk)
    return h.hexdigest()

def verify():
    require((ROOT / 'PUBLIC_MANIFEST.tsv').is_file(), 'Run ./verify from the release root.')
    entries = table('PUBLIC_MANIFEST.tsv')
    expected = {r['path'] for r in entries}
    require(len(expected) == len(entries), 'Duplicate public manifest path')
    observed = set()
    for base, dirs, files in os.walk(ROOT, followlinks=False):
        # An owner may later initialize Git; its metadata are not release payloads.
        if Path(base) == ROOT and '.git' in dirs:
            dirs.remove('.git')
        for name in dirs:
            d = Path(base) / name
            require(not d.is_symlink(), 'Directory symlink: ' + str(d.relative_to(ROOT)))
            require(not name.startswith('.'), 'Unexpected hidden directory: ' + name)
        for name in files:
            p = Path(base) / name
            rel = str(p.relative_to(ROOT))
            st = p.lstat()
            require(stat.S_ISREG(st.st_mode) and st.st_nlink == 1, 'Not an independent regular file: ' + rel)
            observed.add(rel)
    require(observed == expected, 'Physical file inventory differs: ' + repr(sorted(observed ^ expected)))
    for r in entries:
        p = Path(r['path'])
        require(not p.is_absolute() and '..' not in p.parts, 'Unsafe manifest path')
        q = ROOT / p
        require(q.stat().st_size == int(r['bytes']), 'Size mismatch: ' + str(p))
        if str(p) == 'PUBLIC_MANIFEST.tsv':
            require(r['sha256'] == 'SELF_HASH_EXCLUDED', 'Manifest self-hash convention')
        else:
            require(sha(q) == r['sha256'], 'SHA-256 mismatch: ' + str(p))
    orig = table('results/integration/original/hypothesis_coverage.tsv')
    a = table('results/supplementary/external_bulk/all_inferential_results.tsv')
    b = table('results/supplementary/epithelial/all_results.tsv')
    cd = table('results/integration/supplementary_tissue_protein_hypotheses.tsv')
    for name, values, n in [('original',orig,146),('A',a,17),('B',b,12),('CD',cd,13)]:
        require(len(values) == n, name + ' hypothesis coverage')
        ids = [r['hypothesis_id'] for r in values]
        require(len(ids) == len(set(ids)) and all(ids), name + ' hypothesis identities')
    fam = table('results/integration/family_reconciliation.tsv')
    require(len(fam) == 41 and len({(r['namespace'],r['family']) for r in fam}) == 41, 'Complete family inventory')
    for r in fam:
        require(r['frozen_size'] == r['observed_slots'] and r['pending'] == '0', 'Family slot/closure mismatch')
    summary = table('results/integration/coverage_summary.tsv')
    require(sum(int(r['hypotheses']) for r in summary) == 188, 'Total coverage')
    for r in table('publication/SOURCE_INDEX.tsv'):
        require(r['canonical_path'] in expected, 'Unresolved figure/table source: ' + r['canonical_path'])
    forbidden_suffixes = {'.pdf','.svg','.png','.jpg','.jpeg','.gif','.docx','.xlsx','.pptx','.ipynb','.zip','.tar','.rds'}
    require(not any(Path(p).suffix.lower() in forbidden_suffixes for p in expected),
            'Finished/embedded/container payload is outside this working-tree scope')
    require(all(p.endswith('.tsv.gz') for p in expected if p.endswith('.gz')), 'Unexpected compressed container')
    figures = [p for p in expected if p.startswith('publication/Figure') and '/source/' in p]
    from verify_current import check
    current = check(ROOT, table, require, expected)
    for r in table('workflow/config/bundled_input_aliases.tsv'):
        require(r['bundled_path'] in expected, 'Unresolved bundled input alias')
    for r in table('workflow/config/reconstruction_inputs.tsv'):
        require(not r['bundled_path'] or r['bundled_path'] in expected, 'Unresolved reconstruction input')
    for p in ROOT.rglob('*.py'):
        if '.git' not in p.relative_to(ROOT).parts:
            compile(p.read_text(), str(p), 'exec')  # No bytecode/cache files.
    return {**current, 'status':'PASS','files':len(entries),'bytes':sum(int(r['bytes']) for r in entries),
            'original_hypotheses':146,'supplementary_hypotheses':42,
            'original_families':29,'supplementary_families':12,
            'network_access':False,'scientific_recomputation':False,
            'scope':'File integrity, independent files, coverage/closure, source navigation, Python syntax; not rights certification or raw-input reconstruction.'}

if __name__ == '__main__':
    try:
        print(json.dumps(verify(),indent=2))
    except (OSError, ValueError, KeyError) as e:
        print('Verification failed: ' + str(e) + '. See REPRODUCIBILITY.md; no source-project fallback or download is attempted.',file=sys.stderr)
        sys.exit(1)
