"""Regenerate the bundled Arabidopsis MCScanX subsets with default parameters.

Requires Python 3 and an installed MCScanX executable.
Run from the repository root: python3 inst/examples/prepare-synteny-data.py
  --mcscanx /path/to/MCScanX
"""
import argparse
import csv
import pathlib
import re
import subprocess
import urllib.request

parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument('--mcscanx', default='MCScanX')
parser.add_argument('--work', default='work/synteny-source')
args = parser.parse_args()
work = pathlib.Path(args.work).resolve()
work.mkdir(parents=True, exist_ok=True)
source = 'https://raw.githubusercontent.com/wyp1125/MCScanX/0956cb8f900c152e2be8ba3196829e50ab454b94/data/'
for name in ('at.gff', 'at.blast'):
    path = work / name
    if not path.exists():
        with urllib.request.urlopen(source + name, timeout=60) as response:
            path.write_bytes(response.read())
subprocess.run([args.mcscanx, '-a', str(work / 'at')], check=True)
genes = {}
for line in (work / 'at.gff').read_text().splitlines():
    chr, gene, start, end = line.split()
    genes[gene] = (chr, int(start), int(end))
blocks = []
for line in (work / 'at.collinearity').read_text().splitlines():
    if line.startswith('## Alignment'):
        match = re.search(r'Alignment (\d+): score=(\S+) e_value=(\S+) N=(\d+) (\S+)&(\S+) (plus|minus)', line)
        block = dict(zip(('Block', 'Score', 'Evalue', 'N', 'Chr1', 'Chr2', 'Orientation'), match.groups()))
        block['pairs'] = []
        blocks.append(block)
    elif re.match(r'\s*\d+-\s*\d+:', line):
        block['pairs'].append(line.split(':', 1)[1].split()[:2])
selected = [b for b in blocks if b['Block'] in {'15', '74', '75', '76', '77', '78', '79', '80', '81', '82'}]
out = pathlib.Path('inst/extdata')
out.mkdir(parents=True, exist_ok=True)
with (out / 'arabidopsis-synteny-blocks.tsv').open('w') as handle:
    writer = csv.writer(handle, delimiter='\t', lineterminator='\n')
    writer.writerow(('Block', 'Chr1', 'Start1', 'End1', 'Chr2', 'Start2', 'End2', 'Orientation', 'N', 'Score'))
    for b in selected:
        left = [genes[a] for a, _ in b['pairs']]
        right = [genes[z] for _, z in b['pairs']]
        writer.writerow((b['Block'], b['Chr1'].removeprefix('at'), min(v[1] for v in left), max(v[2] for v in left),
            b['Chr2'].removeprefix('at'), min(v[1] for v in right), max(v[2] for v in right),
            '+' if b['Orientation'] == 'plus' else '-', b['N'], b['Score']))
with (out / 'arabidopsis-synteny-pairs.tsv').open('w') as handle:
    writer = csv.writer(handle, delimiter='\t', lineterminator='\n')
    writer.writerow(('Block', 'Gene1', 'Chr1', 'Start1', 'End1', 'Gene2', 'Chr2', 'Start2', 'End2', 'Orientation'))
    for b in selected:
        for a, z in b['pairs']:
            x, y = genes[a], genes[z]
            writer.writerow((b['Block'], a, x[0].removeprefix('at'), x[1], x[2],
                z, y[0].removeprefix('at'), y[1], y[2], '+' if b['Orientation'] == 'plus' else '-'))
