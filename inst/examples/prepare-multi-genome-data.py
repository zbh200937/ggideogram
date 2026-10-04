"""Regenerate the multi-genome example tables from public sequence records.

Python 3.9+, standard library only. Run from the repository root.
FASTA is streamed to retain only chromosome headers. Ensembl gene coordinates
are 1-based closed; LASTZ_NET records retain the two returned genomic intervals.
Versioned downloads are cached in work/multi-genome-source. REST assembly
records and selected alignment endpoints are replayed from the bundled snapshot.
"""
import concurrent.futures
import csv
import gzip
import hashlib
import io
import json
import pathlib
import re
import urllib.request

work = pathlib.Path('work/multi-genome-source')
work.mkdir(parents=True, exist_ok=True)
out = pathlib.Path('inst/extdata')
out.mkdir(parents=True, exist_ok=True)
snapshot = json.loads((out / 'multi-genome-rest-snapshot.json').read_text())


def download(name, url):
    path = work / name
    if not path.exists():
        with urllib.request.urlopen(url, timeout=60) as response:
            path.write_bytes(response.read())
    return path


def table(name, fields, rows):
    handle = io.StringIO()
    writer = csv.writer(handle, delimiter='\t', lineterminator='\n')
    writer.writerow(fields)
    writer.writerows(rows)
    content = handle.getvalue().encode('utf-8')
    if hashlib.sha256(content).hexdigest() != snapshot['tables_sha256'][name]:
        raise ValueError(f'Regenerated {name} differs from the documented source snapshot')
    (out / name).write_bytes(content)


# Populus trichocarpa, two haplotype assemblies from Gao et al. (2025).
poplar = []
for hap, accession in [('A', 'GWHERCL00000000'), ('B', 'GWHERCM00000000')]:
    path = work / (accession + '-headers.txt')
    url = f'https://download.cncb.ac.cn/gwh/Plants/Populus_trichocarpa_Pt_hap_{accession}/{accession}.genome.fasta.gz'
    if not path.exists():
        headers = []
        with urllib.request.urlopen(url, timeout=60) as response:
            with gzip.GzipFile(fileobj=response) as stream:
                for line in stream:
                    if line.startswith(b'>'):
                        headers.append(line.decode().strip())
        path.write_text('\n'.join(headers) + '\n')
    for line in path.read_text().splitlines():
        original = re.search(r'OriSeqID=(\S+)', line).group(1)
        chr = re.fullmatch(r'Chr(\d+)' + hap, original).group(1)
        length = int(re.search(r'Len=(\d+)', line).group(1))
        seq = line.split()[0].lstrip('>')
        poplar.append([f'Ptr_{hap}', accession, chr, 0, length, hap, chr,
                       chr + hap, seq, original])
table('poplar-haplotypes.tsv', ['Genome', 'Assembly', 'Chr', 'Start', 'End',
      'Haplotype', 'Homolog', 'Label', 'Accession', 'Original'], poplar)

# Triticum aestivum Chinese Spring: all 21 nuclear chromosomes.
url = 'https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/018/294/505/GCF_018294505.1_IWGSC_CS_RefSeq_v2.1/GCF_018294505.1_IWGSC_CS_RefSeq_v2.1_assembly_report.txt'
wheat = []
for line in download('wheat-assembly-report.txt', url).read_text().splitlines():
    if line.startswith('#'):
        continue
    values = line.split('\t')
    if values[3] != 'Chromosome' or values[7] != 'Primary Assembly':
        continue
    name = values[2]
    wheat.append([name[-1], 'CSv2.1', name[:-1], 0, int(values[8]),
                  name[-1], name[:-1], name, values[6]])
table('wheat-subgenomes.tsv', ['Genome', 'Assembly', 'Chr', 'Start', 'End',
      'Subgenome', 'Homolog', 'Label', 'Accession'], wheat)


def rest(name, path):
    record = snapshot['responses'][name]
    if record['request'] != path:
        raise ValueError(f'Request for {name} differs from the bundled REST snapshot')
    return record['data']


species = [('rice', 'oryza_sativa', 'IRGSP-1.0'),
           ('wild', 'oryza_rufipogon', 'OR_W1943')]
assemblies = {}
for genome, species_name, expected in species:
    assembly = rest(genome + '-assembly', '/info/assembly/' + species_name)
    if assembly['assembly_name'] != expected:
        raise ValueError(f'Assembly changed: {assembly["assembly_name"]}; expected {expected}')
    assemblies[genome] = assembly
rice = []
for genome, species_name, expected in species:
    chromosome = next(v for v in assemblies[genome]['top_level_region'] if v['name'] == '1')
    rice.append([genome, expected, '1', 0, chromosome['length'], '1',
                 'O. sativa' if genome == 'rice' else 'O. rufipogon', species_name])
table('rice-comparison-karyotype.tsv', ['Genome', 'Assembly', 'Chr', 'Start', 'End',
      'Homolog', 'Label', 'Species'], rice)

genes = []
for genome, species_name, expected in species:
    path = work / (genome + '-chr1-gff-genes.tsv')
    fields = ['Genome', 'Assembly', 'Chr', 'Start', 'End', 'Strand', 'Gene']
    if not path.exists():
        name = 'Oryza_sativa' if genome == 'rice' else 'Oryza_rufipogon'
        url = f'https://ftp.ensemblgenomes.ebi.ac.uk/pub/plants/release-62/gff3/{species_name}/{name}.{expected}.62.gff3.gz'
        rows = []
        with urllib.request.urlopen(url, timeout=60) as response:
            with gzip.GzipFile(fileobj=response) as stream:
                for line in stream:
                    if line.startswith(b'#'):
                        continue
                    values = line.decode().strip().split('\t')
                    if len(values) != 9 or values[0] != '1' or values[2] != 'gene':
                        continue
                    attrs = dict(part.split('=', 1) for part in values[8].split(';') if '=' in part)
                    if attrs.get('biotype') == 'protein_coding':
                        rows.append([genome, expected, '1', values[3], values[4], values[6],
                                     attrs['ID'].removeprefix('gene:')])
        with path.open('w') as handle:
            writer = csv.writer(handle, delimiter='\t', lineterminator='\n')
            writer.writerow(fields)
            writer.writerows(rows)
    with path.open() as handle:
        for row in csv.DictReader(handle, delimiter='\t'):
            genes.append([row['Genome'], row['Assembly'], row['Chr'], int(row['Start']),
                          int(row['End']), row['Strand'], row['Gene']])
table('rice-comparison-genes.tsv', fields, sorted(genes, key=lambda v: (v[0], v[3])))

jobs = []
queries = [
    ('rice-alignment-1000000', 1000000, 1999999),
    ('rice-alignment-7000000', 7000000, 7999999),
    ('rice-alignment-14000000-14000000-14199999', 14000000, 14199999),
    ('rice-alignment-22000000-22200000-22399999', 22200000, 22399999)]
for name, lo, hi in queries:
    jobs.append((name,
                 f'/alignment/region/oryza_sativa/1:{lo}-{hi}?compara=plants;method=LASTZ_NET;species_set=oryza_sativa;species_set=oryza_rufipogon', f'{lo}-{hi}'))


def fetch(job):
    name, path, genome = job
    data = rest(name, path)
    print(name, len(data), flush=True)
    return genome, name, data


blocks = []
alignments = {}
with concurrent.futures.ThreadPoolExecutor(max_workers=2) as pool:
    for genome, name, data in pool.map(fetch, jobs):
        candidates = alignments.setdefault(genome, [])
        for block in data:
            a = next(v for v in block['alignments'] if v['species'] == 'oryza_sativa')
            b = next(v for v in block['alignments'] if v['species'] == 'oryza_rufipogon')
            if a['seq_region'] == b['seq_region'] == '1':
                candidates.append((a, b))
for query, candidates in alignments.items():
    candidates.sort(key=lambda pair: pair[0]['end'] - pair[0]['start'], reverse=True)
    for a, b in candidates[:2]:
        blocks.append(['rice', 'IRGSP-1.0', '1', a['start'], a['end'],
            'wild', 'OR_W1943', '1', b['start'], b['end'],
            '+' if a['strand'] == b['strand'] else '-', query])
table('rice-comparison-blocks.tsv', ['Genome1', 'Assembly1', 'Chr1', 'Start1', 'End1',
      'Genome2', 'Assembly2', 'Chr2', 'Start2', 'End2', 'Orientation', 'Query'], blocks)
print('Saved', len(poplar), 'poplar chromosomes,', len(wheat), 'wheat chromosomes,',
      len(genes), 'rice genes and', len(blocks), 'alignment blocks')
