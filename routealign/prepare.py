"""Prepare a self-contained workspace from the official AlignDockBench dataset."""
import argparse
import csv
import hashlib
import json
from pathlib import Path
import shutil
import tarfile
import urllib.request
from routealign import settings

DATA_URL = 'https://zenodo.org/api/records/15395813/files/AlignDockBench.tar.gz/content'
DATA_MD5 = 'ca8d9076efa4b4661c74c47decee2ff1'


def download(url, destination):
    destination.parent.mkdir(parents=True, exist_ok=True)
    partial = destination.with_suffix(destination.suffix + '.part')
    start = partial.stat().st_size if partial.exists() else 0
    request = urllib.request.Request(url, headers={'User-Agent': 'RouteAlign/1.0', **({'Range': f'bytes={start}-'} if start else {})})
    with urllib.request.urlopen(request, timeout=60) as response:
        mode = 'ab' if start and response.status == 206 else 'wb'
        with partial.open(mode) as out:
            shutil.copyfileobj(response, out, 1024 * 1024)
    partial.replace(destination)


def main():
    p = argparse.ArgumentParser(description=__doc__)
    g = p.add_mutually_exclusive_group(required=True)
    g.add_argument('--dataset', type=Path, help='Existing extracted AlignDockBench directory')
    g.add_argument('--archive', type=Path, help='Original AlignDockBench.tar.gz')
    g.add_argument('--download', action='store_true')
    p.add_argument('--workspace', type=Path, default=settings.workspace())
    args = p.parse_args(); root = args.workspace.resolve()
    for f in ['cohort.csv', 'protocol.json', 'roshambo_protocol.md']:
        src = settings.PACKAGE_ROOT / 'config' / f; dest = root / 'config' / f
        dest.parent.mkdir(parents=True, exist_ok=True)
        if src.resolve() != dest.resolve():
            shutil.copyfile(src, dest)
    rows = list(csv.DictReader((root / 'config/cohort.csv').open()))
    needed = set()
    for row in rows:
        if row['final_primary_eligible'] != 'True': continue
        needed.add(f"{row['family']}/{row['template_id']}/ligand.sdf")
        needed.add(f"{row['family']}/{row['template_id']}/protein.pdb")
        needed.add(f"{row['family']}/{row['query_id']}/ligand.sdf")
    destination = root / 'data/AlignDockBench'
    archive = args.archive
    if args.download:
        archive = root / 'cache/AlignDockBench.tar.gz'
        if not archive.exists(): download(DATA_URL, archive)
    if archive:
        digest = hashlib.md5()
        with archive.open('rb') as f:
            for block in iter(lambda: f.read(1024 * 1024), b''): digest.update(block)
        if digest.hexdigest() != DATA_MD5:
            raise RuntimeError('Dataset archive checksum differs from the official record')
        found = set()
        with tarfile.open(archive) as tar:
            for member in tar:
                if not member.isfile(): continue
                parts = Path(member.name).parts
                if 'AlignDockBench' not in parts: continue
                relative = '/'.join(parts[parts.index('AlignDockBench') + 1:])
                if relative not in needed: continue
                target = destination / relative; target.parent.mkdir(parents=True, exist_ok=True)
                with tar.extractfile(member) as source, target.open('wb') as output:
                    shutil.copyfileobj(source, output)
                found.add(relative)
        if found != needed: raise RuntimeError('Archive lacks required cohort structures')
    else:
        source = args.dataset.resolve()
        for relative in needed:
            target = destination / relative; target.parent.mkdir(parents=True, exist_ok=True)
            if (source / relative).resolve() != target.resolve(): shutil.copyfile(source / relative, target)
    for row in rows:
        if row['final_primary_eligible'] != 'True': continue
        for role in ['query', 'template']:
            file = destination / row['family'] / row[role + '_id'] / 'ligand.sdf'
            assert hashlib.sha256(file.read_bytes()).hexdigest() == row[role + '_sdf_sha256'], (row['family'], role)
    print(json.dumps({'status': 'PASS', 'cohort_pairs': 44, 'structure_files': len(needed)}))


if __name__ == '__main__':
    main()
