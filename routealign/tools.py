"""Install operator-provided binaries and fetch the pinned TEMPL source snapshot."""
import argparse
import hashlib
import json
from pathlib import Path
import shutil
from routealign import settings
from routealign.prepare import download

TEMPL_COMMIT = '727fd306e4066aa39b05a6e51e7833c86d90df98'
TEMPL_FILES = {
    'templ_pipeline/core/mcs.py': 'caf216b39a87ce928cf61921269720ac62103bd7791fc3f2ec4aead1e41ee109',
    'templ_pipeline/core/scoring.py': 'd538ecc50fc4d471a3388eb63e32ccc547ca1a0ff4e1827a5748af6877ce3d76',
    'templ_pipeline/core/chemistry.py': 'df9f6bd0c3e295349984fa78bb79f92c5510d776b251ba70096163b08f034dd4',
}


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--workspace', type=Path, default=settings.workspace())
    p.add_argument('--smina', type=Path)
    p.add_argument('--lsalign', type=Path)
    p.add_argument('--templ', action='store_true')
    p.add_argument('--smina-download', action='store_true')
    p.add_argument('--published-baselines', action='store_true')
    args = p.parse_args(); root = args.workspace.resolve(); tools = root / 'tools'; tools.mkdir(parents=True, exist_ok=True)
    for source, name, expected in [(args.smina, 'smina', 'ffe5e1e78c947f76b0df8805e2c54383d0bbaf2e827a633b643a708cf682a958'),
                                    (args.lsalign, 'LSalign', '598afc0cd5f8994d9aba159b45315de12fa5ecdff4a226de07b281959c9bf1bc')]:
        if source:
            if hashlib.sha256(source.read_bytes()).hexdigest() != expected:
                raise RuntimeError(f'{name} differs from the binary used for the reported experiments')
            target = tools / name
            if source.resolve() != target.resolve(): shutil.copyfile(source, target)
            target.chmod(target.stat().st_mode | 0o111)
    if args.smina_download:
        file = tools / 'smina'
        download('https://sourceforge.net/projects/smina/files/smina.static/download', file)
        expected = 'ffe5e1e78c947f76b0df8805e2c54383d0bbaf2e827a633b643a708cf682a958'
        if hashlib.sha256(file.read_bytes()).hexdigest() != expected:
            raise RuntimeError('Downloaded Smina differs from the reported binary; use --smina with the matching binary')
        file.chmod(file.stat().st_mode | 0o111)
    if args.templ:
        contents = {}
        for relative in [*TEMPL_FILES, 'LICENSES/MIT.txt']:
            target = tools / 'templ' / relative
            download(f'https://raw.githubusercontent.com/fulopjoz/templ-pipeline/{TEMPL_COMMIT}/{relative}', target)
            text = target.read_text(encoding='utf-8')
            if relative in TEMPL_FILES and hashlib.sha256(text.encode()).hexdigest() != TEMPL_FILES[relative]:
                raise RuntimeError('TEMPL source checksum mismatch')
            contents[relative] = text
        (tools / 'templ_source.json').write_text(json.dumps({'commit': TEMPL_COMMIT, 'contents': contents}), encoding='utf-8')
    if args.published_baselines:
        base = 'https://raw.githubusercontent.com/iktos/engine/747d842/'
        for relative, expected in [
            ('3d-align/results/3DAlign_benchmark_results.csv', '1062dbf414525b0d0093ce59c4a85ddeefbf338ad8ead8c7c0c1e12d99064be6'),
            ('3d-dock/results/3DDock_benchmark_results.csv', 'fcd080d8cfebe1801919f0ec5a4fda188cc08dfeb3bfeeaa0d29136e2b2ec0b3')]:
            file = root / 'data/published_baselines' / Path(relative).name
            download(base + relative, file)
            if hashlib.sha256(file.read_bytes()).hexdigest() != expected:
                raise RuntimeError('Published baseline CSV checksum mismatch')
        download(base + 'LICENSE', root / 'data/published_baselines/LICENSE')
    print(json.dumps({'status': 'PASS'}))


if __name__ == '__main__':
    main()
