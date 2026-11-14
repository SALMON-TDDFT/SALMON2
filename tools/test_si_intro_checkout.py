"""Verify starter provenance survives Git's Windows-style text conversion."""
import json
import hashlib
from pathlib import Path
import shutil
import subprocess
import tempfile

root = Path(__file__).resolve().parents[1]
data = Path('salmonpy/examples/si_intro_data')
with tempfile.TemporaryDirectory(prefix='si-checkout-') as directory:
    temporary = Path(directory)
    repo = temporary/'repo'
    repo.mkdir()
    shutil.copytree(root/data, repo/data)
    shutil.copyfile(root/'.gitattributes', repo/'.gitattributes')
    def git(*args):
        subprocess.run(['git', '-C', str(repo), *args], check=True,
                       stdout=subprocess.DEVNULL, stderr=subprocess.PIPE)
    git('init', '-q')
    git('config', 'core.autocrlf', 'true')
    git('add', '.gitattributes', 'salmonpy')
    exported = temporary/'exported'
    exported.mkdir()
    git('checkout-index', '-a', '--prefix='+exported.as_posix()+'/')
    provenance = json.loads((exported/data/'provenance.json').read_bytes())
    for name, expected in provenance['files'].items():
        actual = hashlib.sha256((exported/data/name).read_bytes()).hexdigest()
        assert actual == expected, f'checkout changed immutable starter file: {name}'
print('Si input hashes preserved with core.autocrlf=true: PASS')
