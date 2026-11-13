"""Bundle compiler runtimes, then refresh native checksum and wheel RECORD."""
import base64
import csv
import hashlib
import io
from pathlib import Path
import subprocess
import sys
import tempfile
import zipfile


def refresh_hashes(wheel):
    wheel = Path(wheel)
    with zipfile.ZipFile(wheel) as z:
        infos = z.infolist()
        contents = {i.filename: z.read(i) for i in infos}
    libs = [n for n in contents if n.startswith('salmonpy/_native/libsalmonpy_ld.')
            and n.endswith(('.dll', '.so', '.dylib'))]
    if len(libs) != 1:
        raise ValueError(f'Expected one LD library, found {libs}')
    contents[libs[0]+'.sha256'] = (hashlib.sha256(contents[libs[0]]).hexdigest()+'\n').encode()
    record, = [n for n in contents if n.endswith('.dist-info/RECORD')]
    rows = []
    for name, data in contents.items():
        if name != record:
            digest = base64.urlsafe_b64encode(hashlib.sha256(data).digest()).rstrip(b'=').decode()
            rows.append((name, 'sha256='+digest, str(len(data))))
    stream = io.StringIO(newline='')
    csv.writer(stream).writerows(rows + [(record, '', '')])
    contents[record] = stream.getvalue().encode()
    temporary = wheel.with_suffix('.tmp')
    with zipfile.ZipFile(temporary, 'w', zipfile.ZIP_DEFLATED) as z:
        for name, data in contents.items():
            z.writestr(name, data)
    temporary.replace(wheel)


def main():
    wheel = Path(sys.argv[1]).resolve(); dest = Path(sys.argv[2]).resolve()
    dest.mkdir(parents=True, exist_ok=True)
    with tempfile.TemporaryDirectory() as tmp:
        if sys.platform == 'darwin':
            cmd = ['delocate-wheel', '-w', tmp, str(wheel)]
        elif sys.platform == 'win32':
            cmd = [sys.executable, '-m', 'delvewheel', 'repair', '--analyze-existing', '-w', tmp, str(wheel)]
        else:
            cmd = ['auditwheel', 'repair', '-w', tmp, str(wheel)]
        subprocess.run(cmd, check=True)
        repaired = list(Path(tmp).glob('*.whl'))
        if len(repaired) != 1:
            raise RuntimeError('Expected one repaired wheel')
        refresh_hashes(repaired[0])
        import shutil
        shutil.copy2(repaired[0], dest / repaired[0].name)

if __name__ == '__main__':
    main()
