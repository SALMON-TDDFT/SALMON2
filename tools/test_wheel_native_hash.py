"""Check native hashes in a wheel, before or after dependency repair."""
import hashlib
from pathlib import Path
import sys
import zipfile


def check(wheel):
    with zipfile.ZipFile(wheel) as archive:
        names = [n for n in archive.namelist()
                 if n.startswith('salmonpy/_native/') and n.endswith(('.so', '.dylib', '.dll'))]
        if not names:
            raise AssertionError('wheel contains no native library')
        for name in names:
            expected = archive.read(name + '.sha256').decode().strip()
            actual = hashlib.sha256(archive.read(name)).hexdigest()
            if expected != actual:
                raise AssertionError('installed native checksum mismatch: ' + name)
    print(Path(wheel).name, 'native hashes PASS')


if __name__ == '__main__':
    for wheel in sys.argv[1:]:
        check(wheel)
