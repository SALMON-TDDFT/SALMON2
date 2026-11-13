import base64
import csv
import hashlib
import io
from pathlib import Path
import tempfile
import unittest
import zipfile
from repair_ld_wheel import refresh_hashes


class RepairTests(unittest.TestCase):
    def test_relocated_binary_and_record(self):
        for suffix in ('.dll','.so','.dylib'):
            with self.subTest(suffix=suffix), tempfile.TemporaryDirectory() as tmp:
                path=Path(tmp)/'sample.whl'
                name='salmonpy/_native/libsalmonpy_ld'+suffix
                with zipfile.ZipFile(path,'w') as z:
                    z.writestr(name,b'relocated native binary')
                    z.writestr(name+'.sha256',b'obsolete hash')
                    z.writestr('salmonpy-0.1.dist-info/RECORD','')
                refresh_hashes(path)
                with zipfile.ZipFile(path) as z:
                    self.assertEqual(z.read(name+'.sha256').decode().strip(),hashlib.sha256(z.read(name)).hexdigest())
                    rows=csv.reader(io.StringIO(z.read('salmonpy-0.1.dist-info/RECORD').decode()))
                    for file,digest,size in rows:
                        if digest:
                            self.assertEqual(int(size),len(z.read(file)))
                            actual=base64.urlsafe_b64encode(hashlib.sha256(z.read(file)).digest()).rstrip(b'=').decode()
                            self.assertEqual(digest,'sha256='+actual)

if __name__=='__main__': unittest.main()
