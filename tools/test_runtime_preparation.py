"""Check pinned bundled-source preparation; no build, download or GS job."""
from pathlib import Path
import subprocess
import tempfile
import unittest
from prepare_salmon_runtime import prepare


class PreparationTests(unittest.TestCase):
    def test_nested_unrelated_git_repository(self):
        with tempfile.TemporaryDirectory() as tmp:
            root=Path(tmp)
            subprocess.run(['git','init','-q',str(root)],check=True)
            dest=prepare(root/'runtime')
            self.assertTrue((dest/'src/rt/tddft_worker.f90').is_file())
            self.assertTrue((dest/'src/maxwell/lorentz_drude.f90').is_file())
            self.assertFalse((root/'src').exists())
            self.assertFalse((dest/'.git').exists())

    def test_existing_destination_preserved(self):
        with tempfile.TemporaryDirectory() as tmp:
            dest=Path(tmp)/'runtime';dest.mkdir();(dest/'keep').write_text('existing data')
            with self.assertRaises(FileExistsError):prepare(dest)
            self.assertEqual((dest/'keep').read_text(),'existing data')


if __name__=='__main__':unittest.main()
