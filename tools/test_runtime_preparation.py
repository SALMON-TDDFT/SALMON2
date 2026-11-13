"""Integration checks for external source preparation; requires pinned upstream Git checkout."""
import argparse
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
            dest=prepare(UPSTREAM,root/'runtime')
            self.assertTrue((dest/'src/rt/tddft_worker.f90').is_file())
            self.assertTrue((dest/'src/maxwell/lorentz_drude.f90').is_file())
            self.assertFalse((root/'src').exists())
            self.assertFalse((dest/'.git').exists())

    def test_existing_destination_preserved(self):
        with tempfile.TemporaryDirectory() as tmp:
            dest=Path(tmp)/'runtime';dest.mkdir();(dest/'keep').write_text('existing data')
            with self.assertRaises(FileExistsError):prepare(UPSTREAM,dest)
            self.assertEqual((dest/'keep').read_text(),'existing data')


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('upstream')
    UPSTREAM=Path(parser.parse_args().upstream).resolve()
    unittest.main(argv=['test_runtime_preparation'])
