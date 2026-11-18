from pathlib import Path
import unittest
from importlib.resources import files

class ParserOwnershipTests(unittest.TestCase):
    def test_owned_parser_entry_is_packaged(self):
        root=Path(__file__).resolve().parents[2]/'runtime/sources/standard/src/io'
        if not root.is_dir():root=Path(str(files('salmonpy').joinpath('_runtime_source/runtime/sources/standard/src/io')))
        self.assertIn('type, public :: s_input_parser',(root/'inputoutput.f90').read_text())
