"""The actual Fortran model record must own independently allocated input values."""
from pathlib import Path
import re
import unittest

class ModelInputOwnershipTests(unittest.TestCase):

    def test_supported_initialization_sections_have_owned_fields(self):
        from importlib.resources import files
        root=Path(__file__).resolve().parents[2]/'runtime/sources/solver/src/io'
        if not root.is_dir():
            root=Path(str(files('salmonpy').joinpath('_runtime_source/runtime/sources/solver/src/io')))
        self.assertIn('type, public :: s_input_parser',(root/'inputoutput.f90').read_text())
        sections={'calculation','control','units','parallel','system','pseudo','functional',
                  'rgrid','kgrid','tgrid','propagation','scf','emfield','analysis','poisson',
                  'ewald','atomic_coor','atomic_red_coor','code','jellium','unfolding'}
        expected={'iflag_atom_coor'};section=None
        for line in (root/'salmon_global.f90').read_text().splitlines():
            # Audit input declarations, not the derived conversion view below.
            if re.match(r'\s*end\s+module\s+salmon_global\b',line,re.I):break
            marker=re.match(r'\s*!!\s*&([\w]+)',line)
            if marker:section=marker.group(1)
            clean=line.split('!')[0]
            if section in sections and '::' in clean:
                expected.update(re.findall(r'(\w+)\s*(?:\([^)]*\))?\s*(?:,|$)',clean.split('::',1)[1]))
        text=(root/'model_inputs.f90').read_text()
        for name in expected:
            self.assertRegex(text,r'(?i)self%'+name+r'\s*=\s*'+name+r'\b')
            self.assertRegex(text,r'(?i)\b'+name+r'\s*=\s*self%'+name+r'\b')
