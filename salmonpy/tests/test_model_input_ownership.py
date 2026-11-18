"""The actual Fortran model record must own independently allocated input values."""
from pathlib import Path
import re
import unittest

class ModelInputOwnershipTests(unittest.TestCase):
    def test_ownership_source_is_packaged(self):
        root=Path(__file__).resolve().parents[2]
        # Native runtime validation compiles the Fortran round-trip driver.
        source=root/'runtime/sources/standard/src/io/model_inputs.f90'
        if not source.is_file():
            import importlib.resources
            source=Path(str(importlib.resources.files('salmonpy').joinpath('_runtime_source/runtime/sources/standard/src/io/model_inputs.f90')))
        self.assertTrue(source.is_file(),'model-owned initialization input record missing')

    def test_supported_initialization_sections_have_owned_fields(self):
        from importlib.resources import files
        root=Path(__file__).resolve().parents[2]/'runtime/sources/standard/src/io'
        if not root.is_dir():
            root=Path(str(files('salmonpy').joinpath('_runtime_source/runtime/sources/standard/src/io')))
        sections={'calculation','control','units','parallel','system','pseudo','functional',
                  'rgrid','kgrid','tgrid','propagation','scf','emfield','analysis','poisson',
                  'ewald','atomic_coor','atomic_red_coor','code','jellium','unfolding'}
        expected={'iflag_atom_coor'};section=None
        for line in (root/'salmon_global.f90').read_text().splitlines():
            marker=re.match(r'\s*!!\s*&([\w]+)',line)
            if marker:section=marker.group(1)
            clean=line.split('!')[0]
            if section in sections and '::' in clean:
                expected.update(re.findall(r'(\w+)\s*(?:\([^)]*\))?\s*(?:,|$)',clean.split('::',1)[1]))
        text=(root/'model_inputs.f90').read_text()
        for name in expected:
            self.assertRegex(text,r'(?i)self%'+name+r'\s*=\s*'+name+r'\b')
            self.assertRegex(text,r'(?i)\b'+name+r'\s*=\s*self%'+name+r'\b')
