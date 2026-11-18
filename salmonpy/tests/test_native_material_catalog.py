import json
from pathlib import Path
from tempfile import TemporaryDirectory
import unittest

class NativeMaterialCatalogTests(unittest.TestCase):
    def load(self,folder,materials):
        from salmonpy.examples._tddft_materials import load_materials
        path=Path(folder)/'materials.json'
        path.write_text(json.dumps({'schema':'salmonpy.native-tddft-materials/1','materials':materials}))
        return load_materials(path)

    def test_material_local_settings_and_paths(self):
        with TemporaryDirectory() as d:
            values=self.load(d,{'oxide':{'bundle':'oxide-gs','input_file':'response.inp',
                'decomposition':[1,2,1,1,1],'points_per_group':2},'metal':{'bundle':'metal-gs'}})
            self.assertEqual(values['oxide'].bundle,Path(d).resolve()/'oxide-gs')
            self.assertEqual(values['oxide'].ranks,2)
            self.assertEqual(values['oxide'].points_per_group,2)
            self.assertEqual(values['metal'].input_file,'rt.inp')
            self.assertEqual(values['metal'].decomposition,(1,1,1,1,1))
            with self.assertRaises(AttributeError):values['oxide'].points_per_group=3

    def test_invalid_entries_and_duplicate_ids_fail_before_native_work(self):
        with TemporaryDirectory() as d:
            for entry in ({'bundle':'gs','typo':2},{'bundle':'gs','decomposition':[1,True,1,1,1]},
                          {'bundle':'gs','input_file':'../rt.inp'},{'bundle':'gs','points_per_group':0}):
                with self.assertRaises(ValueError):self.load(d,{'oxide':entry})
            with self.assertRaises(ValueError):self.load(d,{'sbe-point':{'bundle':'gs'}})
            p=Path(d)/'materials.json';p.write_text('{"schema":"salmonpy.native-tddft-materials/1","materials":{"oxide":{"bundle":"a"},"oxide":{"bundle":"b"}}}')
            from salmonpy.examples._tddft_materials import load_materials
            with self.assertRaisesRegex(ValueError,'duplicate'):load_materials(p)
