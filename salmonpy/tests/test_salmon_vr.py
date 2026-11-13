import json
from pathlib import Path
import tempfile
import unittest
import numpy as np
from salmonpy.outputs import read_cube
from salmonpy.salmon_vr import export_cubes


def cube(path, value=1., origin=0.):
    path.write_text(f'test\natomic units\n1 {origin} 0 0\n2 1 0 0\n2 0 1 0\n2 0 0 1\n14 4 0 0 0\n'+' '.join([str(value)]*8)+'\n')


class VRTests(unittest.TestCase):
    def test_zero_placeholder_and_timeline(self):
        with tempfile.TemporaryDirectory() as tmp:
            root=Path(tmp);a=root/'a.cube';b=root/'b.cube';cube(a);cube(b,2.)
            result=export_cubes([a,b],root/'view',steps=[0,8],dt_au=.02,missing_ac='zero')
            config=json.loads((result.directory/'config.json').read_text())
            self.assertIn('placeholder',config['InfomationText'])
            self.assertEqual(config['AcFileFormat'],'*Ac_x.cube')
            np.testing.assert_array_equal(read_cube(result.directory/'Ac_Je_drho_CUBE/it000008_drho.cube').values,2.)
            np.testing.assert_array_equal(read_cube(result.directory/'Ac_Je_drho_CUBE/it000008_Ac_x.cube').values,0.)
            record=json.loads((result.directory/'salmonpy-vr.json').read_text())
            self.assertEqual(record['ac_source'],'zero_placeholder')
            self.assertEqual(record['frames'][1]['time_au'],.16)
            with self.assertRaises(FileExistsError):
                export_cubes([a,b],result.directory,steps=[0,8],dt_au=.02,missing_ac='zero')

    def test_actual_ac_and_reject_mismatch(self):
        with tempfile.TemporaryDirectory() as tmp:
            root=Path(tmp);a=root/'a.cube';b=root/'b.cube';cube(a);cube(b,.3)
            result=export_cubes([a],root/'valid',steps=[1],dt_au=.02,ac_files=[b])
            self.assertEqual(read_cube(result.directory/'Ac_Je_drho_CUBE/it000001_Ac_x.cube').values[0,0,0],.3)
            cube(b,.3,origin=1.)
            with self.assertRaisesRegex(ValueError,'grid'):
                export_cubes([a],root/'bad',steps=[1],dt_au=.02,ac_files=[b])
            self.assertFalse((root/'bad').exists())
            with self.assertRaises(ValueError): export_cubes([a],root/'missing',steps=[1],dt_au=.02)
            with self.assertRaises(ValueError): export_cubes([a,a],root/'duplicate',steps=[1,1],dt_au=.02,missing_ac='zero')

    def test_application_is_optional_and_never_launched_on_export(self):
        from unittest.mock import patch
        with tempfile.TemporaryDirectory() as tmp:
            root=Path(tmp);a=root/'a.cube';cube(a)
            app=root/'SALMON VR.app';(app/'Contents/MacOS').mkdir(parents=True)
            (app/'Contents/MacOS/SALMON VR').write_text('fixture')
            with patch('subprocess.run',side_effect=AssertionError('must not launch')):
                result=export_cubes([a],root/'view',steps=[0],dt_au=.02,missing_ac='zero',application=app)
            self.assertTrue((result.directory/app.name/'Contents/MacOS/SALMON VR').is_file())

    def test_reject_skew_grid(self):
        with tempfile.TemporaryDirectory() as tmp:
            root=Path(tmp);a=root/'a.cube';cube(a)
            a.write_text(a.read_text().replace('2 1 0 0','2 1 .1 0'))
            with self.assertRaisesRegex(ValueError,'axis-aligned'):
                export_cubes([a],root/'view',steps=[0],dt_au=.02,missing_ac='zero')
            self.assertFalse((root/'view').exists())

    def test_launch_preflight(self):
        import os
        import plistlib
        from unittest.mock import patch
        from salmonpy.salmon_vr import VRExport
        with tempfile.TemporaryDirectory() as tmp:
            root=Path(tmp);app=root/'SALMON VR.app'
            binary=app/'Contents/MacOS/SALMON VR';binary.parent.mkdir(parents=True)
            info=app/'Contents/Info.plist'
            info.write_bytes(plistlib.dumps({'CFBundleExecutable':'SALMON VR'}))
            result=VRExport(root,app)
            with patch('salmonpy.salmon_vr.sys.platform','darwin'),patch('salmonpy.salmon_vr.subprocess.run') as launch:
                with self.assertRaisesRegex(ValueError,'missing'):
                    result.launch()
                binary.write_text('fixture');binary.chmod(0o644)
                with patch('salmonpy.salmon_vr.os.access',return_value=False):
                    with self.assertRaisesRegex(PermissionError,'chmod'):
                        result.launch()
                self.assertEqual(binary.stat().st_mode & 0o777,0o644)
                launch.assert_not_called()
                with patch('salmonpy.salmon_vr.os.access',return_value=True):
                    result.launch()
                launch.assert_called_once_with(['/usr/bin/open','-n',str(app)],check=True)
                info.write_bytes(plistlib.dumps({'CFBundleExecutable':'../outside'}))
                with self.assertRaisesRegex(ValueError,'CFBundleExecutable'):
                    result.launch()
