"""Run source tests with an explicit, checksum-verified LD/SBE library.

Usage: PYTHONPATH=. python tools/test_source_materials.py LIBRARY [TEST_PATTERN]
No fixed binary is replaced; both LD and SBE use the specified library.
"""
from pathlib import Path
import sys,unittest
from salmonpy import ld,sbe,_material_handles
from salmonpy._material_handles import LDAPI as API
from salmonpy._material_handles import SBEAPI
api=API(sys.argv[1]);ld._handles_api=lambda:api;ld._native=lambda:api.legacy
sbeapi=SBEAPI(sys.argv[1]);sbe._handles_api=lambda:sbeapi;sbe._native=lambda:sbeapi.legacy
_material_handles.default_ld_api=lambda:api
_material_handles.default_sbe_api=lambda:sbeapi
root=Path(__file__).resolve().parents[1]
suite=unittest.defaultTestLoader.discover(str(root/'salmonpy/tests'),pattern=sys.argv[2] if len(sys.argv)>2 else 'test_*.py')
result=unittest.TextTestRunner(verbosity=2).run(suite)
raise SystemExit(not result.wasSuccessful())
