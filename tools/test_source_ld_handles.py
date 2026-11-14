"""Source suite with explicit new LD/SBE library, without replacing fixed binaries."""
from pathlib import Path
import sys,unittest
from salmonpy import ld,sbe
from salmonpy._ld_handles import API
from salmonpy._sbe_handles import API as SBEAPI
api=API(sys.argv[1]);ld._handles_api=lambda:api;ld._native=lambda:api.legacy
sbeapi=SBEAPI(sys.argv[1]);sbe._handles_api=lambda:sbeapi;sbe._native=lambda:sbeapi.legacy
root=Path(__file__).resolve().parents[1]
suite=unittest.defaultTestLoader.discover(str(root/'salmonpy/tests'),pattern=sys.argv[2] if len(sys.argv)>2 else 'test_*.py')
result=unittest.TextTestRunner(verbosity=2).run(suite)
raise SystemExit(not result.wasSuccessful())
