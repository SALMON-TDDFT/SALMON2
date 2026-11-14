"""Run source tests against an explicitly built, checksum-verified new library.

Does not overwrite the source checkout's existing native binary or build on import.
Usage: PYTHONPATH=. python tools/test_source_sbe_handles.py LIBRARY [TEST_PATTERN]
"""
from pathlib import Path
import sys,unittest
from salmonpy._sbe_handles import API
from salmonpy import sbe
api=API(sys.argv[1]);sbe._handles_api=lambda:api;sbe._native=lambda:api.legacy
root=Path(__file__).resolve().parents[1]
suite=unittest.defaultTestLoader.discover(str(root/'salmonpy/tests'),pattern=sys.argv[2] if len(sys.argv)>2 else 'test_*.py')
result=unittest.TextTestRunner(verbosity=2).run(suite)
raise SystemExit(not result.wasSuccessful())
