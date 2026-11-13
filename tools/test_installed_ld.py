"""Run installed-package tests without allowing checkout imports."""
from pathlib import Path
import sys
import unittest
import salmonpy

root = Path(sys.argv[1]).resolve()
loaded = Path(salmonpy.__file__).resolve()
if loaded.is_relative_to(root):
    raise RuntimeError(f'Test imported checkout instead of wheel: {loaded}')
print(f'Testing installed package: {loaded}', flush=True)
suite = unittest.defaultTestLoader.discover(str(root / 'salmonpy/tests'))
result = unittest.TextTestRunner(verbosity=2).run(suite)
raise SystemExit(not result.wasSuccessful())
