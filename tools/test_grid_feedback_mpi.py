"""Spatial molecular/FDTD validation; physical reproduction stays manual."""
from pathlib import Path
import sys

# Direct source invocation must use this checkout, not an older installed wheel.
root=Path(__file__).resolve().parents[1]
if (root/'salmonpy/__init__.py').is_file():sys.path.insert(0,str(root))
from salmonpy.examples._grid_coupling import entry

if __name__=='__main__':entry()
