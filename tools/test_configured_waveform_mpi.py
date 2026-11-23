"""Three-rank bounded-run worker: different valid waveform lengths are rejected."""
import sys
from pathlib import Path
import numpy as np
from mpi4py import MPI
from salmonpy.examples import mpi_configured_native as example

# Arguments: catalog normal_wave shorter_wave fresh_output
catalog, normal, shorter, output = sys.argv[1:]
original = np.load

def altered(path, *args, **kwargs):
    if MPI.COMM_WORLD.rank == 1 and Path(path) == Path(normal):
        return original(shorter, *args, **kwargs)
    return original(path, *args, **kwargs)

example.np.load = altered
sys.argv = ['sample', catalog, 'metal-A', normal, output, '--mode', 'ld-endpoints']
rejected = False
try:
    example.main()
except ValueError as exc:
    assert 'waveform differs across MPI ranks' in str(exc), str(exc)
    rejected = True
finally:
    np.load = original
assert MPI.COMM_WORLD.allreduce(int(rejected)) == 3
if MPI.COMM_WORLD.rank == 0: print('PASS waveform mismatch before collective stepping')
