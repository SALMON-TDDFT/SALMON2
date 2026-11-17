"""Four host ranks: warning-only budgets, unknown probes and rank agreement."""
import warnings
from unittest.mock import patch
from mpi4py import MPI
from salmonpy.memory import MemoryEstimate,MemoryBudgetWarning,memory_advice,process_memory

comm=MPI.COMM_WORLD
assert comm.size==4
MiB=2**20
estimate=MemoryEstimate(100*MiB,200*MiB,32*MiB)
with warnings.catch_warnings(record=True) as caught:
    warnings.simplefilter('always')
    small=memory_advice(estimate,budget_bytes=1024*MiB,comm=comm)
assert small.known_peak_bytes==4*332*MiB and small.exceeds_budget
assert len(caught)==(1 if comm.rank==0 else 0)
large=memory_advice(estimate,budget_bytes=4*2**30,comm=comm)
assert large.exceeds_budget is False
mixed=memory_advice(estimate,budget_bytes=(4 if comm.rank else 3)*2**30,comm=comm,warn=False)
assert mixed.exceeds_budget is None and mixed.unknown
with patch('salmonpy.memory._peak_rss',side_effect=OSError('measurement unavailable')):
    rss,basis=process_memory()
unknown=memory_advice(MemoryEstimate(rss,20*MiB,basis=basis),budget_bytes=4*2**30,comm=comm,warn=False)
assert unknown.estimated_peak_bytes is None and unknown.exceeds_budget is None
with warnings.catch_warnings():
    warnings.simplefilter('error',MemoryBudgetWarning)
    filtered=memory_advice(estimate,budget_bytes=1,comm=comm)
assert filtered.exceeds_budget and filtered.warning_error
assert comm.allreduce(1)==4
if comm.rank==0:print('memory advice MPI PASS: aggregate4GiB, warning-only excess, unknowns, policy mismatch, warning-filter errors; computation continues')
