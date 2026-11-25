"""One launch: LD/FDTD spatial E/H drives a molecule without current feedback.

mpiexec -n 3 python -m salmonpy.examples.molecular_fdtd_response FIELD_BUILD
WORKER VALIDATION EXISTING_C2H2_BUNDLE NEW_OUTPUT --field-ranks 2 --ld

Small fixed-geometry example, not a converged optical spectrum. Reuses a GS;
never builds, downloads or starts SCF. --stop-after saves field/LD/molecule.
"""
import sys
from ._grid_coupling import entry as _run


def main():
    _run(sys.argv[1:],mode='one-way')


if __name__=='__main__':main()
