"""Typical collective GS/2→retained GS→transverse charge-current workflow.

Call current_from_archive with matching arguments on all session.comm ranks.
Keep one caller-owned NativeGSSession for repeated calculations. The explicit
converter/paired serial executable is prepared beforehand; no implicit build,
download or SCF. A_steps contains A/c at dt,2dt,... The new converted bundle
must be visible to subgroup ranks and remains caller-owned after return.
"""
import numpy as np


def current_from_archive(session,archive,directory,*,converter,executable,
                         pseudopotentials,A_steps,dt,material_id=None,reexport_to=None):
    """Return j[n,3], releasing this calculation's GS/model/state owners."""
    gs=session.load_archive(archive,directory,converter=converter,executable=executable,
                            pseudopotentials=pseudopotentials,material_id=material_id)
    try:
        if reexport_to is not None:gs.export_archive(reexport_to)
        model=gs.tddft_model(dt=dt,max_steps=max(1,len(A_steps)))
    finally:
        gs.close()
    try:
        state=model.new_state('response')
        try:
            return np.asarray([state.step(A) for A in A_steps]).reshape(-1,3)
        finally:
            state.close()
    finally:
        model.close()
