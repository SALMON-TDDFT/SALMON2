"""Native factory ownership and scoped calls; no loading at import."""
from .errors import RuntimeFailure

class _SessionGroup:
    """Session-local model owner, including any returned but unadmitted states."""
    def __init__(self,model):
        self.model=model
        self.states=[]
    def new_state(self,point_id):
        state=self.model.new_state(point_id)
        self.states.append(state)
        return state
    def close(self):
        errors=[]
        for state in self.states:
            try:state.close()
            except Exception as exc:errors.append(str(exc))
        try:self.model.close()
        except Exception as exc:errors.append(str(exc))
        if errors:raise RuntimeFailure('shared native group cleanup failed: '+str(errors))


def _directory_call(lib,comm,work,fn,*args):
    """One main-thread native I/O scope, with collective admission/cleanup."""
    import ctypes as c
    import os
    from .errors import InputError
    entered=False;error=None
    try:
        try:
            probe=getattr(lib,'sp_io_directory_protocol',None)
            if probe is None:raise InputError('runtime requires directory protocol3; rebuild or install a matching runtime')
            probe.argtypes=[];probe.restype=c.c_int
            if probe()!=3:raise InputError('runtime requires directory protocol3; rebuild or install a matching runtime')
            enter=lib.sp_io_directory_enter;leave=lib.sp_io_directory_leave
            enter.argtypes=[c.c_char_p];enter.restype=c.c_int
            leave.argtypes=[];leave.restype=None
            status=enter(os.fsencode(work));entered=status==0
            if status:raise InputError('native directory scope rejected: '+str(status))
        except Exception as exc:error=str(exc)
        errors=comm.allgather(error)
        if any(errors):raise InputError('native directory scope: '+str(errors))
        return fn(*args)
    finally:
        if entered:leave()


def _bind_performance_directory(lib):
    """Bind once during checked setup; subsequent scopes have no MPI admission."""
    import ctypes as c
    from .errors import InputError
    query=lib.sp_io_directory_protocol;query.argtypes=[];query.restype=c.c_int
    if query()!=3:raise InputError('performance mode requires directory protocol3')
    enter=lib.sp_io_directory_enter;leave=lib.sp_io_directory_leave
    enter.argtypes=[c.c_char_p];enter.restype=c.c_int
    leave.argtypes=[];leave.restype=None
    return enter,leave
