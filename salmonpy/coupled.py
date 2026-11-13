"""Experimental native E/H-to-A/J coupling on separate field/material MPI groups.

E is at native time n*dt. Electronic local time n*dt corresponds to field clock
(n+1/2)*dt; native saved H is centered at E time, not the raw H endpoint.
"""
from pathlib import Path
import ctypes as c
import json
import os
import threading
import numpy as np
from . import embedded_tddft as embedded
from .distributed import NativeParallelContext
from .fdtd import FDTD,FDTDMPI,PlaneWave
from .tddft import _digest


def _agree(comm,error):
    errors=comm.allgather(error)
    if any(errors):raise ValueError(str(errors))


class NativeEHField:
    """Collective field-group context; points map IDs to zero-based grid indices.

    Each macro point represents one uniform field-cell volume. Coupled sites and
    their constitutive neighbors must be vacuum in the background geometry; LD
    may occupy other regions. Only zero-initialized plane-wave runs are supported.
    """
    scheme='salmonpy.eh-centered-current/1'

    def __init__(self,context,calculation,build_record,directory,*,points,grid):
        comm=context.solver_comm;error=None
        try:
            if context.group.role!='field':raise ValueError('field group required')
            if embedded._used:raise ValueError('one native lifetime per process')
            if threading.current_thread() is not threading.main_thread():raise ValueError('main thread required')
            if not isinstance(calculation,FDTD) or not isinstance(calculation.source,PlaneWave):
                raise ValueError('zero-initialized plane-wave FDTD required')
            parallel=FDTDMPI(tuple(grid))
            if parallel.ranks!=comm.size:raise ValueError('field decomposition must match group size')
            text=parallel.input_text(calculation)
            self.point_ids=context.plan.point_ids
            if set(points)!=set(self.point_ids):raise ValueError('one field index per material point required')
            indices=np.asarray([points[p] for p in self.point_ids])
            if indices.shape!=(len(self.point_ids),3) or indices.dtype.kind not in 'iu':raise ValueError('integer grid triples required')
            counts=np.array(calculation.geometry.grid.counts)
            if np.any(indices<10) or np.any(indices>=counts-10):raise ValueError('coupled cells must be at least ten cells from boundaries')
            if len({tuple(p) for p in indices})!=len(indices):raise ValueError('one macro state per field cell required')
            axis=calculation.source.axis
            source_index=int(np.argmin(np.abs(calculation.geometry.grid.coordinates[axis]-calculation.source.position)))
            if np.any(np.abs(indices[:,axis]-source_index)<=2):
                raise ValueError('coupled cells must be clear of the incident source plane and its support')
            ids=calculation.geometry.material_ids()
            for point in indices:
                # Conservative guard covers all constitutive cells touching mapped edges.
                i,j,k=point
                if np.any(ids[i-1:i+2,j-1:j+2,k-1:k+2]):raise ValueError('coupled cells require a vacuum constitutive neighborhood')
            record=json.loads(Path(build_record).read_text())
            if record.get('schema')!='salmonpy.embedded-tddft/1':raise ValueError('embedded host library record required')
            library=Path(record['library']);sha=record['sha256']
            if _digest(library)!=sha:raise ValueError('native library checksum mismatch')
            lib=c.CDLL(str(library))
            ints=np.ctypeslib.ndpointer(dtype=np.int32,ndim=2,flags='C_CONTIGUOUS')
            floats=np.ctypeslib.ndpointer(dtype=np.float64,ndim=2,flags='C_CONTIGUOUS')
            lib.sp_eh_init.argtypes=[c.c_char_p,c.c_int,ints];lib.sp_eh_init.restype=c.c_int
            lib.sp_eh_step.argtypes=[c.c_int,floats,floats];lib.sp_eh_step.restype=c.c_int
            lib.sp_eh_close.argtypes=[];lib.sp_eh_close.restype=c.c_int
            self.directory=Path(directory).expanduser().resolve()
            descriptor=(text,indices.tolist(),sha,str(self.directory))
        except Exception as exc:error=str(exc)
        _agree(comm,error)
        descriptors=comm.allgather(descriptor)
        _agree(comm,'field configuration differs across ranks' if any(x!=descriptor for x in descriptors) else None)
        error=None
        if comm.rank==0:
            try:
                self.directory.mkdir(parents=True,exist_ok=False)
                (self.directory/'inputfile').write_text(text)
            except Exception as exc:error=str(exc)
        _agree(comm,error)
        error=None
        try:
            if (self.directory/'inputfile').read_text()!=text:raise ValueError('shared field directory required')
        except Exception as exc:error=str(exc)
        _agree(comm,error)
        self.context=context;self.lib=lib;self.dt=calculation.dt;self.capacity=calculation.steps
        self.step_index=0;self.A=np.zeros((len(self.point_ids),3));self.E=np.zeros_like(self.A)
        self.closed=False;self.initialized=False;self.failed=False
        self.manifest={'schema':self.scheme,'status':'started','point_ids':self.point_ids,
            'grid_indices':indices.tolist(),'cell_volume':float(np.prod(calculation.geometry.grid.spacing)),
            'dt':self.dt,'steps':self.capacity,'field_ranks':comm.size,'field_grid':list(grid),
            'library_sha256':sha,'sample_map':'half of same-index and previous-axis E edges',
            'deposit_map':'transpose of sample_map, equal cell volumes',
            'native_E_time':'n*dt','native_saved_H_time':'n*dt',
            'electronic_local_time':'n*dt','electronic_field_clock':'(n+0.5)*dt',
            'initial_feedback_current':'zero; not substituted with initial GS current'}
        self._write_manifest()
        self.native=NativeParallelContext(context,library,sha)
        embedded._used=True
        try:
            code=self._call(lib.sp_eh_init,b'inputfile',len(self.point_ids),np.array(indices,dtype=np.int32,order='C'))
            self.initialized=(code==0)
            _agree(comm,'native field initialization rejected' if code else None)
        except BaseException:
            self.failed=True
            self.close();raise

    def _write_manifest(self):
        error=None
        if self.context.is_leader:
            try:(self.directory/'salmonpy-coupled-field.json').write_text(json.dumps(self.manifest,indent=2)+'\n')
            except Exception as exc:error=str(exc)
        _agree(self.context.solver_comm,error)

    def _call(self,fn,*args):
        self.native._live();previous=Path.cwd()
        try:
            os.chdir(self.directory)
            return fn(*args)
        finally:os.chdir(previous)

    def advance(self,J_previous):
        comm=self.context.solver_comm;error=None
        try:
            if self.closed or self.failed or self.step_index>=self.capacity:raise ValueError('closed, failed or exhausted field')
            J=np.array(J_previous,dtype=np.float64,order='C')
            if J.shape!=self.A.shape or not np.isfinite(J).all():raise ValueError('finite point currents required')
            if self.step_index==0 and np.any(J):raise ValueError('initial feedback current must be zero')
        except Exception as exc:error=str(exc)
        _agree(comm,error)
        values=comm.allgather(J.tolist())
        _agree(comm,'field currents differ across ranks' if any(v!=values[0] for v in values) else None)
        E=np.empty_like(self.A)
        try:
            code=self._call(self.lib.sp_eh_step,len(self.point_ids),J,E)
            _agree(comm,'native field advance failed' if code or not np.isfinite(E).all() else None)
            with np.errstate(over="ignore",invalid="ignore"):
                updated_A=self.A-self.dt*E
            _agree(comm,'nonfinite integrated A' if not np.isfinite(updated_A).all() else None)
        except BaseException as exc:
            self.failed=True
            self.manifest['failure']=str(exc)
            raise
        self.E=E;self.A=updated_A;self.step_index+=1
        return self.A.copy()

    def close(self):
        if self.closed:return
        code=0
        if self.initialized:
            code=self._call(self.lib.sp_eh_close)
        _agree(self.context.solver_comm,'native field close failed' if code else None)
        self.native.close();self.closed=True
        self.manifest['completed_steps']=self.step_index
        self.manifest['status']='failed' if self.failed else ('completed' if self.step_index==self.capacity else 'stopped')
        self._write_manifest()

    def __enter__(self):return self
    def __exit__(self,*args):self.close()


class CoupledStepper:
    """Collective outer schedule; material implementations remain replaceable.

    Construct ResponseSession on every rank; pass NativeEHField on field ranks
    and None on material ranks. Call step everywhere. Close the response, field,
    and ParallelContext in that order. A step failure invalidates the response.
    """
    def __init__(self,response,field=None):
        self.response=response;self.field=field;self.step_index=0
        context=response.context;error=None
        if response.closed or response.failed or response.step_index!=0:
            error='coupling requires a fresh response session'
        if context.group.role=='field':
            if not isinstance(field,NativeEHField) or field.context is not context:
                error='matching NativeEHField required on field ranks'
            elif field.closed or field.failed or field.step_index!=0:
                error='coupling requires a fresh field'
            elif field.dt!=response.dt or field.point_ids!=context.plan.point_ids:
                error='field/response dt and point order must match'
        elif field is not None:error='material ranks must not provide a field'
        response._agree(error)
        self.J=np.zeros((len(context.plan.point_ids),3))

    def step(self):
        r=self.response;context=r.context;A=None;error=None
        if r.closed or r.failed:error='response is not usable'
        if r.step_index!=self.step_index or (self.field is not None and self.field.step_index!=self.step_index):
            error='field and response histories are not synchronized'
        try:r._agree(error)
        except Exception:
            r.failed=True
            raise
        if self.field is not None:
            try:A=self.field.advance(self.J)
            except Exception as exc:error=str(exc)
        try:r._agree(error)
        except Exception:
            r.failed=True
            raise
        try:J=r.advance(A if context.group.role=='field' and context.is_leader else None)
        except BaseException:
            r.failed=True
            raise
        self.step_index+=1
        if self.field is not None:
            previous=self.J.copy();self.J=J
            n=self.field.step_index;dt=self.field.dt
            return {'field_time':n*dt,'electronic_local_time':n*dt,
                    'electronic_field_clock':(n+.5)*dt,'A_over_c':A.copy(),
                    'E_at_points':self.field.E.copy(),'J_previous':previous,'J_charge':J.copy()}
        return None
