"""Opt-in same-lifetime native step backups and audited append-output rollback.

No compilation, native probing or MPI initialization occurs on import. Protocol
queries happen only when explicitly constructing targets for a live coupled run.
"""
import ctypes as c
import os
from pathlib import Path
import stat
import numpy as np
from .errors import InputError,RuntimeFailure
from .stepping import StepTarget

__all__=['native_step_targets']


class _AppendJournal:
    def __init__(self,root,files,selection=None):
        self.root=root
        self.files=files
        self.selection=selection

    @staticmethod
    def _inventory(root,selection=None):
        records={}
        paths=root.glob('*.data') if selection is None else (root/name for name in selection)
        for path in paths:
            value=path.lstat()
            if not stat.S_ISREG(value.st_mode) or value.st_nlink!=1:
                raise InputError('native step output must be a private regular file: '+path.name)
            records[path.name]=(value.st_dev,value.st_ino,value.st_size)
        return records

    @classmethod
    def capture(cls,root,*,files=None):
        root=Path(root).resolve()
        selection=None if files is None else tuple(files)
        if selection is not None and any(Path(name).name!=name for name in selection):
            raise InputError("output journal requires local basenames")
        return cls(root,cls._inventory(root,selection),selection)

    def restore(self):
        current=self._inventory(self.root,self.selection)
        if set(current)!=set(self.files):
            raise RuntimeFailure('native step stream inventory changed; no automatic deletion/recreation')
        for name,(device,inode,size) in self.files.items():
            actual=current[name]
            if actual[:2]!=(device,inode) or actual[2]<size:
                raise RuntimeFailure('native step stream replaced or shortened: '+name)
        # Writers must already be closed by the native transaction endpoint.
        for name,(device,inode,size) in self.files.items():
            fd=os.open(self.root/name,os.O_WRONLY|getattr(os,'O_NOFOLLOW',0))
            try:
                value=os.fstat(fd)
                if (value.st_dev,value.st_ino)!=(device,inode) or value.st_size<size:
                    raise RuntimeFailure('native step stream changed before truncation: '+name)
                os.ftruncate(fd,size)
            finally:
                os.close(fd)


def _snapshot_attrs(owner,names):
    result={}
    for name in names:
        value=getattr(owner,name)
        result[name]=(value,value.copy()) if isinstance(value,np.ndarray) else value
    return result


def _restore_attrs(owner,saved):
    for name,value in saved.items():
        if isinstance(value,tuple) and len(value)==2 and isinstance(value[0],np.ndarray):
            original,values=value
            if not np.array_equal(original,values):original[:]=values
            setattr(owner,name,original)
        else:
            setattr(owner,name,value)


class _NativeBackup:
    def __init__(self,owner,comm,symbol,root,attrs):
        self.owner=owner
        self.comm=comm
        self.root=Path(root)
        self.attrs=attrs
        self.active=False
        self.unsafe=False
        probe=getattr(owner.lib,symbol+'_protocol',None)
        if probe is None:raise InputError('native whole-step transaction protocol missing; explicitly rebuild current sources')
        probe.argtypes=[];probe.restype=c.c_int
        if probe()!=1:raise InputError('unsupported native whole-step transaction protocol')
        self.fn=getattr(owner.lib,symbol)
        self.fn.argtypes=[c.c_int];self.fn.restype=c.c_int

    def _check(self,error,*,admission=False):
        errors=self.comm.allgather(error)
        if any(errors):
            if admission:raise InputError(str(errors))
            self.unsafe=True;self.owner.failed=True
            raise RuntimeFailure('native transaction group rejected: '+str(errors))

    def _mode(self,mode):
        error=None
        try:
            folder=getattr(self.owner,'directory',None)
            folder=Path(folder if folder is not None else self.owner.work)
            if not folder.is_dir() or folder.resolve()!=self.root.resolve():
                raise RuntimeFailure('native owner output directory changed')
            native=getattr(self.owner,'native',None) or getattr(self.owner,'_native',None)
            native._live()
        except Exception as exc:error=type(exc).__name__+': '+str(exc)
        self._check(error,admission=mode==0)
        code=self.owner._call(self.fn,mode)
        codes=self.comm.allgather(int(code))
        if any(codes):
            if mode==0 and len(set(codes))==1 and codes[0] in (10,11,12,13,14,15):
                raise InputError('native transaction admission rejected: '+str(codes))
            self.unsafe=True
            self.owner.failed=True
            raise RuntimeFailure('native transaction failed: '+str(codes))

    def capture(self):
        error=None
        if getattr(self.owner,'_step_guard',None) is not None:error='one native step backup per owner'
        if getattr(self.owner,'_restore_candidate',None) is not None:error='checkpoint candidate already owns this native lifetime'
        if self.owner.closed or getattr(self.owner,'failed',False):error='usable native owner required'
        if not self.root.is_dir():error='native output directory missing'
        self._check(error,admission=True)
        self.owner._step_guard=self
        self.unsafe=False
        try:
            self._mode(0)
            self.active=True
        except Exception:
            self.owner._step_guard=None
            raise
        saved=_snapshot_attrs(self.owner,self.attrs)
        # Root failure returns to the common capture vote. release(None) still
        # retires the collective native backup on all subgroup ranks.
        journal=None
        if self.comm.rank==0:
            selected=getattr(self.owner,'output_files',None)
            journal=(_AppendJournal.capture(self.root) if selected is None else
                     _AppendJournal.capture(self.root,files=selected))
        return saved,journal

    def rollback(self,token):
        error=None
        if not self.active or token is None:error='live native backup required'
        if self.owner.closed or getattr(self.owner,'failed',False):error='unknown native failure forbids transaction reuse'
        self._check(error)
        self._mode(1)
        saved,journal=token
        error=None
        if self.comm.rank==0:
            try:journal.restore()
            except Exception as exc:error=type(exc).__name__+': '+str(exc)
        errors=self.comm.allgather(error)
        if any(errors):
            self.unsafe=True;self.owner.failed=True
            raise RuntimeFailure('native output rollback failed: '+str(errors))
        self._mode(3)
        _restore_attrs(self.owner,saved)

    def release(self,token):
        error='unsafe native transaction requires shutdown' if self.unsafe or getattr(self.owner,'failed',False) else None
        self._check(error)
        active=self.comm.allgather(self.active)
        self._check('inconsistent native backup activity' if any(v!=active[0] for v in active) else None)
        if self.active:
            self._mode(2)
            self.active=False
        if getattr(self.owner,'_step_guard',None) is self:self.owner._step_guard=None


def native_step_targets(coupled,*,custom_targets=()):
    """Collectively construct explicit audited targets; no capture or stepping.

    Covers exact built-in MultiscaleField, NativeEHField, native semilocal TDDFT, fast-PR HSE and
    SBE endpoint wrapper types. Other owners require explicit custom StepTargets.
    Native libraries must declare the whole-step protocol, not just per-component
    numerical recovery. Native begin also rejects unaudited output profiles.
    Streams must be private, regular, append-only and have no concurrent reader/
    writer contract. Rollback is not durable/crash-atomic or externally invisible.
    """
    from .coupled import NativeEHField,MultiscaleField
    from .embedded_tddft import _NativeState
    from ._shared_tddft import _SharedState as SharedTDDFTState
    from .fast_hse_mpi import _State as FastHSEState
    from .materials import _EndpointState
    from .sbe import SBEState
    r=coupled.response
    comm=r.context.solver_comm
    error=None;targets=[]
    try:
        supplied=tuple(custom_targets)
        if any(not isinstance(t,StepTarget) for t in supplied):raise InputError('explicit custom StepTargets required')
        if len({id(t.owner) for t in supplied})!=len(supplied):raise InputError('duplicate custom owner')
        owners=([coupled.field] if coupled.field is not None else [])+list(r.states.values())
        if any(not any(t.owner is owner for owner in owners) for t in supplied):
            raise InputError('custom target does not match a local owner')
        for index,owner in enumerate(owners):
            override=next((t for t in supplied if t.owner is owner),None)
            if override is not None:
                targets.append(override);continue
            if type(owner) is MultiscaleField:
                def capture(f=owner):return f.state.snapshot(),f.closed,f.failed
                def rollback(token,f=owner):
                    f.state.restore(token[0]);f.closed=token[1];f.failed=token[2]
                targets.append(StepTarget('field-snapshot',owner,capture,rollback))
                continue
            if type(owner) is NativeEHField:
                backup=_NativeBackup(owner,comm,'sp_eh_step_transaction',owner.directory,
                    ('A','E','step_index','failed','closed'))
            elif type(owner) in (_NativeState,SharedTDDFTState):
                backup=_NativeBackup(owner,comm,'sp_tddft_step_transaction',owner.work,
                    ('A','steps','failed','closed')+(('current',) if type(owner) is SharedTDDFTState else ()))
            elif type(owner) is FastHSEState:
                backup=_NativeBackup(owner,comm,'sp_fast_hse_step_transaction',owner.work,
                    ('A','current','energy','steps','failed','closed'))
            elif type(owner) is _EndpointState and type(owner._raw) is SBEState:
                def capture(s=owner):return s._raw.snapshot(),s._steps,s._failed,s._closed
                def rollback(token,s=owner):
                    snapshot,s._steps,s._failed,s._closed=token
                    s._raw.restore(snapshot)
                targets.append(StepTarget('sbe-'+str(index),owner,capture,rollback))
                continue
            else:
                raise InputError('unsupported whole-step owner; supply an explicit complete custom target')
            targets.append(StepTarget('native-'+str(index),owner,backup.capture,backup.rollback,backup.release))
    except Exception as exc:
        error=type(exc).__name__+': '+str(exc)
    r.agree(error)
    return tuple(targets)
