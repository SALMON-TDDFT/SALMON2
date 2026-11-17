"""Private same-lifetime native candidate protocol; not a reentrant solver ABI.

Hosts pause evolution, stage on every solver group, vote globally, then publish.
Only numerical payloads are detached; initialized resource handles stay live.
"""
import ctypes as c
import os
from .errors import InputError, RuntimeFailure


class NativeRestoreCandidate:
    """One detached numerical candidate for an initialized native owner.

    Construct before staging so discard() can release partial preparation too.
    No checkpoint file is accessed by publish(). Python-side counters/arrays are
    owned by the coordinator and must be prepared before the global ready vote.
    """
    def __init__(self, owner, symbol, path):
        self.owner = owner
        self.path = os.fsencode(path)
        self.ready = False
        self.finished = False
        self.reserved = False
        probe = getattr(owner.lib, symbol+'_protocol', None)
        if probe is None:
            raise InputError('native runtime lacks isolated checkpoint candidates; rebuild explicitly')
        probe.argtypes = []
        probe.restype = c.c_int
        if probe() != 1:
            raise InputError('unsupported native checkpoint candidate protocol')
        self.fn = getattr(owner.lib, symbol)
        self.fn.argtypes = [c.c_char_p, c.c_int]
        self.fn.restype = c.c_int

    def stage(self):
        if self.ready or self.finished:
            raise InputError('candidate already staged or retired')
        if getattr(self.owner, '_step_guard', None) is not None:
            raise InputError('cannot stage checkpoint during a protected native step')
        if getattr(self.owner, "_restore_candidate", None) is not None:
            raise InputError("one pending candidate per native lifetime")
        self.owner._restore_candidate = self
        self.reserved = True
        code = self.owner._call(self.fn, self.path, 2)
        if code:
            raise RuntimeFailure('native candidate stage rejected: '+str(code))
        self.ready = True
        return self

    def discard(self):
        if self.finished or not self.reserved:
            return
        code = self.owner._call(self.fn, b'', 4)
        if code:
            raise RuntimeFailure('native candidate discard failed: '+str(code))
        self.finished = True
        self.owner._restore_candidate = None

    def publish(self):
        if not self.ready or self.finished:
            raise InputError('prepared native candidate required')
        code = self.owner._call(self.fn, b'', 3)
        if code:
            raise RuntimeFailure('native candidate publication violated prepared-state contract: '+str(code))
        self.finished = True
        self.owner._restore_candidate = None
