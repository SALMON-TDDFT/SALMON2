"""Advisory aggregate memory estimates; never a scheduler or hard RSS limit."""
from dataclasses import dataclass, replace
import operator
import sys
import warnings
from .errors import InputError

__all__=['MemoryEstimate','MemoryReport','MemoryBudgetWarning','process_memory','memory_advice']


def _bytes(value, name, *, optional=False):
    if value is None and optional:
        return None
    try:
        if isinstance(value,bool):
            raise TypeError()
        result=operator.index(value)
    except TypeError as exc:
        raise InputError(name+' must be a nonnegative integer byte count') from exc
    if result<0:
        raise InputError(name+' must be nonnegative')
    return result


@dataclass(frozen=True,slots=True)
class MemoryEstimate:
    """One rank's live storage, extra candidate/copy storage and reserve.

    None means unavailable, not zero. Components are additive; candidate_bytes
    must be *extra* storage, not a second estimate of the whole process.
    """
    live_bytes: int | None
    candidate_bytes: int | None = 0
    reserve_bytes: int = 0
    basis: str = 'caller supplied'
    unknown: tuple = ()

    def __post_init__(self):
        for name in ('live_bytes','candidate_bytes'):
            object.__setattr__(self,name,_bytes(getattr(self,name),name,optional=True))
        object.__setattr__(self,'reserve_bytes',_bytes(self.reserve_bytes,'reserve_bytes'))
        if not isinstance(self.basis,str) or not self.basis:
            raise InputError('memory estimate basis must be a nonempty string')
        if isinstance(self.unknown,str):
            raise InputError('unknown must be a sequence of descriptions')
        values=tuple(self.unknown)
        if any(not isinstance(value,str) or not value for value in values):
            raise InputError('unknown entries must be nonempty strings')
        object.__setattr__(self,'unknown',values)


class MemoryBudgetWarning(UserWarning):
    """Advisory only: exceeding a budget never refuses scientific work."""


@dataclass(frozen=True,slots=True)
class MemoryReport:
    label: str
    budget_bytes: int | None
    ranks: tuple
    known_peak_bytes: int
    estimated_peak_bytes: int | None
    unknown: tuple
    warning_error: str | None = None

    @property
    def exceeds_budget(self):
        if self.budget_bytes is None:
            return None
        if self.known_peak_bytes>self.budget_bytes:
            return True
        return None if self.unknown else False

    def to_dict(self):
        return {'schema':'salmonpy.memory-advice/1','label':self.label,
            'budget_bytes':self.budget_bytes,'rank_count':len(self.ranks),
            'known_peak_bytes':self.known_peak_bytes,'estimated_peak_bytes':self.estimated_peak_bytes,
            'exceeds_budget':self.exceeds_budget,
            'excess_bytes':(max(0,self.known_peak_bytes-self.budget_bytes) if self.exceeds_budget is not None else None),
            'unknown':list(self.unknown),
            'warning_error':self.warning_error,'hard_limit':False,
            'ranks':[{'live_bytes':r.live_bytes,'candidate_bytes':r.candidate_bytes,
                      'reserve_bytes':r.reserve_bytes,'basis':r.basis,'unknown':list(r.unknown)} for r in self.ranks]}


def _peak_rss():
    if sys.platform not in ('darwin','linux'):
        raise OSError('process peak RSS probe unavailable on '+sys.platform)
    import resource
    value=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
    return int(value)*(1 if sys.platform=='darwin' else 1024)


def process_memory():
    """Return (peak_RSS_bytes, basis), or (None, reason) on probe failure.

    Linux/macOS process high-water RSS, not current use or an aggregate peak.
    Child processes are not included. Imports no MPI and performs no allocation
    probe. Other platforms can supply their own MemoryEstimate.
    """
    try:
        value=_bytes(_peak_rss(),'peak RSS')
        if value==0:
            raise OSError('peak RSS probe returned zero')
        return value,'observed process peak RSS (children excluded)'
    except Exception as exc:
        return None,type(exc).__name__+': '+str(exc)


def memory_advice(estimate, *, budget_bytes=None, comm=None, label='memory', warn=True):
    """Collect supplied rank estimates and optionally warn on rank0 only.

    Budget covers the sum across this communicator, not each node separately.
    No job admission, MPI initialization or scientific callback. Unknown data
    stays unknown. Warning-filter errors are captured, never computation errors.
    All participating ranks must call; communicator failure is not recoverable.
    """
    rank=0 if comm is None else comm.rank
    gather=(lambda value:[value]) if comm is None else comm.allgather
    error=None
    try:
        if estimate is None:
            estimate=MemoryEstimate(None,None,basis='estimate unavailable')
        if not isinstance(estimate,MemoryEstimate):
            raise InputError('MemoryEstimate required')
        budget_bytes=_bytes(budget_bytes,'budget_bytes',optional=True)
        if not isinstance(label,str) or not label or type(warn) is not bool:
            raise InputError('nonempty label and boolean warn required')
    except Exception as exc:
        error=type(exc).__name__+': '+str(exc)
    rows=gather((error,estimate if error is None else None,budget_bytes if error is None else None,
                 label if error is None else None,warn if error is None else None))
    if any(row[0] is not None for row in rows):
        raise InputError('memory advice input error: '+str([r[0] for r in rows]))
    estimates=tuple(row[1] for row in rows)
    unknown=[]
    for i,r in enumerate(estimates):
        unknown.extend('rank'+str(i)+': '+value for value in r.unknown)
        for name in ('live_bytes','candidate_bytes'):
            if getattr(r,name) is None:
                unknown.append('rank'+str(i)+': '+name+' unavailable')
    if any(row[2]!=rows[0][2] for row in rows):
        budget_bytes=None
        unknown.append('memory budgets differ across ranks')
    else:
        budget_bytes=rows[0][2]
    if any(row[3:]!=rows[0][3:] for row in rows):
        unknown.append('memory advice labels/warning settings differ across ranks')
    total=sum((r.live_bytes or 0)+(r.candidate_bytes or 0)+r.reserve_bytes for r in estimates)
    report=MemoryReport(rows[0][3],budget_bytes,estimates,total,None if unknown else total,tuple(unknown))
    message=None
    if report.exceeds_budget:
        message=(report.label+': estimated aggregate memory '+format(total/2**30,'.3f')+' GiB exceeds budget '
                 +format(budget_bytes/2**30,'.3f')+' GiB by '+format((total-budget_bytes)/2**30,'.3f')
                 +' GiB; continuing (advisory, not a hard RSS limit)')
    elif any(row[2] is not None for row in rows) and unknown:
        message=report.label+': incomplete memory estimate; '+ '; '.join(unknown)+'; continuing (advisory only)'
    warning_error=None
    if rank==0 and rows[0][4] and message:
        try:
            warnings.warn(message,MemoryBudgetWarning,stacklevel=2)
        except Exception as exc:
            warning_error=type(exc).__name__+': '+str(exc)
    # Every rank receives the same warning-delivery diagnostic.
    errors=gather(warning_error)
    return replace(report,warning_error=next((e for e in errors if e),None))
