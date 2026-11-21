"""Small pinned CC0 numerical optical library. Import/listing is offline.

Raw source tables and fitted scalar nonmagnetic models are separate records.
Loading never fits, downloads, compiles, initializes MPI or runs a solver.
"""
from dataclasses import dataclass,field
from importlib.resources import files
from pathlib import Path
import ctypes as c
import hashlib
import json
import re
import warnings
import numpy as np
from .errors import InputError,RuntimeFailure
from .ld import LorentzDrude,_array,_handles_api,_isotropic_tensor
from .fields import _readonly

__all__=['OpticalMaterial','list_materials','load_material','wavelength_to_frequency']


def _root():return files('salmonpy').joinpath('optical_data')


def _index():
    data=json.loads(_root().joinpath('catalog.json').read_text(encoding='utf-8'))
    if data.get('schema')!='salmonpy.optical-materials/1' or data.get('license')!='CC0-1.0':
        raise InputError('unsupported optical data inventory')
    if hashlib.sha256(_root().joinpath('CC0-1.0.txt').read_bytes()).hexdigest()!=data['license_sha256']:
        raise InputError('optical data license checksum mismatch')
    return data


def list_materials():
    """Independent metadata summaries; no native call or network access."""
    return tuple({k:v for k,v in m.items() if k not in ('fits',)} |
                 {'fits':[{'poles':r['poles'],'status':r['report']['status'],
                           'relative_rms_error':r['report']['relative_rms_error'],
                           'optimizer_converged':r['report']['optimizer_converged'],
                           'eligible_default':r.get('curation',{}).get('eligible_default',False)}  for r in m['fits']]}
                 for m in _index()['materials'])


def wavelength_to_frequency(wavelength_um):
    """Vacuum wavelength in micrometers -> angular omega(AU), photon energy(eV).

    Conversion runs in Fortran using SALMON's CODATA-2018 constants.
    """
    wavelength=_array(wavelength_um,'wavelength_um')
    omega=np.empty_like(wavelength);energy=np.empty_like(wavelength)
    try:fn=_handles_api().library.sp_ld_wavelength
    except AttributeError as exc:raise RuntimeFailure('native wavelength ABI missing; install matching wheel or explicitly rebuild kernels') from exc
    fn.argtypes=[c.c_int]+[c.c_void_p]*3;fn.restype=c.c_int
    if wavelength.size>2**31-1 or fn(wavelength.size,*(x.ctypes.data for x in (wavelength,omega,energy))):
        raise InputError('wavelength must be nonempty, positive and finite with representable frequency')
    return _readonly(omega),_readonly(energy)


def _table(raw):
    # Deliberately narrow parser for the bundled single 'tabulated nk' YAML block.
    # No general YAML construction or tagged objects are executed.
    text=raw.decode('utf-8');kind=False;active=False;rows=[];blocks=0
    for line in text.splitlines():
        stripped=line.strip()
        if stripped.startswith('- type:'):
            blocks+=1;kind=stripped=='- type: tabulated nk';active=False
        elif stripped=='data: |':
            if not kind:raise InputError('only a single tabulated nk block is supported')
            active=True
        elif active:
            if not stripped:continue
            if not line.startswith('        '):active=False;continue
            fields=stripped.split()
            if len(fields)!=3:raise InputError('expected wavelength,n,k columns')
            rows.append([float(x) for x in fields])
    if blocks!=1 or len(rows)<2:raise InputError('single tabulated nk block required')
    table=_array(rows,'tabulated nk')
    if np.any(table[:,0]<=0) or np.any(np.diff(table[:,0])<=0) or np.any(table[:,1:]<0):
        raise InputError('invalid tabulated wavelengths or n,k')
    return table


@dataclass(frozen=True,slots=True)
class OpticalMaterial:
    omega: np.ndarray
    energy_eV: np.ndarray
    nk: np.ndarray
    _metadata_json: str=field(repr=False)
    _fits_json: str=field(repr=False)

    @property
    def metadata(self):return json.loads(self._metadata_json)
    @property
    def fits(self):return json.loads(self._fits_json)
    @property
    def source_path(self):return Path(str(_root().joinpath(self.metadata['raw_file'])))

    @property
    def epsilon(self):
        """Measured n,k converted to scalar epsilon with the stated mu=1 assumption."""
        n=np.ascontiguousarray(self.nk[:,0]);k=np.ascontiguousarray(self.nk[:,1]);eps=np.empty(len(n),dtype=np.complex128)
        fn=_handles_api().library.sp_ld_nk;fn.argtypes=[c.c_int]+[c.c_void_p]*3;fn.restype=c.c_int
        if fn(len(n),n.ctypes.data,k.ctypes.data,eps.ctypes.data):raise InputError('nonrepresentable measured epsilon')
        return np.frombuffer(eps.tobytes(),dtype=np.complex128)

    @property
    def epsilon_tensor(self):
        """Measured relative epsilon(N,3,3), explicitly isotropic xyz embedding.

        The source provides scalar n,k, not independent tensor measurements.
        This is a model assumption; it does not infer missing crystal axes.
        """
        return _isotropic_tensor(self.epsilon)

    def model(self,*,poles=None):
        """Load stored coefficients without refitting; inaccurate candidates warn.

        Default uses the curated residual/resolution choice, or best passive candidate
        with a warning when no fit reaches the declared tolerance.
        """
        fits=[r for r in self.fits if r['parameters'] is not None]
        if poles is not None:fits=[r for r in fits if r['poles']==poles]
        if not fits:raise InputError('no stored finite model for the requested pole count')
        preferred=self.metadata.get('recommended_fit')
        chosen=next((r for r in fits if r['poles']==preferred),None) if poles is None else fits[0]
        if chosen is None:chosen=min(fits,key=lambda r:r['report']['relative_rms_error'])
        digest=hashlib.sha256(json.dumps(chosen['parameters'],sort_keys=True,separators=(',',':'),allow_nan=False).encode()).hexdigest()
        if digest!=chosen.get('parameters_sha256'):raise InputError('optical fit parameters checksum mismatch')
        report=chosen['report']
        report['catalog_curation']=chosen.get('curation',{})
        report['catalog_quality']=self.metadata.get('quality')
        notes=list(report.get('warnings',[]))
        if not chosen.get('curation',{}).get('eligible_default',False):notes.append('not a curated time-domain default; inspect residual, pointwise errors and fast internal scales')
        if notes:warnings.warn('; '.join(notes),RuntimeWarning,stacklevel=2)
        params={k:v for k,v in chosen['parameters'].items() if k!='units'}
        return LorentzDrude(**params),report


def load_material(material_id):
    """Verify original numerical source bytes, convert units and read fit records."""
    index=_index();entry=next((x for x in index['materials'] if x['id']==material_id),None)
    if entry is None:raise KeyError(material_id)
    name=entry['raw_file']
    if not re.fullmatch(r'raw/[A-Za-z0-9_.-]+\.yml',name):raise InputError('unsafe optical data resource path')
    raw=_root().joinpath(name).read_bytes()
    if hashlib.sha256(raw).hexdigest()!=entry['source_sha256']:raise InputError('optical source checksum mismatch')
    table=_table(raw);omega,energy=wavelength_to_frequency(table[:,0])
    metadata={k:v for k,v in entry.items() if k!='fits'}
    metadata.update(source_commit=index['source_commit'],database_citation=index['database_citation'])
    return OpticalMaterial(_readonly(omega[::-1]),_readonly(energy[::-1]),_readonly(table[::-1,1:]),
                           json.dumps(metadata,ensure_ascii=False),json.dumps(entry['fits']))
