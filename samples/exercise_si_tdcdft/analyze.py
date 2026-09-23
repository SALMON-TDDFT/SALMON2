"""Finite-window transverse response from SALMON atomic-unit electron-number currents.

Use identical GS, grids, dt, field and xc settings for pump and pump+probe.
The probe amplitude is the signed step in A/c, not the electric-field integral.
"""
import argparse
import json
from pathlib import Path
import numpy as np

HARTREE_EV=27.211386245988

def difference(time, current, reference_time, reference_current):
    if current.shape!=time.shape or reference_current.shape!=reference_time.shape:
        raise ValueError("Current and time arrays must have identical shapes")
    if time.shape!=reference_time.shape or not np.allclose(time,reference_time,atol=1e-10,rtol=0):
        raise ValueError('Pump and probe time grids must match exactly')
    return current-reference_current

def response(time, current, impulse, energy_eV, probe_time=0):
    time=np.asarray(time); current=np.asarray(current); energy_eV=np.asarray(energy_eV)
    if len(time)<3 or current.shape!=time.shape or not np.all(np.isfinite(time+current)):
        raise ValueError('Need matching finite time/current arrays')
    dt=time[1]-time[0]
    if dt<=0 or not np.allclose(np.diff(time),dt,atol=1e-10,rtol=1e-8):
        raise ValueError('Need a uniform increasing time grid')
    if not np.isfinite(impulse) or impulse==0 or not np.isfinite(probe_time):
        raise ValueError('Need a finite nonzero signed impulse and finite probe time')
    if np.any(energy_eV<=0) or not np.all(np.isfinite(energy_eV)):
        raise ValueError('Energy must be positive and finite')
    keep=time>probe_time+1e-10
    t=time[keep]-probe_time
    if len(t)<3: raise ValueError('Insufficient post-probe time samples')
    if abs(t[0]-dt)>dt+1e-10: raise ValueError('Current trace does not start at the probe')
    window=1-3*(t/t[-1])**2+2*(t/t[-1])**3
    weighted=current[keep]*window*dt/impulse
    # Chunk energies to keep memory bounded for long propagation traces.
    sigma=np.empty(energy_eV.size,dtype=complex)
    omega=energy_eV/HARTREE_EV
    for start in range(0,len(omega),64):
        sigma[start:start+64]=np.exp(1j*omega[start:start+64,None]*t)@weighted
    return 1+4*np.pi*1j*sigma/omega

def peak_metrics(energy, absorption, low, high):
    mask=(energy>=low)&(energy<=high)
    if np.count_nonzero(mask)<3: raise ValueError('Peak interval must contain at least three points')
    x=energy[mask]; y=absorption[mask]; imax=np.argmax(y)
    return dict(peak_eV=float(x[imax]),height=float(y[imax]),area_eV=float(np.trapezoid(y,x)),
                boundary_peak=bool(imax==0 or imax==len(x)-1))

def read_current(path):
    with path.open() as stream:
        header=''.join(next(stream,'') for _ in range(7))
    if 'Time[a.u.]' not in header or 'Jm_z[a.u.]' not in header:
        raise ValueError('Expected SALMON rt.data in atomic units')
    data=np.loadtxt(path)
    if data.ndim!=2 or data.shape[1]<16:
        raise ValueError('Expected time, fields and current columns')
    return data

def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('current',type=Path)
    p.add_argument('--subtract',type=Path,help='Pump-only rt.data on the identical time grid')
    p.add_argument('--impulse',type=float,required=True,help='Signed A/c step in atomic units')
    p.add_argument('--probe-time',type=float,default=0,help='Probe threshold in atomic time units')
    p.add_argument('--axis',choices=['x','y','z'],default='z')
    p.add_argument('--emax',type=float,default=8)
    p.add_argument('--de',type=float,default=.01)
    p.add_argument('--output',type=Path,required=True)
    p.add_argument('--peak-range',nargs=2,type=float,help='Report peak and integrated Im(epsilon) in eV interval')
    a=p.parse_args()
    if a.de<=0 or a.emax<=a.de: p.error('Need emax > de > 0')
    raw=read_current(a.current); col=13+'xyz'.index(a.axis)
    time=raw[:,0]; current=raw[:,col]
    if a.subtract:
        ref=read_current(a.subtract)
        current=difference(time,current,ref[:,0],ref[:,col])
    energy=np.arange(a.de,a.emax+a.de*.1,a.de)
    eps=response(time,current,a.impulse,energy,a.probe_time)
    metadata=dict(current=str(a.current),subtract=str(a.subtract) if a.subtract else None,
                  impulse_au=a.impulse,probe_time_au=a.probe_time,axis=a.axis,
                  window='1-3(t/T)^2+2(t/T)^3',end_time_au=float(time[-1]),
                  unit_system='a.u.',energy_shift_eV=0,additional_broadening_eV=0)
    if a.peak_range: metadata['peak']=peak_metrics(energy,eps.imag,*a.peak_range)
    np.savetxt(a.output,np.column_stack([energy,eps.real,eps.imag]),header='energy_eV Re_epsilon Im_epsilon')
    a.output.with_suffix(a.output.suffix+'.json').write_text(json.dumps(metadata,indent=2)+'\n')
    print(json.dumps(metadata,indent=2))

if __name__=='__main__': main()
