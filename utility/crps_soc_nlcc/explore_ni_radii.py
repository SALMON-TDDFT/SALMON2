#!/usr/bin/env python3
"""Preserved larger-Rps exploration: atomic tables, shape metrics and radial screening.

No production SALMON or material GS calculations are launched. Reuses an audited
unchanged utility executable only after verifying its recorded SHA256.
"""
import argparse
import hashlib
import importlib.util
import json
import os
from pathlib import Path
import subprocess
import sys
import numpy as np


def shape(path):
    lines=path.read_text().splitlines();nc=int(lines[3].split()[2]);ch=np.array([list(map(float,x.split())) for x in lines[4:4+nc]])
    a=np.array([list(map(float,x.split())) for x in lines[4+nc:-1]]);r=a[:,0];h=r[1];out=[]
    for k,c in enumerate(ch):
        n=int(round(c[2]/h));v=a[:n+1,5+2*k];d=np.diff(v)/h;dd=np.diff(v,2)/h**2
        # Ignore forced origin duplication; count resolved stationary points.
        ds=d[1:];ds=ds[abs(ds)>1e-3]
        out.append(dict(l=int(c[0]),two_j=int(c[1]),rps_angstrom=c[2],minimum_eV=float(min(v)),maximum_eV=float(max(v)),
            interior_stationary_points=int(np.sum(ds[1:]*ds[:-1]<0)),total_variation_eV=float(sum(abs(np.diff(v)))),
            excess_variation_eV=float(sum(abs(np.diff(v)))-abs(v[-1]-v[0])),
            curvature_rms_eV_angstrom2=float(np.sqrt(np.mean(dd[2:]**2))),curvature_p99_eV_angstrom2=float(np.quantile(abs(dd[2:]),.99))))
    return out


def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('executable',type=Path);p.add_argument('recorded_sha256');p.add_argument('destination',type=Path);a=p.parse_args()
    assert hashlib.sha256(a.executable.read_bytes()).hexdigest()==a.recorded_sha256,'unchanged validated executable required'
    a.destination.mkdir(exist_ok=False);source=Path(__file__).resolve().parent;base=(source/'inputs/Ni.nml').read_text()
    (a.destination/'baseline-input.nml').write_text(base)
    candidates=[(s,2.,2.) for s in [2.6,2.8,3.,3.2]]+[(2.2,p,2.) for p in [2.4,2.8,3.2]]+[(2.2,2.,d) for d in [2.4,2.8,3.2]]+[(2.8,2.8,2.8),(3.2,3.2,3.2),(2.8,3.2,3.2)]
    records=[]
    for s,pp,d in candidates:
        name=f's{s:.1f}-p{pp:.1f}-d{d:.1f}-local0';run=a.destination/name;run.mkdir()
        text=base.replace('cutoff_radius_bohr=2.2,2.0,2.0,2.0,2.0',f'cutoff_radius_bohr={s},{pp},{pp},{d},{d}')
        assert text!=base
        (run/'Ni.nml').write_text(text)
        with (run/'stdout.txt').open('w') as out,(run/'stderr.txt').open('w') as err:
            q=subprocess.run([str(a.executable),'Ni.nml'],cwd=run,stdout=out,stderr=err,timeout=120)
        log=(run/'stdout.txt').read_text();r=dict(name=name,radii_bohr=[s,pp,d],returncode=q.returncode,written='CRPS_SOC_NLCC_WRITTEN' in log,
            rejections=[x.strip() for x in log.splitlines() if 'REJECT' in x])
        if r['written']:
            r['shape']=shape(run/'Ni_soc_nlcc_rps.dat')
            q=subprocess.run([sys.executable,str(source/'screen_ni_radial.py'),str(run/'Ni_soc_nlcc_rps.dat'),'--output',str(run/'radial-screen.json'),'--grids','8:.04,8:.02'],check=True,capture_output=True,text=True)
            screen=json.loads((run/'radial-screen.json').read_text());r['deep_extra_state_flags']=sum(x['deep_extra_state_flag'] for x in screen['checks'])
            r['screened_additional_bound_state_flags']=sum(x['kb_bound_count']>x['semilocal_bound_count'] for x in screen['checks'] if x['screened'])
            r['max_reference_error_eV_h02']=max(abs(x['semilocal_nearest_reference_eV']-x['reference_eV']) for x in screen['checks'] if x['screened'] and x['step_angstrom']==.02)
        records.append(r);print(json.dumps({k:v for k,v in r.items() if k!='shape'}),flush=True)
        (a.destination/'exploration.json').write_text(json.dumps(records,indent=2)+'\n')

if __name__=='__main__':main()
