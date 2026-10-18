#!/usr/bin/env python3
"""Compare SALMON 13-column *_response.data, using sheet response Lz*(eps-1)."""
import argparse
import cmath
import json
import math
from pathlib import Path

def read(path):
    rows=[]
    for line in path.read_text().splitlines():
        if not line.strip() or line.lstrip().startswith('#'):continue
        row=[float(x) for x in line.split()]
        assert len(row)>=13 and all(math.isfinite(x) for x in row),'invalid response table'
        rows.append(row)
    return rows

def main():
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument('first',type=Path);ap.add_argument('second',type=Path)
    ap.add_argument('--axis',choices=['x','y'],default='x')
    ap.add_argument('--second-axis',choices=['x','y'])
    ap.add_argument('--lz',type=float,required=True);ap.add_argument('--second-lz',type=float,required=True)
    ap.add_argument('--emin',type=float,default=1.0);ap.add_argument('--emax',type=float,default=6.0)
    args=ap.parse_args()
    a,b=read(args.first),read(args.second);assert len(a)==len(b),'use identical output energy meshes'
    i='xy'.index(args.axis);j='xy'.index(args.second_axis or args.axis)
    values=[]
    for r,s in zip(a,b):
        assert abs(r[0]-s[0])<1e-9,'energy mesh differs'
        e=r[0]*27.211386245988
        if args.emin<=e<=args.emax:
            values.append((e,args.lz*(complex(r[7+i],r[10+i])-1),args.second_lz*(complex(s[7+j],s[10+j])-1)))
    assert values,'no data in energy window'
    scale=max(abs(x[1]) for x in values);assert scale>0
    relative=max(abs(x[1]-x[2]) for x in values)/scale
    peak1=max(values,key=lambda x:x[1].imag)[0];peak2=max(values,key=lambda x:x[2].imag)[0]
    print(json.dumps({'sheet_response_max_relative_difference':relative,
       'dominant_peak_shift_eV':abs(peak1-peak2),'window_eV':[args.emin,args.emax],
       'candidate_acceptance':relative<.02 and abs(peak1-peak2)<.02,
       'note':'Candidate criterion only; inspect individual A/B and higher peaks, integrated weights and shoulders.'},indent=2))

if __name__=='__main__':main()
