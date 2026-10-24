#!/usr/bin/env python3
"""Compare SALMON *_rt_spin.data files without resampling or changing spin units."""
import argparse
import csv
import json
import math
from pathlib import Path
import re


def read_spin(path):
    """Read native spin output; preserve its spin convention (no mu_B conversion)."""
    times, vectors, unit = [], [], None
    for line in Path(path).read_text().splitlines():
        line=line.strip()
        if not line:
            continue
        if line.startswith('#'):
            match=re.search(r'1:Time\[([^]]+)\]',line)
            if match:
                unit=match.group(1)
            continue
        fields=line.replace('D','E').replace('d','e').split()
        if len(fields)<4:
            raise ValueError(f'{path}: expected time and three magnetization components')
        row=[float(value) for value in fields[:4]]
        if not all(math.isfinite(value) for value in row):
            raise ValueError(f'{path}: nonfinite data')
        if times and row[0]<=times[-1]:
            raise ValueError(f'{path}: times must be strictly increasing')
        times.append(row[0]); vectors.append(row[1:])
    if not times or unit is None:
        raise ValueError(f'{path}: missing data or SALMON Time unit header')
    return {'time':times,'vectors':vectors,'time_unit':unit,'source':str(Path(path).resolve())}


def compare(runs):
    if not runs:
        raise ValueError('at least one run is required')
    ref=next(iter(runs.values()))
    result={'time':ref['time'],'time_unit':ref['time_unit'],'reference_time':ref['time'][0],
            'magnetization_units':'native SALMON spin output; no magnetic-moment conversion', 'runs':{}}
    for name,run in runs.items():
        if run['time_unit']!=ref['time_unit'] or len(run['time'])!=len(ref['time']):
            raise ValueError(f'{name}: incompatible time unit or number of samples')
        if any(not math.isclose(a,b,rel_tol=1e-12,abs_tol=1e-12) for a,b in zip(run['time'],ref['time'])):
            raise ValueError(f'{name}: time grids differ; no implicit interpolation is performed')
        magnitudes=[math.hypot(*v) for v in run['vectors']]
        initial=magnitudes[0]
        result['runs'][name]={'source':run['source'],'components':run['vectors'],'magnitude':magnitudes,
            'delta_magnitude':[x-initial for x in magnitudes],
            'normalized_magnitude':[x/initial for x in magnitudes] if initial>1e-14 else None}
    return result


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--run',action='append',required=True,metavar='LABEL=FILE')
    parser.add_argument('--output',type=Path,required=True,help='new output directory in canonical calc2026')
    parser.add_argument('--plot',action='store_true',help='also create PNG plots (requires matplotlib)')
    args=parser.parse_args()
    paths={}
    for item in args.run:
        label,sep,path=item.partition('=')
        if not sep or not label or label in paths:
            parser.error('--run requires unique LABEL=FILE entries')
        paths[label]=read_spin(path)
    result=compare(paths)
    if args.plot:
        import matplotlib
        matplotlib.use('Agg')
        import matplotlib.pyplot as plt
    args.output.mkdir(parents=True,exist_ok=False)
    (args.output/'comparison.json').write_text(json.dumps(result,indent=2,allow_nan=False)+'\n')
    with (args.output/'comparison.csv').open('w',newline='') as stream:
        writer=csv.writer(stream)
        writer.writerow(['time','label','Mx','My','Mz','magnitude','delta_magnitude','normalized_magnitude'])
        for label,run in result['runs'].items():
            for i,time in enumerate(result['time']):
                norm=run['normalized_magnitude']
                writer.writerow([time,label,*run['components'][i],run['magnitude'][i],
                                 run['delta_magnitude'][i],norm[i] if norm is not None else ''])
    if args.plot:
        fig,axes=plt.subplots(2,2,figsize=(10,7),constrained_layout=True)
        for label,run in result['runs'].items():
            for component,style in enumerate(['-', '--', ':']):
                axes[0,0].plot(result['time'],[v[component] for v in run['components']],style,
                               label=f'{label} M{"xyz"[component]}')
            axes[0,1].plot(result['time'],run['magnitude'],label=label)
            axes[1,0].plot(result['time'],run['delta_magnitude'],label=label)
            if run['normalized_magnitude'] is not None:
                axes[1,1].plot(result['time'],run['normalized_magnitude'],label=label)
        for ax,title in zip(axes.flat,['Components','Magnitude','Change from reference magnitude','Normalized to reference magnitude']):
            ax.set_title(title); ax.set_xlabel(f'Time [{result["time_unit"]}]'); ax.grid(alpha=.25)
            if ax.lines: ax.legend(fontsize=8)
        fig.suptitle(f'Reference time: {result["reference_time"]:g} {result["time_unit"]}')
        fig.savefig(args.output/'magnetization.png',dpi=180)
        plt.close(fig)
    print(args.output.resolve())

if __name__=='__main__':
    main()
