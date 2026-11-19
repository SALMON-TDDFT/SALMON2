"""Offline grid/dt planning only; illustrative radius/hash, not a real Si setup."""
import argparse,json
from pathlib import Path
from salmonpy.convergence import plan_grid_scan,plan_dt_scan

def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('output',type=Path);a=p.parse_args()
    spatial=plan_grid_scan([[5.43,0,0],[0,5.43,0],[0,0,5.43]],material_id='example-only',
        radii={'example':{'radius_angstrom':1.2,'definition':'illustrative explicit cutoff; replace with actual PP channel','pp_sha256':'0'*64}},
        metadata={'illustrative_only':True,'xc':'pz','kpoints':[1,1,1]})
    temporal=plan_dt_scan(material_id='example-only',impulse=[1e-4,0,0])
    with a.output.open('x') as stream:json.dump({'grid':spatial,'dt':temporal},stream,indent=2);stream.write('\n')
    print(a.output)

if __name__=='__main__':main()
