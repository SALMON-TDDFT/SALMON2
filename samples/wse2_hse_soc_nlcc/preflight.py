#!/usr/bin/env python3
"""Static package check; never runs SALMON or submits jobs. Python 3.6 compatible."""
import argparse
import hashlib
import json
import math
from pathlib import Path
import re

def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def check_factor_mode(text):
    # Accept the documented omitted/default value; inspect only &functional.
    text = re.sub(r'!.*', '', text)
    block = re.search(r'&functional\b(.*?)/', text, re.I | re.S)
    assert block, 'functional namelist missing'
    assignments = re.findall(r"\bexx_factor_mode\s*=\s*([^,\n]+)", block.group(1), re.I)
    mode = 'off'
    for value in assignments:
        token = re.fullmatch(r"(['\"])(.*?)\1\s*", value.strip())
        assert token, 'exx_factor_mode must be a quoted value'
        mode = token.group(2).lower()
    assert mode == 'off', 'WSe2 SOC input requires factor history off'

def main():
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--case',type=Path)
    ap.add_argument('--binary',type=Path)
    args=ap.parse_args();root=Path(__file__).resolve().parent
    manifest=json.loads((root/'manifest.json').read_text())
    for name,p in manifest['potentials'].items():
        f=root/'pseudo'/(name+'_PBE19.vps');s=f.read_text()
        assert sha(f)==p['sha256'],'pseudo hash mismatch'
        assert re.search(r'j.dependent.pseudo.potentials\s+on',s),'SOC channels missing'
        assert re.search(r'charge.pcc.calc\s+on',s),'NLCC disabled'
        vals=s.split('<density.PCC',1)[1].split('density.PCC>',1)[0].splitlines()
        density=[float(v.split()[2]) for v in vals if len(v.split())>=3]
        assert density and max(density)>0,'empty NLCC density'
        assert int(float(re.search(r'valence.electron\s+([0-9.]+)',s).group(1)))==p['valence']
    cases=[json.loads((args.case.resolve()/'case.json').read_text())] if args.case else manifest['cases']
    for c in cases:
        case=root/c['name'];assert case.is_dir()
        nk=c['k'][0]*c['k'][1]*c['k'][2]
        assert c['k'][0]%6==0 and c['k'][1]%2==0 and c['k'][2]==1
        assert all(x>=8 and x%2==0 for x in c['grid'])
        # Same uniform-index construction as lattice.f90; includes Gamma and folded K.
        xs=[i/c['k'][0]-.5 for i in range(1,c['k'][0]+1)]
        ys=[i/c['k'][1]-.5 for i in range(1,c['k'][1]+1)]
        assert any(abs(x-round(x))<1e-12 for x in xs)
        assert any(abs(y-round(y))<1e-12 for y in ys)
        assert all(any(abs(x-v-round(x-v))<1e-12 for x in xs) for v in (-1/3,1/3))
        for phase in ('gs','gs-x','rt-x','rt-y'):
            s=(case/phase/'inputfile').read_text()
            assert "spin='noncollinear'" in s and "yn_spinorbit='y'" in s
            mask='nyy' if phase in ('gs-x','rt-x') else ('yny' if phase=='rt-y' else 'yyy')
            assert "yn_symmetry='%s'"%mask in s and (case/phase/'sym.dat').is_file()
            assert 'nelec=48,nstate=48' in s and "yn_dc='n'" in s
            assert 'hse_block_rows=1' in s
            check_factor_mode(s)
            for path in re.findall(r"file_pseudo\(\d\)='([^']+)'",s):
                assert (case/phase/path).resolve().is_file()
            assert (case/phase/'job.sh').is_file()
            ranks=c['stage_ranks'][phase]
            assert ranks<=c['irreducible_kpoints'][phase]
            assert 'nproc_k=%d,'%ranks in s
            ops=[]
            rows=[[float(x) for x in row.split()] for row in (case/phase/'sym.dat').read_text().splitlines()]
            assert len(rows)%3==0
            for i in range(0,len(rows),3):
                r=rows[i:i+3]
                diag=[r[j][j] for j in range(3)];t=[r[j][3] for j in range(3)]
                assert all(abs(abs(x)-1)<1e-12 for x in diag)
                assert all(abs(r[a][b])<1e-12 for a in range(3) for b in range(3) if a!=b)
                assert all(abs(t[j]*c['grid'][j]-round(t[j]*c['grid'][j]))<1e-10 for j in range(3))
                if phase in ('rt-x','gs-x'):assert diag[0]==1
                if phase=='rt-y':assert diag[1]==1
                ops.append((diag,t))
            coords=[(line.split()[0],list(map(float,[x.replace('d','e') for x in line.split()[1:4]]))) for line in s.split('&atomic_red_coor')[1].split('/')[0].splitlines() if line.strip()]
            for diag,t in ops:
                for symbol,xyz in coords:
                    new=[diag[j]*xyz[j]+t[j] for j in range(3)]
                    assert any(other==symbol and all(abs(new[j]-q[j]-round(new[j]-q[j]))<1e-10 for j in range(3)) for other,q in coords),'operation breaks atomic structure'
                for d2,t2 in ops:
                    dc=[diag[j]*d2[j] for j in range(3)]
                    tc=[diag[j]*t2[j]+t[j] for j in range(3)]
                    assert any(dd==dc and all(abs(tc[j]-tt[j]-round(tc[j]-tt[j]))<1e-10 for j in range(3)) for dd,tt in ops),'group is not closed'

        assert all(v['planning_estimate_24_matrices_plus_tiles']<27 for v in c['stage_memory_model_GiB'].values()),'memory estimate exceeds 28GiB; increase k ranks'
    if args.binary:
        pin=root/'binary.json'
        assert pin.is_file(),'New binary must be pinned in binary.json; old Si binary is prohibited'
        b=json.loads(pin.read_text());assert sha(args.binary.resolve())==b['sha256'],'binary mismatch'
        assert b['source_commit'].startswith(manifest['engine_commit']),'engine version mismatch'
    print('PASS: %d cases, original SOC/NLCC potentials, 48 occupied spinors, uniform Gamma/K meshes'%len(cases))

if __name__=='__main__':main()
