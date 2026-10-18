#!/usr/bin/env python3
"""Require actual convergence/completion and every rank exit; do not infer from restart files."""
import argparse
import json
import math
from pathlib import Path
import re

def main():
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument('folder',type=Path)
    ap.add_argument('--phase',choices=['gs','rt'],required=True)
    args=ap.parse_args();folder=args.folder.resolve()
    case=json.loads((folder.parent/'case.json').read_text())
    assert (folder/'launcher-exit.txt').read_text().strip()=='0','launcher did not finish successfully'
    for rank in range(case['stage_ranks'][folder.name]):
        assert (folder/('exit.rank-%d'%rank)).read_text().strip()=='0','rank failure'
    log=(folder/'stdout.rank-0').read_text()
    assert 'end SALMON' in log and 'ERROR STOP' not in log,'incomplete SALMON log'
    assert not re.search(r'(?i)(?<![a-z])(nan|infinity)(?![a-z])',log),'nonfinite log'
    if args.phase=='gs':
        rs=re.findall(r'iter and int_x.*?=\s*\d+\s+([0-9.Ee+-]+)',log)
        assert rs and math.isfinite(float(rs[-1])) and float(rs[-1])<1e-8,'unconverged GS'
        assert '#GS converged at' in log,'missing GS convergence'
        ns=re.findall(r'Ne=\s*([0-9.Ee+-]+)',log)
        assert ns and abs(float(ns[-1])-48)<1e-6,'incorrect electron count'
        assert (folder/'data_for_restart').is_dir(),'missing GS restart data'
        if case['pre_scf']>0:
            assert 'EXX_PRE_SCF switch to target hybrid' in log,'did not reach HSE06'
        core=re.findall(r'Maximal rho_NLCC=\s*([0-9.Ee+-]+)',log)
        assert core and float(core[-1])>0,'NLCC not active'
    else:
        lines=[s.split() for s in (folder/'WSe2_rt.data').read_text().splitlines() if s.strip() and not s.lstrip().startswith('#')]
        assert len(lines)>=round(case['duration_fs']/case['dt_fs']),'incomplete current trajectory'
        for row in lines:
            assert all(math.isfinite(float(x)) for x in row),'nonfinite current'
        # Input/output are atomic units: compare elapsed time after conversion to fs.
        elapsed=float(lines[-1][0])*0.024188843265857
        assert abs(elapsed-case['duration_fs'])<case['dt_fs']*1.1,'incorrect final RT time'
    print('PASS '+args.phase+' '+str(folder))

if __name__=='__main__':main()
