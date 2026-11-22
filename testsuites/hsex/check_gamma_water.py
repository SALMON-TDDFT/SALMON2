#!/usr/bin/env python3
"""Offline assertions for the small periodic-water Gamma PDEP campaign.

SALMON performs all screening, MPI and cache I/O in native code. No GS payload
is copied by this checker. The expected stages are independent new run folders.
"""
from pathlib import Path
import argparse,json,re
import numpy as np


def read_w(path):
    data=path.read_bytes()
    if data[:4] not in (b'HXL1',b'HXT1'):
        raise AssertionError('Unexpected Gamma W format')
    d=np.frombuffer(data,'=i4',9,4)
    ng=int(np.prod(d[:3]));nk=int(np.prod(d[3:6]));na=int(d[6])**3;ni=int(d[8])
    if nk!=1 or d[7]!=4:raise AssertionError('Not the occupied-only water Gamma fixture')
    offset=40+8*(9+3*nk+3*ni+ng+1+3*nk)
    w=np.frombuffer(data,'=f8',na*na,offset).reshape(na,na)
    if len(data)!=offset+8*na*na or not np.isfinite(w).all():raise AssertionError('Invalid Gamma W payload')
    np.testing.assert_allclose(w,w.T,atol=1e-10,rtol=0)
    return w


def check(root):
    result={}
    reference=read_w(root/'reference/hsex-w.bin')
    for name in ('reference','full','warm','reduced','warm-reduced','auto','serial-full'):
        directory=root/name;text=(directory/'run.log').read_text()
        if not re.search(r'^end SALMON$',text,re.M) or 'HSEX_GAMMA_SCREEN_PASS empty_band_sum=0' not in text:
            raise AssertionError(f'Incomplete Gamma screening: {name}')
        if 'number of threads :          12' not in text:
            raise AssertionError('Requested OMP team was not observed')
        actual=read_w(directory/'hsex-w.bin')
        tensor=[]
        for row in range(1,4):
            tensor.append([float(x) for x in re.search(r'HSEX_GAMMA_EPS row='+str(row)+r'\s+([^\n]+)',text)[1].split()])
        dimensions=[int(x) for x in re.search(r'HSEX_GAMMA_PDEP dimensions G\+optical/rank/actions\s+([^\n]+)',text)[1].split()]
        result[name]=dict(W_max_absolute_error=float(np.max(abs(actual-reference))),
                          W_relative_norm_error=float(np.linalg.norm(actual-reference)/np.linalg.norm(reference)),
                          dielectric_tensor=tensor,dimensions_rank_actions=dimensions)
    for name in ('full','warm','serial-full'):
        np.testing.assert_allclose(read_w(root/name/'hsex-w.bin'),reference,atol=2e-9,rtol=2e-9)
    np.testing.assert_allclose(read_w(root/'warm-reduced/hsex-w.bin'),read_w(root/'reduced/hsex-w.bin'),atol=2e-8,rtol=2e-8)
    assert result['auto']['dimensions_rank_actions'][1]==8
    for name in ('warm','warm-reduced'):
        assert 'HSEX_GAMMA_PDEP_CACHE_READ_PASS' in (root/name/'run.log').read_text()
    for name in ('full','reduced'):
        assert 'HSEX_GAMMA_PDEP_CACHE_WRITE_PASS' in (root/name/'run.log').read_text()
    return result


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('work',type=Path);parser.add_argument('--output',type=Path)
    args=parser.parse_args();result=check(args.work);text=json.dumps(result,indent=2)+'\n'
    if args.output:args.output.write_text(text)
    print(text)
