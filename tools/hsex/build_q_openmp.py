#!/usr/bin/env python3
"""Build optional q kernel in a NEW directory; audit source/objects, record hash.

On non-Apple platforms set FC, Q_OMP_FLAGS and Q_BLAS_FLAGS explicitly for
an OpenMP compiler and a sequential, reentrant LP64 BLAS/LAPACK library.
"""
import argparse
import hashlib
import json
import os
from pathlib import Path
import platform
import re
import shlex
import subprocess


def build(directory):
    root=Path(__file__).resolve().parents[2]
    checker=root/'tools/check_forbidden_eigensolvers.py'
    native=root/'tools/hsex/native'
    for source in (root/'src',native):
        subprocess.run(['python3',str(checker),'--source-root',str(source)],check=True)
    directory=directory.resolve();directory.mkdir(parents=True,exist_ok=False)
    apple=platform.system()=='Darwin'
    flags=shlex.split(os.environ.get('Q_OMP_FLAGS','-fopenmp'))
    compiler=shlex.split(os.environ.get('FC','gfortran'))
    libraries=shlex.split(os.environ.get('Q_BLAS_FLAGS','-framework Accelerate' if apple else ''))
    compile_flags=shlex.split(os.environ.get('Q_COMPILE_FLAGS','-O3 -fPIC -cpp'))
    for token in compiler+flags+libraries+compile_flags:
        if '@' in token or re.search(r'\.(?:o|obj)(?:$|[, ])',token,re.I):
            raise ValueError('external application objects and response files are not allowed by this builder')
    if not libraries:raise ValueError('set Q_BLAS_FLAGS to a sequential reentrant LP64 BLAS/LAPACK')
    obj=directory/'q_response.o'
    cmd=compiler+compile_flags+flags+(['-DHSEX_ACCELERATE'] if apple else [])
    subprocess.run(cmd+['-c',str(native/'q_response.f90'),'-o',str(obj)],cwd=directory,check=True)
    subprocess.run(['python3',str(checker),'--object-tree',str(directory)],check=True)
    output=directory/('libhsex_q.dylib' if apple else 'libhsex_q.so')
    subprocess.run(compiler+['-dynamiclib' if apple else '-shared']+flags+[str(obj)]+libraries+['-o',str(output)],check=True)
    record=dict(compiler=compiler,compile=cmd,blas=libraries,sha256=hashlib.sha256(output.read_bytes()).hexdigest(),
                source_sha256=hashlib.sha256((native/'q_response.f90').read_bytes()).hexdigest())
    (directory/'audit.json').write_text(json.dumps(record,indent=2)+'\n')
    print(output)

if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('directory',type=Path)
    build(p.parse_args().directory)
