"""Explicit small Fortran coordinate probe; no MPI, GS or network.
Builds this isolated probe only when the CLI is explicitly invoked.

NEW_OUTPUT [--source-root SOLVER] [--compiler gfortran]
Extracts the exact dependency-free routine; audits sources/objects before linking.
"""
import argparse
import hashlib
import json
from pathlib import Path
import shlex
import subprocess
import sys


def main():
    root=Path(__file__).resolve().parents[1]
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('output',type=Path)
    parser.add_argument('--source-root',type=Path,default=root/'runtime/sources/solver')
    parser.add_argument('--compiler',default='gfortran')
    args=parser.parse_args();source=args.source_root.resolve();out=args.output.resolve()
    if out.is_relative_to(root) or out.is_relative_to(source):
        parser.error('new external output directory required')
    checker=root/'tools/check_forbidden_eigensolvers.py'
    subprocess.run([sys.executable,str(checker),'--source-root',str(source/'src')],check=True)
    origin=source/'src/maxwell/common_maxwell.f90';text=origin.read_text()
    begin=text.index('  subroutine set_coo_em(')
    end=text.index('  end subroutine set_coo_em',begin)+len('  end subroutine set_coo_em')
    out.mkdir(parents=True,exist_ok=False)
    (out/'coordinate_kernel.f90').write_text('module coordinate_kernel\ncontains\n'+text[begin:end]+'\nend module\n')
    (out/'probe.f90').write_bytes((root/'tests/native/test_fdtd_coordinates.f90').read_bytes())
    commands=[]
    with (out/'build.log').open('w') as log:
        def run(cmd):
            commands.append([str(x) for x in cmd])
            subprocess.run(commands[-1],cwd=out,stdout=log,stderr=subprocess.STDOUT,check=True)
        run([sys.executable,checker,'--source-root',out])
        compiler=shlex.split(args.compiler)
        for name in ('coordinate_kernel','probe'):
            run([*compiler,'-c',name+'.f90','-o',name+'.o'])
        run([sys.executable,checker,'--object-tree',out])
        run([*compiler,'coordinate_kernel.o','probe.o','-o','probe'])
        result=subprocess.run([str(out/'probe')],cwd=out,stdout=log,stderr=subprocess.STDOUT)
    (out/'result.json').write_text(json.dumps(dict(status='passed' if result.returncode==0 else 'failed',
        returncode=result.returncode,executable_sha256=hashlib.sha256((out/'probe').read_bytes()).hexdigest(),
        origin=str(origin),origin_sha256=hashlib.sha256(origin.read_bytes()).hexdigest(),
        commands=commands,scope='native coordinate/halo arithmetic and deterministic unused padding'),indent=2)+'\n')
    raise SystemExit(result.returncode)


if __name__=='__main__':main()
