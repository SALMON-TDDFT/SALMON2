"""Compile the actual model input ownership record and round-trip test offline."""
import argparse
import json
from pathlib import Path
import subprocess
import sys


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('source',type=Path)
    parser.add_argument('output',type=Path)
    parser.add_argument('--fc',default='gfortran')
    parser.add_argument('--checker',type=Path,default=Path(__file__).with_name('check_forbidden_eigensolvers.py'))
    args=parser.parse_args()
    source=args.source.resolve();out=args.output.resolve();root=Path(__file__).resolve().parents[1]
    checker=args.checker.resolve()
    subprocess.run([sys.executable,str(checker),'--source-root',str(source/'src')],check=True)
    out.mkdir(parents=True,exist_ok=False)
    files=[source/'src/io/salmon_global.f90',source/'src/io/model_inputs.f90',root/'tests/native/test_model_inputs.f90']
    options=['-cpp']
    if 'SALMONPY_FAST_CORE' in (source/'src/config.h.in').read_text():options.append('-DSALMONPY_FAST_CORE')
    objects=[]
    for i,path in enumerate(files):
        obj=out/('part-'+str(i)+'.o')
        subprocess.run([args.fc,*options,'-J'+str(out),'-I'+str(out),'-c',str(path),'-o',str(obj)],check=True)
        objects.append(str(obj))
    subprocess.run([sys.executable,str(checker),'--object-tree',str(out)],check=True)
    executable=out/'test-model-inputs'
    subprocess.run([args.fc,*objects,'-o',str(executable)],check=True)
    result=subprocess.run([str(executable)],check=True,capture_output=True,text=True)
    (out/'result.json').write_text(json.dumps({'status':'pass','stdout':result.stdout},indent=2)+'\n')
    print(result.stdout.strip())


if __name__=='__main__':main()
