"""Sequential MPI2 reference/streamed seed RSS, including snapshot coefficient I/O."""
from pathlib import Path
import array, hashlib, json, math, os, subprocess, tempfile
here = Path(__file__).resolve().parent
root = here.parents[1]
results = []
with tempfile.TemporaryDirectory() as folder:
    p = Path(folder);(p/'config.h').write_text('')
    subprocess.run(['cc', '-c', str(here/'peak_rss.c'), '-o', str(p/'rss.o')], check=True)
    subprocess.run(['mpifort', '-cpp', '-DUSE_MPI', '-I'+folder, '-O2', '-fexternal-blas', '-fno-tree-loop-vectorize',
                    str(root/'src/xc/exx_wannier_gauge.f90'), str(root/'src/xc/lcfo_dist_rows.f90'),
                    str(root/'src/xc/lcfo_seed.f90'), str(here/'seed_stream_memory_probe.f90'), str(p/'rss.o'),
                    '-L/opt/homebrew/opt/openblas/lib', '-lopenblas', '-o', str(p/'probe')], cwd=p, check=True)
    for n in (256,512):
        reference = None;checksum = None
        for mode in ('reference','streamed'):
            output = subprocess.check_output(['mpirun','-np','2',str(p/'probe'),mode,str(n)], cwd=p,
                universal_newlines=True,env=dict(os.environ,OMP_NUM_THREADS='1',OPENBLAS_NUM_THREADS='1'))
            data = array.array('d');data.frombytes((p/'seed-u.bin').read_bytes())
            digest = hashlib.sha256((p/'seed-coeff.bin').read_bytes()).hexdigest()
            if reference is None:reference=data;checksum=digest
            assert digest==checksum,'Snapshot differs'
            error=max(math.hypot(data[i]-reference[i],data[i+1]-reference[i+1]) for i in range(0,len(data),2))
            assert error<1e-10,error
            for line in output.splitlines():
                fields=line.split()
                results.append(dict(mode=mode,occupied=n,basis=16*n,rank=int(fields[2]),peak_rss_bytes=int(fields[3]),
                                    unitarity_error=float(fields[4]),rotation_difference=error,coefficient_sha256=digest))
print(json.dumps(results,indent=2))
