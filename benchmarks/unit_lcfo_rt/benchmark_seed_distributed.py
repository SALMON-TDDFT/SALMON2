"""Sequential MPI2/MPI4 root/distributed seed RSS, including snapshot coefficient I/O."""
from pathlib import Path
import array, hashlib, json, math, os, subprocess, tempfile
root = Path(__file__).resolve().parents[2]
here = root / "developer_tests" / "652_dc_lcfo/rt"
results = []
with tempfile.TemporaryDirectory() as folder:
    p = Path(folder);(p/'config.h').write_text('')
    subprocess.run(['cc', '-c', str(here/'peak_rss.c'), '-o', str(p/'rss.o')], check=True)
    subprocess.run(['mpifort', '-cpp','-ffree-line-length-none','-fallow-argument-mismatch', '-DUSE_MPI', '-DUSE_SCALAPACK', '-I'+folder, '-O2', '-fexternal-blas', '-fno-tree-loop-vectorize',
                    str(root/'src/misc/nvtx_wrapper.f90'), str(root/'src/parallel/communication.f90'),
                    str(root/'src/xc/exx_wannier_gauge.f90'), str(root/'src/xc/lcfo_dist_rows.f90'),
                    str(root/'src/xc/lcfo_seed.f90'), str(here/'seed_stream_memory_probe.f90'), str(p/'rss.o'),
                    '-L/opt/homebrew/opt/scalapack/lib','-lscalapack',
                    '-L/opt/homebrew/opt/openblas/lib', '-lopenblas', '-o', str(p/'probe')], cwd=p, check=True)
    for ranks,n in ((2,256),(2,512),(4,512)):
        reference = None;checksum = None
        for mode in ('root','distributed'):
            output = subprocess.check_output(['mpirun','-np',str(ranks),str(p/'probe'),'streamed',str(n),mode], cwd=p,
                universal_newlines=True,env=dict(os.environ,OMP_NUM_THREADS='1',OPENBLAS_NUM_THREADS='1'))
            data = array.array('d');data.frombytes((p/'seed-u.bin').read_bytes())
            digest = hashlib.sha256((p/'seed-coeff.bin').read_bytes()).hexdigest()
            if reference is None:reference=data;checksum=digest
            assert digest==checksum,'Snapshot differs'
            error=max(math.hypot(data[i]-reference[i],data[i+1]-reference[i+1]) for i in range(0,len(data),2))
            assert error<1e-10,error
            for line in output.splitlines():
                if line.startswith('LCFO seed:'):continue
                fields=line.split()
                results.append(dict(mode=mode,ranks=ranks,seed_seconds=float(fields[5]),occupied=n,basis=16*n,rank=int(fields[2]),peak_rss_bytes=int(fields[3]),
                                    unitarity_error=float(fields[4]),rotation_difference=error,coefficient_sha256=digest))
print(json.dumps(results,indent=2))
