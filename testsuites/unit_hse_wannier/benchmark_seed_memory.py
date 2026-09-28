"""Sequential root seed RSS at nbasis=16*noccupied; not full RT memory."""
from pathlib import Path
import array, json, math, os, subprocess, tempfile
here = Path(__file__).resolve().parent
root = here.parents[1]
results = []
with tempfile.TemporaryDirectory() as folder:
    p = Path(folder)
    subprocess.run(['cc', '-c', str(root/'testsuites/unit_lcfo_rt/peak_rss.c'), '-o', str(p/'rss.o')], check=True)
    subprocess.run([os.environ.get('FC', 'gfortran'), '-O2', '-fexternal-blas', '-fno-tree-loop-vectorize',
                    str(root/'src/xc/exx_wannier_gauge.f90'), str(here/'seed_memory_probe.f90'),
                    str(p/'rss.o'), '-L/opt/homebrew/opt/openblas/lib', '-lopenblas', '-o', str(p/'probe')], cwd=p, check=True)
    for n in (256, 512):
        reference = None
        for mode in ('reference', 'gamma2d'):
            output = subprocess.check_output([str(p/'probe'), mode, str(n)], cwd=p, universal_newlines=True,
                                              env=dict(os.environ, OMP_NUM_THREADS='1', OPENBLAS_NUM_THREADS='1'))
            fields = output.split()
            data = array.array('d');data.frombytes((p/'seed-u.bin').read_bytes())
            if reference is None:
                reference = data
            error = max(math.hypot(data[i]-reference[i], data[i+1]-reference[i+1]) for i in range(0,len(data),2))
            assert error < 1e-10, error
            results.append(dict(mode=mode, occupied=n, basis=16*n, peak_rss_bytes=int(fields[2]),
                                unitarity_error=float(fields[3]), rotation_difference=error))
print(json.dumps(results, indent=2))
