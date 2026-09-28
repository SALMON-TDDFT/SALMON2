"""Sequential process peak RSS for Gamma initial evaluation (not full RT peak)."""
from pathlib import Path
import json, os, subprocess, tempfile
here = Path(__file__).resolve().parent
root = here.parents[1]
results = []
with tempfile.TemporaryDirectory() as folder:
    p = Path(folder)
    subprocess.run(['cc', '-c', str(root/'testsuites/unit_lcfo_rt/peak_rss.c'), '-o', str(p/'rss.o')], check=True)
    subprocess.run([os.environ.get('FC', 'gfortran'), '-O2', '-fexternal-blas', '-fno-tree-loop-vectorize',
                    str(root/'src/xc/exx_wannier_gauge.f90'), str(here/'gamma_memory_probe.f90'),
                    str(p/'rss.o'), '-L/opt/homebrew/opt/openblas/lib', '-lopenblas', '-o', str(p/'probe')], cwd=p, check=True)
    for n in (512, 1024):
        for mode in ('reference', 'inplace'):
            output = subprocess.check_output([str(p/'probe'), mode, str(n)], cwd=p, text=True,
                                              env=dict(os.environ, OMP_NUM_THREADS='1', OPENBLAS_NUM_THREADS='1'))
            fields = output.split()
            results.append(dict(mode=mode, occupied=n, peak_rss_bytes=int(fields[2]), spread=float(fields[3]), gradient=float(fields[4])))
print(json.dumps(results, indent=2))
