#!/usr/bin/env python3
from pathlib import Path
import os,re,shutil,subprocess,tempfile
root=Path(__file__).resolve().parents[2]
with tempfile.TemporaryDirectory(prefix="ow-construction-") as name:
    build=Path(name);(build/"config.h").write_text("")
    exe=build/"construction"
    subprocess.run([shutil.which("mpifort"),"-cpp","-DUSE_MPI","-I",str(build),"-J",str(build),
      "-fcheck=all","-ffpe-trap=invalid,zero,overflow","-fbacktrace",
      str(root/"src/gs/dc/dg_overlapping_wannier_types.f90"),
      str(root/"src/gs/dc/dg_overlapping_wannier_metric.f90"),
      str(root/"src/gs/dc/dg_overlapping_wannier_construction.f90"),
      str(root/"tests/dg/test_dg_overlapping_wannier_construction_mpi.f90"),
      "-llapack","-lblas","-o",str(exe)],check=True)
    env=os.environ.copy();env.setdefault("OMPI_MCA_rmaps_base_oversubscribe","1")
    signatures=[];inverse_fingerprints=[];spatial_sector_fingerprints=[];spectral_window_fingerprints=[];spectral_density_fingerprints=[];spectral_basin_fingerprints=[];spectral_operator_fingerprints=[];direct_frame_fingerprints=[];center_gauge_receipts=[]
    for n in (1,2,4,8):
      p=subprocess.run([shutil.which("mpiexec"),"-n",str(n),str(exe)],capture_output=True,text=True,env=env)
      assert p.returncode==0,(n,p.stdout,p.stderr)
      assert f"PASS overlapping-Wannier construction on {n} ranks" in p.stdout
      match=re.search(r"CONSTRUCTION ranks=\d+ fingerprint=(-?\d+) centers=([^\n]+)",p.stdout)
      assert match,p.stdout
      signatures.append((int(match.group(1)),tuple(int(x) for x in match.group(2).split())))
      inverse_match=re.search(r"INVERSE_CHARACTER_FINGERPRINT\s+(-?\d+)",p.stdout)
      assert inverse_match,p.stdout
      inverse_fingerprints.append(int(inverse_match.group(1)))
      spatial_match=re.search(r"SPATIAL_SECTOR_FINGERPRINT\s+(-?\d+)",p.stdout)
      assert spatial_match,p.stdout
      spatial_sector_fingerprints.append(int(spatial_match.group(1)))
      spectral_match=re.search(r"SPECTRAL_WINDOW_FINGERPRINT\s+(-?\d+)",p.stdout)
      assert spectral_match,p.stdout
      spectral_window_fingerprints.append(int(spectral_match.group(1)))
      density_match=re.search(r"SPECTRAL_DENSITY_FINGERPRINT\s+(-?\d+)",p.stdout)
      assert density_match,p.stdout
      spectral_density_fingerprints.append(int(density_match.group(1)))
      basin_match=re.search(r"SPECTRAL_BASIN_FINGERPRINT\s+(-?\d+)",p.stdout)
      assert basin_match,p.stdout
      spectral_basin_fingerprints.append(int(basin_match.group(1)))
      operator_match=re.search(r"SPECTRAL_OPERATOR_FINGERPRINT\s+(-?\d+)",p.stdout)
      assert operator_match,p.stdout
      spectral_operator_fingerprints.append(int(operator_match.group(1)))
      direct_match=re.search(r"DIRECT_FRAME_FINGERPRINT\s+(-?\d+)",p.stdout)
      assert direct_match,p.stdout
      direct_frame_fingerprints.append(int(direct_match.group(1)))
      gauge_match=re.search(r"POINT_CENTER_GAUGE\s+([^\n]+)",p.stdout)
      assert gauge_match,p.stdout
      center_gauge_receipts.append(gauge_match.group(1).split())
    assert len(set(signatures))==1,signatures
    assert len(set(inverse_fingerprints))==1,inverse_fingerprints
    assert len(set(spatial_sector_fingerprints))==1,spatial_sector_fingerprints
    assert len(set(spectral_window_fingerprints))==1,spectral_window_fingerprints
    assert len(set(spectral_density_fingerprints))==1,spectral_density_fingerprints
    assert len(set(spectral_basin_fingerprints))==1,spectral_basin_fingerprints
    assert len(set(spectral_operator_fingerprints))==1,spectral_operator_fingerprints
    assert len(set(direct_frame_fingerprints))==1,direct_frame_fingerprints
    assert len({tuple(x[:3]) for x in center_gauge_receipts})==1,center_gauge_receipts
    assert all(int(x[3])>0 for x in center_gauge_receipts),center_gauge_receipts
print("PASS overlapping-Wannier construction fixture on 1, 2, 4, and 8 ranks")
