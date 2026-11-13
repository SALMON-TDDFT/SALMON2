"""Opt-in native output regression. Usage: EXE CONVERGED_GS_BUNDLE NEW_DIRECTORY.
The bundle must contain rt.inp, data_for_restart and relative assets. Audit EXE first.
"""
from pathlib import Path
import hashlib,json,os,re,shutil,subprocess,sys
import numpy as np
from salmonpy.tddft import TDDFTInitialState
from salmonpy.fdtd import SalmonExecutable
from salmonpy.fields import FieldSamples
from salmonpy.crystal import respond
from salmonpy.outputs import OutputArchive
exe,bundle,dest=map(Path,sys.argv[1:]);dest.mkdir(exist_ok=False)
base=(bundle/'rt.inp').read_text();base=re.sub(r'&analysis\b.*?/','',base,flags=re.S|re.I)
base=re.sub(r'&emfield\b.*?/','',base,flags=re.S|re.I)
base=re.sub(r'&tgrid\b.*?/','',base,flags=re.S|re.I)
base+="\n&tgrid\n dt=.05,nt=24\n/\n&emfield\n ae_shape1='input',file_input1='drive.data'\n/\n"
t=np.arange(25)*.05;drive=FieldSamples(t,np.column_stack([np.zeros(25),np.zeros(25),.001*np.sin(t)]),'A')
# Native interpolation excludes its final input knot; provide a tail beyond the run.
tail=np.arange(27)*.05
file_drive=FieldSamples(tail,np.column_stack([np.zeros(27),np.zeros(27),.001*np.sin(tail)]),'A')
outputs="""&analysis
 out_rt_energy_step=8, format_voxel_data='cube'
 yn_out_dns_rt='y',out_dns_rt_step=8
 yn_out_micro_je='y',out_micro_je_step=8
 yn_out_elf_rt='y',out_elf_rt_step=8
 yn_out_estatic_rt='y',out_estatic_rt_step=8
 yn_out_current_decomposed='y',out_current_decomposed_step=8
/
"""
for name,analysis in [('plain',"&analysis\n out_rt_energy_step=8\n/\n"),('selected',outputs)]:
 d=dest/name;shutil.copytree(bundle,d);inp=base+analysis;(d/'rt.inp').write_text(inp);(d/'drive.data').unlink(missing_ok=True);file_drive.write_salmon(d/'drive.data')
 with (d/'stdout.log').open('w') as log:
  p=subprocess.run([str(exe)],input=inp,text=True,cwd=d,env={**os.environ,'OMP_NUM_THREADS':'1','OPENBLAS_NUM_THREADS':'1'},stdout=log,stderr=subprocess.STDOUT,timeout=120)
 assert p.returncode==0 and 'end SALMON' in (d/'stdout.log').read_text()
plain=np.loadtxt(dest/'plain/Si_rt.data');selected=np.loadtxt(dest/'selected/Si_rt.data')
np.testing.assert_allclose(selected[:,1:4],drive.values[1:],rtol=1e-12,atol=1e-15)
np.testing.assert_allclose(selected[:,13:16],plain[:,13:16],rtol=1e-10,atol=1e-13)
# Build a clean initial bundle: generated reference files must not count as input assets.
clean=dest/'bundle';shutil.copytree(bundle,clean);(clean/'rt.inp').write_text(base+outputs);(clean/'drive.data').unlink(missing_ok=True);file_drive.write_salmon(clean/'drive.data')
sha=hashlib.sha256(exe.read_bytes()).hexdigest()
with TDDFTInitialState.from_salmon(clean,'rt.inp',executable=SalmonExecutable(exe,sha),max_steps=24,threads=1) as initial:
 result=respond(initial,drive,output_directory=dest/'archive')
 with initial.new_state(dt=.05) as state:
  state.step(drive.values[1]);state.finish(dest/'partial')
  with np.testing.assert_raises(RuntimeError):state.step(drive.values[2])
archive=OutputArchive.open(dest/'archive');np.testing.assert_allclose(archive.array('J')[1:],selected[:,13:16],rtol=1e-10,atol=1e-13)
assert len(archive.array('time'))==25 and len(OutputArchive.open(dest/'partial').array('time'))==2
files=[]
for p in sorted((dest/'selected').glob('Si*')):
 if p.suffix not in ('.cube','.data') or p.name in ('Si_rps.dat','Si_pulse.data'):continue
 q=archive.directory/'native'/p.name
 assert q.is_file(),p.name
 assert p.read_bytes()==q.read_bytes(),p.name
 if p.suffix=='.cube':
  grid=archive.cube('native/'+p.name)
  assert grid.values.shape==(12,12,12) and np.isfinite(grid.values).all()
 files.append(p.name)
assert any('_dns_' in n or '_dns.' in n for n in files),files
assert any('_je_micro_' in n for n in files)
assert any('_elf_' in n for n in files)
assert any('_Exsta_' in n for n in files)
assert 'Si_current_decomposed.data' in files
assert not (archive.directory/'native/data_for_restart/wfn.bin').exists()
print(json.dumps({'native_files_byte_equal':files,'output_selection_current_max_difference':float(np.max(abs(selected[:,13:16]-plain[:,13:16]))),'worker_current_max_difference':float(np.max(abs(archive.array('J')[1:]-selected[:,13:16]))),'initial_assets_omitted':True,'partial_finish':True},indent=2))
