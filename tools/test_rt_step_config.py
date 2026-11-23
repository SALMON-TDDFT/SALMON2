"""Build/run native RT config + actual field/Taylor kernel ownership probes.

python tools/test_rt_step_config.py SOURCE REFERENCE_SOURCE NEW_OUTPUT --fc gfortran
No solver build, physical GS or existing source modification. hpsi is a linear
synthetic operator only for the Taylor loop; actual field/Taylor source is sliced.
"""
import argparse,hashlib,json,shutil,subprocess,sys
from pathlib import Path


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('source',type=Path);p.add_argument('reference',type=Path);p.add_argument('output',type=Path)
    p.add_argument('--fc',default='/opt/homebrew/bin/gfortran');a=p.parse_args()
    a.source=a.source.resolve();a.reference=a.reference.resolve();a.output=a.output.resolve()
    root=Path(__file__).resolve().parents[1];checker=Path(__file__).with_name('check_forbidden_eigensolvers.py')
    subprocess.run([sys.executable,str(checker),'--source-root',str(a.source/'src')],check=True)
    subprocess.run([sys.executable,str(checker),'--source-root',str(a.reference/'src')],check=True)
    a.output.mkdir(parents=True,exist_ok=False);src=a.output/'src';src.mkdir();(src/'config.h').write_text('')
    shutil.copy2(a.source/'src/common/structures.f90',src/'structures.f90')
    for prefix,tree in [('current',a.source),('reference',a.reference)]:
        text=(tree/'src/rt/em_field.f90').read_text();start=text.index('subroutine calc_emfields(');end=text.index('end subroutine calc_emfields',start)+len('end subroutine calc_emfields')
        (src/(prefix+'_field.f90')).write_text('module '+prefix+'_field\ncontains\n'+text[start:end]+'\nend module\n')
        text=(tree/'src/rt/taylor.f90').read_text().replace('module taylor_sub','module '+prefix+'_taylor')
        (src/(prefix+'_taylor.f90')).write_text(text)
    (src/'stubs.f90').write_text('''module salmon_global
real(8) :: dt,film_thickness,epsilon_em(2)
character(2) :: trans_longi
integer :: n_hamil
end module
module inputoutput
use salmon_global, only: n_hamil
end module
module math_constants
real(8),parameter :: pi=acos(-1d0)
end module
module phys_constants
real(8),parameter :: cspeed_au=137.035999084d0
end module
module nvtx_wrapper
contains
subroutine nvtxStartRange(name,id)
character(*) :: name
integer,optional :: id
end subroutine
subroutine nvtxEndRange
end subroutine
end module
module sendrecv_grid
use structures, only: s_sendrecv_grid
end module
module hamiltonian
use structures
contains
subroutine hpsi(input,output,info,mg,V_local,system,stencil,srg,ppg)
type(s_orbital) :: input,output
type(s_parallel_info) :: info
type(s_rgrid) :: mg
type(s_scalar) :: V_local(*)
type(s_dft_system) :: system
type(s_stencil) :: stencil
type(s_sendrecv_grid) :: srg
type(s_pp_grid) :: ppg
output%zwf=2d0*input%zwf
end subroutine
end module
''')
    shutil.copy2(root/'tests/native/test_rt_step_config.f90',src/'probe.f90')
    commands=[]
    with (a.output/'build.log').open('w') as log:
        def run(cmd):
            commands.append([str(x) for x in cmd]);subprocess.run(cmd,cwd=a.output,stdout=log,stderr=subprocess.STDOUT,check=True)
        run([sys.executable,checker,'--source-root',src])
        names=['structures','stubs','current_field','reference_field','current_taylor','reference_taylor','probe']
        for name in names:run([a.fc,'-cpp','-O2','-fopenmp','-ffree-line-length-none','-I'+str(src),'-c',src/(name+'.f90'),'-o',a.output/(name+'.o')])
        run([sys.executable,checker,'--object-tree',a.output])
        run([a.fc,'-fopenmp',*[a.output/(n+'.o') for n in names],'-o',a.output/'probe'])
        run([a.output/'probe'])
    record={'status':'passed','scope':'native configuration/field and synthetic-Hamiltonian Taylor ownership; not full native contexts',
        'sources':{str(f.relative_to(src)):hashlib.sha256(f.read_bytes()).hexdigest() for f in src.iterdir() if f.is_file()},'commands':commands}
    (a.output/'result.json').write_text(json.dumps(record,indent=2)+'\n');print('Native step-config probe passed')


if __name__=='__main__':main()
