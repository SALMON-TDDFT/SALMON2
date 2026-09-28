"""Same H2 geometry as GS measurements; full-cell DC export then mesh RT."""
import sys,math
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'benchmark_h2_pbeh'))
from run import input_text,SHAPES,LAYOUT

def dc_block(shape,ranks):
    return f"""
&dc
 num_fragment=1,1,1
 num_rgrid_buffer=0,0,0
 nproc_rgrid_tot={','.join(map(str,LAYOUT[ranks]))}
 nstate_frag={math.prod(shape)}
 yn_dc_lcfo='y'
 energy_cut=100d0
 lambda_cut=1d-7
/
"""

def gs_input(shape,ranks):
    return input_text(shape,ranks).replace("yn_dc='n'","yn_dc='y'").replace('&system',"&system\n temperature_k=0d0")+dc_block(shape,ranks)

def rt_input(shape,ranks):
    s=input_text(shape,ranks).replace("theory='dft'","theory='tddft_response'").replace("yn_dc='n'","yn_dc='n'\n yn_conventional_from_dcdft='y'")
    return s+dc_block(shape,ranks)+"""
&tgrid
 dt=0.02d0
 nt=16
/
&emfield
 ae_shape1='impulse'
 e_impulse=0.0001d0
 epdir_re1=1d0,0d0,0d0
/
&analysis
 out_rt_energy_step=1
 nenergy=10
 de=0.01d0
/
"""
