"""Prepare a PZ native-host fixture by referencing an existing GS (no solver).

GS_INPUT RESTART_DIRECTORY PSEUDOPOTENTIAL NEW_OUTPUT
The supplied periodic, single-species PZ GS must match its original GS input
or validated RT template and pseudo.
Inline atomic coordinates are retained; no GS is calculated, copied or downloaded.
This small host-ABI fixture is not a physical convergence workflow.
"""
import argparse
import hashlib
import json
from pathlib import Path
import re
import sys
# Direct checkout CLI and installed SDK both use the same package helpers.
if __name__=='__main__' and (Path(__file__).resolve().parents[1]/'salmonpy/__init__.py').is_file():
    sys.path.insert(0,str(Path(__file__).resolve().parents[1]))
from salmonpy.bands import _input_groups, _scalar


def prepare(gs_input, restart, pseudo, output):
    gs_input=Path(gs_input).resolve(strict=True)
    restart=Path(restart).resolve(strict=True)
    pseudo=Path(pseudo).resolve(strict=True)
    output=Path(output).resolve()
    if output.exists():raise FileExistsError(output)
    for name in ('info.bin','occupation.bin','wfn.bin'):
        if not (restart/name).is_file():raise ValueError(f'existing GS is missing {name}')
    text=gs_input.read_text();groups=_input_groups(text)
    if _scalar(groups.get('functional',''),'xc')!='pz':
        raise ValueError('this host fixture needs an existing PZ GS; no XC conversion is performed')
    system=groups.get('system','')
    if _scalar(system,'yn_periodic')!='y' or _scalar(system,'nelem')!='1':
        raise ValueError('this fixture preparer needs periodic single-species input')
    if not any(name in groups for name in ('atomic_coor','atomic_red_coor')):
        raise ValueError('provide a matching GS input with inline atomic coordinates')
    body,count=re.subn(r"(?i)\bfile_pseudo\s*\(\s*1\s*\)\s*=\s*(?:'(?:''|[^'])*'|\"(?:\"\"|[^\"])*\")",
                       "file_pseudo(1)='pseudo1_rps.dat'",groups.get('pseudo',''))
    if count!=1:raise ValueError('expected one explicit file_pseudo(1) assignment')
    groups['pseudo']=body
    groups['calculation']="\n theory='tddft_pulse'\n"
    groups['control']="\n sysname='rt', directory_read_data='./data_for_restart/', yn_restart='n'\n"
    for name in ('parallel','scf'):groups.pop(name,None)
    groups['tgrid']='\n dt=.003, nt=128\n'
    groups['emfield']="\n ae_shape1='none', trans_longi='tr'\n"
    groups['propagation']="\n propagator='middlepoint', yn_predictor_corrector='n'\n"
    output.mkdir(parents=True,exist_ok=False)
    (output/'data_for_restart').symlink_to(restart,target_is_directory=True)
    (output/'pseudo1_rps.dat').symlink_to(pseudo)
    (output/'rt.inp').write_text('\n'.join('&'+name+body+'/' for name,body in groups.items())+'\n')
    (output/'fixture.json').write_text(json.dumps({
        'schema':'salmonpy.pz-host-fixture/1','gs_input':str(gs_input),
        'gs_input_sha256':hashlib.sha256(gs_input.read_bytes()).hexdigest(),
        'restart_reference':str(restart),'pseudo_reference':str(pseudo),
        'pseudo_sha256':hashlib.sha256(pseudo.read_bytes()).hexdigest(),
        'GS_copied':False,'solver_executed':False,
        'scope':'existing GS host-ABI fixture; matching/physical convergence remain caller-owned'
    },indent=2)+'\n')
    return output


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    for name in ('gs_input','restart','pseudo','output'):parser.add_argument(name,type=Path)
    args=parser.parse_args()
    print(prepare(args.gs_input,args.restart,args.pseudo,args.output))


if __name__=='__main__':main()
