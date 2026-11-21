"""Prepare DFT or supplied TDDFT inputs only; never execute or submit.

python -m salmonpy.examples.prepare_electrons dft Si_rps.dat NEW_DIRECTORY
python -m salmonpy.examples.prepare_electrons tddft rt.inp NEW_DIRECTORY --asset pseudo1_rps.dat Si_rps.dat
The Si model is a small staging example, not a converged physical reference.
"""
import argparse
import json
from pathlib import Path
import numpy as np
from salmonpy.dft import DFT,Pseudopotential
from salmonpy.preparation import prepare_dft,prepare_tddft


def silicon_model(pseudo_file):
    pseudo=Pseudopotential.from_file(pseudo_file,atomic_number=14,local_channel=2)
    positions=((0,0,0),(0,.5,.5),(.5,0,.5),(.5,.5,0),
               (.25,.25,.25),(.25,.75,.75),(.75,.25,.75),(.75,.75,.25))
    electrons=8*pseudo.valence_electrons
    return DFT([14]*8,np.eye(3)*10.26,positions,[pseudo],electrons=electrons,
               states=electrons//2+4,grid=(16,16,16),kmesh=(1,1,1))


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    modes=parser.add_subparsers(dest='kind',required=True)
    dft=modes.add_parser('dft',help='write a small Si DFT input only')
    dft.add_argument('pseudo_file',type=Path);dft.add_argument('directory',type=Path)
    dft.add_argument('--no-save-restart',action='store_true',help='disable restart output in the future DFT run')
    rt=modes.add_parser('tddft',help='copy supplied RT input and explicitly listed assets only')
    rt.add_argument('input_file',type=Path);rt.add_argument('directory',type=Path)
    rt.add_argument('--asset',nargs=2,action='append',default=[],metavar=('NAME','LOCAL_FILE'))
    args=parser.parse_args()
    if args.kind=='dft':
        directory=prepare_dft(silicon_model(args.pseudo_file),args.directory,save_restart=not args.no_save_restart)
    else:
        files=dict(args.asset)
        if len(files)!=len(args.asset):parser.error('duplicate asset name')
        directory=prepare_tddft(args.input_file,args.directory,files=files)
    print(json.dumps({'status':'prepared-not-run','kind':args.kind,'directory':str(directory),
                      'scope':'input files only; no solver, worker, GS discovery or job submission'},indent=2))


if __name__=='__main__':main()
