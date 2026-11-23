"""Prepare a WSe2 SOC input from an explicit structure; never run a solver.

python -m salmonpy.examples.prepare_wse2_soc STRUCTURE.cif POTENTIAL_DIRECTORY NEW_DIRECTORY
Potentials must already exist at W/W-sp_r.upf and Se/Se-d_r.upf below POTENTIAL_DIRECTORY.
The default grid/k sampling are an interface example, not converged optical settings.
"""
import argparse
from pathlib import Path
from salmonpy.pseudopotentials import catalog_reference
from salmonpy.dft import Pseudopotential, SCFSettings
from salmonpy.ase import dft_from_atoms
from salmonpy.preparation import prepare_dft


def prepare(structure, potentials, directory, *, grid_spacing=0.16, kmesh=(2,2,1), states=None):
    from ase.io import read
    atoms=read(structure)
    symbols=atoms.get_chemical_symbols()
    if set(symbols)!={'W','Se'} or symbols.count('Se')!=2*symbols.count('W'):
        raise ValueError('supply an explicit WSe2 structure or supercell')
    pseudos=[]
    for element,z,channel in [('W',74,3),('Se',34,2)]:
        reference=catalog_reference('pseudodojo-pbe-0.4-fr-standard-'+element.lower())
        path=Path(potentials)/element/reference.filename
        pseudos.append(Pseudopotential.from_file(path,atomic_number=z,local_channel=channel,
                                               reference=reference))
    electrons=sum(p.valence_electrons*symbols.count(element) for p,element in zip(pseudos,['W','Se']))
    model=dft_from_atoms(atoms,pseudopotentials=pseudos,states=2*electrons if states is None else states,
                         kmesh=kmesh,grid_spacing=grid_spacing,grid_spacing_unit='angstrom',
                         xc='libxc_pbe',spin_orbit=True,electronic_temperature_k=300)
    return prepare_dft(model,directory,scf=SCFSettings(max_iterations=400,tolerance=1e-8))


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('structure',type=Path)
    parser.add_argument('potentials',type=Path)
    parser.add_argument('directory',type=Path)
    parser.add_argument('--grid-spacing',type=float,default=0.16,help='angstrom; not a convergence guarantee')
    parser.add_argument('--kmesh',type=int,nargs=3,default=(2,2,1))
    parser.add_argument('--states',type=int)
    args=parser.parse_args()
    result=prepare(args.structure,args.potentials,args.directory,grid_spacing=args.grid_spacing,
                   kmesh=args.kmesh,states=args.states)
    print(result)
    print('Input preparation only. No calculation or submission was started.')


if __name__=='__main__':
    main()
