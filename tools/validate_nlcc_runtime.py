"""Bind completed native consistency tests and whole-app regression to a binary.

Explicit validation only: no build, GS collection, or automatic certification on
import. This checks implementation regression, not physical material convergence.
"""
import argparse,copy,hashlib,json,os
from pathlib import Path
import numpy as np
from salmonpy.dft import DFT,Pseudopotential,SCFSettings,_read_status,_read_forces
from salmonpy._runtime_capabilities import validate_nlcc_receipt,reference_sha256,NLCC_XC
from measure_nlcc_force_work import run_scalar_only,write_json
from validate_nlcc_ci_results import validate_native_records


def whole_application(executable,potential,output):
    from salmonpy import _runtime_capabilities
    ref=json.loads(Path(_runtime_capabilities.__file__).with_name('_nlcc_reference.json').read_text())
    executable=Path(executable);raw=Path(potential).read_bytes()
    if hashlib.sha256(raw).hexdigest()!=ref['pseudo_sha256']:raise ValueError('NLCC regression UPF checksum mismatch')
    before=hashlib.sha256(executable.read_bytes()).hexdigest()
    pseudo=Pseudopotential(14,2,'.upf',raw)
    positions=np.array([[0,0,0],[0,.5,.5],[.5,0,.5],[.5,.5,0],[.25,.25,.25],[.25,.75,.75],[.75,.25,.75],[.75,.75,.25]],float)
    positions[0,0]+=ref['initial_atom_x_displacement_bohr']/ref['lattice_bohr']
    model=DFT((14,)*8,np.eye(3)*ref['lattice_bohr'],positions,(pseudo,),ref['electrons'],ref['states'],tuple(ref['grid']),tuple(ref['kmesh']),ref['xc'])
    text=model.input_text(SCFSettings(max_iterations=600,tolerance=ref['scf_tolerance'])).replace("write_gs_restart_data='all'","write_gs_restart_data='no'")
    output=Path(output)
    if output.exists():raise FileExistsError(output)
    env={k:v for k,v in os.environ.items() if k!='SALMONPY_TDDFT_WORKER'};env.update(OMP_NUM_THREADS='1',OPENBLAS_NUM_THREADS='1')
    status=run_scalar_only(executable,text,raw,output,env)
    if status or 'end SALMON' not in (output/'stdout.log').read_text():raise RuntimeError('NLCC whole-application regression failed')
    energy=_read_status(output/'gs_status.dat')['energy_hartree'];force=_read_forces(output/'gs_forces.dat',8)
    after=hashlib.sha256(executable.read_bytes()).hexdigest()
    if after!=before:raise ValueError('runtime changed during validation')
    e=abs(energy-ref['energy_hartree']);f=float(np.max(np.abs(force-np.asarray(ref['forces_hartree_per_bohr']))))
    if not np.isfinite(e) or not np.isfinite(f) or e>ref['energy_atol'] or f>ref['force_atol']:
        raise AssertionError('NLCC whole-application output differs from fixed reference')
    return dict(status='passed',executable_sha256=before,reference_sha256=reference_sha256(),
                energy_error_hartree=e,maximum_force_error_hartree_per_bohr=f)


def certify_relocated(receipt,executable,potential,output):
    """Recheck relocated executable; preserve original object-test build binding."""
    receipt=copy.deepcopy(receipt)
    validate_nlcc_receipt(receipt,executable_sha256=receipt['build_executable_sha256'],source_patch_sha256=receipt['source_patch_sha256'],source_tree_sha256=receipt['source_tree_sha256'])
    receipt['pre_relocation_application_validation']=receipt['application_validation']
    receipt['application_validation']=whole_application(executable,potential,output)
    receipt['executable_sha256']=receipt['application_validation']['executable_sha256']
    validate_nlcc_receipt(receipt,executable_sha256=receipt['executable_sha256'],source_patch_sha256=receipt['source_patch_sha256'],source_tree_sha256=receipt['source_tree_sha256'])
    return receipt


def main():
    p=argparse.ArgumentParser(description=__doc__)
    for name in ('build_record','xc_record','core_record','potential','output'):p.add_argument(name,type=Path)
    p.add_argument('--mpi-record',type=Path,help='Completed2/4-rank force evidence for this MPI-linked build')
    p.add_argument('--application-record',type=Path,help='Explicit completed whole_application scalar comparison, bound to the same executable/reference')
    a=p.parse_args();spec=json.loads((Path(__file__).resolve().parents[1]/'runtime/source.json').read_text())
    native=validate_native_records(a.build_record,a.xc_record,a.core_record,source_spec=spec,
        mpi_build_record=a.build_record if a.mpi_record else None,mpi_record=a.mpi_record)
    build=json.loads(a.build_record.read_text());a.output.mkdir(parents=True,exist_ok=False)
    app=json.loads(a.application_record.read_text()) if a.application_record else whole_application(build['executable'],a.potential,a.output/'application')
    receipt=dict(schema='salmonpy.nlcc-forces/1',status='passed',executable_sha256=build['sha256'],
        build_executable_sha256=build['sha256'],source_patch_sha256=spec['patch_sha256'],
        source_tree_sha256=spec['source_tree_sha256'],
        scope={'backend':'cpu','spin':'unpolarized','potential':'NC-UPF','xc':list(NLCC_XC)},
        native_validation=native,application_validation=app)
    validate_nlcc_receipt(receipt,executable_sha256=build['sha256'],source_patch_sha256=spec['patch_sha256'],source_tree_sha256=spec['source_tree_sha256'])
    write_json(a.output/'validation.json',receipt);print(json.dumps(receipt,indent=2))

if __name__=='__main__':main()
