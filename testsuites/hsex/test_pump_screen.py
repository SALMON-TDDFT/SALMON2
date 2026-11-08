"""Check the finite-q pump response against the complete fractional band sum."""
from pathlib import Path
import sys
import os
from functools import partial
from tempfile import TemporaryDirectory

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[2] / 'tools/hsex'))
from pump_screen import solve_pump_target
from sternheimer import polarizability


def run(comm=None):
    solver=(partial(solve_pump_target,q_kernel=os.environ['HSEX_Q_KERNEL'],q_threads=int(os.environ.get("HSEX_Q_THREADS","2")))
            if os.environ.get('HSEX_Q_KERNEL') else solve_pump_target)
    rng = np.random.default_rng(1301)
    n = np.array([4, 4, 4])
    mesh = np.array([2, 1, 1])
    order = np.array([[0, 0, 0], [1, 0, 0]])
    R = -order
    cell = np.eye(3)*4
    frac = order/mesh
    k = 2*np.pi*frac/4
    shells = [np.array([[0, 0, 0]]), np.array([[0, 0, 0], [1, 0, 0]])]
    eta = 0.017
    grid = np.indices((4, 4, 4)).reshape(3, -1)[::-1].T/4
    occ = np.zeros((2, 64))
    occ[:, :3] = [[1.8, .15, .05], [1.65, .25, .10]]
    energies = np.array([np.arange(64)*.1+i*.023 for i in range(2)])
    states = []
    with TemporaryDirectory() as tmp:
        for j in range(2):
            u = np.linalg.qr(rng.normal(size=(64,64))+1j*rng.normal(size=(64,64)))[0]
            states.append(u)
            with (Path(tmp)/f'hsex-pump-{j+1:06d}.bin').open('wb') as f:
                f.write(b'HSP1')
                np.array([*n, 2, j+1, 3, 80], dtype='=i4').tofile(f)
                for x in [cell,k[j],np.zeros(3),energies[j,:3],occ[j,:3]]:
                    np.asarray(x,dtype='=f8').tofile(f)
                u[:,:3].T.tofile(f)
                ((u*energies[j])@u.conj().T).T.tofile(f)
        jobs = [(tmp,n,k,cell,frac,order,mesh,R,shells,j,eta) for j in range(2)]
        result = polarizability(jobs,R,mesh,shells,1,solver=solver,comm=comm)
        paired_jobs = [(*job[:-1], (eta, eta/2)) for job in jobs]
        paired = polarizability(paired_jobs,R,mesh,shells,1,solver=solver,channels=2,comm=comm)
        owned = polarizability(paired_jobs,R,mesh,shells,1,solver=solver,
                               channels=2,comm=comm,distributed=True)
        pieces = [owned] if comm is None else comm.gather(owned,root=0)
        if comm is not None and comm.rank != 0:
            return
        gathered={}
        for part in pieces:
            assert not (set(gathered) & set(part)), 'q owners must be unique'
            gathered.update(part)
        assert set(gathered)==set(range(len(R)))
        for iq,value in gathered.items():
            np.testing.assert_allclose(value,paired[iq],rtol=1e-12,atol=1e-12)
        for single, dual in zip(result,paired):
            np.testing.assert_allclose(single,dual[0],rtol=1e-12,atol=1e-12)
        for channel, broadening in enumerate((eta,eta/2)):
            expected = np.zeros((2,2),complex)
            intraband = np.zeros_like(expected)
            for j in range(2):
                i = 1-j
                G = shells[1]+np.rint(frac[i]+R[1]/mesh-frac[j]).astype(int)
                phase = np.exp(2j*np.pi*grid@G.T)
                m = np.einsum('rm,rg,rn->mng',states[j].conj(),phase,states[i])
                gap = energies[j,:,None]-energies[i,None,:]
                weights = .5*(occ[i,None,:]-occ[j,:,None])*gap/(gap**2+broadening**2)
                expected += np.einsum('mng,mnh,mn->gh',m.conj(),m,weights)
                intraband += np.einsum('mng,mnh,mn->gh',m.conj(),m,weights*np.eye(64))
            assert np.linalg.norm(intraband) > 1e-3, 'fixture must resolve same-band response'
            np.testing.assert_allclose(paired[1][channel],expected,rtol=1e-11,atol=1e-12)
        print('fractional intraband + interband full-sum / paired eta PASS')



if __name__ == '__main__':
    if "--mpi" in sys.argv:
        from mpi4py import MPI
        run(MPI.COMM_WORLD)
    else:
        run()
