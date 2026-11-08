"""Si probe-only current and dielectric spectra for fixed and pumped screening."""
import argparse
from pathlib import Path

import numpy as np

from compare import read_case


def compare_pump(root, step, output):
    trajectories = [read_case(root/name) for name in ('fixed','fixed0','updated','updated0')]
    time = trajectories[0][:,0]
    if not 0 < step < len(time)-1:
        raise ValueError('probe step outside trajectory')
    dt = time[0]
    np.testing.assert_allclose(time, dt*np.arange(1,len(time)+1),atol=1e-8,rtol=0)
    for a in trajectories:
        np.testing.assert_array_equal(a[:,0],time)
        np.testing.assert_allclose(a[:step,1:],trajectories[0][:step,1:],atol=1e-11,rtol=1e-8)
    delays = time[step:]-time[step-1]
    energies = np.arange(1,6.001,.01)
    omega = energies/27.211386245988
    u = delays/delays[-1]
    window = 1-3*u**2+2*u**3
    kernel = np.exp(1j*omega[:,None]*delays)*window
    currents = []
    spectra = []
    amplitude = trajectories[0][step,1]-trajectories[1][step,1]
    if abs(amplitude)<1e-14: raise ValueError('missing delayed probe')
    for signal,control in zip(trajectories[::2],trajectories[1::2]):
        np.testing.assert_allclose(signal[step:,1]-control[step:,1],amplitude,rtol=0,atol=1e-10)
        delta_j = signal[step:,2]-control[step:,2]
        eps = 1+4j*np.pi*dt/(amplitude*omega)*(kernel@delta_j)
        if not np.isfinite(eps).all(): raise ValueError('nonfinite pump-probe response')
        currents.append(delta_j)
        spectra.append(eps)
    fixed, updated = spectra
    np.savetxt(output.with_suffix('.tsv'),np.column_stack([energies,fixed.real,fixed.imag,updated.real,updated.imag]),
               delimiter='\t',header='energy_eV Re_fixed Im_fixed Re_updated Im_updated')
    np.savetxt(output.with_name(output.name+'-Jz.tsv'),np.column_stack([delays*.0241888433,*currents]),
               delimiter='\t',header='delay_fs Jz_probe_fixed_au Jz_probe_updated_au')
    import matplotlib.pyplot as plt
    fig,axes = plt.subplots(3,1,sharex=True,figsize=(8,8),layout='constrained')
    for ax, baseline, pumped, label in zip(axes[:2], (fixed.real,fixed.imag),
                                           (updated.real,updated.imag), ('Re εzz','Im εzz')):
        ax.plot(energies,baseline,label='Fixed W')
        ax.plot(energies,pumped,'--',label='Updated W')
        ax.set_ylabel(label)
    axes[0].legend()
    axes[2].plot(energies,(updated-fixed).imag,color='tab:purple')
    axes[2].set(xlabel='Probe energy (eV)',ylabel='Δ Im εzz')
    for ax in axes: ax.grid(alpha=.2);ax.set_xlim(1,6)
    axes[0].set_title(f'Si pump–probe · common probe window {delays[-1]*.0241888433:.2f} fs')
    fig.savefig(output.with_suffix('.png'),dpi=160)
    plt.close(fig)


if __name__ == '__main__':
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('root',type=Path)
    p.add_argument('step',type=int,help='last time step before the probe')
    p.add_argument('output',type=Path)
    args=p.parse_args()
    compare_pump(args.root,args.step,args.output)
