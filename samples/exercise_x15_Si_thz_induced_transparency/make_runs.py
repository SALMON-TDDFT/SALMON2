#!/usr/bin/env python3
"""Generate the x15 run tree: one SELF-CONTAINED directory per run.

    python3 make_runs.py [--fields 100,300,1000,3000] [--nk 9] [--nstate 32]
                         [--dt-fs 0.05] [--t-end-fs 600] [--dt-scan] [--outdir runs]

WHY SUBDIRECTORIES.  The EPM ground state writes `Si_prim_k.data`,
`Si_prim_eigen.data` and `Si_prim_tm.data` into the CURRENT directory, and the SBE
step reads them from there and writes `Si_prim_sbe_*.data` back into it.  Two runs
started in the same directory therefore overwrite each other's ground state AND
their outputs.  Every run here gets its own directory holding both inputs, so runs
are independent and can be submitted in parallel.

WHY THERE IS NO FIELD FILE.  The measured DAST single-cycle transient is reproduced
analytically by `ae_shape1 = 'Acos2'` with the carrier-envelope phase set so that
A(t) is a single hump, hence E = -dA/dt is one clean cycle:

    A(t) = -(f0/w) cos^2(pi*tt/tw) sin(w*tt + 2*pi*phi_CEP),   |tt| < tw/2

    phi_cep1 = 0.75  ->  sin(...) = -cos(w*tt)
    omega1   = hbar*pi/tw  ->  the cosine runs 0 -> 1 -> 0 across the whole support

Against `x12/DAST_singlecycle_100kV.txt` this gives corr(A) = +0.9994 (1.8 % rms),
corr(E) = +0.9952 (6.3 % rms), spectral centroid 3.81 vs 3.83 THz, and peak |A|
low by 2.6 %.  Peak |E| is exact by construction: the envelope makes
peak|E| = (2/sqrt(3)) f0, so E_amplitude1 = E0 * sqrt(3)/2.

PHYSICS SETTINGS follow the cluster production configuration, with two corrections
established in this fork (see README sections 3 and 4):
  * `yn_sbe_vg_sumrule = 'y'` is MANDATORY.  Without it the velocity-gauge
    truncation leaves an uncancelled diamagnetic current ~ eta*N_e*A/V which at THz
    dominates the absorbed energy; it cannot be removed afterwards because eta grows
    with field (0.23 % -> 7.35 % along this pulse at nstate = 20).
  * the e-ph detailed-balance split now acts for every material, not only the 2D
    Dirac registry.
"""
import argparse
import os

HBAR_EVFS = 0.6582119569
PEAK_OVER_F0 = 2.0 / 3.0 ** 0.5          # = 1.154700, exact for the cos^2 envelope

GS = """!### Si primitive-cell EPM ground state (step 1 of 2) -- run in THIS directory.
!### Writes Si_prim_k/_eigen/_tm.data, which Si_prim_sbe_rt.inp reads from here.
&calculation
    theory = 'epm'
/
&control
    sysname = 'Si_prim'
/
&units
    unit_system = 'A_eV_fs'
/
&system
    yn_periodic = 'y'
    al(1:3) = 5.42935818d0, 5.42935818d0, 5.42935818d0
    nelec  = 8
    nstate = {nstate}
/
&kgrid
    num_kgrid(1:3) = {nk}, {nk}, {nk}
/
&epm
    epm_material            = 'Si'
    epm_lattice_constant_au = 10.26d0
    epm_pw_cutoff_ry        = 27.0d0      ! raised: the Si Delta-valley camel-back needs it
    epm_cell                = 'primitive'
/
"""

RT = """!### {title}
!### step 2 of 2 -- run AFTER Si_prim_epm_gs.inp, in THIS directory.
!###
!### drive: {drive}
!### peak |E| = {ekv:g} kV/cm; A0 = {a0:.4e} fs*V/Ang = {a0au:.4e} a.u.
&calculation
  theory = 'sbe'
/
&control
  sysname = 'Si_prim'
/
&units
  unit_system = 'A_eV_fs'
/
&system
  yn_periodic  = 'y'
  al_vec1(1:3) = 0.00000000d0, 2.71467909d0, 2.71467909d0
  al_vec2(1:3) = 2.71467909d0, 0.00000000d0, 2.71467909d0
  al_vec3(1:3) = 2.71467909d0, 2.71467909d0, 0.00000000d0
  nelec  = 8
  nstate = {nstate}
/
&kgrid
  num_kgrid(1:3) = {nk}, {nk}, {nk}
/
&tgrid
  dt = {dt}d0
  nt = {nt}
/
&emfield
{emfield}/
&analysis
  out_projection_step      = {proj}
  out_rt_energy_step       = 10
  yn_out_intraband_current = 'y'
/
&epm
  epm_material            = 'Si'
  epm_lattice_constant_au = 10.26d0
/
&sbe
  yn_sbe_full_dressed      = 'y'
  frozen_core_threshold_ev = -10.0d0
  frozen_free_threshold_ev =  10.0d0
  yn_sbe_dressed_ref       = 'y'
  yn_sbe_colmem            = '{ring}'
  yn_sbe_colmem_pop        = '{ring}'
  yn_sbe_superres          = 'y'
  yn_sbe_eph               = '{ring}'
  yn_sbe_impact_ionization = '{ring}'
  yn_sbe_auger             = 'n'
  yn_sbe_eeh               = '{ring}'
  yn_sbe_coulomb           = 'n'
  yn_sbe_bgr_threshold     = 'n'
  sbe_eph_temperature_k    = 300.0d0
  yn_sbe_ii_holes          = 'n'
  yn_sbe_eph_acoustic      = '{ring}'
  sbe_ii_phassist          = 1.0d0
  yn_sbe_vg_sumrule        = 'y'    ! MANDATORY at THz -- see the module docstring
  sbe_checkpoint_step      = 200    ! so a wall-clock kill costs <= 200 steps
/
"""

EM_PULSE = """  ! analytic reproduction of the DAST single-cycle transient (no field file).
  ! phi_cep = 0.75 makes A(t) a single hump, so E = -dA/dt is one clean cycle.
  ae_shape1      = 'Acos2'
  epdir_re1(1:3) = 1.0d0, 0.0d0, 0.0d0
  E_amplitude1   = {amp}   ! = {ekv:g} kV/cm * sqrt(3)/2
  tw1            = {tw:.1f}d0
  omega1         = {om:.8f}d0    ! = hbar*pi/tw1
  phi_cep1       = 0.75d0
"""
EM_DARK = """  ! zero-field control: everything below must stay at zero.
  ae_shape1 = 'none'
"""


def fort(x):
    """Fortran double literal: 8.660258e-03 -> 8.660258d-03 (never both e and d)."""
    return ('%.6e' % x).replace('e', 'd')


def write_run(outdir, tag, title, drive, emfield, ekv, a0, nk, nstate, dt, nt, ring, proj):
    d = os.path.join(outdir, tag)
    os.makedirs(d, exist_ok=True)
    with open(os.path.join(d, 'Si_prim_epm_gs.inp'), 'w') as fh:
        fh.write(GS.format(nk=nk, nstate=nstate))
    body = RT.format(title=title, drive=drive, emfield=emfield, ekv=ekv, a0=a0,
                     a0au=a0 * 0.0241888 / 0.0529177 / 10.0,
                     nk=nk, nstate=nstate, dt=dt, nt=nt, ring=ring, proj=proj)
    with open(os.path.join(d, 'Si_prim_sbe_rt.inp'), 'w') as fh:
        fh.write(body)
    # Resume twin. A wall-clock kill is normal on a cluster, and without this the
    # restart begins from t = 0 (the solver opens the outputs with status='replace')
    # -- the checkpoint is then written but never read, which is the worst of both.
    with open(os.path.join(d, 'Si_prim_sbe_rt_resume.inp'), 'w') as fh:
        fh.write(body.replace("  sbe_checkpoint_step      = 200",
                              "  sbe_checkpoint_step      = 200\n"
                              "  yn_sbe_checkpoint_restart = 'y'   ! continue from the checkpoint"))
    return d


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('--fields', default='100,300,1000,3000', help='peak |E| [kV/cm]')
    ap.add_argument('--nk', type=int, default=9)
    ap.add_argument('--nstate', type=int, default=32)
    ap.add_argument('--dt-fs', type=float, default=0.05)
    ap.add_argument('--t-end-fs', type=float, default=600.0)
    ap.add_argument('--tw-fs', type=float, default=273.0)
    ap.add_argument('--dt-scan', action='store_true',
                    help='also emit the dt-convergence set (mandatory first step, README 4)')
    ap.add_argument('--no-ring', action='store_true', help='dissipators off (coherent)')
    ap.add_argument('--outdir', default='runs')
    a = ap.parse_args()

    om = HBAR_EVFS * 3.141592653589793 / a.tw_fs
    ring = 'n' if a.no_ring else 'y'
    made = []

    def emit(tag, title, ekv, dt, ring_):
        nt = int(round(a.t_end_fs / dt))
        proj = max(1, int(round(20.0 / dt)))
        if ekv > 0:
            amp = (ekv / 1e5) / PEAK_OVER_F0
            a0 = amp * PEAK_OVER_F0 / (om / HBAR_EVFS)
            em = EM_PULSE.format(amp=fort(amp), ekv=ekv, tw=a.tw_fs, om=om)
            drive = 'Acos2 single-cycle, tw = %.0f fs, centroid 3.8 THz' % a.tw_fs
        else:
            a0 = 0.0
            em = EM_DARK
            drive = 'none (dark control)'
        d = write_run(a.outdir, tag, title, drive, em, ekv, a0,
                      a.nk, a.nstate, dt, nt, ring_, proj)
        made.append('%-22s  E0=%-7g kV/cm  dt=%-6g fs  nt=%-6d ring=%s' % (tag, ekv, dt, nt, ring_))
        return d

    for s in a.fields.split(','):
        ekv = float(s)
        emit('E%gkVcm' % ekv, 'Si THz induced transparency, %g kV/cm' % ekv, ekv, a.dt_fs, ring)
    emit('dark', 'Si zero-field control (ring on, no drive)', 0.0, a.dt_fs, ring)
    if a.dt_scan:
        for dt in (0.1, 0.05, 0.025):
            emit('dtscan_%g_E1000' % dt, 'dt convergence at 1000 kV/cm', 1000.0, dt, ring)

    with open(os.path.join(a.outdir, 'MANIFEST.txt'), 'w') as fh:
        fh.write('\n'.join(made) + '\n')
    print('\n'.join(made))
    print('\nwrote %d self-contained run directories under %s/' % (len(made), a.outdir))
    print('each holds Si_prim_epm_gs.inp + Si_prim_sbe_rt.inp and nothing else;')
    print('run the GS first, then the RT, in that directory.')


if __name__ == '__main__':
    main()
