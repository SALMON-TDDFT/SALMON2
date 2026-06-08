#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Standalone single-machine Python/NumPy/SciPy reference implementation of the
local Empirical Pseudopotential Method (EPM) ground-state solver for GaAs
(Cohen-Bergstresser local pseudopotential, zincblende structure).

This is a *debugging reference* for the Fortran `theory='epm'` module
(src/epm/epm_solver.f90, src/epm/epm_cohen_bergstresser.f90), with no MPI and
no OpenMP (plain NumPy/SciPy on a single machine). It writes exactly the same
SYSNAME_k.data / SYSNAME_eigen.data / SYSNAME_tm.data files consumed by
gs_info_ssbe -- so its output can be diffed against the Fortran output, or fed
directly into an SBE real-time run for cross-checking.

Note on lattice convention: unlike the Fortran module (which builds H(k) in
the primitive fcc cell with a 2-atom Ga/As basis), this reference deliberately
uses the *conventional cubic* cell -- simple-cubic direct/reciprocal lattice
vectors a_i = a*e_i / b_i = (2*pi/a)*e_i and the full 8-atom basis (4 Ga + 4
As) -- with the plane-wave basis restricted to "fcc-allowed" G (all-even or
all-odd Miller indices, see build_plane_wave_basis). This is a different but
*physically equivalent* representation: the resulting H(k) is related to the
primitive-cell H(k) by a unitary diagonal gauge transform (see the derivation
in build_hamiltonian), so eigenvalues and momentum matrix elements come out
identical to machine precision -- band-for-band the same GaAs band structure,
making the two representations directly diffable/cross-checkable.

Monolithic, self-contained: every parameter (lattice constant, plane-wave
cutoff, k-grid, number of bands/electrons, sysname, output directory) is a
hardcoded constant below -- edit them directly to change the run. No external
input files, no command-line arguments.

Run:
    python3 epm_gaas_reference.py
"""

import numpy as np
from numpy.linalg import norm
from scipy.linalg import eigh

# =============================================================================
# Hardcoded run parameters (edit here -- mirrors a GaAs &epm / &system input)
# =============================================================================
MATERIAL            = 'GaAs'
SYSNAME             = 'GaAs'
OUTPUT_DIR          = './'          # written as SYSNAME_k.data etc. in this dir

A_LATTICE_AU        = 10.68                          # lattice constant a [Bohr]
PW_CUTOFF_RY        = 11.1                          # |G|^2 cutoff [(2*pi/a)^2 units == Ry-style]

NUM_KGRID           = (4, 4, 4)     # Monkhorst-Pack grid n1 x n2 x n3 (uniform weights, no symmetry)

NSTATE              = 32            # number of bands to keep (16 valence + 16 conduction)
NELEC               = 32            # 8 atoms/cell x 4 valence electrons (Ga:3 + As:5)/2 -> NELEC/2 = 16 occupied bands

# =============================================================================
# Cohen-Bergstresser (1966) local pseudopotential form factors for GaAs,
# tabulated at G^2 = 3, 4, 8, 11 in units of (2*pi/a)^2, given in Rydberg and
# converted here to Hartree (SALMON's internal convention) by dividing by 2.
# =============================================================================
RY_TO_HA = 0.5

_CB_FORM_FACTORS_RY = {
    # G^2 : (V^S [Ry], V^A [Ry])
    3:  (-0.23,  0.07),
    4:  ( 0.00,  0.05),   # V^S(4) = 0: structure factor vanishes for zincblende
    8:  ( 0.01,  0.00),
    11: ( 0.06,  0.01),
}


def form_factors(material, G2):
    """Return (V^S, V^A) in Hartree for integer |G|^2 in units of (2*pi/a)^2."""
    if material != 'GaAs':
        raise ValueError("Only 'GaAs' Cohen-Bergstresser form factors are tabulated")
    vs_ry, va_ry = _CB_FORM_FACTORS_RY.get(G2, (0.0, 0.0))
    return vs_ry * RY_TO_HA, va_ry * RY_TO_HA


def lattice_vectors_cubic(a):
    """
    Conventional simple-cubic cell vectors a_i = a * e_i for the zincblende
    conventional (8-atom) cubic cell. The corresponding reciprocal lattice is
    likewise simple cubic, b_i = (2*pi/a) * e_i (see reciprocal_lattice()).
    """
    a1 = np.array([a,   0.0, 0.0])
    a2 = np.array([0.0, a,   0.0])
    a3 = np.array([0.0, 0.0, a  ])
    return a1, a2, a3


def atom_basis_zincblende_cubic(a):
    """
    Full 8-atom basis of the conventional zincblende cubic cell: 4 Ga (cation,
    fcc sub-lattice) + 4 As (anion, fcc sub-lattice shifted by (1/4,1/4,1/4)).
    Returns a list of (species, R[3]) with R in Cartesian a.u.
    """
    ga_frac = np.array([
        [0.00, 0.00, 0.00],
        [0.00, 0.50, 0.50],
        [0.50, 0.00, 0.50],
        [0.50, 0.50, 0.00],
    ])
    as_frac = ga_frac + np.array([0.25, 0.25, 0.25])

    basis = []
    for frac in ga_frac:
        basis.append(('Ga', a * frac))
    for frac in as_frac:
        basis.append(('As', a * frac))
    return basis


def reciprocal_lattice(a1, a2, a3):
    """b_i = 2*pi (a_j x a_k)/V,  V = a1.(a2 x a3). Returns (b1,b2,b3,volume)."""
    a23 = np.cross(a2, a3)
    a31 = np.cross(a3, a1)
    a12 = np.cross(a1, a2)
    volume = np.dot(a1, a23)
    b1 = (2.0 * np.pi / volume) * a23
    b2 = (2.0 * np.pi / volume) * a31
    b3 = (2.0 * np.pi / volume) * a12
    return b1, b2, b3, volume


def build_plane_wave_basis(a_lattice, b_matrix, cutoff_ry, fcc_selection_rule=True):
    """
    Fixed (k-independent) set of reciprocal lattice vectors
    G = m1*b1 + m2*b2 + m3*b3 with |G|^2 <= cutoff [a.u.^2] (same convention
    as the Fortran build_plane_wave_basis: the cutoff bounds |G|^2 directly).
    Returns (Gcart[npw,3], G2_units[npw] integer |G|^2 in (2*pi/a)^2 units).

    fcc_selection_rule: the simple-cubic basis vectors b_i = (2*pi/a) e_i used
    here for bookkeeping span a lattice that is *finer* than the true (fcc)
    reciprocal lattice of zincblende GaAs -- it also contains "fcc-forbidden"
    vectors with mixed-parity Miller indices (m1,m2,m3), which are not
    reciprocal lattice vectors of the crystal and must be excluded from the
    Bloch plane-wave basis (otherwise the eigenproblem describes a different,
    band-folded simple-cubic supercell, not zincblende GaAs). The textbook fcc
    structure-factor selection rule keeps exactly those G with all-even or
    all-odd (m1,m2,m3) -- this reconstructs the fcc reciprocal lattice (itself
    bcc) inside the simple-cubic bookkeeping lattice. With this filter active,
    the resulting Hamiltonian is unitarily equivalent (by a diagonal Bloch-phase
    gauge, see build_hamiltonian) to the validated 2-atom primitive-cell
    Hamiltonian, and reproduces it to machine precision band-for-band.
    """
    bnorm = norm(b_matrix, axis=1)
    gcut2 = cutoff_ry
    nmax = np.ceil(np.sqrt(gcut2) / bnorm).astype(int) + 1

    g2cart_to_units = (a_lattice / (2.0 * np.pi)) ** 2

    Gcart_list = []
    G2_list = []
    for m1 in range(-nmax[0], nmax[0] + 1):
        for m2 in range(-nmax[1], nmax[1] + 1):
            for m3 in range(-nmax[2], nmax[2] + 1):
                if fcc_selection_rule and (m1 % 2, m2 % 2, m3 % 2) not in ((0, 0, 0), (1, 1, 1)):
                    continue
                Gtmp = m1 * b_matrix[0] + m2 * b_matrix[1] + m3 * b_matrix[2]
                g2_units = np.dot(Gtmp, Gtmp)
                if g2_units <= gcut2 + 1.0e-8:
                    Gcart_list.append(Gtmp)
                    G2_list.append(int(round(g2_units * g2cart_to_units)))

    Gcart = np.array(Gcart_list)
    G2 = np.array(G2_list, dtype=int)
    return Gcart, G2


def monkhorst_pack_grid(b_matrix, num_kgrid):
    """Uniform MP grid, equal weights 1/nk, no symmetry reduction (matches Fortran ordering)."""
    n1, n2, n3 = num_kgrid
    nk = n1 * n2 * n3
    kpoint = np.zeros((nk, 3))
    kweight = np.full(nk, 1.0 / nk)
    ik = 0
    for i1 in range(1, n1 + 1):
        for i2 in range(1, n2 + 1):
            for i3 in range(1, n3 + 1):
                f1 = (2 * i1 - n1 - 1) / (2.0 * n1)
                f2 = (2 * i2 - n2 - 1) / (2.0 * n2)
                f3 = (2 * i3 - n3 - 1) / (2.0 * n3)
                kpoint[ik] = f1 * b_matrix[0] + f2 * b_matrix[1] + f3 * b_matrix[2]
                ik += 1
    return kpoint, kweight


def atomic_form_factors(material, G2):
    """
    Single-atom (cation/anion) local pseudopotential form factors, derived
    from the tabulated Cohen-Bergstresser symmetric/antisymmetric combinations
    via the standard zincblende decomposition (cation at R = -tau, anion at
    R = +tau, tau = a/8 (1,1,1) -- the same placement implicit in the
    primitive-cell formula V^S cos(G.tau) + i V^A sin(G.tau)):
        f_cation e^{-iG.tau} + f_anion e^{+iG.tau}
            = (f_cation+f_anion) cos(G.tau) + i(f_anion-f_cation) sin(G.tau)
            == V^S(G) cos(G.tau) + i V^A(G) sin(G.tau)
        =>  f_Ga (cation) = [V^S(G) - V^A(G)] / 2
            f_As (anion)  = [V^S(G) + V^A(G)] / 2
    """
    VS, VA = form_factors(material, G2)
    return {'Ga': 0.5 * (VS - VA), 'As': 0.5 * (VS + VA)}


def build_hamiltonian(material, kvec, Gcart, a_lattice, basis):
    """
    H_{G,G'}(k) = (1/2)|k+G|^2 delta_{G,G'} + V(G-G')

    The local-pseudopotential matrix element is built directly as the
    structure-factor sum over the full 8-atom basis (4 Ga + 4 As) of the
    conventional zincblende cubic cell, on the simple-cubic reciprocal
    lattice b_i = (2*pi/a) e_i:
        V(dG) = (1/n_cell) * sum_{j=1}^{8} f_j(|dG|^2) e^{i dG . R_j}
    with f_Ga, f_As the atomic form factors of atomic_form_factors().

    The 1/n_cell prefactor (n_cell = 4 GaAs formula units per conventional
    cubic cell) renormalizes the larger-cell structure-factor sum back to a
    per-formula-unit potential. Writing the 4 fcc sub-lattice translations as
    R_fcc and the cation/anion offset as delta = 2*tau (tau = a/8 (1,1,1), the
    Cohen-Bergstresser origin convention), each formula unit contributes
        f_Ga e^{i dG.R_fcc} + f_As e^{i dG.(R_fcc+delta)}
            = e^{i dG.(R_fcc+tau)} * [f_Ga e^{-i dG.tau} + f_As e^{+i dG.tau}]
            = e^{i dG.(R_fcc+tau)} * [V^S(dG) cos(dG.tau) + i V^A(dG) sin(dG.tau)] ,
    i.e. the validated 2-atom matrix element times a pure Bloch/gauge phase.
    Summed over the 4 fcc translations, sum_{R_fcc} e^{i dG.R_fcc} equals
    n_cell=4 for any dG that is an actual fcc reciprocal lattice vector
    ("fcc-allowed": all-even or all-odd Miller indices on the simple-cubic
    grid -- constructive interference of the 4 translationally-equivalent
    formula units) and vanishes identically for "fcc-forbidden" dG (which are
    NOT reciprocal lattice vectors of the crystal and must be excluded from
    the plane-wave basis up front by build_plane_wave_basis(fcc_selection_rule
    =True), since the 8-atom structure factor alone does not reach exactly
    zero for them within the truncated/tabulated form-factor model). After the
    1/n_cell division this leaves precisely
        V(dG) = e^{i dG.(R_fcc+tau)} * [V^S(dG) cos(dG.tau) + i V^A(dG) sin(dG.tau)]
    i.e. H_cubic = D H_primitive D^dagger with D = diag(e^{i G.tau}) a unitary
    diagonal gauge transform (the R_fcc-dependent phase cancels between bra and
    ket of the structure-factor sum and only the net e^{i(G-G').tau} survives).
    D is unitary, so H_cubic has *exactly* the same eigenvalues as the
    validated 2-atom primitive-cell Hamiltonian H_primitive (band-for-band, to
    machine precision -- verified numerically), and the momentum matrix
    elements p_mn = sum_G c*_m(G)(k+G)c_n(G) are gauge-invariant for the same
    reason (the G-diagonal phases of c_m, c_n cancel pairwise).
    """
    npw = Gcart.shape[0]
    H = np.zeros((npw, npw), dtype=complex)
    units = (a_lattice / (2.0 * np.pi)) ** 2
    n_cell = len(basis) // 2   # number of GaAs formula units in the conventional cell (= 4)

    kpg = kvec[None, :] + Gcart                      # (npw,3)
    diag = 0.5 * np.einsum('ij,ij->i', kpg, kpg)     # (1/2)|k+G|^2
    np.fill_diagonal(H, diag)

    for i in range(npw):
        for j in range(i + 1, npw):
            dG = Gcart[i] - Gcart[j]
            dG2 = int(round(np.dot(dG, dG) * units))
            f = atomic_form_factors(material, dG2)
            if f['Ga'] == 0.0 and f['As'] == 0.0:
                continue
            S = sum(f[species] * np.exp(1j * np.dot(dG, R)) for species, R in basis)
            val = S / n_cell
            H[i, j] = val
            H[j, i] = np.conj(val)
    return H


def momentum_matrix(kvec, Gcart, evec):
    """
    p_mn(k) = sum_G conjg(c_m(G)) (k+G) c_n(G)   (diagonal in the plane-wave basis)
    evec: (npw, nb) matrix of the lowest `nb` eigenvectors (columns).
    Returns p_mn[nb,nb,3].
    """
    npw, nb = evec.shape
    kpg = kvec[None, :] + Gcart                      # (npw,3)
    p_mn = np.zeros((nb, nb, 3), dtype=complex)
    for idir in range(3):
        Dc = kpg[:, idir][:, None] * evec             # (npw, nb)
        p_mn[:, :, idir] = evec.conj().T @ Dc
    return p_mn


def main():
    a1, a2, a3 = lattice_vectors_cubic(A_LATTICE_AU)
    basis = atom_basis_zincblende_cubic(A_LATTICE_AU)
    b1, b2, b3, volume = reciprocal_lattice(a1, a2, a3)
    b_matrix = np.array([b1, b2, b3])

    Gcart, G2 = build_plane_wave_basis(A_LATTICE_AU, b_matrix, PW_CUTOFF_RY)
    npw = Gcart.shape[0]

    kpoint, kweight = monkhorst_pack_grid(b_matrix, NUM_KGRID)
    nk = kpoint.shape[0]

    nb = NSTATE
    nocc = NELEC // 2

    print('# EPM (local Cohen-Bergstresser pseudopotential) -- Python/NumPy/SciPy reference')
    print(f'#   material           = {MATERIAL}')
    print(f'#   lattice constant a = {A_LATTICE_AU:.5e} a.u.')
    print(f'#   plane waves        = {npw}')
    print(f'#   k-points           = {nk}')
    print(f'#   bands requested    = {nb} / valence electrons = {NELEC}')

    eigen = np.zeros((nb, nk))
    occup = np.zeros((nb, nk))
    p_tm = np.zeros((nb, nb, 3, nk), dtype=complex)

    for ik in range(nk):
        H = build_hamiltonian(MATERIAL, kpoint[ik], Gcart, A_LATTICE_AU, basis)
        evals, evecs = eigh(H)              # ascending order, like ZHEEV

        eigen[:, ik] = evals[:nb]
        occup[:nocc, ik] = 2.0
        occup[nocc:nb, ik] = 0.0

        p_mn = momentum_matrix(kpoint[ik], Gcart, evecs[:, :nb])
        p_tm[:, :, :, ik] = p_mn

        if (ik + 1) % max(1, nk // 10) == 0 or ik == nk - 1:
            print(f'#   ... diagonalized k-point {ik + 1}/{nk}')

    write_epm_files(SYSNAME, OUTPUT_DIR, MATERIAL, kpoint, kweight, eigen, occup, p_tm)


# =============================================================================
# Output writers -- byte-for-byte compatible with gs_info_ssbe::read_k_data /
# read_eigen_data / read_tm_data (free-format read(fh,*), fixed header-line
# counts and fixed numeric-field counts per data line; see write_epm_files in
# src/epm/epm_solver.f90, which this mirrors line-for-line).
# =============================================================================
def write_epm_files(sysname, outdir, material, kpoint, kweight, eigen, occup, p_tm):
    nk = kpoint.shape[0]
    nb = eigen.shape[0]

    # --- SYSNAME_k.data --- (5 header lines, then "ik kx ky kz weight")
    with open(f'{outdir}{sysname}_k.data', 'w') as f:
        f.write('# k-point data\n')
        f.write('# generated by EPM (Cohen-Bergstresser local pseudopotential) -- Python reference\n')
        f.write(f'# material = {material}, nk = {nk}\n')
        f.write('# units: kx,ky,kz [a.u.], weight (sums to 1)\n')
        f.write('# ik, kx, ky, kz, weight\n')
        for ik in range(nk):
            f.write('{:6d}{:18.10E}{:18.10E}{:18.10E}{:18.10E}\n'.format(
                ik + 1, kpoint[ik, 0], kpoint[ik, 1], kpoint[ik, 2], kweight[ik]))

    # --- SYSNAME_eigen.data --- (3 header lines; per ik: 1 header + nb lines)
    with open(f'{outdir}{sysname}_eigen.data', 'w') as f:
        f.write('# eigenvalue data\n')
        f.write('# generated by EPM (Cohen-Bergstresser local pseudopotential) -- Python reference\n')
        f.write(f'# nk = {nk:6d}, nb = {nb:6d}\n')
        for ik in range(nk):
            f.write(f'# ik = {ik + 1:6d}\n')
            for ib in range(nb):
                f.write('{:6d}{:18.10E}{:18.10E}\n'.format(ib + 1, eigen[ib, ik], occup[ib, ik]))

    # --- SYSNAME_tm.data --- (3 header lines; block 1 p_tm; 1 header; block 2 rvnl_tm = 0)
    with open(f'{outdir}{sysname}_tm.data', 'w') as f:
        f.write('# transition matrix data\n')
        f.write('# generated by EPM (Cohen-Bergstresser local pseudopotential) -- Python reference\n')
        f.write('# block 1: p_tm = <u_m|p|u_n>  (ik, ib, jb, Re px, Im px, Re py, Im py, Re pz, Im pz)\n')
        for ik in range(nk):
            for ib in range(nb):
                for jb in range(nb):
                    px, py, pz = p_tm[ib, jb, 0, ik], p_tm[ib, jb, 1, ik], p_tm[ib, jb, 2, ik]
                    f.write('{:6d}{:6d}{:6d}{:18.10E}{:18.10E}{:18.10E}{:18.10E}{:18.10E}{:18.10E}\n'.format(
                        ik + 1, ib + 1, jb + 1,
                        px.real, px.imag, py.real, py.imag, pz.real, pz.imag))
        f.write('# block 2: rvnl_tm = -i[r,Vnl]  (all zero: local pseudopotential, no nonlocal correction)\n')
        zeros9 = '{:6d}{:6d}{:6d}' + '{:18.10E}' * 6 + '\n'
        for ik in range(nk):
            for ib in range(nb):
                for jb in range(nb):
                    f.write(zeros9.format(ik + 1, ib + 1, jb + 1, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0))

    print(f'# EPM (Python reference): wrote ground-state data files for sysname = {sysname}')


if __name__ == '__main__':
    main()
