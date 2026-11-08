"""Guarded finite-q extrapolation for a smooth, statically screened metal.

Fit the FULL box-projected W (including local fields), not epsilon^-1/q^2.
No fitted Drude mass, relaxation time, isotropic tail or positivity clipping.
This closure is invalid for unresolved screening or a directionally nonanalytic
limit. Two stencils and a fit residual detect many, not all, such failures;
independent k/q refinement remains mandatory.
"""
import numpy as np
from q_mesh import signed_indices


def population_head(occupations, energies, valence, occupation_tolerance=1e-7):
    """Allow an insulating head only for a common, filled, gapped manifold.

    Integer electron/hole pockets and snapshots without an empty level cannot
    establish that condition. Returning insulating=False requires the guarded
    finite-q treatment; it is not a claim that a smooth metallic limit exists.
    """
    occ, levels = np.asarray(occupations), np.asarray(energies)
    if (occ.ndim != 2 or levels.shape != occ.shape or not np.isfinite(occ).all()
            or not np.isfinite(levels).all() or not 0 < valence <= occ.shape[1]
            or np.any(np.diff(levels, axis=1) < -1e-10)):
        raise ValueError('invalid occupations or ordered projection spectrum')
    fractional = bool(np.any((occ > occupation_tolerance) & (occ < 2-occupation_tolerance)))
    filled = (np.all(abs(occ[:, :valence]-2) <= occupation_tolerance)
              and np.all(abs(occ[:, valence:]) <= occupation_tolerance))
    gap = (float(levels[:, valence:].min()-levels[:, :valence].max())
           if valence < levels.shape[1] else None)
    insulating = bool(filled and gap is not None and gap > 0)
    reason = ('common filled manifold with positive indirect gap' if insulating else
              'non-common or partially filled occupied manifold' if not filled else
              'empty projection levels unavailable' if gap is None else
              'no positive indirect gap')
    return dict(insulating=insulating, fractional_populations=fractional,
                indirect_gap_Ha=gap, head_classification_reason=reason)


def relative_channel_errors(first, second, aggregate=False):
    """Resolve uniform charge and both charge-neutral wings independently.

    Equal-volume box densities make the normalized uniform vector the charge
    mode. Roundoff floors are relative to the full matrix norm, not physical
    tolerances; a large neutral block must not mask a significant charge error.
    Aggregate=False compares each q separately (eta); True pools stencil data.
    """
    a, b = np.asarray(first), np.asarray(second)
    if a.shape != b.shape or a.ndim < 2 or a.shape[-1] != a.shape[-2]:
        raise ValueError('screening channel comparison shape mismatch')
    size = a.shape[-1]
    def channels(matrix):
        right = np.sum(matrix, axis=-1)/np.sqrt(size)
        left = np.sum(matrix, axis=-2)/np.sqrt(size)
        head = np.sum(right, axis=-1)/np.sqrt(size)
        return dict(matrix=matrix.reshape(*matrix.shape[:-2],size*size),
                    charge=head[...,None],
                    right_wing=right-head[...,None]/np.sqrt(size),
                    left_wing=left-head[...,None]/np.sqrt(size))
    ca, cb = channels(a), channels(b)
    def norm(values):
        return np.linalg.norm(values.ravel()) if aggregate else np.linalg.norm(values, axis=-1)
    floor = np.maximum(256*np.finfo(float).eps*np.maximum(norm(ca['matrix']),norm(cb['matrix'])),1e-14)
    errors = {}
    for name in ca:
        difference = norm(ca[name]-cb[name])
        # A nominally zero wing carries cancellation noise from the full matrix.
        # Do not reinterpret error below that roundoff bound as a physical
        # relative error and compare it with the user's convergence tolerance.
        relative = np.where(difference <= floor, 0., difference/np.maximum(norm(cb[name]),floor))
        errors[name] = float(np.max(relative))
    return errors


def check_eta(coarse, fine, tolerance):
    """Compare W(eta) and W(eta/2), independently of a physical scattering rate."""
    if not 0 < tolerance < 1 or not all(np.isfinite(x).all() for x in (coarse, fine)):
        raise ValueError('invalid eta comparison')
    errors = relative_channel_errors(coarse, fine)
    relative = max(errors.values())
    if relative > tolerance:
        raise ValueError(f'pump W is not eta-converged: {errors}; decrease --eta')
    return relative


def _metallic_stencils(mesh, reciprocal):
    sites = signed_indices(mesh)
    q = sites / mesh @ reciprocal
    r2 = np.sum(q*q, axis=1)
    valid = np.flatnonzero((r2 > 0) & np.all(2*abs(sites) < mesh, axis=1))
    ranked = valid[np.argsort(r2[valid])]
    if len(ranked) < 52:
        raise ValueError('too few inversion-paired q samples for metallic head')
    stencils=[]
    for count in (26,52):
        radius=r2[ranked[count-1]]
        stencils.append(valid[r2[valid] <= radius*(1+1e-10)])
    return q,r2,stencils


def metallic_sample_indices(mesh, reciprocal):
    """Only these q matrices are needed on the head-fitting rank."""
    _,_,stencils=_metallic_stencils(np.asarray(mesh),np.asarray(reciprocal))
    return np.unique(np.concatenate(stencils))


def metallic_cell(values, mesh, reciprocal, output_mesh, tolerance=.05, indices=None):
    """Return W(0), cell-average W, diagnostics using Cartesian quadratic fits.

    Input order is x-fastest, with optical q=0 at index zero (excluded).
    The integration cell is [-1/2,1/2]^3 / output_mesh in reciprocal coordinates.
    Its covariance contracts the quadratic coefficients analytically.
    """
    values = np.asarray(values)
    mesh, output_mesh = np.asarray(mesh), np.asarray(output_mesh)
    reciprocal = np.asarray(reciprocal)
    sample_indices = np.arange(int(np.prod(mesh))) if indices is None else np.asarray(indices)
    if (mesh.shape != (3,) or output_mesh.shape != (3,) or
            np.any(mesh < 6) or np.any(output_mesh < mesh) or
            reciprocal.shape != (3, 3) or not np.isfinite(reciprocal).all() or
            abs(np.linalg.det(reciprocal)) < 1e-20 or
            not 0 < tolerance < 1 or values.ndim != 3 or
            sample_indices.ndim != 1 or values.shape[0] != len(sample_indices) or
            not np.issubdtype(sample_indices.dtype,np.integer) or len(np.unique(sample_indices)) != len(sample_indices) or
            np.any(sample_indices < 0) or np.any(sample_indices >= np.prod(mesh)) or
            values.shape[1] != values.shape[2] or
            not np.isfinite(values).all()):
        raise ValueError('metallic head needs finite matrices and a >=6^3 q mesh')
    q,r2,stencils=_metallic_stencils(mesh,reciprocal)
    lookup={int(i):j for j,i in enumerate(sample_indices)}
    if any(int(i) not in lookup for stencil in stencils for i in stencil):
        raise ValueError('incomplete metallic head samples')
    scale = np.sqrt(max(r2[stencils[-1]]))
    u = q / scale
    design = np.column_stack([np.ones(len(q)), u, u*u,
                              u[:, 0]*u[:, 1], u[:, 0]*u[:, 2], u[:, 1]*u[:, 2]])
    covariance = (reciprocal/output_mesh[:, None]).T @ (reciprocal/output_mesh[:, None]) / (12*scale**2)
    moment = np.r_[1., np.zeros(3), np.diag(covariance),
                   covariance[0, 1], covariance[0, 2], covariance[1, 2]]
    fits, averages, errors, channel_errors = [], [], [], []
    flat = values.reshape(len(values), -1)
    for idx in stencils:
        rows=np.array([lookup[int(i)] for i in idx])
        c, _, rank, singular = np.linalg.lstsq(design[idx], flat[rows], rcond=None)
        if rank != 10 or singular[0]/singular[-1] > 1e6:
            raise ValueError('ill-conditioned metallic head stencil')
        comparison = relative_channel_errors((design[idx]@c).reshape(values[rows].shape), values[rows], aggregate=True)
        error = max(comparison.values())
        channel_errors.append(comparison)
        fits.append(c)
        averages.append(moment@c)
        errors.append(float(error))
    nested_channels = relative_channel_errors(averages[1].reshape(values.shape[1:]),
                                              averages[0].reshape(values.shape[1:]), aggregate=True)
    nested = max(nested_channels.values())
    if max(*errors, nested) > tolerance:
        raise ValueError(f'unresolved metallic q limit: fit={errors}, nested={nested:.3e}; refine k/q mesh')
    shape = values.shape[1:]
    point, average = fits[0][0].reshape(shape), averages[0].reshape(shape)
    # The static real-space box kernel is real at q=0 and reciprocal Hermitian.
    for matrix in (point, average):
        size = max(float(np.linalg.norm(matrix)), 1e-14)
        if (np.linalg.norm(matrix.imag) > 1e-8*size or
                np.linalg.norm(matrix-matrix.conj().T) > 1e-8*size or
                np.linalg.eigvalsh(matrix).min() < -1e-8*size):
            raise ValueError('non-real, non-Hermitian or nonpositive metallic head')
    return point.real, average.real, dict(method='finite-q full-W quadratic cell average',
        samples=[len(i) for i in stencils], radii_au=[float(np.sqrt(max(r2[i]))) for i in stencils],
        fit_relative_errors=errors, nested_relative_error=nested,
        fit_channel_relative_errors=channel_errors, nested_channel_relative_errors=nested_channels)
