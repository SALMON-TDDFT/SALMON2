"""Fail-closed capability checks for trusted, checksum-verified runtime records.

These attest implementation regression tests, not material convergence or
cryptographic publisher identity. A feature name alone grants no capability.
"""
import hashlib
import math
import re
from pathlib import Path
from .errors import InputError

NLCC_PATCH_SHA256 = 'a25f4cf5e71700546836e7310ca40f4458e8320dfab40a550f90595259d182a8'
NLCC_SOURCE_TREE_SHA256 = 'd8419ef56f116480129ba0d49e4369825ed1ce8ecb8f8dbd895f0f2083c3ef78'
NLCC_XC = ('pz', 'libxc_pz', 'libxc_pbe')


def reference_sha256():
    return hashlib.sha256(Path(__file__).with_name('_nlcc_reference.json').read_bytes()).hexdigest()


def _bounded(value, limit):
    return type(value) in (int, float) and math.isfinite(value) and 0 <= value <= limit


def validate_nlcc_receipt(receipt, *, executable_sha256, source_patch_sha256, source_tree_sha256=None):
    """Validate evidence and its binary/source/scope binding; return XC scope."""
    try:
        if (receipt['schema'] != 'salmonpy.nlcc-forces/1' or receipt['status'] != 'passed'
                or receipt['executable_sha256'] != executable_sha256
                or source_patch_sha256 != NLCC_PATCH_SHA256
                or source_tree_sha256 != NLCC_SOURCE_TREE_SHA256
                or receipt['source_patch_sha256'] != source_patch_sha256
                or receipt['source_tree_sha256'] != source_tree_sha256):
            raise InputError('NLCC validation binary/source mismatch')
        if not re.fullmatch('[0-9a-f]{64}', receipt['build_executable_sha256']):
            raise InputError('invalid NLCC build binding')
        scope=receipt['scope']
        if scope != {'backend':'cpu','spin':'unpolarized','potential':'NC-UPF','xc':list(NLCC_XC)}:
            raise InputError('unsupported NLCC validation scope')
        native=receipt['native_validation']
        if (native['schema'] != 'salmonpy.nlcc-native-validation/1' or native['status'] != 'passed'
                or native['build_executable_sha256'] != receipt['build_executable_sha256']
                or native['source_patch_sha256'] != source_patch_sha256
                or native['source_tree_sha256'] != source_tree_sha256):
            raise InputError('NLCC native evidence mismatch')
        for name, test, limit in (('xc','test_nlcc_xc_variation',5e-9),('core','test_nlcc_core_force',1e-7)):
            evidence=native['tests'][name]
            if (evidence['test'] != test or type(evidence['cases']) is not int or evidence['cases'] != 2
                    or evidence['rank_thread_pairs'] != [[1,1],[1,2]]
                    or not re.fullmatch('[0-9a-f]{64}', evidence['record_sha256'])
                    or not re.fullmatch('[0-9a-f]{64}', evidence['test_executable_sha256'])
                    or not _bounded(evidence['maximum_reported_error'],limit)):
                raise InputError('NLCC native validation incomplete or outside tolerance')
        if 'mpi' in native['tests']:
            evidence=native['tests']['mpi']
            if (evidence['test']!='test_nlcc_mpi_force' or evidence['cases']!=4
                    or evidence['rank_thread_pairs']!=[[2,1],[2,2],[4,1],[4,2]]
                    or evidence['build_executable_sha256']!=receipt['build_executable_sha256']
                    or not re.fullmatch('[0-9a-f]{64}',evidence['record_sha256'])
                    or not re.fullmatch('[0-9a-f]{64}',evidence['test_executable_sha256'])
                    or not _bounded(evidence['maximum_reported_error'],1e-13)):
                raise InputError('NLCC MPI validation incomplete or mismatched')
        app=receipt['application_validation']
        if (app['status'] != 'passed' or app['executable_sha256'] != executable_sha256
                or app['reference_sha256'] != reference_sha256()
                or not _bounded(app['energy_error_hartree'],1e-8)
                or not _bounded(app['maximum_force_error_hartree_per_bohr'],1e-8)):
            raise InputError('NLCC whole-application validation mismatch')
    except (KeyError,TypeError,ValueError,OverflowError) as exc:
        if isinstance(exc,InputError): raise
        raise InputError('invalid NLCC validation record') from exc
    return NLCC_XC
