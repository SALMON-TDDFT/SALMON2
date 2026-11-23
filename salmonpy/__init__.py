"""Composable optical/electronic calculation components, in atomic units.

Import named modules for their APIs; see salmonpy.api.PUBLIC_API. Importing the
package performs no native-library loading, compilation, MPI or network work.
"""
from ._version import __version__
from .errors import (SalmonpyError, InputError, InputTypeError, NumericalError,
                     RuntimeFailure, SnapshotError, SCFReviewRequired)
from .ld import LorentzDrude

__all__ = ['__version__', 'LorentzDrude', 'SalmonpyError', 'InputError',
           'InputTypeError', 'NumericalError', 'RuntimeFailure', 'SnapshotError', 'SCFReviewRequired']
