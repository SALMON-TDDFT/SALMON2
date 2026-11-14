"""Package-authored errors, retaining the previous built-in catch families.

Operating-system, dependency and user-callback exceptions are not normalized.
An exception type alone does not guarantee that a whole coupled run rolled back.
"""

__all__ = ['SalmonpyError', 'InputError', 'InputTypeError', 'NumericalError',
           'RuntimeFailure', 'SnapshotError']


class SalmonpyError(Exception):
    """Base for errors explicitly raised by SALMONpy's library APIs."""


class InputError(SalmonpyError, ValueError):
    """An input value or requested state transition is unsupported."""


class InputTypeError(InputError, TypeError):
    """An input has an unsupported Python type (also catchable as TypeError)."""


class NumericalError(SalmonpyError, FloatingPointError):
    """A numerical operation produced an invalid/nonfinite result."""


class RuntimeFailure(SalmonpyError, RuntimeError):
    """A runtime is unavailable, failed, or did not establish a valid result."""


class SnapshotError(InputError):
    """A snapshot is incompatible or invalid; restore leaves the state intact."""
