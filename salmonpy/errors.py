"""Package-authored errors, retaining the previous built-in catch families.

Operating-system, dependency and user-callback exceptions are not normalized.
An exception type alone does not guarantee that a whole coupled run rolled back.
"""

__all__ = ['SalmonpyError', 'InputError', 'InputTypeError', 'NumericalError',
           'RuntimeFailure', 'SnapshotError', 'SCFReviewRequired']


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


class SCFReviewRequired(SalmonpyError, RuntimeError):
    """SCF ended unconverged; inspect retained data before continuing.

    result is the saved GroundStateResult or live NativeGSState (not closed).
    report_directory contains diagnostics, or None if report writing failed.
    """

    def __init__(self, result, report_directory):
        self.result = result
        self.report_directory = report_directory
        location = str(report_directory) if report_directory is not None else str(result.directory)
        super().__init__('SCF not converged: automatic continuation paused. '
                         'Inspect residual history and threshold in ' + location +
                         '. Retained results/state remain available as exception.result; '
                         'do not bypass review automatically.')
