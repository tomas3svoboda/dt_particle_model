"""Stable checkpoint codec for the existing declared Faner numerical controls."""

from dataclasses import dataclass, fields
import math
from dtdc_simulator.core2.particle import cut_integrator as ci


@dataclass(frozen=True)
class FanerStudyControls(ci.CutSolverControls):
    """Declared wider numerical chart; original acceptance fields validated.

    The native constructor additionally caps a solver coordinate at 35.
    The owner's RULE 2 authorizes this separately named study range through
    70. The validation probe advances nothing and is never used in a solve.
    """

    def __post_init__(self):
        if (
            not math.isfinite(self.maximum_logit_magnitude)
            or not 10 <= self.maximum_logit_magnitude <= 70
        ):
            raise ValueError("declared Faner coordinate magnitude must lie in [10,70]")
        if (
            not math.isfinite(self.nonlinear_residual_tolerance)
            or not 0 < self.nonlinear_residual_tolerance <= 1e-9
        ):
            raise ValueError("declared local numerical residual tolerance cannot exceed 1e-9")
        kwargs = {f.name: getattr(self, f.name) for f in fields(ci.CutSolverControls) if f.init}
        kwargs["maximum_logit_magnitude"] = min(self.maximum_logit_magnitude, 35.0)
        kwargs["nonlinear_residual_tolerance"] = min(self.nonlinear_residual_tolerance, 1e-10)
        ci.CutSolverControls(**kwargs)
