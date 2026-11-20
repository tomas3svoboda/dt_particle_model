"""Fatal per-advance numerical-work exhaustion, distinct from an invalid trial."""


class StepBudgetExceeded(RuntimeError):
    """The current physical advance has exhausted its declared work allowance.

    Assemblers must propagate this exception rather than reclassify it as an
    inadmissible phase trial. The public integrator still owns exact rollback;
    no continuation, seed restart or line search receives a fresh allowance.
    """
