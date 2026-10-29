"""Declared absolute-clock step caps for whole-history particle qualification.

The DT film-resolved family follows the measured 2026-10-23 onset/depletion
windows. Its landmarks are fixed physical clocks, shared across meshes; they
are not shifted to each trajectory's events. Native advances and all rejection,
proof, assembly and conservation accounting remain with the existing driver.

The last band follows the film (ruling packet 7, 2026-10-27): after the 570 s
landmark it requests the film-window cap (base 0.625 s) while the external water
film exists on the particle and the declared base 5 s after the depletion event.
The caller states the film's presence on every request (``film_present``); the
depletion event itself is still located natively. Every request before 570 s is
the same either way.
"""

from dataclasses import dataclass
import math


@dataclass(frozen=True)
class HistoryTimeControl:
    """One common refinement multiplier for feed, front, endpoint and dry caps."""

    study: str
    level: int
    dry_schedule: str = "uniform"

    def __post_init__(self):
        if self.study not in ("DT", "Faner"):
            raise ValueError("unknown particle study")
        if isinstance(self.level, bool) or self.level not in (0, 1, 2):
            raise ValueError("time level must be an integer in 0, 1, 2")
        if not isinstance(self.level, int):
            raise ValueError("time level must be an integer in 0, 1, 2")
        if self.dry_schedule not in ("uniform", "dt-film-resolved"):
            raise ValueError("unknown dry time schedule")
        if self.dry_schedule == "dt-film-resolved" and self.study != "DT":
            raise ValueError("dt-film-resolved is declared only for DT")

    @property
    def multiplier(self):
        return 2**self.level

    def dry_bands(self, *, film_present):
        """Pairs of exclusive upper physical clock and maximum duration, seconds.

        ``film_present`` (required, keyword-only): whether the external water film
        exists on the particle at the request's start clock. For dt-film-resolved the
        last band follows the film: the film-window cap (base 0.625 s) while the film
        exists, the declared base 5 s after the depletion event. The two landmarks and
        the caps before 570 s do not depend on it; the uniform schedule ignores it.
        """
        if self.dry_schedule == "dt-film-resolved":
            # WAS (None, 5.0); 2026-10-27, ruling packet 7, record GT_PS2_PAPER1_N80_RUNG_COMPLETED_ON_FILM_AWARE_SCHEDULE_2026-10-27.md
            base = ((515.0, 0.3125), (570.0, 0.625), (None, 0.625 if film_present else 5.0))
        else:
            base = ((None, 5.0 if self.study == "DT" else 1.0),)
        return tuple((end, cap / self.multiplier) for end, cap in base)

    def dry_request(self, clock_s, horizon_s, *, film_present):
        """Cap the next request and land at the next fixed landmark or horizon.

        ``film_present`` (required, keyword-only) selects the last band (dry_bands).
        The caller may shorten a request by the existing event/halving machinery.
        It must call again from the actual accepted clock. No state is snapped to
        a landmark, and this function neither retries nor spends solver work.
        """
        if not all(math.isfinite(t) and t >= 0 for t in (clock_s, horizon_s)):
            raise ValueError("time controls require finite nonnegative clocks")
        if clock_s >= horizon_s:
            raise ValueError("next request requires clock before horizon")
        for end, cap in self.dry_bands(film_present=film_present):
            if end is None or clock_s < end:
                target = horizon_s if end is None else min(horizon_s, end)
                requested = min(cap, target - clock_s)
                if not requested > 0 or clock_s + requested <= clock_s:
                    raise ValueError("time request does not advance the physical clock")
                if clock_s + requested > target:
                    raise ValueError("time request overshoots its declared landmark")
                return requested
        raise AssertionError("dry bands must include an unbounded final band")

    def describe(self):
        factor = self.multiplier
        return {
            "dry_schedule": self.dry_schedule,
            "time_level": self.level,
            "common_refinement_multiplier": factor,
            "feed_dt_s": (0.1 if self.study == "DT" else 0.25) / factor,
            "DT_steam_feed_dt_s": 0.02 / factor if self.study == "DT" else None,
            "birth_dt_s": 0.01 / factor,
            "front_maximum_dt_s": (0.1 if self.study == "DT" else 0.4) / factor,
            "core_remaining_time_ceiling_s": 0.002 / factor,
            # The bands without the external film (the former listing), then those with it.
            "dry_bands": [
                {"exclusive_upper_clock_s": end, "maximum_dt_s": cap}
                for end, cap in self.dry_bands(film_present=False)
            ],
            "dry_bands_while_external_film_present": [
                {"exclusive_upper_clock_s": end, "maximum_dt_s": cap}
                for end, cap in self.dry_bands(film_present=True)
            ],
            "landmarks": "fixed absolute physical clocks; no event alignment or state reset",
            "last_band_follows_external_film": self.dry_schedule == "dt-film-resolved",
            "depletion_advance_allowance_factor": (
                "2 on an advance with film whose predicted depletion duration is in (0, 2 requested],"
                " else 1 (ruling packet 7, 2026-10-27)"
            ),
            "dry_physical_advance_limit": 1500,
            "changed_native_tolerances_or_work_allowances": False,
        }


def depletion_advance_allowance_factor(
    film_loading_kg_m2, predicted_depletion_duration_s, requested_s
):
    """Multiplier of one dry advance's shared residual-assembly allowance: 2 on the advance
    predicted to deplete the external film, else 1.

    2 when the film exists (loading > 0) and the runner's predicted depletion duration is known
    and lies in (0, 2 requested_s] (the runner's own event-first seeding window); 1 otherwise.
    The per-solve Newton budget is not touched. Measured need (2026-10-27, follow-up 4, mirror):
    one advance of the 80-cell history used 4,121 assemblies against the declared 2,910 (the
    depletion event solve 2,422, then the halved step); the 40-cell event advance used 856 of
    1,470. Provenance: ruling packet 7, record
    docs/GT_PS2_PAPER1_N80_RUNG_COMPLETED_ON_FILM_AWARE_SCHEDULE_2026-10-27.md.
    """
    if (
        film_loading_kg_m2 > 0
        and predicted_depletion_duration_s is not None
        and 0 < predicted_depletion_duration_s <= 2 * requested_s
    ):
        return 2
    return 1
