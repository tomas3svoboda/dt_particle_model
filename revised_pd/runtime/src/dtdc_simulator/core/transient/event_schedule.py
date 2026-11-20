"""Shared exact event alignment for offline transient integrations."""

from __future__ import annotations

from collections.abc import Iterator
from typing import Protocol, TypeVar


class TimedEvent(Protocol):
    time_s: float


EventT = TypeVar("EventT", bound=TimedEvent)


def event_aligned_steps(
    *,
    end_time_s: float,
    maximum_step_s: float,
    events: tuple[EventT, ...],
    event_name: str,
    horizon_name: str = "the horizon",
    positive_error: str = "end time and maximum step must be positive",
) -> Iterator[tuple[float, float, EventT]]:
    """Yield ``(step_end, dt, active_event)`` without straddling an event."""

    if end_time_s <= 0.0 or maximum_step_s <= 0.0:
        raise ValueError(positive_error)
    if not events or events[0].time_s != 0.0:
        raise ValueError(f"the first {event_name} event must start at t=0")
    if any(
        later.time_s <= earlier.time_s
        for earlier, later in zip(events, events[1:])
    ):
        raise ValueError(f"{event_name} events must be strictly ordered")
    if events[-1].time_s > end_time_s:
        raise ValueError(f"{event_name} event lies beyond {horizon_name}")

    time_s = 0.0
    event_index = 0
    epsilon = 1.0e-12 * max(1.0, end_time_s)
    while time_s < end_time_s - epsilon:
        next_event = (
            events[event_index + 1].time_s
            if event_index + 1 < len(events)
            else end_time_s
        )
        step_end = min(time_s + maximum_step_s, next_event, end_time_s)
        dt_s = step_end - time_s
        if dt_s <= epsilon:
            event_index += 1
            continue
        yield step_end, dt_s, events[event_index]
        time_s = step_end
        if (
            event_index + 1 < len(events)
            and abs(time_s - next_event) <= epsilon
        ):
            event_index += 1
