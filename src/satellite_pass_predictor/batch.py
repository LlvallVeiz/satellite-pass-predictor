"""Batch prediction: many satellites x many stations, with summaries."""

from __future__ import annotations

from dataclasses import dataclass, field
from datetime import datetime

from .ground_station import GroundStation
from .orbit import Propagator
from .passes import Pass, predict_passes
from .tle import TLE


@dataclass
class BatchResult:
    """All passes found for a satellite-station grid."""

    passes: list[Pass] = field(default_factory=list)

    @property
    def total_passes(self) -> int:
        """Total number of passes."""
        return len(self.passes)

    def passes_for(self, satellite: str, station: str) -> list[Pass]:
        """Passes for one satellite-station pair."""
        return [p for p in self.passes
                if p.satellite == satellite and p.station == station]

    def best_pass(self, satellite: str, station: str) -> Pass | None:
        """Highest-elevation pass for one pair, or None."""
        found = self.passes_for(satellite, station)
        if not found:
            return None
        return max(found, key=lambda p: p.max_elevation_deg)

    def total_visible_minutes(self, satellite: str, station: str) -> float:
        """Total time above the mask for one pair, minutes."""
        return sum(p.duration_min for p in self.passes_for(satellite, station))


def predict_batch(tles: list[TLE], stations: list[GroundStation],
                  start: datetime, end: datetime,
                  min_elevation_deg: float = 10.0,
                  scan_step_s: float = 30.0) -> BatchResult:
    """Predict passes for every satellite-station combination."""
    result = BatchResult()
    for tle in tles:
        propagator = Propagator(tle)
        for station in stations:
            result.passes.extend(predict_passes(
                propagator, station, start, end,
                min_elevation_deg=min_elevation_deg,
                scan_step_s=scan_step_s,
            ))
    result.passes.sort(key=lambda p: p.aos)
    return result
