"""Pass prediction: find satellite passes over a ground station.

Algorithm
---------
1. Coarse scan: sample the elevation every ``scan_step_s`` over the
   requested window.  A pass is a contiguous run of samples above
   the elevation mask.
2. Refine AOS/LOS: bisection on (elevation - mask) between the last
   below-mask and first above-mask samples, to ``refine_tol_s``.
3. Maximum elevation: golden-section search on the elevation curve
   inside the refined pass (elevation is unimodal within one pass
   for the near-circular LEO orbits handled here).
4. Track: azimuth/elevation/range sampled across the pass.

A pass is reported only if its refined maximum elevation exceeds the
mask.  All times are naive UTC datetimes.
"""

from __future__ import annotations

import math
from dataclasses import dataclass, field
from datetime import datetime, timedelta

from .ground_station import GroundStation
from .orbit import Propagator


@dataclass
class Pass:
    """One satellite pass over a ground station."""

    satellite: str
    station: str
    aos: datetime  # acquisition of signal
    los: datetime  # loss of signal
    max_elevation_deg: float
    max_elevation_time: datetime
    max_elevation_azimuth_deg: float
    aos_azimuth_deg: float
    los_azimuth_deg: float
    min_range_km: float
    track: list[tuple[datetime, float, float, float]] = field(default_factory=list)
    # track entries: (time, azimuth_deg, elevation_deg, range_km)

    @property
    def duration_s(self) -> float:
        """Pass duration in seconds."""
        return (self.los - self.aos).total_seconds()

    @property
    def duration_min(self) -> float:
        """Pass duration in minutes."""
        return self.duration_s / 60.0

    def summary(self) -> str:
        """One-line human-readable summary."""
        return (f"{self.satellite} over {self.station}: "
                f"AOS {self.aos:%Y-%m-%d %H:%M:%S} az {self.aos_azimuth_deg:.1f} deg, "
                f"LOS {self.los:%Y-%m-%d %H:%M:%S} az {self.los_azimuth_deg:.1f} deg, "
                f"max el {self.max_elevation_deg:.1f} deg at "
                f"{self.max_elevation_time:%H:%M:%S}, "
                f"duration {self.duration_min:.1f} min")


def _elevation_at(propagator: Propagator, station: GroundStation,
                  time: datetime) -> float:
    return station.elevation_deg(propagator.propagate(time))


def _refine_crossing(propagator: Propagator, station: GroundStation,
                     t_low: datetime, t_high: datetime,
                     min_elevation_deg: float,
                     tol_s: float = 0.5) -> datetime:
    """Bisection on elevation(t) - mask between t_low (below) and t_high (above)."""
    f_low = _elevation_at(propagator, station, t_low) - min_elevation_deg
    f_high = _elevation_at(propagator, station, t_high) - min_elevation_deg
    if f_low * f_high > 0:
        # No sign change (should not happen from the scan); return midpoint.
        return t_low + (t_high - t_low) / 2
    while (t_high - t_low).total_seconds() > tol_s:
        mid = t_low + (t_high - t_low) / 2
        f_mid = _elevation_at(propagator, station, mid) - min_elevation_deg
        if f_mid * f_low <= 0:
            t_high, f_high = mid, f_mid
        else:
            t_low, f_low = mid, f_mid
    return t_low + (t_high - t_low) / 2


def _golden_max_elevation(propagator: Propagator, station: GroundStation,
                          t0: datetime, t1: datetime,
                          tol_s: float = 1.0) -> tuple[datetime, float]:
    """Golden-section maximisation of elevation on [t0, t1]."""
    gr = (math.sqrt(5.0) - 1.0) / 2.0
    a = t0
    b = t1
    span = (b - a).total_seconds()
    c = b - timedelta(seconds=gr * span)
    d = a + timedelta(seconds=gr * span)
    fc = _elevation_at(propagator, station, c)
    fd = _elevation_at(propagator, station, d)
    while (b - a).total_seconds() > tol_s:
        if fc > fd:
            b, d, fd = d, c, fc
            span = (b - a).total_seconds()
            c = b - timedelta(seconds=gr * span)
            fc = _elevation_at(propagator, station, c)
        else:
            a, c, fc = c, d, fd
            span = (b - a).total_seconds()
            d = a + timedelta(seconds=gr * span)
            fd = _elevation_at(propagator, station, d)
    t_best = a + (b - a) / 2
    return t_best, _elevation_at(propagator, station, t_best)


def predict_passes(propagator: Propagator, station: GroundStation,
                   start: datetime, end: datetime,
                   min_elevation_deg: float = 10.0,
                   scan_step_s: float = 30.0,
                   track_step_s: float = 30.0,
                   refine_tol_s: float = 0.5) -> list[Pass]:
    """Predict all passes of one satellite over one station in [start, end].

    ``min_elevation_deg`` is the elevation mask (default 10 deg, a
    typical amateur ground-station mask).  ``scan_step_s`` must be
    small compared with a pass duration (minutes) so no pass is
    missed; 30 s is conservative for LEO.
    """
    if end <= start:
        raise ValueError("end must be after start")
    if scan_step_s <= 0:
        raise ValueError("scan_step_s must be positive")

    # Coarse scan.
    samples: list[tuple[datetime, float]] = []
    t = start
    while t <= end:
        samples.append((t, _elevation_at(propagator, station, t)))
        t += timedelta(seconds=scan_step_s)
    if samples[-1][0] < end:
        samples.append((end, _elevation_at(propagator, station, end)))

    passes: list[Pass] = []
    i = 0
    n = len(samples)
    while i < n:
        if samples[i][1] >= min_elevation_deg:
            # Start of a pass run; find its end.
            j = i
            while j + 1 < n and samples[j + 1][1] >= min_elevation_deg:
                j += 1
            # Refine AOS: between sample i-1 (below) and i (above).
            if i > 0:
                aos = _refine_crossing(propagator, station, samples[i - 1][0],
                                       samples[i][0], min_elevation_deg, refine_tol_s)
            else:
                aos = samples[i][0]
            # Refine LOS: between sample j (above) and j+1 (below).
            if j + 1 < n:
                los = _refine_crossing(propagator, station, samples[j][0],
                                       samples[j + 1][0], min_elevation_deg, refine_tol_s)
            else:
                los = samples[j][0]
            # Maximum elevation inside the pass.
            t_max, el_max = _golden_max_elevation(propagator, station, aos, los)
            if el_max >= min_elevation_deg:
                state_max = propagator.propagate(t_max)
                look_max = station.look_angles(state_max)
                look_aos = station.look_angles(propagator.propagate(aos))
                look_los = station.look_angles(propagator.propagate(los))
                # Track samples.
                track: list[tuple[datetime, float, float, float]] = []
                tt = aos
                while tt <= los:
                    look = station.look_angles(propagator.propagate(tt))
                    track.append((tt, look.azimuth_deg, look.elevation_deg, look.range_km))
                    tt += timedelta(seconds=track_step_s)
                if track[-1][0] < los:
                    look = station.look_angles(propagator.propagate(los))
                    track.append((los, look.azimuth_deg, look.elevation_deg, look.range_km))
                min_range = min(entry[3] for entry in track)
                passes.append(Pass(
                    satellite=propagator.tle.name or f"NORAD {propagator.tle.sat_number}",
                    station=station.name,
                    aos=aos, los=los,
                    max_elevation_deg=el_max,
                    max_elevation_time=t_max,
                    max_elevation_azimuth_deg=look_max.azimuth_deg,
                    aos_azimuth_deg=look_aos.azimuth_deg,
                    los_azimuth_deg=look_los.azimuth_deg,
                    min_range_km=min_range,
                    track=track,
                ))
            i = j + 1
        else:
            i += 1
    return passes
