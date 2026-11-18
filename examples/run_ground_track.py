"""Example 5: ISS ground track over one orbit from the epoch.

Run from the repository root::

    python3 examples/run_ground_track.py

Prints the sub-satellite point (latitude/longitude/height) every
5 minutes for one orbital period, plus the J2 nodal regression
accumulated over that orbit.
"""

from __future__ import annotations

import math
import sys
from datetime import timedelta
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "src"))

from satellite_pass_predictor.coordinates import ecef_to_geodetic, eci_to_ecef  # noqa: E402
from satellite_pass_predictor.orbit import Propagator  # noqa: E402
from satellite_pass_predictor.tle import load_tle_file  # noqa: E402


def main() -> None:
    data = Path(__file__).resolve().parent.parent / "data" / "satellites.tle"
    iss = next(t for t in load_tle_file(str(data)) if t.sat_number == 25544)
    prop = Propagator(iss)
    t0 = iss.epoch
    period = iss.period_seconds
    print(f"ISS (ZARYA) ground track from epoch {t0:%Y-%m-%d %H:%M:%S} UTC")
    print(f"Period {iss.period_minutes:.2f} min; "
          f"J2 nodal regression {prop.rates.d_raan * 86400 * 180 / math.pi:.2f} deg/day")
    print(f"{'Time (UTC)':<12} {'Lat (deg)':>10} {'Lon (deg)':>10} {'Height (km)':>12}")
    t = t0
    while t <= t0 + timedelta(seconds=period):
        state = prop.propagate(t)
        geo = ecef_to_geodetic(eci_to_ecef(state.position_km, t))
        print(f"{t:%H:%M:%S}     {geo.latitude_deg:10.2f} {geo.longitude_deg:10.2f} "
              f"{geo.height_m / 1000:12.1f}")
        t += timedelta(minutes=5)


if __name__ == "__main__":
    main()
