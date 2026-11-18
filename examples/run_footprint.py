"""Example 4: footprint / coverage geometry for the four satellites.

Run from the repository root::

    python3 examples/run_footprint.py

Prints the horizon footprint and the footprint at 10- and 30-degree
elevation masks, plus the ground-track speed, for each satellite.
"""

from __future__ import annotations

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "src"))

from satellite_pass_predictor.footprint import (  # noqa: E402
    compute_footprint,
    ground_track_speed_km_s,
)
from satellite_pass_predictor.tle import load_tle_file  # noqa: E402


def main() -> None:
    data = Path(__file__).resolve().parent.parent / "data" / "satellites.tle"
    tles = load_tle_file(str(data))
    for tle in tles:
        print(f"{tle.name} (NORAD {tle.sat_number}): "
              f"altitude {tle.altitude_km:.1f} km, "
              f"ground-track speed {ground_track_speed_km_s(tle.altitude_km):.2f} km/s")
        for mask in (0.0, 10.0, 30.0):
            fp = compute_footprint(tle.altitude_km, mask)
            print(f"  mask {mask:4.0f} deg: central angle {fp.central_angle_deg:5.2f} deg, "
                  f"radius {fp.footprint_radius_km:7.0f} km, "
                  f"area {fp.footprint_area_km2:12,.0f} km^2")
        print()


if __name__ == "__main__":
    main()
