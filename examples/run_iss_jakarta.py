"""Example 1: ISS passes over Jakarta for 24 hours.

Run from the repository root::

    python3 examples/run_iss_jakarta.py

Prints every pass above a 10-degree elevation mask, the best pass,
and the Doppler shift of a 437 MHz downlink along that pass.
"""

from __future__ import annotations

import sys
from datetime import datetime, timedelta
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "src"))

from satellite_pass_predictor.doppler import doppler_shift_hz  # noqa: E402
from satellite_pass_predictor.ground_station import JAKARTA  # noqa: E402
from satellite_pass_predictor.orbit import Propagator  # noqa: E402
from satellite_pass_predictor.passes import predict_passes  # noqa: E402
from satellite_pass_predictor.tle import load_tle_file  # noqa: E402

FREQUENCY_HZ = 437e6  # typical amateur-satellite UHF downlink


def main() -> None:
    data = Path(__file__).resolve().parent.parent / "data" / "satellites.tle"
    iss = next(t for t in load_tle_file(str(data)) if t.sat_number == 25544)
    prop = Propagator(iss)
    start = datetime(2026, 11, 16, 0, 0, 0)
    end = start + timedelta(hours=24)
    passes = predict_passes(prop, JAKARTA, start, end, min_elevation_deg=10.0)

    print(f"ISS (ZARYA) over Jakarta, {start:%Y-%m-%d} UTC, mask 10 deg")
    print(f"Orbital period {iss.period_minutes:.2f} min, "
          f"mean altitude {iss.altitude_km:.1f} km")
    print(f"Passes in 24 h: {len(passes)}")
    print()
    for p in passes:
        print(p.summary())
    best = max(passes, key=lambda p: p.max_elevation_deg)
    print()
    print(f"Best pass: max elevation {best.max_elevation_deg:.1f} deg, "
          f"duration {best.duration_min:.1f} min")
    print(f"Doppler shift at {FREQUENCY_HZ / 1e6:.0f} MHz along the best pass:")
    for (t, _az, el, _rng) in best.track:
        shift = doppler_shift_hz(prop.propagate(t), JAKARTA, FREQUENCY_HZ)
        print(f"  {t:%H:%M:%S}  el {el:6.2f} deg  {shift / 1000:+8.2f} kHz")


if __name__ == "__main__":
    main()
