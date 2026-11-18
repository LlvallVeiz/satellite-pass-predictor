"""Example 3: Doppler shift of a 437 MHz downlink during the best ISS pass.

Run from the repository root::

    python3 examples/run_doppler.py

Prints the Doppler curve (one sample per minute) and the classic
S-shape: positive shift while approaching, ~0 at maximum elevation,
negative while receding.
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

FREQUENCY_HZ = 437e6


def main() -> None:
    data = Path(__file__).resolve().parent.parent / "data" / "satellites.tle"
    iss = next(t for t in load_tle_file(str(data)) if t.sat_number == 25544)
    prop = Propagator(iss)
    start = datetime(2026, 11, 16, 0, 0, 0)
    passes = predict_passes(prop, JAKARTA, start, start + timedelta(hours=24),
                            min_elevation_deg=10.0)
    best = max(passes, key=lambda p: p.max_elevation_deg)
    print(best.summary())
    print(f"Downlink {FREQUENCY_HZ / 1e6:.0f} MHz — Doppler shift:")
    print(f"{'Time (UTC)':<12} {'El (deg)':>9} {'Range (km)':>11} {'Shift (kHz)':>12}")
    t = best.aos
    while t <= best.los:
        state = prop.propagate(t)
        look = JAKARTA.look_angles(state)
        shift = doppler_shift_hz(state, JAKARTA, FREQUENCY_HZ)
        print(f"{t:%H:%M:%S}     {look.elevation_deg:9.2f} {look.range_km:11.0f} "
              f"{shift / 1000:12.2f}")
        t += timedelta(seconds=60)
    print()
    print(f"Max shift: +{doppler_shift_hz(prop.propagate(best.aos), JAKARTA, FREQUENCY_HZ) / 1000:.2f} kHz at AOS")
    print(f"Min shift: {doppler_shift_hz(prop.propagate(best.los), JAKARTA, FREQUENCY_HZ) / 1000:.2f} kHz at LOS")
    print(f"Shift at max elevation: "
          f"{doppler_shift_hz(prop.propagate(best.max_elevation_time), JAKARTA, FREQUENCY_HZ):.0f} Hz")


if __name__ == "__main__":
    main()
