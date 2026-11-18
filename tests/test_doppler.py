"""Tests for Doppler shift (doppler.py).

Physics checks
--------------
* Doppler shift is ~0 at maximum elevation (closest approach, where
  the range rate vanishes) — verified against the analytic range rate
  and an independent finite-difference range rate.
* Shift is positive (blueshift) while the satellite approaches and
  negative (redshift) while it recedes.
* Magnitude: ISS range rate up to ~7 km/s at 437 MHz gives a shift
  up to ~10 kHz (f * v / c = 437e6 * 7e3 / 3e8 ~ 10.2 kHz).
"""

from __future__ import annotations

from datetime import timedelta

import pytest

from satellite_pass_predictor.doppler import (
    doppler_shift_hz,
    ecef_velocity_km_s,
    range_rate_km_s,
    range_rate_numerical_km_s,
)
from satellite_pass_predictor.passes import predict_passes


def _best_iss_pass(iss_prop, jakarta, day_start):
    passes = predict_passes(iss_prop, jakarta, day_start,
                            day_start + timedelta(hours=24),
                            min_elevation_deg=10.0)
    assert passes, "expected ISS passes over Jakarta in 24 h"
    return max(passes, key=lambda p: p.max_elevation_deg)


def test_doppler_zero_at_max_elevation(iss_prop, jakarta, day_start):
    best = _best_iss_pass(iss_prop, jakarta, day_start)
    shift = doppler_shift_hz(
        iss_prop.propagate(best.max_elevation_time), jakarta, 437e6)
    assert abs(shift) < 500.0  # Hz; range rate ~0 at closest approach


def test_doppler_sign_approach_vs_recede(iss_prop, jakarta, day_start):
    best = _best_iss_pass(iss_prop, jakarta, day_start)
    shift_aos = doppler_shift_hz(iss_prop.propagate(best.aos), jakarta, 437e6)
    shift_los = doppler_shift_hz(iss_prop.propagate(best.los), jakarta, 437e6)
    assert shift_aos > 5000.0  # approaching: positive shift
    assert shift_los < -5000.0  # receding: negative shift


def test_doppler_magnitude_bounded(iss_prop, jakarta, day_start):
    best = _best_iss_pass(iss_prop, jakarta, day_start)
    for (t, _az, _el, _rng) in best.track:
        shift = doppler_shift_hz(iss_prop.propagate(t), jakarta, 437e6)
        assert abs(shift) < 12000.0  # < f * 8.2 km/s / c


def test_range_rate_analytic_vs_numerical(iss_prop, jakarta, day_start):
    best = _best_iss_pass(iss_prop, jakarta, day_start)
    for t in (best.aos, best.max_elevation_time, best.los):
        state = iss_prop.propagate(t)
        analytic = range_rate_km_s(state, jakarta)
        numeric = range_rate_numerical_km_s(iss_prop, t, jakarta)
        assert analytic == pytest.approx(numeric, abs=0.01)


def test_range_rate_sign(iss_prop, jakarta, day_start):
    best = _best_iss_pass(iss_prop, jakarta, day_start)
    assert range_rate_km_s(iss_prop.propagate(best.aos), jakarta) < 0.0
    assert range_rate_km_s(iss_prop.propagate(best.los), jakarta) > 0.0


def test_doppler_scales_with_frequency(iss_prop, jakarta, day_start):
    best = _best_iss_pass(iss_prop, jakarta, day_start)
    state = iss_prop.propagate(best.aos)
    s1 = doppler_shift_hz(state, jakarta, 437e6)
    s2 = doppler_shift_hz(state, jakarta, 874e6)
    assert s2 == pytest.approx(2.0 * s1, rel=1e-9)


def test_ecef_velocity_magnitude(iss_prop):
    # ECEF speed differs from ECI speed by at most Earth rotation
    # (~0.46 km/s at the equator).
    state = iss_prop.propagate(iss_prop.tle.epoch)
    v_ecef = ecef_velocity_km_s(state)
    import numpy as np
    assert abs(float(np.linalg.norm(v_ecef)) - state.speed_km_s) < 0.6


def test_doppler_noaa_437mhz(noaa_tle, jakarta, day_start):
    from satellite_pass_predictor.orbit import Propagator
    prop = Propagator(noaa_tle)
    passes = predict_passes(prop, jakarta, day_start,
                            day_start + timedelta(hours=24),
                            min_elevation_deg=10.0)
    assert passes
    best = max(passes, key=lambda p: p.max_elevation_deg)
    shift = doppler_shift_hz(prop.propagate(best.aos), jakarta, 437e6)
    assert 5000.0 < shift < 12000.0
