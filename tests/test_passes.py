"""Tests for pass prediction (passes.py).

Numbers quoted here come from real runs of this code on the real
TLEs in data/satellites.tle (see examples/RESULTS.md); the tests
assert the structural properties (AOS/LOS at the mask, max elevation
inside the pass, durations in the LEO range) plus the observed
counts, so a regression in the propagator or the search would fail.
"""

from __future__ import annotations

from datetime import datetime, timedelta

import pytest

from satellite_pass_predictor.ground_station import GroundStation
from satellite_pass_predictor.passes import predict_passes


def test_iss_passes_jakarta_24h_count(iss_prop, jakarta, day_start):
    passes = predict_passes(iss_prop, jakarta, day_start,
                            day_start + timedelta(hours=24),
                            min_elevation_deg=10.0)
    assert len(passes) == 3


def test_iss_best_pass_jakarta(iss_prop, jakarta, day_start):
    passes = predict_passes(iss_prop, jakarta, day_start,
                            day_start + timedelta(hours=24),
                            min_elevation_deg=10.0)
    best = max(passes, key=lambda p: p.max_elevation_deg)
    assert best.max_elevation_deg == pytest.approx(50.0, abs=0.5)
    assert best.duration_min == pytest.approx(6.6, abs=0.3)
    assert best.aos.hour == 18


def test_aos_los_elevation_equals_mask(iss_prop, jakarta, day_start):
    passes = predict_passes(iss_prop, jakarta, day_start,
                            day_start + timedelta(hours=24),
                            min_elevation_deg=10.0)
    for p in passes:
        assert jakarta.elevation_deg(iss_prop.propagate(p.aos)) == pytest.approx(10.0, abs=0.1)
        assert jakarta.elevation_deg(iss_prop.propagate(p.los)) == pytest.approx(10.0, abs=0.1)


def test_max_elevation_inside_pass(iss_prop, jakarta, day_start):
    passes = predict_passes(iss_prop, jakarta, day_start,
                            day_start + timedelta(hours=24),
                            min_elevation_deg=10.0)
    for p in passes:
        assert p.aos < p.max_elevation_time < p.los
        assert p.max_elevation_deg > 10.0
        assert p.duration_s > 0


def test_pass_durations_leo_range(iss_prop, jakarta, day_start):
    passes = predict_passes(iss_prop, jakarta, day_start,
                            day_start + timedelta(hours=24),
                            min_elevation_deg=10.0)
    for p in passes:
        assert 30.0 < p.duration_s < 15 * 60.0


def test_pass_track_elevations_above_mask(iss_prop, jakarta, day_start):
    passes = predict_passes(iss_prop, jakarta, day_start,
                            day_start + timedelta(hours=24),
                            min_elevation_deg=10.0)
    for p in passes:
        assert p.track, "track must not be empty"
        assert p.track[0][0] == p.aos
        for (_t, _az, el, _rng) in p.track:
            assert el >= 10.0 - 0.05


def test_pass_summary_string(iss_prop, jakarta, day_start):
    passes = predict_passes(iss_prop, jakarta, day_start,
                            day_start + timedelta(hours=24),
                            min_elevation_deg=10.0)
    summary = passes[0].summary()
    assert "ISS (ZARYA)" in summary
    assert "Jakarta" in summary
    assert "AOS" in summary and "LOS" in summary


def test_higher_mask_fewer_passes(iss_prop, jakarta, day_start):
    p10 = predict_passes(iss_prop, jakarta, day_start,
                         day_start + timedelta(hours=24), min_elevation_deg=10.0)
    p40 = predict_passes(iss_prop, jakarta, day_start,
                         day_start + timedelta(hours=24), min_elevation_deg=40.0)
    assert len(p40) <= len(p10)
    assert len(p40) == 1  # only the 50-deg pass survives a 40-deg mask
    assert p40[0].max_elevation_deg == pytest.approx(50.0, abs=0.5)


def test_no_passes_over_north_pole_station(iss_prop, day_start):
    pole = GroundStation("North Pole", 89.0, 0.0, 0.0)
    passes = predict_passes(iss_prop, pole, day_start,
                            day_start + timedelta(hours=24),
                            min_elevation_deg=10.0)
    assert passes == []  # ISS inclination 51.6 deg never reaches lat 89


def test_noaa_passes_exist(noaa_tle, jakarta, day_start):
    from satellite_pass_predictor.orbit import Propagator
    prop = Propagator(noaa_tle)
    passes = predict_passes(prop, jakarta, day_start,
                            day_start + timedelta(hours=24),
                            min_elevation_deg=10.0)
    assert len(passes) >= 2
    for p in passes:
        assert p.max_elevation_deg > 10.0


def test_empty_window_rejected(iss_prop, jakarta, day_start):
    with pytest.raises(ValueError):
        predict_passes(iss_prop, jakarta, day_start, day_start)


def test_short_window_finds_nothing(iss_prop, jakarta):
    # A 10-minute window far from any pass.
    start = datetime(2026, 11, 16, 0, 0, 0)
    passes = predict_passes(iss_prop, jakarta, start, start + timedelta(minutes=10))
    assert passes == []


def test_min_range_is_smallest_track_range(iss_prop, jakarta, day_start):
    passes = predict_passes(iss_prop, jakarta, day_start,
                            day_start + timedelta(hours=24),
                            min_elevation_deg=10.0)
    for p in passes:
        assert p.min_range_km == pytest.approx(
            min(entry[3] for entry in p.track), rel=1e-12)


def test_passes_sorted_by_aos(iss_prop, jakarta, day_start):
    passes = predict_passes(iss_prop, jakarta, day_start,
                            day_start + timedelta(hours=24),
                            min_elevation_deg=10.0)
    aos_list = [p.aos for p in passes]
    assert aos_list == sorted(aos_list)
