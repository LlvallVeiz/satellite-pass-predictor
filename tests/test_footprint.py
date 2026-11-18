"""Tests for footprint / coverage geometry (footprint.py).

Reference values (spherical Earth, Re = 6378.137 km)
-----------------------------------------------------
* ISS altitude 418.4 km, mask 0 deg:
  central angle = acos(Re/(Re+h)) = 20.20 deg,
  radius = Re * lambda = 2248 km,
  nadir angle = asin(Re/(Re+h)) = 69.79 deg.
* Footprint area = 2 pi Re^2 (1 - cos lambda).
"""

from __future__ import annotations

import math

import pytest

from satellite_pass_predictor.constants import EARTH_RADIUS_KM
from satellite_pass_predictor.footprint import (
    compute_footprint,
    ground_track_speed_km_s,
    max_slant_range_km,
)


def test_iss_footprint_horizon(iss_tle):
    fp = compute_footprint(iss_tle.altitude_km, 0.0)
    lam = math.acos(EARTH_RADIUS_KM / (EARTH_RADIUS_KM + iss_tle.altitude_km))
    assert fp.central_angle_deg == pytest.approx(math.degrees(lam), abs=0.01)
    assert fp.central_angle_deg == pytest.approx(20.20, abs=0.05)
    assert fp.footprint_radius_km == pytest.approx(2248.0, abs=5.0)
    assert fp.nadir_angle_deg == pytest.approx(69.79, abs=0.05)
    assert fp.diameter_km == pytest.approx(2 * fp.footprint_radius_km)


def test_footprint_area_formula():
    fp = compute_footprint(500.0, 0.0)
    expected = 2 * math.pi * EARTH_RADIUS_KM ** 2 * (
        1 - math.cos(math.radians(fp.central_angle_deg)))
    assert fp.footprint_area_km2 == pytest.approx(expected, rel=1e-9)


def test_mask_shrinks_footprint():
    fp0 = compute_footprint(418.4, 0.0)
    fp10 = compute_footprint(418.4, 10.0)
    fp30 = compute_footprint(418.4, 30.0)
    assert fp10.footprint_radius_km < fp0.footprint_radius_km
    assert fp30.footprint_radius_km < fp10.footprint_radius_km
    assert fp10.central_angle_deg == pytest.approx(12.46, abs=0.05)


def test_higher_altitude_bigger_footprint():
    assert (compute_footprint(850.0, 0.0).footprint_radius_km
            > compute_footprint(400.0, 0.0).footprint_radius_km)


def test_max_slant_range_matches_geometry():
    # At the mask, range from slant_range_km must equal the max range.
    from satellite_pass_predictor.coordinates import slant_range_km
    assert max_slant_range_km(418.4, 10.0) == pytest.approx(
        slant_range_km(418.4, 10.0), rel=1e-12)


def test_invalid_altitude_rejected():
    with pytest.raises(ValueError):
        compute_footprint(0.0)
    with pytest.raises(ValueError):
        compute_footprint(-100.0)


def test_invalid_mask_rejected():
    with pytest.raises(ValueError):
        compute_footprint(400.0, 90.0)
    with pytest.raises(ValueError):
        compute_footprint(400.0, -5.0)


def test_ground_track_speed_iss(iss_tle):
    # v_orbit * Re/(Re+h): ~7.66 * 6378/6796 = 7.19 km/s.
    v = ground_track_speed_km_s(iss_tle.altitude_km)
    assert v == pytest.approx(7.19, abs=0.05)


def test_ground_track_speed_decreases_with_altitude():
    assert (ground_track_speed_km_s(850.0)
            < ground_track_speed_km_s(400.0))
