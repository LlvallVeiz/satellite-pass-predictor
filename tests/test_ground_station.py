"""Tests for ground stations (ground_station.py)."""

from __future__ import annotations

import pytest

from satellite_pass_predictor.ground_station import (
    BANDUNG,
    JAKARTA,
    SINGAPORE,
    GroundStation,
)


def test_jakarta_coordinates():
    assert JAKARTA.latitude_deg == pytest.approx(-6.2)
    assert JAKARTA.longitude_deg == pytest.approx(106.8)


def test_station_ecef_magnitude_sea_level():
    # Sea-level station: |ECEF| between polar and equatorial radius.
    station = GroundStation("Test", -6.2, 106.8, 0.0)
    import numpy as np
    r = float(np.linalg.norm(station.ecef_km))
    assert 6356.0 < r < 6378.2


def test_invalid_latitude_rejected():
    with pytest.raises(ValueError):
        GroundStation("Bad", 91.0, 0.0)


def test_invalid_longitude_rejected():
    with pytest.raises(ValueError):
        GroundStation("Bad", 0.0, 181.0)


def test_look_angles_iss_overhead_near_epoch(iss_prop, jakarta):
    # At the epoch the ISS sub-satellite point is far from Jakarta;
    # elevation must simply be a valid angle.
    state = iss_prop.propagate(iss_prop.tle.epoch)
    look = jakarta.look_angles(state)
    assert -90.0 <= look.elevation_deg <= 90.0
    assert 0.0 <= look.azimuth_deg < 360.0
    assert look.range_km > 400.0


def test_elevation_consistency(jakarta, iss_prop):
    state = iss_prop.propagate(iss_prop.tle.epoch)
    assert jakarta.elevation_deg(state) == pytest.approx(
        jakarta.look_angles(state).elevation_deg)


def test_known_stations():
    assert BANDUNG.height_m == pytest.approx(768.0)
    assert SINGAPORE.latitude_deg == pytest.approx(1.35)
