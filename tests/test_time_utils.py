"""Tests for time utilities: Julian dates and GMST (time_utils.py).

Reference values
----------------
* JD 2451545.0 = 2000-01-01 12:00 UTC (J2000.0, by definition).
* JD 2440587.5 = 1970-01-01 00:00 UTC (Unix epoch).
* GMST at J2000.0 = 280.46061837 deg (Meeus eq. 12.4, exact by construction).
* GMST on 2025-06-21 20:00 UTC = 210.2666 deg (Meeus eq. 12.4;
  cross-checked against an independent online Meeus-based calculator).
"""

from __future__ import annotations

from datetime import datetime, timezone

import pytest

from satellite_pass_predictor.time_utils import (
    datetime_to_julian_date,
    gmst_degrees,
    julian_date_to_datetime,
    tle_epoch_to_datetime,
)


def test_julian_date_j2000():
    assert datetime_to_julian_date(datetime(2000, 1, 1, 12, 0, 0)) == pytest.approx(2451545.0)


def test_julian_date_unix_epoch():
    assert datetime_to_julian_date(datetime(1970, 1, 1, 0, 0, 0)) == pytest.approx(2440587.5)


def test_julian_date_fractional_day():
    # 2025-06-21 20:00 UTC -> JD 2460848.33333 (hand-checked day count).
    assert datetime_to_julian_date(datetime(2025, 6, 21, 20, 0, 0)) == pytest.approx(2460848.33333, abs=1e-4)


def test_julian_date_round_trip():
    for dt in (datetime(2000, 1, 1, 12, 0, 0),
               datetime(2026, 11, 16, 13, 27, 14),
               datetime(1999, 12, 31, 23, 59, 59)):
        back = julian_date_to_datetime(datetime_to_julian_date(dt))
        assert abs((back - dt).total_seconds()) < 1.0


def test_julian_date_aware_input():
    aware = datetime(2000, 1, 1, 12, 0, 0, tzinfo=timezone.utc)
    assert datetime_to_julian_date(aware) == pytest.approx(2451545.0)


def test_gmst_j2000():
    assert gmst_degrees(datetime(2000, 1, 1, 12, 0, 0)) == pytest.approx(280.46061837, abs=1e-6)


def test_gmst_known_value():
    # Independent Meeus-based calculator gives 210.2666 deg.
    assert gmst_degrees(datetime(2025, 6, 21, 20, 0, 0)) == pytest.approx(210.2666, abs=0.01)


def test_gmst_advances_one_sidereal_day():
    # GMST advances 360.9856... deg per solar day (one extra rotation
    # per year is spread over the year: 360 deg + 360/365.25 deg/day).
    t0 = datetime(2026, 11, 16, 0, 0, 0)
    t1 = datetime(2026, 11, 17, 0, 0, 0)
    delta = (gmst_degrees(t1) - gmst_degrees(t0)) % 360.0
    assert delta == pytest.approx(0.98564736629, abs=1e-6)


def test_gmst_range():
    for hour in range(0, 24, 3):
        g = gmst_degrees(datetime(2026, 11, 16, hour, 0, 0))
        assert 0.0 <= g < 360.0


def test_tle_epoch_year_2000s():
    dt = tle_epoch_to_datetime(26, 320.56058678)
    assert (dt.year, dt.month, dt.day) == (2026, 11, 16)
    assert dt.hour == 13 and dt.minute == 27


def test_tle_epoch_year_1900s():
    # Year 98 -> 1998 (57-99 map to 1957-1999).
    dt = tle_epoch_to_datetime(98, 1.0)
    assert dt.year == 1998


def test_tle_epoch_year_boundary():
    assert tle_epoch_to_datetime(56, 1.0).year == 2056
    assert tle_epoch_to_datetime(57, 1.0).year == 1957
