"""Tests for TLE parsing and checksum validation (tle.py)."""

from __future__ import annotations

import pytest

from satellite_pass_predictor.tle import (
    TLEError,
    compute_checksum,
    load_tle_file,
    parse_tle,
    parse_tle_text,
)

ISS_L1 = "1 25544U 98067A   26320.56058678  .00009493  00000+0  17760-3 0  9996"
ISS_L2 = "2 25544  51.6320 266.0553 0007391  32.4819 327.6623 15.49415609595375"


def test_iss_line_lengths():
    assert len(ISS_L1) == 69
    assert len(ISS_L2) == 69


def test_checksum_known_values():
    assert compute_checksum(ISS_L1) == 6
    assert compute_checksum(ISS_L2) == 5


def test_checksum_counts_minus_as_one():
    # A line of 68 '-' would sum to 68 -> checksum 8.
    assert compute_checksum("-" * 68) == 8
    # Letters, blanks, periods and plus signs count as zero.
    assert compute_checksum("A . + " + "0" * 62) == 0


def test_parse_iss_fields(iss_tle):
    assert iss_tle.name == "ISS (ZARYA)"
    assert iss_tle.sat_number == 25544
    assert iss_tle.classification == "U"
    assert iss_tle.intl_designator == "98067A"
    assert iss_tle.epoch_year == 26
    assert iss_tle.epoch.day == 16 and iss_tle.epoch.month == 11
    assert iss_tle.epoch.year == 2026


def test_parse_iss_elements(iss_tle):
    assert iss_tle.inclination_deg == pytest.approx(51.6320)
    assert iss_tle.raan_deg == pytest.approx(266.0553)
    assert iss_tle.eccentricity == pytest.approx(0.0007391)
    assert iss_tle.arg_perigee_deg == pytest.approx(32.4819)
    assert iss_tle.mean_anomaly_deg == pytest.approx(327.6623)
    assert iss_tle.mean_motion_rev_per_day == pytest.approx(15.49415609)


def test_parse_bstar_implied_decimal(iss_tle):
    # ' 17760-3' means 0.17760e-3.
    assert iss_tle.bstar == pytest.approx(0.17760e-3)
    assert iss_tle.ndot == pytest.approx(0.00009493)
    assert iss_tle.nddot == pytest.approx(0.0)


def test_period_from_mean_motion(iss_tle):
    # Hand check: 86400 s / 15.49415609 rev/day = 5576.30 s = 92.938 min.
    assert iss_tle.period_seconds == pytest.approx(5576.30, abs=0.05)
    assert iss_tle.period_minutes == pytest.approx(92.9383, abs=0.001)


def test_semi_major_axis_iss(iss_tle):
    # Kepler: a = (mu / n^2)^(1/3); ISS a ~ 6796.5 km, altitude ~418 km.
    assert iss_tle.semi_major_axis_km == pytest.approx(6796.5, abs=1.0)
    assert iss_tle.altitude_km == pytest.approx(418.4, abs=1.0)


def test_all_four_satellites_load(tles):
    assert len(tles) == 4
    numbers = sorted(t.sat_number for t in tles)
    assert numbers == [25544, 27607, 33591, 48274]


def test_noaa_is_sun_synchronous(noaa_tle):
    assert noaa_tle.inclination_deg == pytest.approx(98.9372)
    assert noaa_tle.period_minutes == pytest.approx(101.875, abs=0.01)
    assert noaa_tle.altitude_km == pytest.approx(847.4, abs=1.0)


def test_css_and_so50(css_tle, so50_tle):
    assert css_tle.period_minutes == pytest.approx(92.269, abs=0.01)
    assert so50_tle.period_minutes == pytest.approx(97.078, abs=0.01)


def test_bad_checksum_rejected():
    bad = ISS_L1[:-1] + ("7" if ISS_L1[-1] != "7" else "8")
    with pytest.raises(TLEError, match="checksum"):
        parse_tle("ISS", bad, ISS_L2)


def test_wrong_line_length_rejected():
    with pytest.raises(TLEError, match="69 characters"):
        parse_tle("ISS", ISS_L1[:-1], ISS_L2)


def test_wrong_line_number_rejected():
    with pytest.raises(TLEError, match="start with '1'"):
        parse_tle("ISS", "2" + ISS_L1[1:], ISS_L2)


def test_mismatched_sat_numbers_rejected():
    # Same line 1 with the satellite number changed (checksum recomputed).
    line = "1 25545U" + ISS_L1[8:68]
    line = line + str(compute_checksum(line))
    with pytest.raises(TLEError, match="differ"):
        parse_tle("ISS", line, ISS_L2)


def test_parse_tle_text_multiple():
    text = f"ISS (ZARYA)\n{ISS_L1}\n{ISS_L2}\n"
    parsed = parse_tle_text(text)
    assert len(parsed) == 1
    assert parsed[0].sat_number == 25544


def test_parse_tle_text_without_name():
    text = f"{ISS_L1}\n{ISS_L2}\n"
    parsed = parse_tle_text(text)
    assert len(parsed) == 1
    assert parsed[0].name == ""


def test_parse_tle_text_garbage_rejected():
    with pytest.raises(TLEError):
        parse_tle_text("not a tle\nat all\n")


def test_epoch_year_mapping():
    # 57-99 -> 19xx, 00-56 -> 20xx (CelesTrak convention).
    tle = parse_tle("X", ISS_L1, ISS_L2)
    assert tle.epoch.year == 2026
