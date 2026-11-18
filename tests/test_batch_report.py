"""Tests for batch prediction and reports (batch.py, report.py)."""

from __future__ import annotations

from datetime import timedelta

import pytest

from satellite_pass_predictor.batch import BatchResult, predict_batch
from satellite_pass_predictor.report import (
    pass_table_markdown,
    render_html_report,
    render_markdown_report,
)


def test_batch_total_is_sum_of_pairs(tles, stations, day_start):
    result = predict_batch(tles, stations, day_start,
                           day_start + timedelta(hours=24),
                           min_elevation_deg=10.0)
    total = 0
    for tle in tles:
        for st in stations:
            total += len(result.passes_for(
                tle.name or f"NORAD {tle.sat_number}", st.name))
    assert result.total_passes == total
    assert result.total_passes > 10


def test_batch_sorted_by_aos(tles, stations, day_start):
    result = predict_batch(tles, stations, day_start,
                           day_start + timedelta(hours=24),
                           min_elevation_deg=10.0)
    aos_list = [p.aos for p in result.passes]
    assert aos_list == sorted(aos_list)


def test_batch_best_pass(tles, stations, day_start):
    result = predict_batch(tles, stations, day_start,
                           day_start + timedelta(hours=24),
                           min_elevation_deg=10.0)
    best = result.best_pass("ISS (ZARYA)", "Jakarta")
    assert best is not None
    assert best.max_elevation_deg == pytest.approx(50.0, abs=0.5)


def test_batch_best_pass_missing_pair(tles, stations, day_start):
    result = predict_batch(tles, stations, day_start,
                           day_start + timedelta(hours=24),
                           min_elevation_deg=10.0)
    assert result.best_pass("ISS (ZARYA)", "Nowhere") is None


def test_batch_total_visible_minutes(tles, stations, day_start):
    result = predict_batch(tles, stations, day_start,
                           day_start + timedelta(hours=24),
                           min_elevation_deg=10.0)
    total = result.total_visible_minutes("ISS (ZARYA)", "Jakarta")
    # 3 passes: ~1.0 + 5.6 + 6.6 min.
    assert total == pytest.approx(13.2, abs=0.5)


def test_markdown_report_contents(tles, stations, day_start):
    result = predict_batch(tles, stations, day_start,
                           day_start + timedelta(hours=24),
                           min_elevation_deg=10.0)
    report = render_markdown_report(result, tles, stations, day_start,
                                    day_start + timedelta(hours=24),
                                    10.0, frequency_hz=437e6)
    assert "# Satellite Pass Report" in report
    assert "ISS (ZARYA)" in report
    assert "Jakarta" in report
    assert "437.000 MHz" in report or "437" in report
    assert "| # | Satellite |" in report


def test_markdown_report_no_doppler_section_without_frequency(
        tles, stations, day_start):
    result = predict_batch(tles, stations, day_start,
                           day_start + timedelta(hours=24),
                           min_elevation_deg=10.0)
    report = render_markdown_report(result, tles, stations, day_start,
                                    day_start + timedelta(hours=24), 10.0)
    assert "Doppler" not in report


def test_markdown_report_empty_result(stations, day_start):
    from satellite_pass_predictor.ground_station import GroundStation
    result = BatchResult(passes=[])
    report = render_markdown_report(result, [], [GroundStation("X", 89.0, 0.0)],
                                    day_start, day_start + timedelta(hours=24), 10.0)
    assert "No passes above the mask" in report


def test_pass_table_markdown(tles, stations, day_start):
    result = predict_batch(tles, stations, day_start,
                           day_start + timedelta(hours=24),
                           min_elevation_deg=10.0)
    table = pass_table_markdown(result.passes)
    assert table.count("\n") == len(result.passes) + 1  # header + separator + rows
    assert "ISS (ZARYA)" in table


def test_html_report(tles, stations, day_start):
    result = predict_batch(tles, stations, day_start,
                           day_start + timedelta(hours=24),
                           min_elevation_deg=10.0)
    html = render_html_report("Test Report", result)
    assert "<table>" in html
    assert "ISS (ZARYA)" in html
    assert f"Total passes: {result.total_passes}" in html


def test_html_report_empty():
    html = render_html_report("Empty", BatchResult(passes=[]))
    assert "No passes." in html
