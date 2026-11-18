"""Tests for the command-line interface (cli.py)."""

from __future__ import annotations

from pathlib import Path

import pytest

from satellite_pass_predictor.cli import main

DATA = Path(__file__).resolve().parent.parent / "data" / "satellites.tle"


def test_cli_tles(capsys):
    assert main(["tles", str(DATA)]) == 0
    out = capsys.readouterr().out
    assert "ISS (ZARYA)" in out
    assert "92.94 min" in out
    assert "NOAA 19" in out


def test_cli_propagate(capsys):
    assert main(["propagate", str(DATA), "--sat", "25544",
                 "--time", "2026-11-16T13:27:14"]) == 0
    out = capsys.readouterr().out
    assert "ISS (ZARYA)" in out
    assert "ECI position" in out
    assert "sub-satellite point" in out


def test_cli_passes(capsys):
    assert main(["passes", str(DATA), "--sat", "25544",
                 "--station-name", "Jakarta", "--lat", "-6.2", "--lon", "106.8",
                 "--start", "2026-11-16T00:00:00", "--hours", "24"]) == 0
    out = capsys.readouterr().out
    assert "ISS (ZARYA) over Jakarta" in out
    assert out.count("AOS") == 3


def test_cli_passes_no_pass(capsys):
    assert main(["passes", str(DATA), "--sat", "25544",
                 "--station-name", "Pole", "--lat", "89", "--lon", "0",
                 "--start", "2026-11-16T00:00:00", "--hours", "24"]) == 0
    assert "No passes" in capsys.readouterr().out


def test_cli_footprint(capsys):
    assert main(["footprint", "--altitude", "418.4", "--mask", "10"]) == 0
    out = capsys.readouterr().out
    assert "footprint radius" in out
    assert "1386" in out  # ~1386.5 km at 10-deg mask


def test_cli_doppler(capsys):
    assert main(["doppler", str(DATA), "--sat", "25544",
                 "--station-name", "Jakarta", "--lat", "-6.2", "--lon", "106.8",
                 "--start", "2026-11-16T00:00:00", "--hours", "24",
                 "--frequency", "437000000"]) == 0
    out = capsys.readouterr().out
    assert "doppler" in out
    assert "kHz" in out


def test_cli_batch_markdown_and_html(tmp_path, capsys):
    md = tmp_path / "report.md"
    html = tmp_path / "report.html"
    assert main(["batch", str(DATA),
                 "--stations", "Jakarta,-6.2,106.8,8", "Bandung,-6.9,107.6,768",
                 "--start", "2026-11-16T00:00:00", "--hours", "24",
                 "--mask", "10", "--frequency", "437000000",
                 "--output", str(md), "--html", str(html)]) == 0
    assert md.exists() and html.exists()
    assert "ISS (ZARYA)" in md.read_text()
    assert "<table>" in html.read_text()
    assert "passes" in capsys.readouterr().out


def test_cli_batch_stdout(capsys):
    assert main(["batch", str(DATA),
                 "--stations", "Jakarta,-6.2,106.8,8",
                 "--start", "2026-11-16T00:00:00", "--hours", "24"]) == 0
    assert "ISS (ZARYA)" in capsys.readouterr().out


def test_cli_bad_time_rejected():
    with pytest.raises(SystemExit):
        main(["propagate", str(DATA), "--time", "not-a-time"])


def test_cli_missing_sat_rejected():
    with pytest.raises(SystemExit):
        main(["propagate", str(DATA), "--sat", "99999",
              "--time", "2026-11-16T00:00:00"])
