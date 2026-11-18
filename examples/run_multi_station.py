"""Example 2: batch prediction — 4 satellites x 3 stations, 24 hours.

Run from the repository root::

    python3 examples/run_multi_station.py

Writes a Markdown report and a simple HTML summary to
``examples/output/`` and prints the per-pair summary table.
"""

from __future__ import annotations

import sys
from datetime import datetime, timedelta
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "src"))

from satellite_pass_predictor.batch import predict_batch  # noqa: E402
from satellite_pass_predictor.ground_station import BANDUNG, JAKARTA, SINGAPORE  # noqa: E402
from satellite_pass_predictor.report import render_html_report, render_markdown_report  # noqa: E402
from satellite_pass_predictor.tle import load_tle_file  # noqa: E402


def main() -> None:
    root = Path(__file__).resolve().parent.parent
    tles = load_tle_file(str(root / "data" / "satellites.tle"))
    stations = [JAKARTA, BANDUNG, SINGAPORE]
    start = datetime(2026, 11, 16, 0, 0, 0)
    end = start + timedelta(hours=24)
    result = predict_batch(tles, stations, start, end, min_elevation_deg=10.0)

    print(f"Window: {start:%Y-%m-%d %H:%M} to {end:%Y-%m-%d %H:%M} UTC, mask 10 deg")
    print(f"Total passes: {result.total_passes}")
    print()
    print(f"{'Satellite':<20} {'Station':<10} {'Passes':>6} "
          f"{'Best el':>8} {'Longest':>8} {'Visible':>8}")
    for tle in tles:
        name = tle.name
        for st in stations:
            found = result.passes_for(name, st.name)
            if not found:
                print(f"{name:<20} {st.name:<10} {0:>6} {'-':>8} {'-':>8} {'-':>8}")
                continue
            best = max(found, key=lambda p: p.max_elevation_deg)
            longest = max(found, key=lambda p: p.duration_min)
            total = sum(p.duration_min for p in found)
            print(f"{name:<20} {st.name:<10} {len(found):>6} "
                  f"{best.max_elevation_deg:>7.1f}° {longest.duration_min:>7.1f}m "
                  f"{total:>7.1f}m")

    out_dir = root / "examples" / "output"
    out_dir.mkdir(exist_ok=True)
    report = render_markdown_report(result, tles, stations, start, end, 10.0,
                                    frequency_hz=437e6)
    (out_dir / "report.md").write_text(report, encoding="utf-8")
    (out_dir / "report.html").write_text(
        render_html_report("Satellite Pass Report — 2026-11-16", result),
        encoding="utf-8")
    print()
    print(f"Wrote {out_dir / 'report.md'}")
    print(f"Wrote {out_dir / 'report.html'}")


if __name__ == "__main__":
    main()
