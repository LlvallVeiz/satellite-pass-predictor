"""Report generation: Markdown and simple HTML summaries of passes."""

from __future__ import annotations

import html
from datetime import datetime

from .batch import BatchResult
from .doppler import doppler_shift_hz
from .footprint import compute_footprint
from .ground_station import GroundStation
from .orbit import Propagator
from .passes import Pass
from .tle import TLE


def pass_table_markdown(passes: list[Pass]) -> str:
    """Markdown table of passes."""
    lines = [
        "| # | Satellite | Station | AOS (UTC) | LOS (UTC) | Max el (deg) | "
        "Duration (min) | AOS az (deg) | LOS az (deg) | Min range (km) |",
        "|---|---|---|---|---|---|---|---|---|---|",
    ]
    for idx, p in enumerate(passes, start=1):
        lines.append(
            f"| {idx} | {p.satellite} | {p.station} | {p.aos:%Y-%m-%d %H:%M:%S} | "
            f"{p.los:%Y-%m-%d %H:%M:%S} | {p.max_elevation_deg:.1f} | "
            f"{p.duration_min:.1f} | {p.aos_azimuth_deg:.1f} | "
            f"{p.los_azimuth_deg:.1f} | {p.min_range_km:.0f} |")
    return "\n".join(lines)


def render_markdown_report(result: BatchResult, tles: list[TLE],
                           stations: list[GroundStation],
                           start: datetime, end: datetime,
                           min_elevation_deg: float,
                           frequency_hz: float | None = None,
                           title: str = "Satellite Pass Report") -> str:
    """Build a Markdown pass report for a batch result."""
    lines = [
        f"# {title}",
        "",
        f"Window: {start:%Y-%m-%d %H:%M:%S} UTC to {end:%Y-%m-%d %H:%M:%S} UTC.  ",
        f"Elevation mask: {min_elevation_deg:.0f} deg.  ",
        f"Total passes: {result.total_passes}.",
        "",
        "## Satellites",
        "",
    ]
    for tle in tles:
        prop = Propagator(tle)
        fp = compute_footprint(tle.altitude_km, min_elevation_deg)
        lines.append(
            f"- **{tle.name}** (NORAD {tle.sat_number}): period "
            f"{tle.period_minutes:.2f} min, mean altitude "
            f"{tle.altitude_km:.1f} km, inclination "
            f"{tle.inclination_deg:.2f} deg, footprint radius "
            f"{fp.footprint_radius_km:.0f} km at {min_elevation_deg:.0f} deg mask.")
    lines += ["", "## Stations", ""]
    for st in stations:
        lines.append(
            f"- **{st.name}**: {st.latitude_deg:.2f} deg, "
            f"{st.longitude_deg:.2f} deg, {st.height_m:.0f} m.")
    lines += ["", "## Passes", ""]
    if result.passes:
        lines.append(pass_table_markdown(result.passes))
    else:
        lines.append("No passes above the mask in this window.")
    lines += ["", "## Per-pair summary", ""]
    lines.append("| Satellite | Station | Passes | Best max el (deg) | "
                 "Longest (min) | Total visible (min) |")
    lines.append("|---|---|---|---|---|---|")
    for tle in tles:
        name = tle.name or f"NORAD {tle.sat_number}"
        for st in stations:
            found = result.passes_for(name, st.name)
            if not found:
                lines.append(f"| {name} | {st.name} | 0 | - | - | - |")
                continue
            best = max(found, key=lambda p: p.max_elevation_deg)
            longest = max(found, key=lambda p: p.duration_min)
            total = sum(p.duration_min for p in found)
            lines.append(
                f"| {name} | {st.name} | {len(found)} | "
                f"{best.max_elevation_deg:.1f} | {longest.duration_min:.1f} | "
                f"{total:.1f} |")
    if frequency_hz is not None:
        lines += ["", "## Doppler (best pass per pair)", "",
                  f"Downlink frequency: {frequency_hz / 1e6:.3f} MHz.  "
                  "Shift is positive while the satellite approaches.",
                  "",
                  "| Satellite | Station | Max shift (kHz) | Min shift (kHz) | "
                  "Shift at max elevation (Hz) |",
                  "|---|---|---|---|"]
        for tle in tles:
            name = tle.name or f"NORAD {tle.sat_number}"
            prop = Propagator(tle)
            for st in stations:
                best = result.best_pass(name, st.name)
                if best is None:
                    lines.append(f"| {name} | {st.name} | - | - | - |")
                    continue
                shifts = []
                for (t, _az, _el, _rng) in best.track:
                    shifts.append(doppler_shift_hz(prop.propagate(t), st, frequency_hz))
                at_max = doppler_shift_hz(
                    prop.propagate(best.max_elevation_time), st, frequency_hz)
                lines.append(
                    f"| {name} | {st.name} | {max(shifts) / 1000:.2f} | "
                    f"{min(shifts) / 1000:.2f} | {at_max:.0f} |")
    lines.append("")
    return "\n".join(lines)


def render_html_report(markdown_title: str, result: BatchResult) -> str:
    """Minimal standalone HTML summary (table of passes)."""
    rows = []
    for idx, p in enumerate(result.passes, start=1):
        rows.append(
            f"<tr><td>{idx}</td><td>{html.escape(p.satellite)}</td>"
            f"<td>{html.escape(p.station)}</td>"
            f"<td>{p.aos:%Y-%m-%d %H:%M:%S}</td>"
            f"<td>{p.los:%Y-%m-%d %H:%M:%S}</td>"
            f"<td>{p.max_elevation_deg:.1f}</td>"
            f"<td>{p.duration_min:.1f}</td></tr>")
    body = "\n".join(rows) if rows else "<tr><td colspan='7'>No passes.</td></tr>"
    return f"""<!DOCTYPE html>
<html lang="en">
<head><meta charset="utf-8">
<title>{html.escape(markdown_title)}</title>
<style>
body {{ font-family: system-ui, sans-serif; margin: 2em; }}
table {{ border-collapse: collapse; }}
th, td {{ border: 1px solid #999; padding: 4px 10px; text-align: left; }}
th {{ background: #eee; }}
</style></head>
<body>
<h1>{html.escape(markdown_title)}</h1>
<p>Total passes: {result.total_passes}</p>
<table>
<tr><th>#</th><th>Satellite</th><th>Station</th><th>AOS (UTC)</th>
<th>LOS (UTC)</th><th>Max el (deg)</th><th>Duration (min)</th></tr>
{body}
</table>
</body>
</html>
"""
