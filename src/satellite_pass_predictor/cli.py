"""Command-line interface for satellite-pass-predictor.

Subcommands
-----------
tles        parse and summarise a TLE file
propagate   ECI position/velocity of a satellite at a time
passes      predict passes of one satellite over one station
batch       predict passes for all satellites x stations in a TLE file
footprint   coverage geometry for an altitude
doppler     Doppler shift of a downlink along one pass
"""

from __future__ import annotations

import argparse
import sys
from datetime import datetime, timedelta

from .batch import predict_batch
from .coordinates import ecef_to_geodetic, eci_to_ecef
from .doppler import doppler_shift_hz
from .footprint import compute_footprint
from .ground_station import GroundStation
from .orbit import Propagator
from .passes import predict_passes
from .report import render_html_report, render_markdown_report
from .tle import load_tle_file


def _parse_time(text: str) -> datetime:
    try:
        return datetime.strptime(text, "%Y-%m-%dT%H:%M:%S")
    except ValueError:
        raise argparse.ArgumentTypeError(
            f"bad time {text!r}, expected YYYY-MM-DDTHH:MM:SS (UTC)")


def _station(args: argparse.Namespace) -> GroundStation:
    return GroundStation(args.station_name, args.lat, args.lon, args.height)


def _load_single_tle(path: str, sat_number: int | None):
    tles = load_tle_file(path)
    if sat_number is None:
        if len(tles) != 1:
            raise SystemExit(
                f"{path} contains {len(tles)} satellites; pass --sat to choose one")
        return tles[0]
    for tle in tles:
        if tle.sat_number == sat_number:
            return tle
    raise SystemExit(f"satellite {sat_number} not found in {path}")


def cmd_tles(args: argparse.Namespace) -> int:
    tles = load_tle_file(args.tle_file)
    for tle in tles:
        print(f"{tle.name} (NORAD {tle.sat_number}, {tle.intl_designator})")
        print(f"  epoch        : {tle.epoch:%Y-%m-%d %H:%M:%S} UTC")
        print(f"  period       : {tle.period_minutes:.2f} min "
              f"({tle.mean_motion_rev_per_day:.5f} rev/day)")
        print(f"  semi-major a : {tle.semi_major_axis_km:.2f} km, "
              f"altitude {tle.altitude_km:.1f} km")
        print(f"  inclination  : {tle.inclination_deg:.4f} deg, "
              f"e = {tle.eccentricity:.7f}")
        print(f"  RAAN {tle.raan_deg:.4f} deg, arg perigee "
              f"{tle.arg_perigee_deg:.4f} deg, M0 {tle.mean_anomaly_deg:.4f} deg")
        print(f"  B*           : {tle.bstar:.3e} /Earth radii "
              f"(parsed, not used by the Kepler+J2 propagator)")
    return 0


def cmd_propagate(args: argparse.Namespace) -> int:
    tle = _load_single_tle(args.tle_file, args.sat)
    prop = Propagator(tle)
    state = prop.propagate(args.time)
    ecef = eci_to_ecef(state.position_km, args.time)
    geo = ecef_to_geodetic(ecef)
    print(f"{tle.name} at {args.time:%Y-%m-%d %H:%M:%S} UTC")
    print(f"  ECI position (km) : {state.position_km[0]:.3f} "
          f"{state.position_km[1]:.3f} {state.position_km[2]:.3f}")
    print(f"  ECI velocity (km/s): {state.velocity_km_s[0]:.4f} "
          f"{state.velocity_km_s[1]:.4f} {state.velocity_km_s[2]:.4f}")
    print(f"  radius {state.radius_km:.3f} km, speed {state.speed_km_s:.4f} km/s")
    print(f"  sub-satellite point: lat {geo.latitude_deg:.4f} deg, "
          f"lon {geo.longitude_deg:.4f} deg, height {geo.height_m / 1000:.2f} km")
    return 0


def cmd_passes(args: argparse.Namespace) -> int:
    tle = _load_single_tle(args.tle_file, args.sat)
    station = _station(args)
    prop = Propagator(tle)
    end = args.start + timedelta(hours=args.hours)
    passes = predict_passes(prop, station, args.start, end,
                            min_elevation_deg=args.mask,
                            scan_step_s=args.scan_step)
    if not passes:
        print("No passes above the mask in this window.")
        return 0
    for p in passes:
        print(p.summary())
    return 0


def cmd_batch(args: argparse.Namespace) -> int:
    tles = load_tle_file(args.tle_file)
    stations = [
        GroundStation(name, lat, lon, height)
        for name, lat, lon, height in args.stations
    ]
    end = args.start + timedelta(hours=args.hours)
    result = predict_batch(tles, stations, args.start, end,
                           min_elevation_deg=args.mask,
                           scan_step_s=args.scan_step)
    report = render_markdown_report(result, tles, stations, args.start, end,
                                    args.mask, frequency_hz=args.frequency)
    if args.output:
        with open(args.output, "w", encoding="utf-8") as fh:
            fh.write(report)
        print(f"Wrote {args.output} ({result.total_passes} passes)")
    else:
        print(report)
    if args.html:
        with open(args.html, "w", encoding="utf-8") as fh:
            fh.write(render_html_report("Satellite Pass Report", result))
        print(f"Wrote {args.html}")
    return 0


def cmd_footprint(args: argparse.Namespace) -> int:
    fp = compute_footprint(args.altitude, args.mask)
    print(f"altitude {fp.altitude_km:.1f} km, mask {fp.min_elevation_deg:.1f} deg")
    print(f"  central angle : {fp.central_angle_deg:.3f} deg")
    print(f"  footprint radius: {fp.footprint_radius_km:.1f} km "
          f"(diameter {fp.diameter_km:.1f} km)")
    print(f"  nadir angle   : {fp.nadir_angle_deg:.3f} deg")
    print(f"  footprint area: {fp.footprint_area_km2:,.0f} km^2")
    return 0


def cmd_doppler(args: argparse.Namespace) -> int:
    tle = _load_single_tle(args.tle_file, args.sat)
    station = _station(args)
    prop = Propagator(tle)
    end = args.start + timedelta(hours=args.hours)
    passes = predict_passes(prop, station, args.start, end,
                            min_elevation_deg=args.mask)
    if not passes:
        print("No passes above the mask in this window.")
        return 0
    best = max(passes, key=lambda p: p.max_elevation_deg)
    print(best.summary())
    print(f"  downlink {args.frequency / 1e6:.3f} MHz")
    for (t, _az, el, _rng) in best.track:
        shift = doppler_shift_hz(prop.propagate(t), station, args.frequency)
        print(f"  {t:%H:%M:%S}  el {el:6.2f} deg  doppler {shift / 1000:+8.2f} kHz")
    return 0


def _station_arg(text: str) -> tuple[str, float, float, float]:
    parts = text.split(",")
    if len(parts) not in (3, 4):
        raise argparse.ArgumentTypeError(
            "station must be NAME,LAT,LON[,HEIGHT_M]")
    try:
        height = float(parts[3]) if len(parts) == 4 else 0.0
        return parts[0], float(parts[1]), float(parts[2]), height
    except ValueError:
        raise argparse.ArgumentTypeError(f"bad station {text!r}")


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        prog="satellite-pass-predictor",
        description="Predict satellite passes over ground stations "
                    "(Kepler + secular J2 propagation from TLEs).")
    sub = parser.add_subparsers(dest="command", required=True)

    p = sub.add_parser("tles", help="summarise a TLE file")
    p.add_argument("tle_file")
    p.set_defaults(func=cmd_tles)

    p = sub.add_parser("propagate", help="ECI state at a time")
    p.add_argument("tle_file")
    p.add_argument("--sat", type=int, default=None)
    p.add_argument("--time", type=_parse_time, required=True)
    p.set_defaults(func=cmd_propagate)

    p = sub.add_parser("passes", help="predict passes over one station")
    p.add_argument("tle_file")
    p.add_argument("--sat", type=int, default=None)
    p.add_argument("--station-name", default="Station")
    p.add_argument("--lat", type=float, required=True)
    p.add_argument("--lon", type=float, required=True)
    p.add_argument("--height", type=float, default=0.0)
    p.add_argument("--start", type=_parse_time, required=True)
    p.add_argument("--hours", type=float, default=24.0)
    p.add_argument("--mask", type=float, default=10.0)
    p.add_argument("--scan-step", type=float, default=30.0)
    p.set_defaults(func=cmd_passes)

    p = sub.add_parser("batch", help="predict passes for all sats x stations")
    p.add_argument("tle_file")
    p.add_argument("--stations", type=_station_arg, nargs="+", required=True,
                   metavar="NAME,LAT,LON[,HEIGHT_M]")
    p.add_argument("--start", type=_parse_time, required=True)
    p.add_argument("--hours", type=float, default=24.0)
    p.add_argument("--mask", type=float, default=10.0)
    p.add_argument("--scan-step", type=float, default=30.0)
    p.add_argument("--frequency", type=float, default=None,
                   help="downlink frequency in Hz (adds a Doppler section)")
    p.add_argument("--output", default=None, help="Markdown report path")
    p.add_argument("--html", default=None, help="HTML report path")
    p.set_defaults(func=cmd_batch)

    p = sub.add_parser("footprint", help="coverage geometry for an altitude")
    p.add_argument("--altitude", type=float, required=True)
    p.add_argument("--mask", type=float, default=0.0)
    p.set_defaults(func=cmd_footprint)

    p = sub.add_parser("doppler", help="Doppler shift along the best pass")
    p.add_argument("tle_file")
    p.add_argument("--sat", type=int, default=None)
    p.add_argument("--station-name", default="Station")
    p.add_argument("--lat", type=float, required=True)
    p.add_argument("--lon", type=float, required=True)
    p.add_argument("--height", type=float, default=0.0)
    p.add_argument("--start", type=_parse_time, required=True)
    p.add_argument("--hours", type=float, default=24.0)
    p.add_argument("--mask", type=float, default=10.0)
    p.add_argument("--frequency", type=float, default=437e6)
    p.set_defaults(func=cmd_doppler)

    return parser


def main(argv: list[str] | None = None) -> int:
    parser = build_parser()
    args = parser.parse_args(argv)
    return args.func(args)


if __name__ == "__main__":
    sys.exit(main())
