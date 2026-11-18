"""Two-Line Element set (TLE) parsing and validation.

Format reference: Dr. T. S. Kelso, "Frequently Asked Questions: Two-Line
Element Set Format", CelesTrak (January 1998), mirrored at
https://www.scss.tcd.ie/Stephen.Farrell/ipn/background/computers-and-satellites/v04n03.html

Each data line is exactly 69 characters; column 69 is a modulo-10
checksum: sum of all digits in columns 1-68, counting each minus sign
as 1, ignoring letters, blanks, periods and plus signs; the checksum
is the last digit of that sum.

TLE elements are *mean* elements fitted for the SGP4/SDP4 models
(Hoots & Roehrich, Spacetrack Report No. 3, 1980).  This package does
NOT implement SGP4: it propagates the elements with a Kepler solver
plus secular J2 rates (see :mod:`orbit`).  That is an approximation,
accurate to a few kilometres over a day or two near the epoch for
near-circular LEO, and it is documented as such in the README.
"""

from __future__ import annotations

import re
from dataclasses import dataclass
from datetime import datetime

from .constants import DEG2RAD
from .time_utils import tle_epoch_to_datetime


class TLEError(ValueError):
    """Raised when a TLE is malformed or fails its checksum."""


def compute_checksum(line: str) -> int:
    """Modulo-10 checksum of a TLE line (columns 1-68, i.e. first 68 chars)."""
    total = 0
    for ch in line[:68]:
        if ch.isdigit():
            total += int(ch)
        elif ch == "-":
            total += 1
    return total % 10


def _parse_implied_decimal(raw: str) -> float:
    """Parse a TLE field with an implied decimal point, e.g. ' 17760-3'.

    The field is a mantissa with an implied leading decimal point
    followed by a signed exponent, i.e. 0.17760e-3.
    """
    raw = raw.strip()
    if not raw:
        return 0.0
    # Split mantissa / exponent: exponent is the last two chars (+/- digit).
    match = re.fullmatch(r"([+-]?\d+)([+-]\d)", raw)
    if match is None:
        return float(raw)
    mantissa = float("0." + match.group(1).lstrip("+-"))
    if match.group(1).startswith("-"):
        mantissa = -mantissa
    exponent = int(match.group(2))
    return mantissa * (10.0 ** exponent)


def _parse_float(raw: str) -> float:
    raw = raw.strip()
    return float(raw) if raw else 0.0


@dataclass(frozen=True)
class TLE:
    """One parsed two-line element set.

    Angles are stored in degrees (as in the TLE); helper properties
    expose radians.  ``mean_motion`` is in revolutions per day.
    ``bstar`` is the SGP4 drag term (1 / Earth radii); it is parsed
    and recorded but NOT used by the Kepler+J2 propagator.
    """

    name: str
    sat_number: int
    classification: str
    intl_designator: str
    epoch: datetime
    epoch_year: int
    epoch_day: float
    ndot: float  # 1st derivative of mean motion, rev/day^2
    nddot: float  # 2nd derivative of mean motion, rev/day^3
    bstar: float  # drag term, 1/Earth radii
    ephemeris_type: int
    element_set_number: int
    inclination_deg: float
    raan_deg: float
    eccentricity: float
    arg_perigee_deg: float
    mean_anomaly_deg: float
    mean_motion_rev_per_day: float
    rev_number: int

    # -- derived quantities -------------------------------------------------
    @property
    def inclination(self) -> float:
        """Inclination in radians."""
        return self.inclination_deg * DEG2RAD

    @property
    def raan(self) -> float:
        """Right ascension of the ascending node in radians."""
        return self.raan_deg * DEG2RAD

    @property
    def arg_perigee(self) -> float:
        """Argument of perigee in radians."""
        return self.arg_perigee_deg * DEG2RAD

    @property
    def mean_anomaly(self) -> float:
        """Mean anomaly at epoch in radians."""
        return self.mean_anomaly_deg * DEG2RAD

    @property
    def mean_motion_rad_s(self) -> float:
        """Mean motion in rad/s."""
        return self.mean_motion_rev_per_day * 2.0 * 3.141592653589793 / 86400.0

    @property
    def period_seconds(self) -> float:
        """Keplerian orbital period in seconds, from mean motion."""
        return 86400.0 / self.mean_motion_rev_per_day

    @property
    def period_minutes(self) -> float:
        """Keplerian orbital period in minutes."""
        return self.period_seconds / 60.0

    @property
    def semi_major_axis_km(self) -> float:
        """Semi-major axis in km, from mean motion (two-body relation)."""
        from .constants import MU_EARTH_KM3_S2
        n = self.mean_motion_rad_s
        return (MU_EARTH_KM3_S2 / (n * n)) ** (1.0 / 3.0)

    @property
    def altitude_km(self) -> float:
        """Mean altitude above a spherical Earth of radius 6378.137 km."""
        from .constants import EARTH_RADIUS_KM
        return self.semi_major_axis_km - EARTH_RADIUS_KM


def parse_tle(name: str, line1: str, line2: str) -> TLE:
    """Parse a three-line element set (name + two data lines).

    Raises :class:`TLEError` on wrong line lengths, wrong line numbers,
    mismatched satellite numbers, or checksum failure.
    """
    line1 = line1.rstrip("\n")
    line2 = line2.rstrip("\n")
    if len(line1) != 69:
        raise TLEError(f"line 1 must be 69 characters, got {len(line1)}")
    if len(line2) != 69:
        raise TLEError(f"line 2 must be 69 characters, got {len(line2)}")
    if not line1.startswith("1"):
        raise TLEError("line 1 must start with '1'")
    if not line2.startswith("2"):
        raise TLEError("line 2 must start with '2'")
    if compute_checksum(line1) != int(line1[68]):
        raise TLEError(
            f"line 1 checksum mismatch: computed {compute_checksum(line1)}, "
            f"stored {line1[68]}")
    if compute_checksum(line2) != int(line2[68]):
        raise TLEError(
            f"line 2 checksum mismatch: computed {compute_checksum(line2)}, "
            f"stored {line2[68]}")
    sat1 = int(line1[2:7])
    sat2 = int(line2[2:7])
    if sat1 != sat2:
        raise TLEError(f"satellite numbers differ: {sat1} vs {sat2}")

    epoch_year = int(line1[18:20])
    epoch_day = float(line1[20:32])
    epoch = tle_epoch_to_datetime(epoch_year, epoch_day)

    return TLE(
        name=name.strip(),
        sat_number=sat1,
        classification=line1[7:8],
        intl_designator=line1[9:17].strip(),
        epoch=epoch,
        epoch_year=epoch_year,
        epoch_day=epoch_day,
        ndot=_parse_float(line1[33:43]),
        nddot=_parse_implied_decimal(line1[44:52]),
        bstar=_parse_implied_decimal(line1[53:61]),
        ephemeris_type=int(line1[62:63]),
        element_set_number=int(line1[64:68]),
        inclination_deg=_parse_float(line2[8:16]),
        raan_deg=_parse_float(line2[17:25]),
        eccentricity=float("0." + line2[26:33].strip()),
        arg_perigee_deg=_parse_float(line2[34:42]),
        mean_anomaly_deg=_parse_float(line2[43:51]),
        mean_motion_rev_per_day=_parse_float(line2[52:63]),
        rev_number=int(line2[63:68]),
    )


def parse_tle_text(text: str) -> list[TLE]:
    """Parse a block of text containing one or more 3-line element sets.

    Blank lines are ignored.  A set is a name line followed by lines
    starting with '1' and '2'.
    """
    lines = [ln.rstrip("\n") for ln in text.splitlines() if ln.strip()]
    result: list[TLE] = []
    i = 0
    while i < len(lines):
        if lines[i].startswith("1 ") and i + 1 < len(lines) and lines[i + 1].startswith("2 "):
            result.append(parse_tle("", lines[i], lines[i + 1]))
            i += 2
        elif i + 2 < len(lines) and lines[i + 1].startswith("1 ") and lines[i + 2].startswith("2 "):
            result.append(parse_tle(lines[i], lines[i + 1], lines[i + 2]))
            i += 3
        else:
            raise TLEError(f"unrecognised TLE block starting at line {i + 1}: {lines[i]!r}")
    return result


def load_tle_file(path: str) -> list[TLE]:
    """Load TLEs from a text file (3-line format)."""
    with open(path, "r", encoding="utf-8") as fh:
        return parse_tle_text(fh.read())
