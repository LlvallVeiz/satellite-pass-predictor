"""Time utilities: Julian dates and Greenwich Mean Sidereal Time.

The GMST formula implemented here is Meeus, "Astronomical Algorithms",
2nd ed. (1998), eq. (12.4):

    GMST = 280.46061837 + 360.98564736629 * (JD - 2451545.0)
           + 0.000387933 * T^2 - T^3 / 38710000          [degrees]

with T = (JD - 2451545.0) / 36525, JD the Julian Date (UTC).  It is
accurate to well under a second of time for dates within a century of
J2000.0 (Meeus, ch. 12; also the IAU-adopted expression).  UT1 vs UTC
is ignored (difference < 0.9 s), as is nutation (GMST, not GAST).
"""

from __future__ import annotations

import math
from datetime import datetime, timezone

from .constants import J2000_JD, wrap_degrees


def datetime_to_julian_date(dt: datetime) -> float:
    """Convert an aware (or naive-as-UTC) datetime to a Julian Date.

    Uses the standard algorithm (Meeus ch. 7 / Fliegel & Van Flandern).
    """
    if dt.tzinfo is not None:
        dt = dt.astimezone(timezone.utc).replace(tzinfo=None)
    year, month = dt.year, dt.month
    day = dt.day + (dt.hour + dt.minute / 60.0 + dt.second / 3600.0
                    + dt.microsecond / 3.6e9) / 24.0
    if month <= 2:
        year -= 1
        month += 12
    a = math.floor(year / 100.0)
    b = 2 - a + math.floor(a / 4.0)
    jd = (math.floor(365.25 * (year + 4716))
          + math.floor(30.6001 * (month + 1))
          + day + b - 1524.5)
    return jd


def julian_date_to_datetime(jd: float) -> datetime:
    """Convert a Julian Date to a UTC datetime (naive, UTC)."""
    unix_epoch_jd = 2440587.5  # JD of 1970-01-01 00:00 UTC
    seconds = (jd - unix_epoch_jd) * 86400.0
    return datetime.fromtimestamp(seconds, tz=timezone.utc).replace(tzinfo=None)


def gmst_degrees(dt: datetime) -> float:
    """Greenwich Mean Sidereal Time in degrees, [0, 360)."""
    jd = datetime_to_julian_date(dt)
    t = (jd - J2000_JD) / 36525.0
    gmst = (280.46061837
            + 360.98564736629 * (jd - J2000_JD)
            + 0.000387933 * t * t
            - t ** 3 / 38710000.0)
    return wrap_degrees(gmst)


def gmst_radians(dt: datetime) -> float:
    """Greenwich Mean Sidereal Time in radians, [0, 2*pi)."""
    return math.radians(gmst_degrees(dt))


def tle_epoch_to_datetime(epoch_year: int, epoch_day: float) -> datetime:
    """Convert a TLE epoch (two-digit year, day-of-year fraction) to UTC.

    Two-digit years 57-99 map to 1957-1999, 00-56 to 2000-2056
    (CelesTrak TLE format FAQ, Kelso 1998).
    """
    full_year = 1900 + epoch_year if epoch_year >= 57 else 2000 + epoch_year
    base = datetime(full_year, 1, 1, tzinfo=timezone.utc)
    from datetime import timedelta
    dt = base + timedelta(days=epoch_day - 1.0)
    return dt.replace(tzinfo=None)
