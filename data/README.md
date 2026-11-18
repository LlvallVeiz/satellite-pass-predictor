# TLE data sources
#
# data/satellites.tle — real TLEs downloaded from CelesTrak
# (https://celestrak.org/NORAD/elements/gp.php?CATNR=<number>&FORMAT=tle)
# on 2026-11-17, all with epoch 2026-11-16 (day 320 of 2026):
#   ISS (ZARYA)        NORAD 25544
#   NOAA 19            NORAD 33591
#   CSS (TIANHE)       NORAD 48274
#   SAUDISAT 1C (SO-50) NORAD 27607
# TLE data is public orbital data published by the U.S. Space Force
# 18th Space Defense Squadron and redistributed by CelesTrak.
# The parser validates the modulo-10 checksum of every line, so any
# transcription error is caught on load.
