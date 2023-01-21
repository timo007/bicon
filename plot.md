# Plotting BinCon data

```
Usage: plot.jl -i data [--cpt palette] [--inc increment] [--reg region] [--cnt contour]
[-o outfile] [-h]
```

## -i data
data is either the name of a local data file, or a URL to one which can be
downloaded.

## --cpt palette
palette is either the name of one of the pre-defined GMT palettes, or is the
name of a file containing a custom colour palette. The default is the pre-defined
batlow colour palette.

## --cnt contour
Check what this actually does.

## --inc increment
When converting the contours to a grid, use this increment for the grid spacing. Format:
loninc/latinc, where loninc is the longitudinal grid spacing in degrees, and latinc is
the latitudinal grid spacing. Append an "m" to loninc and latinc if the grid spacings
are in arc-minutes. The default is "15m/15m" (0.25°/0.25°).
