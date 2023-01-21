# Plotting BinCon data

```
Usage: plot.jl -i data [--cpt palette] [--inc increment] [--reg region]
    [--cnt contour] [-o outfile] [-h]
```

## -i data
data is either the name of a local data file, or a URL to one which can be
downloaded.

## --cpt palette
palette is either the name of one of the pre-defined GMT palettes, or is the
name of a file containing a custom colour palette. The default is the pre-defined
batlow colour palette.

## --cnt contour
This option might be redundant. Look into it.

## --inc increment
When converting the contours to a grid, use this increment for the grid spacing. Format:
loninc/latinc, where loninc is the longitudinal grid spacing in degrees, and latinc is
the latitudinal grid spacing. Append an "m" to loninc and latinc if the grid spacings
are in arc-minutes. The default is 15m/15m (0.25°/0.25°). It makes sense to set
loninc and latinc to the same values as used by the raw data which the contours were
derived from.

## -o outfile
Write the final plot to outfile. The format of the file will be determined by
the suffix of outfile. The format must
be one supported by the GMT
[psconvert](https://docs.generic-mapping-tools.org/dev/psconvert.html) command:

- BMP
- EPS
- JPEG
- PDF
- PNG
- PPM
- SVG
- TIFF
