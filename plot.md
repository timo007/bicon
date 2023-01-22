# Plotting BinCon data

## Usage
```
Usage: plot.jl -i data [--cpt palette] [--inc increment] [--reg region]
    [-o outfile] [-h]
```

### -i data
data is either the name of a local data file, or a URL to one which can be
downloaded. The default action is to download an [example data
file](https://github.com/timo007/bincon/blob/main/test_data/GFS_NZ_000-003-001_MSL_2023010412_024.bin?raw=true)
stored on GitHub.

### --cpt palette
palette is either the name of one of the pre-defined GMT palettes, or is the
name of a file containing a custom colour palette. The default is the
pre-defined batlow colour palette.

### --inc increment
When converting the contours to a grid, use this increment for the grid
spacing. Format: loninc/latinc, where loninc is the longitudinal grid spacing
in degrees, and latinc is the latitudinal grid spacing. Append an "m" to loninc
and latinc if the grid spacings are in arc-minutes. The default is 15m/15m
(0.25°/0.25°). It makes sense to set loninc and latinc to the same values as
used by the raw data which the contours were derived from.

### --reg region
region is the region to plot the map for. The regions are defined
[here](./regions.md)

### -o outfile
Write the final plot to outfile. The format of the file will be determined by
the suffix of outfile. The format must be one supported by the GMT
[psconvert](https://docs.generic-mapping-tools.org/dev/psconvert.html) command:

- BMP
- EPS
- JPEG
- PDF
- PNG
- PPM
- SVG
- TIFF

## Examples

```
julia plot.jl -i https://github.com/timo007/bincon/blob/main/test_data/GFS_NZ_000-003-001_MSL_2023010412_024.bin?raw=true -o example.pdf
```
This downloads and plots the example dataset stored on GitHub, and saves the
plot to example.pdf.


```
julia plot.jl -i data.bin --cpt turbo -o localdata.png
```
This plots the local data file, data.bin, using the GMT "turbo" colour palette and
saving the plot to localdata.png.
