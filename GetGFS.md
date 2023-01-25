# Downloading GFS data and compressing it with BinCon

This Julia script downloads recent GFS forecasts, and compresses them using
BinCon. The user specifies an NWP base-time, parameter, and the tolerance and
contour interval compression parameters.

The script currently collects the 24, 48, 72, 96 and 120 hour lead time
forecasts. It can easily be modified to collect other forecast lead times.

## Usage
```
Usage: julia GetGFS.jl -t yyyymmddHH [-v var] [-p lev] [--tol tol] [--cnt cint] [--reg region] [-h]
```

### -t yyyymmddHH
Is the base date and time of the GFS data to download. yyyy is year, mm is
month, dd is day, HH is hour. Times and dates are in UTC.

### -v var
Is the GFS parameter to download. Default is "prsmslmsl". Currently recognised
parameters (defined in [BinaryContour.jl](./BinaryContour.jl)) are:

| Var             | Description                                | GRIB codes   | Level     |
|-----------------|--------------------------------------------|--------------|-----------|
| tmax2m          | Maximum temperature at 2 m above surface   | 0, 0, 4      | 2m        |
| tmin2m          | Minimum temperature at 2 m above surface   | 0, 0, 5      | 2m        |
| tmpsfc          | Surface temperature                        | 0, 0, 0      | Surface   |
| tmp2m           | Temperature 2 m above surface              | 0, 0, 0      | 2m        |
| tmpprs          | Temperature on pressure level              | 0, 0, 0      | -p lev    |
| rh2m            | Relative humidity 2 m above surface        | 0, 1, 1      | 2m        |
| apcpsfc         | Accumulated precipitation at the surface   | 0, 1, 8      | Surface   |
| prateavesfc     | Averarge precipitation rate at the surface | 0, 1, 52     | Surface   |
| gustsfc         | Wind gust speed at the surface             | 0, 2, 22     | Surface   |
| prmslmsl        | Mean sea level pressure                    | 0, 3, 1      | MSL       |
| hgtprs          | Geopotential heigh of pressure level       | 0, 3, 5      | -v lev    |
| tozneclm        | Total ozone content in atmospheric column  | 0, 14, 0     | Atmos     |

The GRIB codes refer to the WMO discipline, category and parameter codes, which
are defined in [GRIB2](https://library.wmo.int/?lvl=notice_display&id=10684#.Y9D2OBNBzWQ) 
itables [0.0](http://codes.wmo.int/grib2/codeflag/_0.0),
[4.1](http://codes.wmo.int/grib2/codeflag/_4.1) and
[4.2](http://codes.wmo.int/grib2/codeflag/_4.2) respectively. Where level is "-p
lev" it means the pressure level is specified by the user using the -p lev
option described below.

Adding more parameters is simply a matter of editing the dictionary of
parameters in [BinaryContour.jl](./BinaryContour.jl). Units of the parameters
are as specified in WMO GRIB2 table
[4.2](http://codes.wmo.int/grib2/codeflag/_4.2).

### -p lev
When a pressure level parameter is being downloaded and compressed, specify the
level with this option. The level is in hPa (e.g. 500).

### --tol tol
tol is the tolerance (in degrees) to use when compressing the data. Default:
0.25°.

### --cnt cint
cint is the contour interval to use when compressing. Default: 200. If instead
of a number, the name of a GMT colour palette is provided, then the contour
levels are read from the colour palette file. An example is the file
[rain.cpt](./rain.cpt):

```
1  255/255/102          2     177.24/223.62/104 L
2  177.24/223.62/104    5     118.88/169.88/117 L
5  118.88/169.88/117    10    83.5/133.5/127 L
10 83.5/133.5/127       20    59.119/109.12/145.88 L
20 59.119/109.12/145.88 50    39.381/77.381/165.62 L
50 39.381/77.381/165.62 100   26/51/179   B
B  white
F  26/51/179
N  128
```

The contour levels in this file are 1, 2, 5, 10, 20, 50 and 100.

### --reg region
region is the region to extract the data for. The regions are defined
[here](./regions.md). Some of the regions do not cover rectangular portions of
a longitude-latitude grid. In this case the script will download enough data to
cover the entire non-rectangular region.

## Examples

```
julia GetGFS.jl -t 2023012412 -v prmslmsl --tol 0.125 --cnt 200 --reg AUS
```

Collects mean sea level pressure data from the GFS run on 24 January 2023, 1200
UTC. The data are compressed using a tolerance of 0.125° and contour spacing of
200 Pa. The script creates a series of files containing the raw NetCDF data,
and the compressed contour data for lead times 24, 48, 72, 96 and 120 hours:

```
GFS_000-003-001_MSL_2023012412_024.nc
GFS_000-003-001_MSL_2023012412_048.nc
GFS_000-003-001_MSL_2023012412_072.nc
GFS_000-003-001_MSL_2023012412_096.nc
GFS_000-003-001_MSL_2023012412_120.nc

GFS_AUS_000-003-001_MSL_2023012412_024.bin
GFS_AUS_000-003-001_MSL_2023012412_048.bin
GFS_AUS_000-003-001_MSL_2023012412_072.bin
GFS_AUS_000-003-001_MSL_2023012412_096.bin
GFS_AUS_000-003-001_MSL_2023012412_120.bin
```
