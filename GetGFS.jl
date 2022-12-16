module GetGFS

using NCDatasets
using Dates
using GMT

export opendap_to_gmt

"""
Read a variable from the NCEP OpenDAP server (nomads.ncep.noaa.gov)
and convert it to a GMT grid type. This defaults to reading global data.
If no level is specified as an argument, read a 2D field, otherwise read
a 3D field at the level specified by the level argument.

Arguments:

url:		The OpenDAP URL
var:		The variable to read (e.g. ugrdprs)
south: 	The southern limit of the domain to read.
north:	The northern limit of the domain to read.
west:		The western limit of the domain to read.
east:		The eastern limit of the domain to read.
level:	The vertical level to read (e.g. 925 for 925 hPa). If ommitted, read a 2D field.
fcst:		The forecast lead time (in units of time since the the first/base time. e.g. 24)
"""
function opendap_to_gmt(
    url::String,
    var::String;
    south::Float32 = -90.0,
    north::Float32 = 90.0,
    west::Float32 = 0.0,
    east::Float32 = 360.0,
    level::Float32 = NaN32,
    fcst::Float32 = 0,
)
    ds = NCDataset(url, "r")

    lon = ds["lon"][:]
    lat = ds["lat"][:]
    time = ds["time"][:]
    lead_time = time - time[1]

    time_idx = findall(x -> x == Dates.Hour(fcst), lead_time)[1]
    lon_idx = findall(x -> x >= west && x <= east, lon)
    lat_idx = findall(x -> x >= south && x <= north, lat)

    if (!isnan(level))
        # Extract a slice of 3D data (e.g. a level from pressure level data)
        lev = ds["lev"][:]
        lev_idx = findall(x -> x ≈ level, lev)[1]
        data_grid = mat2grid(
            permutedims(
                nomissing(ds[var][lon_idx, lat_idx, lev_idx, time_idx], NaN),
                (2, 1),
            ),
            x = lon[lon_idx],
            y = lat[lat_idx],
        )
    else
        # 2D data (e.g. surface fields)
        data_grid = mat2grid(
            permutedims(nomissing(ds[var][lon_idx, lat_idx, time_idx], NaN), (2, 1)),
            x = lon[lon_idx],
            y = lat[lat_idx],
        )
    end

    close(ds)

    return data_grid
end

end
