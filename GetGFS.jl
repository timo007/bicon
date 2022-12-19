module GetGFS

using NCDatasets
using Dates
using GMT
using Printf

export opendap_to_gmt, download_var

function download_var(
	 url::String,
	 var::String,
	 outdir::String;
    level::Float32 = NaN32,
)
	 # Collect the data from NCEP's OpenDAP server.
    ds = NCDataset(url, "r")

    lon = ds["lon"][:]
    lat = ds["lat"][:]
    valtime = ds["time"][:]

    if (!isnan(level))
        # Extract a slice of 3D data (e.g. a level from pressure level data)
        lev = ds["lev"][:]
        lev_idx = findall(x -> x ≈ level, lev)[1]
		  data = ds[var][:, :, lev_idx, :]
    else
        # 2D data (e.g. surface fields)
		  data = ds[var][:, :, :]
    end

    close(ds)

	 #
	 # Convert the data to GMT grids, and save to files.
	 #
	 for t = 1:length(valtime)
		 base_time = Dates.format(valtime[1], "yyyymmddHH")
		 lead_time = Dates.value(valtime[t] - valtime[1])/3600000
		 outfile=@sprintf("gfs_%s_%s_%03d.nc", var, base_time, lead_time)

        data_grid = mat2grid(
			   permutedims(nomissing(data[:, :, t], NaN), (2, 1),),
            x = lon[:],
            y = lat[:],
        )
		  println("Writing ", outfile)
		  gmtwrite(outfile, data_grid)
	  end
end


function opendap_to_gmt(
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
