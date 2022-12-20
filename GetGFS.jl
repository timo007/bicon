module GetGFS

include("BinaryContour.jl")

using ArgParse
using .BinaryContour
using NCDatasets
using Dates
using Format
using GMT
using Printf

export opendap_to_gmt, download_var

function download_var(
	 url::String,
	 var::String,
    level::Float32 = NaN32,
)
	 # Collect the data from NCEP's OpenDAP server.
    ds = NCDataset(url, "r")

    lon = ds["lon"][:]
    lat = ds["lat"][:]
    valtime = ds["time"][1:41]

    if (!isnan(level))
        # Extract a slice of 3D data (e.g. a level from pressure level data)
        lev = ds["lev"][:]
        lev_idx = findall(x -> x ≈ level, lev)[1]
		  data = ds[var][:, :, lev_idx, 1:41]
		  levstr = replace(@sprintf("%09.3f", level), "." => "p")
    else
        # 2D data (e.g. surface fields)
		  data = ds[var][:, :, 1:41]
		  levstr = "2dsurface"
    end

    close(ds)

	 #
	 # Convert the data to GMT grids, and save to files.
	 #
	 file_list = String[]
	 for t = 1:length(valtime)
		 base_time = Dates.format(valtime[1], "yyyymmddHH")
		 lead_time = Dates.value(valtime[t] - valtime[1])/3600000 # ms to hours.
		 GRIBparam = NCEPvar_to_GRIBparam(var)
		 outfile=@sprintf("gfs_%03d-%03d-%03d_%s_%s_%03d.nc", GRIBparam[1], GRIBparam[2], GRIBparam[3], levstr, base_time, lead_time)
		 push!(file_list, outfile)

        data_grid = mat2grid(
			   permutedims(nomissing(data[:, :, t], NaN), (2, 1),),
            x = lon[:],
            y = lat[:],
        )
		  println("Writing ", outfile)
		  gmtwrite(outfile, data_grid)
	  end
	  return file_list
end


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


function map_params(region_name::Symbol)
    """
    Convert a region name (e.g. NZ) to GMT map projection parameters.
    """
    proj = Dict(
        :NZ => Dict(
            :proj =>
                (name = :lambertConic, center = [170, -40], parallels = [-35, -45]),
            :mapRegion => "142/-52/-170/-28+r",
            :dataRegion => (140.0f0, 200.0f0, -55.0f0, -25.0f0),
        ),
        :SWP => Dict(
            :proj => (name = :Mercator, center = [175, 0]),
            :mapRegion => "150/240/-35/0",
            :dataRegion => (150.0f0, 240.0f0, -35.0f0, 0.0f0),
        ),
        :UK => Dict(
            :proj => (name = :conicEquidistant, center = [0, 50], parallels = [45, 55]),
            :mapRegion => "-30/40/15/65+r",
            :dataRegion => (0f0, 360f0, 40f0, 70f0),
        ),
        :Russia => Dict(
            :proj => (name = :conicEquidistant, center = [100, 65], parallels = [60, 70]),
            :mapRegion => "50/0/190/50+r",
            :dataRegion => (0.0f0, 200.0f0, 0.0f0, 90.0f0),
        ),
        :World => Dict(
            :proj => (name = :Robinson, center = 175),
            :mapRegion => "0/360/-90/90",
            :dataRegion => (0.0f0, 360.0f0, -90.0f0, 90.0f0),
        ),
    )
    return proj[region_name]
end

function parse_commandline()
    s = ArgParseSettings()
    @add_arg_table! s begin
        "-t"
        help = "NWP base time (yyyymmddHH)"
        default = "mslp.bin"
		  "-v"
		  help = "NCEP variable to collect"
		  default = "prmslmsl"
		  "-p"
		  help = "Vertical level (when required)"
        arg_type = Float32
		  default = NaN32
		  "-s"
		  help = "Scale factor (when required)"
        arg_type = Float32
		  default = 1
        "--tol"
        help = "Tolerance"
        arg_type = Float32
        default = convert(Float32, 0.25)
        "--cnt"
        help = "Contour interval"
        arg_type = Float32
        default = convert(Float32, 2)
        "--reg"
        help = "Region to plot"
        default = :NZ
    end

    return parse_args(s)
end

function main()
    #
    # Read the command line.
    #
    parsed_args = parse_commandline()
	 reg = Symbol(parsed_args["reg"])

	 #
	 # Collect data from NCEP.
	 #
	 ncep_url = string(
        "http://nomads.ncep.noaa.gov:80/dods/gfs_0p25/gfs",
		  parsed_args["t"][1:8], "/gfs_0p25_", parsed_args["t"][9:10], "z",
    )
	 file_list = download_var(ncep_url, parsed_args["v"],)

	 #
	 # Process the files listed on the command line
	 #
	 for file in file_list
		 # Try and extract variable, base time and lead time from the 
		 # file name.
		 cntfile = replace(file, ".nc" => ".bin")
		 GRIBparam = NCEPvar_to_GRIBparam(var)
		 fcst = match(r"^gfs_\d{3}-\d{3}-\d{3}_.{9}_\d{10}_(\d{3}).nc", file)[1]
		 
		 grid = gmtread(file, grid = true, region = map_params(reg)[:dataRegion])
		 header = ContourHeader(
										GRIBparam[1],
										GRIBparam[2],
										GRIBparam[3],
										datetime2unix(DateTime(parsed_args["t"], dateformat"yyyymmddHH")),
										parse(Float32, fcst),
										parsed_args["p"],
										map_params(reg)[:dataRegion][1],
										map_params(reg)[:dataRegion][2],
										map_params(reg)[:dataRegion][3],
										map_params(reg)[:dataRegion][4],
		)
		println("Contouring ", file)
		grid = grid * parsed_args["s"]
      grid_to_contour(grid, header, parsed_args["cnt"], parsed_args["tol"], cntfile)

	 end

end

main()

end
