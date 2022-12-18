include("BinaryContour.jl")
include("GetGFS.jl")

using ArgParse
using .BinaryContour
using Dates
using .GetGFS
using Format
using GMT
using Printf
using Statistics

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
            :mapRegion => "155/220/-35/-5",
            :dataRegion => (155.0f0, 220.0f0, -35.0f0, -5.0f0),
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

function get_data(
    dtstr::String,
    fcst::Float32,
    region::Tuple{Float32,Float32,Float32,Float32},
)
    """
    Collect GFS MSLP data from NCEP NOMADS server.

    Arguments:
    	dtstr:		Date and time string (YYYYMMDDHH)
    	fcst:			Forecast lead time (hours)
    	region:		Region tuple (west, east, south, north)

    Returns:			GMTgrid with MSLP and data header struct.
    """
    nwp_time = DateTime(dtstr, dateformat"yyyymmddHH")
    datedir = "gfs." * dtstr[1:8]
    hourdir = dtstr[9:10]

    println(@sprintf("Collecting +%03d hour forecast MSLP data", fcst))
    mslp_header = ContourHeader(
        UInt8(0),
        UInt8(0),
        UInt8(4),
        datetime2unix(nwp_time),
        fcst,
        NaN32,
        region[1],
        region[2],
        region[3],
        region[4],
    )

    infile = string(
        "http://nomads.ncep.noaa.gov:80/dods/gfs_0p25/gfs",
        Dates.format(nwp_time, "yyyymmdd"),
        "/gfs_0p25_",
        Dates.format(nwp_time, "HH"),
        "z",
    )

    #
    # Collect the MSLP data.
    #
    mslpgrd = opendap_to_gmt(
        infile,
        "prmslmsl",
        south = region[3],
        north = region[4],
        west = region[1],
        east = region[2],
        fcst = fcst,
    )

    return mslpgrd, mslp_header
end

function contour_data(
    mslp_grd::GMTgrid,
    mslp_header::ContourHeader,
    cint::Float32,
    tol::Float32,
    outfile::String,
)
    """
    Contour a MSLP field, then simplify the contours and write the simplified
    contours to a compressed binary file.

    Arguments:
    	mslp_grd:		GMT grid containting MSLP
    	mslp_header:	MSLP header information
    	cint:				Contour interval (hPa)
    	tol:				Contour tolerance (degrees)
  outfile:			Name of binary data file to write.

    Returns: Nothing (data is written to file)
    """

    grdcontour(mslp_grd, cont = cint, dump = "mslpcnt.gmt")
    mslpcnt = gmtread("mslpcnt.gmt", table = true)
    smslp = gmtsimplify(mslpcnt, tol = tol)
    contour_to_bin(smslp, mslp_header, outfile, zval = NaN32)
end

function contours_to_grid(contours, inc, region)
    println("Converting contours to grid")
    mean_contour = blockmean(contours, inc = inc, region = region)
    grid = surface(mean_contour, inc = inc, region = region, tension = 0, A = "m")
    return grid
end

function make_plot(mslp_grid, header::ContourHeader, region, titlestr, cpt)
    #
    # Plot the data on a map.
    #
    grdimage(
        mslp_grid,
        color = cpt,
        proj = map_params(region)[:proj],
        region = map_params(region)[:mapRegion],
        frame = (axes = :wsen, ticks = 360, grid = 360, title = titlestr),
        par = (
            FONT_TITLE = "10,AvantGarde-Book,black",
            MAP_TITLE_OFFSET = "-9p",
            MAP_FRAME_TYPE = "plain",
        ),
    )
    coast!(area = (0, 0, 1), res=:low, shore = "thinnest,white")
    grdcontour!(
        mslp_grid,
        annot = (int = 4, labels = (font = (6, "AvantGarde-Book"),)),
        cont = 2,
        pen = "thin, black",
        labels = (dist = 4,),
    )
end

function parse_commandline()
    s = ArgParseSettings()
    @add_arg_table! s begin
        "-t"
        help = "NWP base time (yyyymmddHH)"
        default = "mslp.bin"
        "-f"
        help = "Forecast lead time"
        arg_type = Float32
        default = convert(Float32, 72)
        "--tol"
        help = "Tolerance"
        arg_type = Float32
        default = convert(Float32, 0.25)
        "--cnt"
        help = "Contour interval"
        arg_type = Float32
        default = convert(Float32, 2)
        "--inc"
        help = "MSLP grid spacing"
        default = "15m/15m"
        "--reg"
        help = "Region to plot"
        default = :NZ
        "-o"
        help = "Name of map file to produce"
        default = "mslp.png"
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
    # Get the raw MSLP data from NCEP.
    #
    raw_grid, mslp_header = get_data(
        parsed_args["t"],
        parsed_args["f"],
        map_params(reg)[:dataRegion],
    )
    raw_grid = raw_grid / 100	# Convert Pa to hPa.
	 #gmtwrite("raw.nc", raw_grid)

    #
    # Contour the data, and save to file.
    #
    outfile = @sprintf(
        "mslp_%s_t%03dc%03d_%s_%03d.bin",
		  parsed_args["reg"],
        parsed_args["tol"] * 100,
        parsed_args["cnt"] * 100,
        parsed_args["t"],
        mslp_header.lead_time
    )
    contour_data(raw_grid, mslp_header, parsed_args["cnt"], parsed_args["tol"], outfile)

    #
    # Read the contours straight back from the file.
    #
    mslp, mslp_header = bin_to_contour(outfile)

    #
    # Grid the MSLP.
    #
    mslp_grid = contours_to_grid(
        mslp,
        parsed_args["inc"],
        (mslp_header.west, mslp_header.east, mslp_header.south, mslp_header.north),
    )
    raw_grid = grdedit(raw_grid, coltype = "g")

    #
    # Make the plot.
    #
    valid_time = Dates.format(
        DateTime(parsed_args["t"], dateformat"yyyymmddHH") + Dates.Hour(parsed_args["f"]),
        "HH:MMZ e d u YYYY",
    )
	 mslp_cpt = grd2cpt(raw_grid, cmap = :batlow, bg=:i,
							  continuous=true, nlevels=true)

    subplot(
        grid = "2x1",
        panels_size = 12,
        region = map_params(reg)[:mapRegion],
        proj = (map_params(reg)[:proj]),
        savefig = parsed_args["o"],
        margin = -0.4,
        title = @sprintf(
            "MSLP +%dh forecast valid at %s",
            round(parsed_args["f"]),
            valid_time
        ),
        par = (FONT_HEADING = "12,AvantGarde-Demi,black",),
    )
    subplot(:set, panel = (1, 1))
    raw_title = @sprintf(
        "\"Gridded data as 32-bit floats: %s bytes\"",
        replace(format(length(raw_grid) * sizeof(Float32), commas = true), "," => " ")
    )
    make_plot(
        raw_grid,
        mslp_header::ContourHeader,
        reg,
        raw_title,
        mslp_cpt,
    )
    subplot(:set, panel = "next")
    cnt_title = @sprintf(
        "\"Contoured data: %s bytes\"",
        replace(format(filesize(contour_file), commas = true), "," => " ")
    )
    make_plot(
        mslp_grid,
        mslp_header::ContourHeader,
        reg,
        cnt_title,
        mslp_cpt,
    )
    subplot(:show)

    #
    # Print some statistics.
    #
    error = mslp_grid - raw_grid
    bias = mean(error)
    rmse = sqrt(mean(error .^ 2))
    println("Bias: ", bias)
    println("RMSE: ", rmse)

end

main()
