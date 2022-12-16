module GetMSLP

include("BinaryContour.jl")
include("GetGFS.jl")

using Dates
using GMT
using Printf
using .BinaryContour
using .GetGFS
using ArgParse

export get_data, contour_data

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
)
    """
    Contour a MSLP field, then simplify the contours and write the simplified
    contours to a compressed binary file.

    Arguments:
    	mslp_grd:		GMT grid containting MSLP
    	mslp_header:	MSLP header information
    	cint:				Contour interval (hPa)
    	tol:				Contour tolerance (degrees)

    Returns: Nothing (data is written to file)
    """

    grdcontour(mslp_grd, cont = cint, dump = "mslpcnt.gmt")
    mslpcnt = gmtread("mslpcnt.gmt", table = true)
    smslp = gmtsimplify(mslpcnt, tol = tol)
    contour_to_bin(
        smslp,
        mslp_header,
        @sprintf(
            "mslp_t%03dc%03d_%s_%03d.bin",
            tol * 100,
            cint * 100,
            Dates.format(unix2datetime(mslp_header.base_time), "yyyymmddHH"),
            mslp_header.lead_time
        ),
        zval = NaN32,
    )
end

end
