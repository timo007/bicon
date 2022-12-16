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

function get_data(dtstr, fcst)

    west::Float32 = 140
    east::Float32 = 200
    south::Float32 = -55
    north::Float32 = -25

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
        west,
        east,
        south,
        north,
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
        south = south,
        north = north,
        west = west,
        east = east,
        fcst = convert(Float32, fcst),
    )

    #
    # Write the raw data to file.
    #
    dtstr = @sprintf("%s_%03d", Dates.format(nwp_time, "yyyymmddHH"), fcst)
    raw_grd = grdedit(mslpgrd, coltype = "g")
    gmtwrite("mslpraw_" * dtstr * ".nc", grdedit(mslpgrd, coltype = "g") / 100)

    return mslpgrd, mslp_header

end

function contour_data(mslp_grd, mslp_header::ContourHeader, cint, tol)
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
        zval = nothing,
    )
end

end
