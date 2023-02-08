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

function download_var(url::String, var::String; level::Number = NaN)
    """
    Collect data from NCEP's OpenDAP server.

    Arguments:
    url:    The URL of the dataset containing the data to collect.
    var:    The name of teh variable to collect.
    level:  The pressure level to collect (not required for 2D fields).

    Returns:
    Writes the downloaded data to files on the local computer, and returns
    a list of these file names.
    """
    ds = NCDataset(url, "r")

    lon = ds["lon"][:]
    lat = ds["lat"][:]
    basetime = ds["time"][1]
    valtime = ds["time"][9:4:41]

    if (!isnan(level))
        # Extract a slice of 3D data (e.g. a level from pressure level data)
        lev = ds["lev"][:]
        lev_idx = findall(x -> x ≈ level, lev)[1]
        data = ds[var][:, :, lev_idx, 9:8:41]
    else
        # 2D data (e.g. surface fields)
        data = ds[var][:, :, 9:8:41]
    end

    close(ds)

    #
    # Convert the data to GMT grids, and save to files.
    #
    file_list = String[]
    for t = 1:length(valtime)
        basetime_str = Dates.format(basetime, "yyyymmddHH")
        lead_time = Dates.value(valtime[t] - basetime) / 3600000 # ms to hours.
        GRIBparam = NCEPvar_to_GRIBparam(var)
        if isnan(level)
            vlevstr = GRIBparam[4]
        else
            vlevstr = @sprintf("%dPa", round(level * 100))
        end
        outfile = @sprintf(
            "GFS_%03d-%03d-%03d_%s_%s_%03d.nc",
            GRIBparam[1],
            GRIBparam[2],
            GRIBparam[3],
            strip(vlevstr, ' '),
            basetime_str,
            lead_time
        )
        push!(file_list, outfile)

        data_grid = permutedims(nomissing(data[:, :, t], NaN), (2, 1))
        data_grid = hcat(data_grid, data_grid[:, 1])

        data_grid = mat2grid(data_grid, x = vcat(lon, 360), y = lat)
        println("Writing ", outfile)
        gmtwrite(outfile, data_grid)
    end
    return file_list
end


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
    """
    Read a variable from the NCEP OpenDAP server (nomads.ncep.noaa.gov)
    and convert it to a GMT grid type. This defaults to reading global data.
    If no level is specified as an argument, read a 2D field, otherwise read
    a 3D field at the level specified by the level argument.

    Arguments:

    url:    The OpenDAP URL
    var:    The variable to read (e.g. ugrdprs)
    south:  The southern limit of the domain to read.
    north:  The northern limit of the domain to read.
    west:   The western limit of the domain to read.
    east:   The eastern limit of the domain to read.
    level:  The vertical level to read (e.g. 925 for 925 hPa). If ommitted, read a 2D field.
    fcst:   The forecast lead time (in units of time since the the first/base time. e.g. 24)
    """
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


function parse_commandline()
    s = ArgParseSettings()
    @add_arg_table! s begin
        "-t"
        help = "NWP base time (yyyymmddHH)"
        "-v"
        help = "NCEP variable to collect"
        default = "prmslmsl"
        "-p"
        help = "Vertical level (when required)"
        arg_type = Float32
        default = NaN32
        "--tol"
        help = "Tolerance"
        arg_type = Float32
        default = Float32(0.25)
        "--cnt"
        help = "Contour interval"
        default = 200
        "--reg"
        help = "Region to plot"
        default = "NZ"
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
        parsed_args["t"][1:8],
        "/gfs_0p25_",
        parsed_args["t"][9:10],
        "z",
    )
    file_list = download_var(ncep_url, parsed_args["v"], level = parsed_args["p"])

    #
    # Compress the NetCDF files which have been downloaded.
    #
    for file in file_list
        # The name of the compressed file is derived from the NetCDF file.
        cntfile = replace(file, "GFS" => "GFS_" * parsed_args["reg"], ".nc" => ".bin")

        GRIBparam = NCEPvar_to_GRIBparam(parsed_args["v"])
        if isnan(parsed_args["p"])
            vlevstr = GRIBparam[4]
        else
            vlevstr = @sprintf("%dPa", round(parsed_args["p"] * 100))
        end
        # Extract the lead time from the file name.
        fcst = match(r"^.*_(\d{3}).nc", file)[1]

        map_region = data_region(reg)
        if map_region[1] < 0 || map_region[2] > 359.75
            # Deal with regions which cross the east/west border of the global
            # NWP data (e.g. :UK).
            grid = grdedit(file, region = (-180, 180, -90, 90), wrap = true, f = "ig")
            grid = grdcut(grid, region = map_region)
        else
            grid = gmtread(file, grid = true, region = map_region)
        end

        header = ContourHeader(
            GRIBparam[1],
            GRIBparam[2],
            GRIBparam[3],
            datetime2unix(DateTime(parsed_args["t"], dateformat"yyyymmddHH")),
            parse(Float32, fcst),
            lpad(vlevstr, 8, ' '),
            data_region(reg)[1],
            data_region(reg)[2],
            data_region(reg)[3],
            data_region(reg)[4],
        )
        println("Contouring ", cntfile)
        grid_to_contour(grid, header, parsed_args["cnt"], parsed_args["tol"], cntfile)

    end

end

main()

end
