include("BinaryContour.jl")
include("get_mslp.jl")

using ArgParse
using .BinaryContour
using .GetMSLP
using Format
using GMT
using Printf
using Statistics

function contours_to_grid(contours, inc, region)
    println("Converting contours to grid")
    mean_contour = blockmean(contours, inc = inc, region = region)
    grid = surface(mean_contour, inc = inc, region = region, tension = 0, A = "m")
    return grid
end

function make_plot(mslp_grid, header::ContourHeader, titlestr, cpt)
    #
    # Plot the data on a map.
    #
    grdimage(
        mslp_grid,
        color = cpt,
        proj = (name = :lambertConic, center = [170, -40], parallels = [-35, -45]),
        R = "142/-52/-170/-28+r",
        frame = (axes = :wsen, ticks = 360, grid = 360, title = titlestr),
        par = (FONT_TITLE = "10,AvantGarde-Book,black", MAP_TITLE_OFFSET = "-9p"),
    )
    coast!(shore = "thinnest,black", land = :darkgreen)
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

    #
    # Get the raw MSLP data from NCEP.
    #
    raw_grid, mslp_header = get_data(
        parsed_args["t"],
        parsed_args["f"],
        (Float32(140), Float32(200), Float32(-55), Float32(-25)),
    )
    raw_grid = raw_grid / 100# Convert Pa to hPa.

    #
    # Contour the data, and save to file.
    #
    contour_data(raw_grid, mslp_header, parsed_args["cnt"], parsed_args["tol"])

    #
    # Read the contours from file.
    #
    contour_file = @sprintf(
        "mslp_t%03dc%03d_%s_%03d.bin",
        parsed_args["tol"] * 100,
        parsed_args["cnt"] * 100,
        parsed_args["t"],
        parsed_args["f"]
    )

    mslp, mslp_header = bin_to_contour(contour_file)

    #
    # Grid the MSLP.
    #
    mslp_grid = contours_to_grid(
        mslp,
        parsed_args["inc"],
        (mslp_header.west, mslp_header.east, mslp_header.south, mslp_header.north),
    )
    raw_grid = grdedit(raw_grid, coltype = "g")
    #gmtwrite(mslp_contour, mslp1_grid)

    #
    # Make the plot.
    #
    mslp_cpt =
        makecpt(cmap = :roma, range = (980, 1030, 10), inverse = true, continuous = true)

    subplot(
        grid = "2x1",
        dims = (panels = (12, 8),),
        margins = 0.15,
        savefig = parsed_args["o"],
    )
    subplot(:set, panel = (1, 1))
    raw_title = @sprintf(
        "\"Gridded data as 32-bit floats: %s bytes\"",
        replace(format(length(raw_grid) * sizeof(Float32), commas = true), "," => " ")
    )
    make_plot(raw_grid, mslp_header::ContourHeader, raw_title, mslp_cpt)
    subplot(:set, panel = "next")
    cnt_title = @sprintf(
        "\"Contoured data: %s bytes\"",
        replace(format(filesize(contour_file), commas = true), "," => " ")
    )
    make_plot(mslp_grid, mslp_header::ContourHeader, cnt_title, mslp_cpt)
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
