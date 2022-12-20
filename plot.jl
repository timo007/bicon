include("BinaryContour.jl")

using ArgParse
using .BinaryContour
using Dates
using Downloads
using Format
using GMT
using Printf

function map_params(region_name::Symbol)
    """
    Convert a region name (e.g. NZ) to GMT map projection parameters.
    """
    proj = Dict(
        :NZ => Dict(
            :dataRegion => (140.0f0, 200.0f0, -55.0f0, -25.0f0),
            :proj =>
                (name = :lambertConic, center = [170, -40], parallels = [-35, -45]),
            :mapRegion => "142/-52/-170/-28+r",
            :frame => (axes = :WSen, ticks = 1, grid = 10, annot = 10),
        ),
        :SWP => Dict(
            :dataRegion => (150.0f0, 240.0f0, -35.0f0, 0.0f0),
            :proj => (name = :Mercator, center = [175, 0]),
            :mapRegion => "150/240/-35/0",
            :frame => (axes = :wsen, ticks = 360, grid = 360),
        ),
        :UK => Dict(
            :dataRegion => (0.0f0, 360.0f0, 40.0f0, 70.0f0),
            :proj => (name = :conicEquidistant, center = [0, 50], parallels = [45, 55]),
            :mapRegion => "-30/40/15/65+r",
            :frame => (axes = :wsen, ticks = 360, grid = 360),
        ),
        :Russia => Dict(
            :dataRegion => (0.0f0, 200.0f0, 0.0f0, 90.0f0),
            :proj => (name = :conicEquidistant, center = [100, 65], parallels = [60, 70]),
            :mapRegion => "50/0/190/50+r",
            :frame => (axes = :wsen, ticks = 360, grid = 360),
        ),
        :World => Dict(
            :dataRegion => (0.0f0, 360.0f0, -90.0f0, 90.0f0),
            :proj => (name = :Robinson, center = 175),
            :mapRegion => "0/360/-90/90",
            :frame => (axes = :wsen, ticks = 360, grid = 360),
        ),
    )
    return proj[region_name]
end

function contours_to_grid(contours, inc, region)
    println("Converting contours to grid")
    mean_contour = blockmean(contours, inc = inc, region = region)
    grid = surface(mean_contour, inc = inc, region = region, tension = 0, A = "m")
    return grid
end

function make_plot(mslp_grid, header::ContourHeader, region, outfile)
    valid_time = Dates.format(
        unix2datetime(header.base_time) + Dates.Hour(header.lead_time),
        "HH:MMZ e d u YYYY",
    )
    if header.lead_time == 0
        fcst_type = "Analysis"
    else
        fcst_type = @sprintf("%dh forecast", header.lead_time)
    end
    var, unit = GRIBparam(header.discipline, header.category, header.parameter)
    title = @sprintf("\"%s\\072 %s valid at %s\"", var, fcst_type, valid_time)

    cpt = grd2cpt(mslp_grid, cmap = :batlow, bg = :i, continuous = true, nlevels = true)

    #
    # Plot the data on a map.
    #
    grdimage(
        mslp_grid,
        color = cpt,
        proj = map_params(region)[:proj],
        region = map_params(region)[:mapRegion],
        frame = (map_params(region)[:frame]..., (title = title)),
        par = (
            FONT_TITLE = "14,AvantGarde-Book,black",
            MAP_TITLE_OFFSET = "-6p",
            MAP_FRAME_TYPE = "plain",
            MAP_GRID_PEN_PRIMARY = "thinnest,158",
        ),
        figsize = 20,
    )
    coast!(area = (0, 0, 1), shore = "thinnest,white")
    grdcontour!(
        mslp_grid,
        annot = (int = 4, labels = (font = (8, "AvantGarde-Book"),)),
        cont = 2,
        pen = "thin, black",
        labels = (dist = 4,),
        savefig = outfile,
        show = true,
    )
end

function parse_commandline()
    s = ArgParseSettings()
    @add_arg_table! s begin
        "-i"
        arg_type = AbstractString
        help = "Name or URL of the data file"
        default = "http://nomuka.com/data/mslp_NZ_t025c200_2022121818_000.bin"
        "-f"
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
    # Download (if required) the data file.
    #
    if isfile(parsed_args["i"])
        infile = parsed_args["i"]
    else
        Downloads.download(parsed_args["i"], "./data.bin")
        println(
            @sprintf("Downloaded %s: %d bytes", parsed_args["i"], filesize("./data.bin"))
        )
        infile = "./data.bin"
    end

    #
    # Read the contours from the file.
    #
    mslp, mslp_header = bin_to_contour(infile)

    #
    # Grid the MSLP.
    #
    mslp_grid = contours_to_grid(
        mslp,
        parsed_args["inc"],
        (mslp_header.west, mslp_header.east, mslp_header.south, mslp_header.north),
    )

    #
    # Make the plot.
    #
    make_plot(mslp_grid, mslp_header::ContourHeader, reg, parsed_args["o"])

end

main()
