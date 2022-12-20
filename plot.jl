include("BinaryContour.jl")

using ArgParse
using .BinaryContour
using Dates
using Downloads
using Format
using GMT
using Printf


function make_plot(mslp_grid, header::ContourHeader, region, contint, outfile)
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
        annot = (int = contint*2, labels = (font = (8, "AvantGarde-Book"),)),
        cont = contint,
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
        "--cnt"
		  help = "Contour spacing"
		  arg_type = Float32
		  default = Float32(2)
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
    mslp_grid = contour_to_grid(
        mslp,
        parsed_args["inc"],
        (mslp_header.west, mslp_header.east, mslp_header.south, mslp_header.north),
    )

    #
    # Make the plot.
    #
	 make_plot(mslp_grid, mslp_header::ContourHeader, reg, parsed_args["cnt"], parsed_args["o"])

end

main()
