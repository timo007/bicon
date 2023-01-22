include("BinaryContour.jl")

using ArgParse
using .BinaryContour
using Dates
using Downloads
using Format
using GMT
using Printf


function make_plot(grid, contour, header::ContourHeader, region, cpt, outfile;)
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
    gmtunit = replace(unit, r"([-\d]+)" => s"@+\1@+", "_" => " ")
    title = @sprintf(
        "\"%s %s (%s)\\072 %s valid at %s\"",
        strip(header.level, ' '),
        lowercasefirst(var),
        gmtunit,
        fcst_type,
        valid_time
    )

    #
    # Make the colour palette.
    #
    cpt = grd2cpt(grid, cmap = cpt, bg = :i, continuous = true, nlevels = true)

    #
    # Plot the data on a map.
    #
    grdimage(
        grid,
        color = cpt,
        proj = map_params(region)[:proj],
        region = map_params(region)[:mapRegion],
        frame = (map_params(region)[:frame]..., (title = title)),
        par = (
            FONT_TITLE = "14,AvantGarde-Book,black",
            FONT_ANNOT = "12,AvantGarde-Book,black",
            MAP_TITLE_OFFSET = "-6p",
            MAP_FRAME_TYPE = "plain",
            MAP_GRID_PEN_PRIMARY = "thinnest,158",
        ),
        figsize = 20,
    )
    coast!(area = (0, 0, 1), shore = "thinnest,purple")
    plot!(
        contour,
        pen = "thin, black",
        decorated = (
            quoted = true,
            n_labels = 1,
            label = :header,
            font = (8, "AvantGarde-Book"),
            n_data = 4,
        ),
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
        default = "https://github.com/timo007/bincon/blob/main/test_data/GFS_NZ_000-003-001_MSL_2023010412_024.bin?raw=true"
        "--cpt"
        help = "Colour palette"
        default = "batlow"
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
        println(@sprintf("Local file %s: %d bytes", parsed_args["i"], filesize(infile)))
    else
        try
            Downloads.download(parsed_args["i"], "./data.bin")
        catch
            println("Unable to download the data file.")
            exit()
        end
        println(
            @sprintf("Downloaded %s: %d bytes", parsed_args["i"], filesize("./data.bin"))
        )
        infile = "./data.bin"
    end

    #
    # Read the contours from the file, if the input is a contour file.
    #
    contour, header = bin_to_contour(infile)

    #
    # Grid the MSLP.
    #
    grid = contour_to_grid(
        contour,
        parsed_args["inc"],
        (header.west, header.east, header.south, header.north),
    )

    #
    # Make the plot.
    #
    make_plot(
        grid,
        contour,
        header::ContourHeader,
        reg,
        parsed_args["cpt"],
        parsed_args["o"],
    )

end

main()
