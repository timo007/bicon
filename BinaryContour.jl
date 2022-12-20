module BinaryContour

using CSV
using GMT
using Dates
using Printf

export contour_to_bin,
    bin_to_contour, ContourHeader, GRIBparam, NCEPvar_to_GRIBparam,
	grid_to_contour, contour_to_grid, map_params

struct ContourHeader
    discipline::UInt8
    category::UInt8
    parameter::UInt8
    base_time::Float64
    lead_time::Float32
    prs::Float32
    west::Float32
    east::Float32
    south::Float32
    north::Float32
end


function contour_to_bin(
    contour::Vector{GMTdataset{Float64,2}},
    header::ContourHeader,
    outfile::String;
    zval::Float32 = 1,
)
	"""
	Convert contours (in GMT data sets) to binary encoded contours.
	"""
    #
    # Write the contours in binary, network endianess.
    #
    open(outfile, "w") do file
        write(file, hton(header.discipline))
        write(file, hton(header.category))
        write(file, hton(header.parameter))
        write(file, hton(header.base_time))
        write(file, hton(header.lead_time))
        write(file, hton(header.prs))
        write(file, hton(header.west))
        write(file, hton(header.east))
        write(file, hton(header.south))
        write(file, hton(header.north))

        for segment in contour
            if isnan(zval)
                clev = segment.data[1, 3]
            else
                clev = zval
            end
            nrows = size(segment.data, 1)
            inc =
                maximum(abs.(segment.data[2:nrows, 1:2] - segment.data[1:nrows-1, 1:2])) /
                127
            write(file, hton(convert(Float32, clev)))
            write(file, hton(convert(Int16, nrows)))
            write(file, hton(convert(Float32, inc)))
            write(file, hton.(convert(Array{Float32}, segment.data[1, 1:2])))
            locs =
                round.(
                    Int8,
                    ((segment.data[2:nrows, 1:2] - segment.data[1:nrows-1, 1:2]) / inc),
                )
            write(file, hton.(locs))
        end
    end
end

function bin_to_contour(infile::String)
    open(infile, "r") do file
        var = Array{UInt8}(undef, 1, 3)
        data = Array{Float32}(undef, 1, 6)
        read!(file, var)
        bt = ntoh(read(file, Float64))
        read!(file, data)
        data = ntoh.(data)
        header = ContourHeader(
            var[1],
            var[2],
            var[3],
            bt,
            data[1],
            data[2],
            data[3],
            data[4],
            data[5],
            data[6],
        )

        gmtds = Vector{GMTdataset{Float32,2}}(undef, 0)
        while !eof(file)
            zval = ntoh(read(file, Float32))
            npts = ntoh(read(file, UInt16))
            inc = ntoh(read(file, Float32))
            locs = Array{Float32}(undef, npts, 3)
            locs[1, 1] = ntoh(read(file, Float32))
            locs[1, 2] = ntoh(read(file, Float32))
            locs[1, 3] = zval
            offsets = Array{Int8}(undef, npts - 1, 2)
            read!(file, offsets)
            for row = 2:npts
                locs[row, 1:2] = locs[row-1, 1:2] + offsets[row-1, 1:2] * inc
                locs[row, 3] = zval
            end
            push!(gmtds, mat2ds(locs, hdr = @sprintf("%g", zval)))
        end
        return gmtds, header
    end
end

function GRIBparam(discipline::Integer, category::Integer, parameter::Integer)
    """
    Read the parameter description (name) and units from the WMO GRIB-2 parameter
    file (downloaded from WMO in CSV format).

	 Data are available here: http://codes.wmo.int/grib2/codeflag/4.2?_format=csv

    Arguments:
    	discipline:	Discipline - GRIB octet ...
    	category:   Category - GRIB octet ...
    	parameter:  Parameter - GRIB octet ..."

    Returns:
    """
    name = "Unknown"
    unit = "Unknown"
    param = @sprintf("%d-%d-%d", discipline, category, parameter)
    data = CSV.File("4.2.csv")
    for row in data
        if row[8] == param
            name = match(r"\'(.*)\'", row[10])[1]
            unit = match(r"unit\/(.*)>", row[2])[1]
        end
    end
    return name, unit
end

function NCEPvar_to_GRIBparam(ncep_var::String)
    """"
    Convert NCEP parameter names to a tuple: (discipline, category, parameter, level)
    as used in GRIB-2. This table needs serious work to be a bit more complete.
    """
    ncep_table = Dict(
		 "tmax2m"   => (0, 0, 4, "2m"),
		 "tmin2m"   => (0, 0, 5, "2m"),
		 "tmpsfc"   => (0, 0, 0, "sfc"),
		 "tmp2m"    => (0, 0, 0, "2m"),
		 "prmslmsl" => (0, 3, 1, "msl"),
	 )

    if haskey(ncep_table, ncep_var)
        GRIB_var = ncep_table[ncep_var]
    else
        GRIB_var = (255, 255, 255)
    end

    return GRIB_var
end

function grid_to_contour(
    grid::GMTgrid,
    header::ContourHeader,
    cint::Float32,
    tol::Float32,
    cntfile::String,
)
    """
    Contour a field, then simplify the contours and write the simplified
    contours to a binary contour file.

    Arguments:
    	grid:		GMT grid containting data to be contoured.
    	header:	Header information
    	cint:		Contour interval
    	tol:		Contour tolerance (degrees)
      cntfile:	Name of binary contour file to write.

    Returns: Nothing (data is written to file)
    """
    grdcontour(grid, cont = cint, dump = "mslpcnt.gmt")
    contour_data = gmtread("mslpcnt.gmt", table = true)
    simplified_contour = gmtsimplify(contour_data, tol = tol)
    contour_to_bin(simplified_contour, header, cntfile, zval = NaN32)
end

function contour_to_grid(contour, inc, region)
    println("Converting contours to grid")
    mean_contour = blockmean(contour, inc = inc, region = region)
    grid = surface(mean_contour, inc = inc, region = region, tension = 0, A = "m")
    return grid
end

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
            :dataRegion => (0f0, 360f0, 40f0, 70f0),
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

end
