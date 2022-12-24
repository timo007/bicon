module BinaryContour

using CSV
using GMT
using Dates
using Printf

export contour_to_bin,
    bin_to_contour,
    ContourHeader,
    GRIBparam,
    NCEPvar_to_GRIBparam,
    grid_to_contour,
    contour_to_grid,
    map_params,
	 data_region

struct ContourHeader
    discipline::UInt8
    category::UInt8
    parameter::UInt8
    base_time::Float64
    lead_time::Float32
    level::String
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
        write(file, header.level)
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
        var = Vector{UInt8}(undef, 3)
        level = Vector{UInt8}(undef, 8)
        region = Vector{Float32}(undef, 4)
        read!(file, var)
        bt = ntoh(read(file, Float64))
        lt = ntoh(read(file, Float32))
        read!(file, level)
        read!(file, region)
        region = ntoh.(region)
        header = ContourHeader(
            var[1],
            var[2],
            var[3],
            bt,
            lt,
            String(level),
            region[1],
            region[2],
            region[3],
            region[4],
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


function grid_to_contour(
    grid::GMTgrid,
    header::ContourHeader,
    cint::Union{Number,String},
    tol::Float32,
    cntfile::String,
)
    """
    Contour a field, then simplify the contours and write the simplified
    contours to a binary contour file.

    Contours are initially written to a temporary file, which is deleted
    after it has served its purpose.

    Arguments:
    	grid:		GMT grid containting data to be contoured.
    	header:	Header information
    	cint:		Contour interval or the name of a colour palette file to get contours from.
    	tol:		Contour tolerance (degrees)
      cntfile:	Name of binary contour file to write.

    Returns: Nothing (data is written to file)
    """
    contour_file = tempname()
    grdcontour(grid, cont = cint, dump = contour_file)
    simplified_contour = gmtsimplify(contour_file, tol = tol)
    rm(contour_file)
    contour_to_bin(simplified_contour, header, cntfile, zval = NaN32)
end

function contour_to_grid(
    contour::Vector{GMTdataset{Float32,2}},
    inc::String,
    region::NTuple{4,Float32},
)
    """
    Convert contour lines to a regular grid.

    Arguments:
    	contour:		The contour lines in GMTdataset format.
    	inc:			The grid resolution (e.g. 15m/15m for 0.25°x0.25°)
    	region:		The region (west, east, south, north)

    Returns:
    	A GMTgrid containing the data on a regular grid.
    """
    println("Converting contours to grid")
    mean_contour = blockmean(contour, inc = inc, region = region, center = true)
    grid = surface(mean_contour, inc = inc, region = region, tension = 0, A = "m")
    return grid
end

function data_region(region_name::Symbol)
	"""
	Work out what rectangular region of data we need to cut out to cover
	the domain covered by the map.

	Arguments:
		region_name:	Symbol representing the region (:NZ etc)
	"""
	region_ds = mapproject(region = map_params(region_name)[:mapRegion],
									proj = map_params(region_name)[:proj],
									map_size="r",)
	west = floor(region_ds.data[1])
	east = ceil(region_ds.data[2])
	south = floor(region_ds.data[3])
	north = ceil(region_ds.data[4])
	south = max(-90, south)
	north = min(90, north)
	return (west, east, south, north)
end

function map_params(region_name::Symbol)
    """
    Convert a region name (e.g. NZ) to GMT map projection parameters.

  Argument:
  	region_name		A symbol representing the region.
  					:NZ = New Zealand
  					:SWP = South West Pacific
  					:TON = Tonga
  					:AUS = Australia
  					:UK = United Kingdom
  					:WORLD = World

  Return:
  	Returns a dictionary with the various map related parameters for the
  	region of interest. More information on these parameters can be found
  	in the GMT documentation.

    """
    proj = Dict(
        :NZ => Dict(
            :proj =>
                (name = :lambertConic, center = [170, -40], parallels = [-35, -45]),
            :mapRegion => "142/-52/-170/-28+r",
            :frame => (axes = :WSen, ticks = 1, grid = 10, annot = 10),
        ),
        :SWP => Dict(
            :proj => (name = :Mercator, center = [175, 0]),
            :mapRegion => "150/240/-35/0",
            :frame => (axes = :WSen, ticks = 2, grid = 10, annot = 10),
        ),
        :TON => Dict(
            :dataRegion => (175.0f0, 195.0f0, -30.0f0, -10.0f0),
            :proj => (name = :Mercator, center = [182, 0]),
            :mapRegion => "175/195/-30/-10",
            :frame => (axes = :WSen, ticks = 1, grid = 5, annot = 5),
        ),
        :AUS => Dict(
            :proj =>
                (name = :lambertConic, center = [130, -30], parallels = [-20, -40]),
            :mapRegion => "80/-50/165/-3+r",
            :frame => (axes = :WSen, ticks = 2, grid = 10, annot = 10),
        ),
        :UK => Dict(
            :proj => (name = :conicEquidistant, center = [0, 50], parallels = [45, 55]),
            :mapRegion => "-30/40/15/65+r",
            :frame => (axes = :WSen, ticks = 1, grid = 10, annot = 10),
        ),
        :WORLD => Dict(
            :proj => (name = :Robinson, center = 175),
            :mapRegion => "0/360/-90/90",
            :frame => (axes = :wsen, ticks = 360, grid = 360),
        ),
    )
    return proj[region_name]
end

function NCEPvar_to_GRIBparam(ncep_var::String)
    """"
    Convert NCEP parameter names to a tuple: (discipline, category, parameter, level)
    as used in GRIB-2. This table needs serious work to be a bit more complete.
    """
    ncep_table = Dict(
        "tmax2m" => (0, 0, 4, "2m"),
        "tmin2m" => (0, 0, 5, "2m"),
        "tmpsfc" => (0, 0, 0, "Surface"),
        "tmp2m" => (0, 0, 0, "2m"),
		  "tmpprs" => (0, 0, 0, "Pressure"),
        "rh2m" => (0, 1, 1, "2m"),
        "apcpsfc" => (0, 1, 8, "Surface"),
		  "prateavesfc" => (0, 1, 52, "Surface"),
        "gustsfc" => (0, 2, 22, "Surface"),
        "prmslmsl" => (0, 3, 1, "MSL"),
		  "hgtprs" => (0, 3, 5, "Pressure"),
        "tozneclm" => (0, 14, 0, "Atmos"),
    )

    if haskey(ncep_table, ncep_var)
        GRIB_var = ncep_table[ncep_var]
    else
        GRIB_var = (255, 255, 255, "Unknown")
    end

    return GRIB_var
end

end
