module BinaryContour

using CSV
using GMT
using Dates
using Printf

export contour_to_bin, bin_to_contour, ContourHeader, GRIBparam

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
	Read the parameter description (name) and units from the WMO GRIB-2 paramete
	file.

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

end
