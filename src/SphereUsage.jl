using DelimitedFiles

struct MatchingSphere
    x::String
    y::String
    z::String
    ms_id::Int64
end

function add_sphere!(ms_dict, receptor, coords)
    if !haskey(ms_dict, receptor)
        ms_dict[receptor] = Dict()
    end

    if !haskey(ms_dict[receptor], coords)
        ms_id = length(ms_dict[receptor]) + 1
        ms_dict[receptor][coords] = ms_id
    end


    sphere = MatchingSphere(
        coords[1], coords[2], coords[3],
        ms_dict[receptor][coords]
    )

    return sphere
end

function add_spheres!(ms_dict, values, spheres)
    receptor = values[1]
    coords_1 = Tuple(values[9:11])
    coords_2 = Tuple(values[12:14])
    coords_3 = Tuple(values[15:17])
    coords_4 = Tuple(values[18:20])

    for (index, coords) in enumerate([coords_1, coords_2, coords_3, coords_4])
        spheres[index] = add_sphere!(ms_dict, receptor, coords)
    end
end

function update_usage!(sph_usage, values, s)
    index = Tuple(values[1:3])
    if !haskey(sph_usage, index)
        sph_usage[index] = Dict()
    end

    if !haskey(sph_usage[index], s)
        sph_usage[index][s] = Dict(
            :Usage => 0,
            :Ligand_Usage => 0,
            :Decoy_Usage => 0
        )
    end

    sph_usage[index][s][:Usage] += 1
    if values[end] == "ligand"
        sph_usage[index][s][:Ligand_Usage] += 1
    else
        sph_usage[index][s][:Decoy_Usage] += 1
    end
end

function update_cooccurrence!(co_occurrence, values, spheres)
    index = Tuple(values[1:3])

    if !haskey(co_occurrence, index)
        co_occurrence[index] = zeros((45,45))
    end

    for s_i in spheres
        for s_j in spheres
            if s_i.ms_id != s_j.ms_id
                co_occurrence[index][s_i.ms_id,s_j.ms_id] += 1
            end
        end
    end
end

function write_usage(sph_usage)
    open("sphere_usage.tab", "w") do io
        header = [
            "receptor", "match_type", "cluster_id",
            "x", "y", "z", "ms_id",
            "Usage", "Ligand_Usage", "Decoy_Usage"
        ]
        writedlm(io, [header])
        for idx in keys(sph_usage)
            for s in keys(sph_usage[idx])
                values = [
                    idx[1], idx[2], idx[3],
                    s.x, s.y, s.z, s.ms_id,
                    sph_usage[idx][s][:Usage],
                    sph_usage[idx][s][:Ligand_Usage],
                    sph_usage[idx][s][:Decoy_Usage],
                    ]
                writedlm(io,[values])
            end
        end

    end
end

function write_cooccurrence(co_occurrence)
    open("co-occurrence.tab", "w") do io
        header = vcat(
            ["receptor", "match_type", "cluster_id"],
            ["sph.$i" for i in 1:45]
            )
        writedlm(io, [header])

        for index in keys(co_occurrence)
            for i in eachrow(co_occurrence[index])
                values = vcat([_ for _ in index], [c for c in i])
                writedlm(io, [values])
            end
        end
    end
end

function main(fn)

    ms_dict = Dict()
    sph_usage = Dict()
    co_occurrence = Dict()

    open(fn) do io

        # preallocated
        spheres = Array{MatchingSphere, 1}(undef, 4)

        while !eof(io)

            # read and split line
            line = readline(io)
            values = split(strip(line))

            # create matching spheres
            add_spheres!(ms_dict, values, spheres)

            # update specific sample usage
            for s in spheres
                update_usage!(sph_usage, values, s)
            end

            # updates sphere co-occurrence matrix
            update_cooccurrence!(co_occurrence, values, spheres)

        end
    end

    write_usage(sph_usage)
    write_cooccurrence(co_occurrence)
end

if length(ARGS) != 1
    println("Error : Provide a coordinate dataframe to process")
    exit()
end

filename = ARGS[1]
main(filename)
