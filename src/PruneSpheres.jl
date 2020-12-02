using Printf

struct Sphere
    sphere_id::Int64
    cluster_id::Int64
    coordinates::Array{Float64, 1}
    line::String
end

function GetClusterHeader(cluster)
    """
    Returns the updated size of cluster
    """

    re_size = r"[0-9]+$"

    num_sph = length(cluster["spheres"])

    return replace(
        cluster["line"],
        re_size => @sprintf "%2i" num_sph
        )
end

function UpdateCluster!(clusters, sph, rm_info)
    """
    Removes a given sphere from its cluster
    """

    cluster_id = sph.cluster_id

    sph_list = clusters[cluster_id]["spheres"]

    rm_idx = 0
    for (idx, s) in enumerate(sph_list)
        if s === sph
            rm_idx = idx
        end
    end

    deleteat!(clusters[cluster_id]["spheres"], rm_idx)
    rm_info[cluster_id] += 1
end

function FindNearest(c, all_spheres)

    d_arr = zeros(length(all_spheres))

    for (idx, s) in enumerate(all_spheres)
        d = sum((s.coordinates - c).^2)
        d_arr[idx] = d
    end

    return all_spheres[argmin(d_arr)]
end

function ParseSphere(line)
    vals = split(line)
    coords = [
        parse(Float64, c) for c in vals[2:4]
    ]
    return coords
end

function ParseMatchingSpheres(fn)

    re_cls = r"^cluster"
    re_sph = r"^ 9"

    headers = []
    clusters = Dict()

    current_cluster = 0
    sphere_id = 1

    all_spheres = []
    open(fn, "r") do io
        while !eof(io)
            line = readline(io)

            if occursin(re_cls, line)
                current_cluster += 1
                clusters[current_cluster] = Dict(
                    "line" => line,
                    "spheres" => []
                )

            elseif occursin(re_sph, line)
                coords = ParseSphere(line)
                sph = Sphere(sphere_id, current_cluster, coords, line)

                push!(all_spheres, sph)
                push!(clusters[current_cluster]["spheres"], sph)
                sphere_id += 1

            else
                push!(headers, line)
            end

        end
    end

    return headers, clusters, all_spheres
end

function ParseCoordinates(fn, rm_set)

    rm_coords = []
    open(fn, "r") do io
        while !eof(io)
            line = readline(io)
            ms_id = split(line)[end]
            if in(ms_id, rm_set)
                coords = ParseSphere(line)
                push!(rm_coords, coords)
            end
        end
    end
    return rm_coords
end

function main(args)

    ms_fn, coord_fn, rm_set = args

    headers, clusters, all_spheres = ParseMatchingSpheres(ms_fn)
    rm_coords = ParseCoordinates(coord_fn, rm_set)


    rm_info = zeros(Int32, length(clusters))
    for c in rm_coords
        closest_sph = FindNearest(c, all_spheres)
        UpdateCluster!(clusters, closest_sph, rm_info)
    end

    # print out headers
    for h in headers
        println(h)
    end

    # print out cluster spheres
    for idx in 1:length(clusters)

        cluster = clusters[idx]
        ch = GetClusterHeader(cluster)

        # cluster header
        println(ch)

        # sphere lines
        for sph in cluster["spheres"]
            println(sph.line)
        end

    end


    total_rm = sum(rm_info)
    println(stderr, "Removed : $total_rm")
    for (idx, rm_size) in enumerate(rm_info)
        println(stderr, "  From Cluster $idx : $rm_size")
    end

end

function print_help()
    println("Error : Expecting at least 3 positional arguments")
    println("julia PruneSpheres.jl <matching_spheres> <usage_coordinates> <indices> ... ")
end

function get_args()
    if length(ARGS) < 3
        print_help()
        exit()
    end

    matching_sphere_fn = ARGS[1]
    coordinate_usage = ARGS[2]
    rm_set = Set(ARGS[3:end])

    return (matching_sphere_fn, coordinate_usage, rm_set)
end

args = get_args()
main(args)
