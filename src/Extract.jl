using ArgParse
using Distributed

@everywhere using DataFrames
@everywhere using Statistics
@everywhere using Printf
@everywhere using GZip
@everywhere using CSV

@everywhere re_MolName = r"^#+ +Name: +"
@everywhere re_Number = r"^#+ +Number: +"
@everywhere re_OXR = r"^#+ +OXR +"
@everywhere re_Total = r"^#+ +Total Energy: +"

@everywhere function build_gzfn(sub_idx, cls_idx)
    # zero padded cls_idx
    zp_idx = @sprintf "%04d" parse(Int64, cls_idx)

    mol2_fn = "$sub_idx/test.$zp_idx.mol2.gz"

    return mol2_fn
end

@everywhere function AddMolecule(line_name, io, mol_set, sub_idx, cls_idx, coord_frame)
    values = [sub_idx, cls_idx]

    name = replace(strip(line_name), re_MolName => "")

    # skip 4
    [readline(io) for i in 1:4]
    number = replace(strip(readline(io)), re_Number => "")

    # exit if number is not found
    if !in(number, mol_set)
        return
    end

    # skip 4
    [readline(io) for i in 1:4]

    # prepare OXR
    oxr = []
    for i in 1:4
        line = split(replace(readline(io), re_OXR => ""))
        append!(oxr, line)
    end

    # skip 14
    [readline(io) for i in 1:14]

    # parse Total
    total = replace(strip(readline(io)), re_Total=>"")


    append!(values, [name, number, total])
    append!(values, oxr)

    push!(coord_frame, values)
end

@everywhere function ParseMol2(sub_idx, cls_idx, mol_frame, coord_frame)

    gz_fn = build_gzfn(sub_idx, cls_idx)
    mol_set = Set(mol_frame.mol_idx)

    GZip.open(gz_fn) do io
        while !eof(io)
            line = readline(io)
            if occursin(re_MolName, line)
                AddMolecule(line, io, mol_set, sub_idx, cls_idx, coord_frame)
            end
        end
    end
end

@everywhere function ParseOUTDOCK(subcluster, score_frame, time_frame)

    # define regex expressions
    re_header = r"mol#"
    re_elapsed = r"elapsed"
    re_values = r"^ +[0-9]"
    re_ignore = r"^ 9|colors|skip_size|poses|<|>|no_match|bump|clashes"
    re_timesub = r"(elapsed.+\(sec\): +| +\(hour\).+)"

    # define file name
    OUTDOCK = "$subcluster/OUTDOCK"

    # start file parse
    current_group = 0
    open(OUTDOCK) do io
        while !eof(io)
            line = readline(io)

            # new group found
            if occursin(re_header, line)
                current_group += 1

            elseif occursin(re_values, line)
                if occursin(re_ignore, line)
                    continue
                end
                values = split(line)
                prepend!(values, [subcluster, string(current_group)])
                push!(score_frame, values)

            # final line found
            elseif occursin(re_elapsed, line)
                time = replace(line, re_timesub => "")
                push!(time_frame, [subcluster, time])
            end

        end
    end
end

@everywhere function ProcessDir(dir_name, q)
    println(dir_name)

    # initialize dataframes
    score_frame = DataFrame(
        sub_idx = String[],
        cls_idx = String[],
        mol_idx = String[],
        mol_name = String[],
        flex_code = String[],
        matched = String[],
        nscored = String[],
        time = String[],
        hac = String[],
        setnum = String[],
        matnum = String[],
        rank = String[],
        cloud = String[],
        elect = String[],
        gist = String[],
        vdW = String[],
        psol = String[],
        asol = String[],
        inter = String[],
        rec_e = String[],
        rec_d = String[],
        r_hyd = String[],
        Total = String[]
    )
    time_frame = DataFrame(
        sub_idx = [],
        time = []
    )
    coord_frame = DataFrame(
        sub_idx = String[],
        cls_idx = String[],
        mol_idx = String[],
        mol_name = String[],
        mol_total = String[],
        OXR_X1 = String[],
        OXR_Y1 = String[],
        OXR_Z1 = String[],
        OXR_X2 = String[],
        OXR_Y2 = String[],
        OXR_Z2 = String[],
        OXR_X3 = String[],
        OXR_Y3 = String[],
        OXR_Z3 = String[],
        OXR_X4 = String[],
        OXR_Y4 = String[],
        OXR_Z4 = String[]
    )

    # generate sdi subcluster list
    subcluster_list = [
        joinpath(dir_name, i) for i in readdir(dir_name) if occursin("subcluster", i)
    ]

    # parse each subcluster OUTDOCK
    for s in subcluster_list
        ParseOUTDOCK(s, score_frame, time_frame)
    end

    # convert types
    time_frame[!, :time] = map(x -> parse(Float64, x), time_frame[!, :time])
    score_frame[!, :Total] = map(x -> parse(Float64, x), score_frame[!, :Total])

    # sort on time and score/molecule
    sort!(time_frame, :time)
    sort!(score_frame, [:Total, :mol_name])

    # take best scoring unique molecules
    unique_score_frame = unique(score_frame, [:mol_name])

    # write dataframes
    CSV.write(
        "$dir_name/elapsed_time.tab",
        time_frame, delim="\t", header=false
        )
    CSV.write(
        "$dir_name/extract_all.sort.txt",
        select(score_frame, Not(:cls_idx)),
        delim="\t", header=false
        )
    CSV.write(
        "$dir_name/extract_all.sort.unique.txt",
        select(unique_score_frame, Not(:cls_idx)),
        delim="\t", header=false
        )


    # calculate percentile
    pc = quantile!(unique_score_frame.Total, q)

    # filter scores for hit molecules
    top_score_frame = filter(x -> x.Total < pc, unique_score_frame)
    sort!(top_score_frame, [:sub_idx, :cls_idx, :mol_idx])

    # parse mol2 gzipped files for coordinates
    gdf = groupby(top_score_frame, [:sub_idx, :cls_idx])
    for (key, mol_frame) in pairs(gdf)
        ParseMol2(key.sub_idx, key.cls_idx, mol_frame, coord_frame)
    end

    # write coordinate dataframe
    CSV.write(
        "$dir_name/coords.tab",
        coord_frame, delim="\t"
    )
end

function get_args()
    s = ArgParseSettings()

    @add_arg_table! s begin
        "--directories", "-i"
            help = "Directories to process"
            required = true
            nargs='+'
        "--quantile", "-q"
            help = "Quantile to filter scores (default = 0.1)"
            required = false
            default = 0.1
            arg_type = Float64
    end

    return parse_args(s)
end

function main()

    parsed_args = get_args()
    directories = parsed_args["directories"]

    arguments = [
        (d, parsed_args["quantile"]) for d in directories
    ]

    # Parallel Map Arguments
    pmap(
        (args)->ProcessDir(args...),
        arguments
        )

end



main()
