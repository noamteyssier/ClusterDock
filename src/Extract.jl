using ArgParse
using Distributed

@everywhere using DataFrames
@everywhere using Statistics
@everywhere using Printf
@everywhere using GZip
@everywhere using CSV

@everywhere const re_MolName = r"^#+ +Name: +"
@everywhere const re_Number = r"^#+ +Number: +"
@everywhere const re_OXR = r"^#+ +OXR +"
@everywhere const re_Total = r"^#+ +Total Energy: +"

@everywhere function build_gzfn(sub_idx, cls_idx)
    # zero padded cls_idx
    zp_idx = @sprintf "%04d" parse(Int64, cls_idx)

    mol2_fn = "$sub_idx/test.$zp_idx.mol2.gz"

    return mol2_fn
end

@everywhere function AUC(x, y)
    """
    Calculates AUC using the trapezoidal rule,
    x = FPR
    y = TPR
    """

    s = length(x) - 1
    areas = zeros(s)
    for i in 1:(s)
        x1 = x[i]
        x2 = x[i+1]
        y1 = y[i]
        y2 = y[i+1]

        dx = x2-x1
        dy = y2-y1

        a_i = dx * y1
        a_j = (dy * dx)/2
        a = a_i + a_j

        areas[i] = a
    end

    area = sum(areas)
    str_area = @sprintf "%.2f" area*100
    return str_area
end

@everywhere function LogAUC(x, y)
    """
    Calculates LogAUC using the adapted trapezoidal rule
    x = FPR
    y = TPR
    """

    # avoid a negative infinite X
    minimum_nonzero = min(x[x.>0]...)

    # set early zeros to earliest non_negative minimum
    x[x.==0] .= minimum_nonzero

    # transform FPR to log10
    log_x = log10.(x)

    # center minimum to zero
    log_x = log_x .+ abs(min(log_x...))

    # scale to maximum
    log_x = log_x ./ max(log_x...)

    s = length(x) - 1
    areas = zeros(s)
    for i in 1:(s)
        x1 = log_x[i]
        x2 = log_x[i+1]
        y1 = y[i]
        y2 = y[i+1]

        dx = x2-x1
        dy = y2-y1

        a_i = dx * y1
        a_j = (dy * dx)/2
        a = a_i + a_j

        areas[i] = a
    end

    area = sum(areas)
    str_area = @sprintf "%.2f" area*100
    return str_area
end

@everywhere function ReadName(name_fn)
    names = []
    open(name_fn) do io
        while !eof(io)
            line = readline(io)
            push!(names, line)
        end
    end
    return Set(names)
end

@everywhere function LabelType(frame, ligand_set, decoy_set)
    """
    Labels the ligand/decoy type for a given dataframe
    """

    frame[!, :Type] = map(
        x -> if in(x, ligand_set) "Ligand" elseif in(x, decoy_set) "Decoy" else missing end,
        frame[!, :mol_name]
    )
end

@everywhere function CalculateAUC(score_frame, ligand_set, decoy_set)

    LabelType(score_frame, ligand_set, decoy_set)

    dropmissing!(score_frame)
    sort!(score_frame, :Total)

    roc = DataFrame(FPR = Float64[], TPR = Float64[])

    type_arr = score_frame[!, :Type]
    total_ligands = sum(map(x -> x == "Ligand", type_arr))
    total_decoys = sum(map(x -> x == "Decoy", type_arr))

    num_ligands = 0
    num_decoys = 0
    total = 0
    for i in 1:size(score_frame)[1]
        values = []
        if type_arr[i] == "Ligand"
            num_ligands += 1
        else
            num_decoys += 1
        end
        total+=1

        TPR = num_ligands / total_ligands
        FPR = num_decoys / total_decoys

        values = [FPR, TPR]
        push!(roc, values)
    end

    auc = AUC(roc[!, :FPR], roc[!, :TPR])
    log_auc = LogAUC(roc[!, :FPR], roc[!, :TPR])

    return [auc, log_auc]
end

@everywhere function AddMolecule(line_name, io, mol_set, sub_idx, cls_idx, coord_frame)
    values = [sub_idx, cls_idx]

    name = replace(strip(line_name), re_MolName => "")

    # skip 4
    [readline(io) for i in 1:4]
    number = replace(strip(readline(io)), re_Number => "")

    # exit if number is not found
    if !in(number, mol_set)
        return false

    else
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

        return true
    end
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
        mol_name = String[],
        mol_idx = String[],
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

    # generate ligand/decoy filenames
    ligand_fn = joinpath(dir_name, "ligands.names")
    decoy_fn = joinpath(dir_name, "decoys.names")

    # generate ligand/decoy molecule name sets
    ligand_set = ReadName(ligand_fn)
    decoy_set = ReadName(decoy_fn)

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
        "$dir_name/extract_all.sort.txt",
        select(score_frame, Not(:cls_idx)),
        delim="\t", header=false
        )
    CSV.write(
        "$dir_name/extract_all.sort.unique.txt",
        select(unique_score_frame, Not(:cls_idx)),
        delim="\t", header=false
        )

    # Calculate AUC and LogAUC
    auc, log_auc = CalculateAUC(score_frame, ligand_set, decoy_set)
    time_frame[!, :AUC] .= auc
    time_frame[!, :LogAUC] .= log_auc

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

    # label type for ligand set
    LabelType(coord_frame, ligand_set, decoy_set)
    coord_frame[!, :mol_total] = map(
        x -> parse(Float64, x),
        coord_frame[!, :mol_total]
        )
    sort!(coord_frame, :mol_total)

    # write timing/enrichment dataframe
    CSV.write(
        "$dir_name/time_and_enrichment.tab",
        time_frame, delim="\t", header=false
        )

    # write coordinate dataframe
    CSV.write(
        "$dir_name/coords.tab",
        coord_frame, delim="\t"
    )
end

@everywhere function parallel_process(d, q)
    # try
    #     ProcessDir(d, q)
    # catch
    #     println("ERROR in {$d}")
    #     return
    # end
    ProcessDir(d, q)
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
        (args)->parallel_process(args...),
        arguments
        )

end



main()
