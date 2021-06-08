#### DirTrees

export DirTree

"""
    DirTree

A struct containing all output directories.
"""
struct DirTree
    """
    the root output directory
    """
    output::String

    """
    the directory for output data
    """
    data::String

    """
    a dictionary, initialized by a Tuple of variable names,
    which contains the output directories for those variables.
    """
    vars::Dict{Symbol, String}
end

Base.getindex(A::DirTree, k::Symbol) = A.vars[k]

path_separator = Sys.iswindows() ? "\\" : "/"

"""
    DirTree(output_folder_name::AbstractString,
            in_vars::Tuple{Vararg{Symbol}})

Gets a DirTree struct so that output directories
can be grabbed using, e.g., dir_tree[:variable_name].
"""
function DirTree(
    output_folder_name::AbstractString,
    in_vars::Tuple{Vararg{Symbol}},
)
    root_dir = "output"
    output_folder_name = split(output_folder_name, ".")[end]
    output_folder_name = replace(output_folder_name, "(" => "")
    output_folder_name = replace(output_folder_name, ")" => "")

    output = joinpath(root_dir, output_folder_name)
    data = joinpath(output, "data")

    vars = Dict{Symbol, String}()
    vars[:initial_conditions] =
        joinpath(data, "InitialConditions") * path_separator
    vars[:processed_initial_conditions] =
        joinpath(data, "ProcessedInitialConditions") * path_separator
    vars[:solution_raw] = joinpath(data, "SolutionRaw") * path_separator
    vars[:solution_processed] =
        joinpath(data, "SolutionProcessed") * path_separator

    mkpath(output)
    mkpath(data)

    for (k, v) in vars
        mkpath(v)
    end

    return DirTree(output, data, vars)
end
