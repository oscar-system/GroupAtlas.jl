
"""
    notify_data_file( url::String, localpath::String )

Tell the system that a JSON format file with an overview of
a part of the Group Atlas can be found at the URL `url`.

If [the Julia package GAP.jl](https://github.com/oscar-system/GAP.jl)
is loaded and if the `pkg` subdirectory of one of the GAP root paths
contains a file with path `localpath` then the same overview file
can be found there, and hence one need not download this file.

The locations of the overview files get collected in the `source_urls`
field of [`GroupAtlasData`](@ref),
where locally avaliable files are preferred to remote files.

(This function is called in the initialisation of the package.)
"""
function notify_data_file(url::String, localpath::String)
    global GroupAtlasData

    if localpath != "" && isdefined(Main, :GAP)
        # Some data may be stored in the local GAP installation.
        for path in GAP.gap_to_julia(Array{String,1}, GAP.Globals.GAPInfo.RootPaths)
            localpath2 = path * "/pkg/" * localpath
            if isfile(localpath2)
                push!(GroupAtlasData.source_urls, localpath2)
                return
            end
        end
    end

    # Take the remote files.
    push!(GroupAtlasData.source_urls, url)
end

# temporary hack:
# The JSON format files of the core part are stored in a different place;
# note that the current public formats that are used for the corresponding
# files from the core database are only GAP or Magma readable.
const source_url_json = "http://www.math.rwth-aachen.de/~Thomas.Breuer/atlasrep/datachar0/"
