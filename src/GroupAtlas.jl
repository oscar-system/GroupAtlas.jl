module GroupAtlas

#############################################################################

# for using LaTex in docstrings processed by Documenter.jl
import Markdown

# because of the "conditional dependency" on GAP (using Browse, see below)
# (see also https://github.com/JuliaLang/Pkg.jl/issues/1285)
import Requires

# the tables of contents of the database are in JSON format,
# as well as the data files for matrix groups over number fields
import JSON

# we want to transfer remote files
import HTTP

import Primes
import AbstractAlgebra
import Nemo

# provide the possibility to read gzipped data files;
# do not use Gzip.jl (errors with libz.so) or Libz.jl (deprecated)
import CodecZlib

#############################################################################

# small auxiliary functions
include("utils.jl")

# an interface between the C-MeatAxe text format of permutations and
# matrices and the corresponding Nemo formats
include("meataxe.jl")

# an interface between the JSON format of matrices in characteristic zero
# and the corresponding Nemo formats
include("charzero.jl")

# administration for the tables of contents
include("toc.jl")

# overview functions
include("overview.jl")

# actual access to remote/local files,
# and interpretation of the contents
include("access.jl")

# function for notifying the table of contents files
include("sources.jl")

function __init__()
    # If GAP is available then use its BrowseAtlasInfo feature in Julia,
    # for selecting available information.
    # This is currently the only feature in GroupAtlas.jl
    # for which GAP is needed.
    Requires.@require GAP = "c863536a-3901-11e9-33e7-d5cd0df7b904" begin
        GAP.Globals.LoadPackage(GAP.julia_to_gap("Browse"), false)
        include("browse.jl")
    end

    # Notify the tables of contents known at runtime;
    # We have the URLs of the JSON files and the local addresses
    # where these files may be stored, relative to a local GAP root path.
    notify_data_file(
        "http://www.math.rwth-aachen.de/~Thomas.Breuer/atlasrep/atlasprm.json",
        "atlasrep/atlasprm.json",
    )

    notify_data_file(
        "http://www.math.rwth-aachen.de/~Thomas.Breuer/atlasrep/datapkg/toc.json",
        "atlasrep/datapkg/toc.json",
    )

    notify_data_file(
        "http://www.math.rwth-aachen.de/~mfer/mfertoc.json",
        "mfer/mfertoc.json",
    )

    notify_data_file(
        "http://www.math.rwth-aachen.de/~Thomas.Breuer/ctblocks/ctblockstoc.json",
        "ctblocks/ctblockstoc.json",
    )

    # Initialize the data overview available at runtime.
    AtlasTOC("all")
end

end # module
