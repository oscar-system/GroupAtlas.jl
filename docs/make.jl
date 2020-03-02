push!(LOAD_PATH, "../src")

using Documenter
using GAP        # then Documenter finds the docstring for `BrowseAtlasInfo`
using GroupAtlas

makedocs(
    sitename = "GroupAtlas",
    format = Documenter.HTML(),
    modules = [GroupAtlas],
    doctest = false,
)

# Documenter can also automatically deploy documentation to gh-pages.
# See "Hosting Documentation" and deploydocs() in the Documenter manual
# for more information.
#=deploydocs(
    repo = "<repository url>"
)=#
