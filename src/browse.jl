
export BrowseAtlasInfo

"""
    BrowseAtlasInfo( groupname::String = "all" )

Open an interactive table for browsing through the Atlas information
that is available for the group `groupname` or for all Atlas groups.
One can scroll, search, select, etc. in this table via key strokes,
and finally return to the Julia prompt by entering `Q`.

This function can be used only in a terminal
and if [the Julia package GAP.jl](https://github.com/oscar-system/GAP.jl)
is loaded and if
[the GAP package Browse](http://www.math.rwth-aachen.de/~Browse) is available.
In this case,
[ncurses](https://invisible-island.net/ncurses/ncurses.faq.html)
based GAP functions are called when the user navigates in the table.
The arrow keys are used to scroll in the table,
entering `?` opens an overview of the supported inputs and their meaning,
and entering `Q` closes the interactive table.
The example shown below lists a series of inputs.

The function returns an array of all those entries that have been selected
in the `Browse` table.
Each such entry is a dictionary with the same format as is returned by
[`OneAtlasInfo`](@ref); such dictionaries can be used as input for [`AtlasGenerators`](@ref).

# Examples
```julia
julia> keys = GAP.Globals.NCurses.keys;
julia> d = keys.DOWN;  r = keys.RIGHT;  c = keys.ENTER;
julia> GAP.Globals.BrowseData.SetReplay(   # simulate keyboard inputs
         GAP.julia_to_gap( 
           [ '/', 'A', '5', # Find the string A5 ...
             d, d, r,       # ... such that just the word matches,
             c,             # start the search,
             c,             # click the table entry A5,
             d, d,          # move down two rows,
             c,             # click the row for this representation,
             'Q',           # quit the second level table,
             d, d,          # move down two rows,
             c,             # click the table entry A6,
             d,             # move down one row,
             c,             # click the first row,
             'Q',           # quit the second level table,
             'Q' ],         # and quit the application.
           Val(true) ) )

julia> tworeps = BrowseAtlasInfo()     # show the overview of Atlas data

julia> GAP.Globals.BrowseData.SetReplay( false )  # reset the replay feature

julia> map( x -> x[ :name ], tworeps )
2-element Array{String,1}:
 "A5"
 "A6"

```
"""
function BrowseAtlasInfo(groupname::String = "all")
    # Check whether the Browse package is available.
    # There was an attempt to load it in the `__init__` function.
    if !GAP.Globals.IsBoundGlobal(GAP.julia_to_gap("BrowseAtlasInfo"))
        println("The GAP function 'BrowseAtlasInfo' is not available.")
        return []
    end

    # Do not distinguish core data from other sources.
    GAP.Globals.SetUserPreference(
        GAP.julia_to_gap("AtlasRep"),
        GAP.julia_to_gap("AtlasRepMarkNonCoreData"),
        GAP.julia_to_gap(""),
    )

    # Show the interactive Browse table.
    # (The condition posed by `Identifier` restricts to representations.)
    if groupname == "all"
        info = GAP.Globals.BrowseAtlasInfo(GAP.Globals.Identifier, GAP.Globals.ReturnTrue)
    else
        info = GAP.Globals.BrowseAtlasInfo(
            GAP.julia_to_gap(groupname),
            GAP.Globals.Identifier,
            GAP.Globals.ReturnTrue,
        )
    end

    # Replace each GAP record by the corresponding Julia dictionary.
    # The keys are different in GAP and Julia,
    # and some results returned by GAP may be not available in Julia.)
    res = []
    for cand in GAP.gap_to_julia(Array{Dict{Symbol,Any},1}, info)
        inpath = "/" * cand[:repname] * "."
        for r in AllAtlasInfos(Dict(:tocID => cand[:contents], :name => cand[:groupname]))
            if !isnothing(findfirst(inpath, r[:path]))
                push!(res, r)
                break
            end
        end
    end

    return res
end
