
##
##  This is just a hack:
##  The file contains a GAP function call,
##  and we cut out the relevant data.
##
# function nemo_matrices_from_GAP_Z_file( filename::String, F = nothing )
#     # Try to read the file.
#     str = contents_textfile( filename )
# 
#     if isnothing( F )
#       F = Nemo.ZZ
#     end
# 
#     gens = []
#     pos = findfirst( "generators", str )
#     pos = findnext( "\n", str, pos[1] )
#     pos2 = pos
#     while ! isnothing( pos )
#       pos2 = findnext( "]", str, pos2[1]+1 )
#       if pos2[1] < pos[1]
#         break
#       end
#       push!( gens, SubString( str, (pos[1]+1):(pos2[1]-1) ) )
#       pos = findnext( "[", str, pos2[1] )
#     end
# 
#     mats = [ split( mat, (',','\n'), keepempty = false ) for mat in gens ]
#     dim = Int( sqrt( length( mats[1] ) ) )
#     return [ Nemo.matrix( F, dim, dim, map( x -> parse( Int, x ), l ) )
#              for l in mats ]
# end
