
##  flaws:
##  Nemo.NmodRing( 7 )  -> error
##  Nemo.NmodRing( UInt64(7) )  -> o.k.


###################################
# some features missing in Nemo

#     quadratic_field( n::Int, nam::String )
# > Return the field extension of the rational number field by the square root
# > of the integer `n`.
# > 
# >   F, z = quadratic_field( 7, "z" )
#
function quadratic_field(n::Int, nam::String)
    R, x = Nemo.PolynomialRing(Nemo.QQ, "x")
    pol = x^2 - n
    if !isirreducible(pol)
        error("$pol is not irreducible")
    end
    return Nemo.NumberField(pol, nam)
end


#T Implement something less simpleminded.
import AbstractAlgebra.Generic.order

# need for fmpz_mat, fq_nmod_mat, ...
function my_order(mat)
    #T prescribe at least a type that forces Nemo objects!
    ord = 1
    pow = mat
    while !Nemo.isone(pow)
        pow = pow * mat
        ord = ord + 1
    end
    return ord
end

order(mat::Nemo.fmpz_mat) = my_order(mat)
order(mat::Nemo.fq_nmod_mat) = my_order(mat)

###################################
