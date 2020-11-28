module fTst

export f

function f(x::AbstractArray)::Real
    f = 2.0*x[1]^2 + 3.0*x[2]
    return f
end

end