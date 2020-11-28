# using Pkg
# pkg"activate ."
# pkg"instantiate"

using ForwardDiff
using Test

function f(x::AbstractArray)::Real
  f = 2.0*x[1]^2 + 3.0*x[2]
  return f
end

function main()
  x::AbstractArray = [2.0, 2.0]
  @show f(x)

  g = ForwardDiff.gradient(f, x)
  @show g
  
  h = ForwardDiff.hessian(f, x)
  @show h
  
end

main()
