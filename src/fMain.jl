using Pkg
pkg"activate ."
pkg"instantiate"

using ForwardDiff
using Test

using fTst

function main()
  x::AbstractArray = [2.0, 2.0]
  @show f(x)

  g = ForwardDiff.gradient(f, x)
  @show g
  
  h = ForwardDiff.hessian(f, x)
  @show h

end

main()
