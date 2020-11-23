using ForwardDiff
using LinearAlgebra

function kAxial(x::AbstractArray)::AbstractArray
  L::Real = 20.0                                  # element length
  A = x[1]                                        # element cross-section area
  E = x[2]                                        # element elastic modulus
  k = A * E / L                                   # axial stiffness
  # ka::AbstractArray = k*[ 1.0 -1.0; -1.0 1.0]
  # ka = Symmetric(ka)
  ka::AbstractArray = Symmetric(k * [1.0 -1.0; -1.0 1.0]) # element stiffness matrix
  # @show ka
  return ka
end

function solver(x::AbstractArray)::AbstractArray
  Up::AbstractArray = [0.0]
  Fu::AbstractArray = [20000.0]
  # compute element stiffness matrix
  # assemble the element stiffness matrix into the global stiffness matrix
  Kg::AbstractArray = kAxial(x)
  # partition global stiffness matrix
  Kpp = Kg[1:1, 1:1]
  Kpu = Kg[1:1, 2:2]
  Kup = Kg[2:2, 1:1]
  Kuu = Kg[2:2, 2:2]
  # solve the displacement eqn
  Uu::AbstractArray = Kuu \ (Fu - Kup * Up)
  # @show Uu
  return Uu
end

function main()
  A::Real = 2.25
  E::Real = 10.0e6
  x::AbstractArray = [A, E]
  # solve the system of eqns
  # chkSolver = isa(solver(x), AbstractArray)
  Uu::AbstractArray = solver(x)
  # chkUu = isa(Uu, AbstractArray)
  @show Uu
  # solve for the gradient
  # cfg = ForwardDiff.JacobianConfig(solver, x, ForwardDiff.Chunk{2}())
  # jacUu = ForwardDiff.jacobian(solver, x, cfg)
  jacUu = ForwardDiff.jacobian(solver, x)
  @show jacUu
end

main()
