module feaTest

  using ForwardDiff
  using LinearAlgebra

  const Ebrass = 15.0e6                                               # brass element elastic modulus
  const Lbrass = 20.0                                                 # brass element length
  const Ealum = 10.0e6                                                # aluminum element elastic modulus
  const Lalum = 20.0                                                  # aluminum element length

  const L  = Lalum                                                    # post length [in]
  const E  = Ealum                                                    # elastic modulus [psi]
  const γ = 0.101                                                     # weight density [lbf/in³]

  const Sy = 40.0*1000.0                                              # yield strength [psi]
  const nSy = 1.11                                                    # factor of safety on yield

  const δAllowAxial = 0.005                                           # axial allowable displacement [in]

  function area(x::AbstractArray)::Real
    # A = Similar(x)
    A::Real = (pi/4.0)*(x[1]^2 - x[2]^2)                              # cross-section area
    return A
  end

  function weight(x::AbstractArray)::Real
    weight::Real = γ*area(x)*L                                        # post weight
    return weight
  end

  function axialDeformConstr(Uu::AbstractArray)::Real
    g1::Real = Uu[2] - δAllowAxial                                    # axial deformation constraint
    return g1
  end

  function axialStress(Uu::AbstractArray)::AbstractArray
    σAllowYield = Sy/nSy                                              # allowable yield stress [psi]
    B = (1.0/L)*[-1 1; -1 1]                                          # strain-displ matrix
    εAxial::AbstractArray = B*Uu                                      # axial strain
    σAxial::AbstractArray = E*[1 0; 0 1]*εAxial                       # axial stress
    return σAxial
  end

  function axialStressConstr(Uu::AbstractArray)::Real
    σAxial::AbstractArray = axialStress(Uu)                           # compute axial stress
    g2::Real = σAxial[1,1] - σAllowYield                              # axial stress constraint
    return g2
  end

  function assemble(Kg, Ke, assyv)::AbstractArray
    n = length(assyv)                                                 # length of the assembly vector
    for eRow = 1:n                                                    # for each row of element stiffness matrix
      gRow = assyv[eRow]                                              #   get the global stiffness matrix assembly index
      for eCol = 1:n                                                  #   for each column of element stiffness matrix
        gCol = assyv[eCol]                                            #     get the global stiffness matrix assembly index
        Kg[gRow, gCol] = Kg[gRow, gCol] + Ke[eRow, eCol]              #     add element stiffness coeff to global stiffness matrix
      end                                                             #   end do for each column
    end                                                               # end dof for each row
    return Kg                                                         # return global stiffness matrix
  end

  function kaxial(x::AbstractArray)::AbstractArray  # compute brass element stiffness matrix
    A = area(x)
    k = A*E/L                                                         # axial stiffness
    ke::AbstractArray = Symmetric(k*[ 1.0 -1.0; -1.0 1.0])            # element stiffness matrix
    return ke                                                         # return brass element stiffness matrix
  end

  function displacement(x::AbstractArray)::AbstractArray
    # Uu = similar(x)
    Kg::AbstractArray = zeros(3, 3)                                   # zero the global stiffness matrix
    # compute the stiffness matrix for the E1
    # A = area(x)                                                       # compute cross-section area
    K1::AbstractArray = kaxial(x)                                     # compute brass element stiffness matrix
    # Kg = assemble(Kg, K1, [3,1])                                      # assemble element stiffness matrices into global stiffness matrix [Kg]
    n = 2
    assyv = [3,1]
    for eRow = 1:n                                                    # for each row of element stiffness matrix
      gRow = assyv[eRow]                                              #   get the global stiffness matrix assembly index
      for eCol = 1:n                                                  #   for each column of element stiffness matrix
        gCol = assyv[eCol]                                            #     get the global stiffness matrix assembly index
        Kg[gRow, gCol] = Kg[gRow, gCol] + K1[eRow, eCol]              #     add element stiffness coeff to global stiffness matrix
      end                                                             #   end do for each column
    end                                                               # end dof for each row

    # compute the "hardcoded" stiffness matrix for E2
    K2::AbstractArray = kaxial(x)                                     # compute aluminum element stiffness matrix
    # Kg = assemble(Kg, K2, [1,2])                                      # assemble element stiffness matrices into global stiffness matrix [Kg]
    assyv = [1,2]
    for eRow = 1:n                                                    # for each row of element stiffness matrix
      gRow = assyv[eRow]                                              #   get the global stiffness matrix assembly index
      for eCol = 1:n                                                  #   for each column of element stiffness matrix
        gCol = assyv[eCol]                                            #     get the global stiffness matrix assembly index
        Kg[gRow, gCol] = Kg[gRow, gCol] + K2[eRow, eCol]              #     add element stiffness coeff to global stiffness matrix
      end                                                             #   end do for each column
    end                                                               # end dof for each row

    # partition the global stiffness matrix for solution
    Kuu::AbstractArray = Kg[1:2, 1:2]                                 # partition [Kg] to get [Kuu]
    Kup::AbstractArray = Kg[1:2, 3:3]                                 # partition [Kg] to get [Kup]
    Kpu::AbstractArray = Kg[3:3, 1:2]                                 # partiiton [Kg] to get [Kpu]
    Kpp::AbstractArray = Kg[3:3, 3:3]                                 # partition [Kg] to get [Kpp]

    # define the displacement BCs
    Up::AbstractArray = [0.0]                                         # prescribed displacement BCs - fixed end

    # define the applied loads
    Fu::AbstractArray = [0.0; 20000.0]                                # applied load at free end

    Uu::AbstractArray = Kuu\(Fu - Kup*Up)                             # solution of the partitioned form of displacement eqn
    return Uu
  end

  # function solver(x::AbstractArray)::AbstractArray              # solve for unknown displacement vector {Uu}
  #   Uu = similar(x)
  #   Uu = displacement(x)                                              # solve the fe displacement
  #   dUdx = ForwardDiff.jacobian(displacement, x)                      # jacobian of displacment
  #
  #   f  = weight(x)
  #   dfdx  = ForwardDiff.gradient(weight, x)                           # gradient of f
  #
  #   g1 = axialDeformConstr(Uu)                                        # axial deformation constraint
  #   dg1dx = ForwardDiff.gradient(axialDeformConstr, Uu)               # gradient of g1
  #
  #   g2 = axialStressConstr(Uu)                                        # axial stress constraint
  #   dg2dx = ForwardDiff.gradient(axialStressConstr, Uu)               # gradient of g2
  #
  #   return Uu                                                         # return unknown displacement solution {Uu}
  # end

  function reaction(Kpu, Kpp, Uu, Up)                                 # solve for unknown reaction force vector {Fr}
    Fr = Kpu*Uu + Kpp*Up                                              # solution of the partitioned form of reaction eqn
    return Fr                                                         # return unknown reaction forces {Fr}
  end

end
