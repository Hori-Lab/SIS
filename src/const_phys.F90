module const_phys

  use const

  real(PREC), parameter :: EPS0 = 8.8541878128             ! Electric constant [F/m]
  real(PREC), parameter :: ELE  = 1.602176634e-19_PREC     ! Elementary charge [C]
  real(PREC), parameter :: BOLTZ_J = 1.380649e-23_PREC     ! Boltzmann constant [J/K]
  real(PREC), parameter :: N_AVO   = 6.02214076e23_PREC    ! Avogadro constant [/mol]
  real(PREC), parameter :: KCAL2JOUL = 4184.0              ! (kcal -> J)  [J/kcal]

  real(PREC), parameter :: JOUL2KCAL = 1.0/KCAL2JOUL   !< (J -> kcal)  [kcal/J]
  real(PREC), parameter :: JOUL2KCAL_MOL  = JOUL2KCAL * N_AVO  !< (J -> kcal/mol)
  real(PREC), parameter :: BOLTZ_KCAL_MOL = BOLTZ_J * JOUL2KCAL_MOL   !< Boltzmann constant [kcal/mol/K]

endmodule const_phys
