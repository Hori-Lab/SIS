module var_potential

   use const

   implicit none
  
   integer, save :: nbond
   integer, allocatable, save :: bond_mp(:,:)
   real(PREC), allocatable, save :: bond_para(:,:)

   integer, save :: nangl
   integer, allocatable, save :: angl_mp(:,:)
   real(PREC), allocatable, save :: angl_para(:,:)

   real(PREC), parameter :: Ubp0 = - 5.0 / 3.0
   real(PREC), parameter :: Ubp_bond_k = 3.0
   real(PREC), parameter :: Ubp_bond_r = 13.8
   real(PREC), parameter :: Ubp_angl_k = 1.5
   real(PREC), parameter :: Ubp_angl_theta1 = 1.8326
   real(PREC), parameter :: Ubp_angl_theta2 = 0.9425
   real(PREC), parameter :: Ubp_dihd_k = 0.5
   real(PREC), parameter :: Ubp_dihd_phi1 = 1.8326
   real(PREC), parameter :: Ubp_dihd_phi2 = 1.1345


end module var_potential
