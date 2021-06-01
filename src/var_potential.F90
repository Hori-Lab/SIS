module var_potential

   use const, only : PREC

   implicit none
  
   integer, save :: nbond
   integer, allocatable, save :: bond_mp(:,:)
   real(PREC), allocatable, save :: bond_para(:,:)
   real(PREC), parameter :: bond_k = 15.0_PREC
   real(PREC), parameter :: bond_r0 = 5.9_PREC

   integer, save :: nangl
   integer, allocatable, save :: angl_mp(:,:)
   real(PREC), allocatable, save :: angl_para(:,:)
   real(PREC), parameter :: angl_k = 10.0_PREC
   real(PREC), parameter :: angl_t0 = 2.618_PREC

   real(PREC), parameter :: Ubp_cutoff = 18.0
   real(PREC), parameter :: Ubp0 = - 5.0 / 3.0
   real(PREC), parameter :: Ubp_bond_k = 3.0
   real(PREC), parameter :: Ubp_bond_r = 13.8
   real(PREC), parameter :: Ubp_angl_k = 1.5
   real(PREC), parameter :: Ubp_angl_theta1 = 1.8326
   real(PREC), parameter :: Ubp_angl_theta2 = 0.9425
   real(PREC), parameter :: Ubp_dihd_k = 0.5
   real(PREC), parameter :: Ubp_dihd_phi1 = 1.8326
   real(PREC), parameter :: Ubp_dihd_phi2 = 1.1345

   integer, save :: nbp
   integer, allocatable, save :: bp_mp(:,:)  ! 1: imp1, 2: imp2, 3: bp type

end module var_potential
