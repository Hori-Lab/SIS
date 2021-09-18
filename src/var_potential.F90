module var_potential

   use const, only : PREC

   implicit none
  
   integer, save :: nbond
   integer, allocatable, save :: bond_mp(:,:)
   real(PREC), allocatable, save :: bond_para(:,:)
   real(PREC), save :: bond_k  ! = 15.0_PREC
   real(PREC), save :: bond_r0  ! = 5.9_PREC

   integer, save :: nangl
   integer, allocatable, save :: angl_mp(:,:)
   real(PREC), allocatable, save :: angl_para(:,:)
   real(PREC), save :: angl_k  ! = 10.0_PREC
   real(PREC), save :: angl_t0  ! = 2.618_PREC

   real(PREC), save :: Ubp_cutoff  ! = 18.0
   real(PREC), save :: Ubp0  ! = - 5.0 / 3.0
   real(PREC), save :: Ubp_bond_k  ! = 3.0
   real(PREC), save :: Ubp_bond_r  ! = 13.8
   real(PREC), save :: Ubp_angl_k  ! = 1.5
   real(PREC), save :: Ubp_angl_theta1  ! = 1.8326
   real(PREC), save :: Ubp_angl_theta2  ! = 0.9425
   real(PREC), save :: Ubp_dihd_k  ! = 0.5
   real(PREC), save :: Ubp_dihd_phi1  ! = 1.8326
   real(PREC), save :: Ubp_dihd_phi2  ! = 1.1345
   integer,    save :: Ubp_min_loop  ! = 4 (in the original CAG work), = 3 (for mRNA)

   integer, save :: nbp
   integer, allocatable, save :: bp_mp(:,:)  ! 1: imp1, 2: imp2, 3: bp type

   real(PREC), save :: wca_sigma  ! = 10.0
   real(PREC), save :: wca_eps    ! = 2.0

   integer, save :: nwca
   integer, allocatable, save :: wca_mp(:,:)  ! 1: imp1, 2: imp2

end module var_potential
