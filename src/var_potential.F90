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

   integer,    save :: bp_min_loop  ! = 4 (in the original CAG work), = 3 (for mRNA)
   real(PREC), save :: bp_cutoff  ! = 18.0
   !real(PREC), save :: bp_U0  ! = - 5.0 / 3.0
   integer,    save :: bp_seqdep
      ! = 0 (Default): No sequence dependence. Only U0_GC, U0_AU, U0_GU are required.
      ! = 1: Sequence dependent parameters. All possible combinations of trinucleotide-dimer are required.
   real(PREC), save :: bp_U0_GC
   real(PREC), save :: bp_U0_AU
   real(PREC), save :: bp_U0_GU
   real(PREC), save :: bp_bond_k  ! = 3.0
   real(PREC), save :: bp_bond_r  ! = 13.8
   real(PREC), save :: bp_angl_k  ! = 1.5
   real(PREC), save :: bp_angl_theta1  ! = 1.8326
   real(PREC), save :: bp_angl_theta2  ! = 0.9425
   real(PREC), save :: bp_dihd_k  ! = 0.5
   real(PREC), save :: bp_dihd_phi1  ! = 1.8326
   real(PREC), save :: bp_dihd_phi2  ! = 1.1345

   integer, save :: nbp
   integer, save :: nbp_max  ! This defines the size of neighbor list
   integer, allocatable, save :: bp_mp(:,:)  ! 1: imp1, 2: imp2, 3: bp type
   real(PREC), allocatable, save :: bp_U0(:)
   real(PREC), save :: bp_nl_cut2

   real(PREC), save :: wca_sigma  ! = 10.0
   real(PREC), save :: wca_eps    ! = 2.0
   real(PREC), save :: wca_nl_cut2 ! = (wca_sigma + nl_margin) ** 2

   integer, save :: nwca
   integer, save :: nwca_max  ! This defines the size of neighbor list
   integer, allocatable, save :: wca_mp(:,:)  ! (2, nwca) or (2, nwca_max), 1:imp1, 2:imp2

   integer, save :: bp_type2nhb(1:3) = (/ 3, 2, 2/)

end module var_potential
