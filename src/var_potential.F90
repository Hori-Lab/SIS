module var_potential

   use const, only : PREC
   use const_idx, only : BPT

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

   ! Basepair
   integer,    save :: max_bp_per_nt
   integer,    save :: bp_min_loop  ! = 4 (in the original CAG work), = 3 (for mRNA)
   real(PREC), save :: bp_cutoff_energy  ! 0.01 kcal/mol
   integer,    save :: bp_seqdep
      ! = 0 (Default): No sequence dependence. Only U0_GC, U0_AU, U0_GU are required.
      ! = 1: Sequence dependent parameters. All possible combinations of trinucleotide-dimer are required.

   type basepair_parameters
      real(PREC) :: cutoff_ddist
      real(PREC) :: U0
      real(PREC) :: bond_k
      real(PREC) :: bond_r
      real(PREC) :: angl_k1
      real(PREC) :: angl_k2
      real(PREC) :: angl_theta1
      real(PREC) :: angl_theta2
      real(PREC) :: dihd_k1
      real(PREC) :: dihd_k2
      real(PREC) :: dihd_phi1
      real(PREC) :: dihd_phi2
   end type basepair_parameters

   type(basepair_parameters), save :: bp_paras(BPT%MAX)

   integer, save :: bp_type2nhb(1:3) = (/ 3, 2, 2/)

   ! Basepair list
   integer, save :: nbp
   integer, save :: nbp_max  ! This defines the size of neighbor list
   integer, allocatable, save :: bp_mp(:,:)  ! 1: imp1, 2: imp2, 3: bp type
   real(PREC), save :: bp_nl_cut2

   ! WCA parameters
   real(PREC), save :: wca_sigma  ! = 10.0
   real(PREC), save :: wca_eps    ! = 2.0
   real(PREC), save :: wca_nl_cut2 ! = (wca_sigma + nl_margin) ** 2

   ! WCA list
   integer, save :: nwca
   integer, save :: nwca_max  ! This defines the size of neighbor list
   integer, allocatable, save :: wca_mp(:,:)  ! (2, nwca) or (2, nwca_max), 1:imp1, 2:imp2

   ! Electrostatic parameters
   logical, save :: flg_ele
   integer, save :: ele_cutoff_type
   real(PREC), save :: ele_cutoff_inp
   real(PREC), save :: ele_cutoff
   real(PREC), save :: ele_coef
   real(PREC), save :: ele_coef_QQ

   ! Electrostatic list
   integer, save :: nele
   integer, save :: nele_max
   integer, allocatable, save :: ele_mp(:, :)

   real(PREC), save :: ele_nl_cut2


end module var_potential
