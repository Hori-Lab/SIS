module var_potential

   use const, only : PREC
   use const_idx, only : BPT

   implicit none
  
   integer, save :: nbond
   integer, allocatable, save :: bond_mp(:,:)
   real(PREC), allocatable, save :: bond_para(:,:)
   real(PREC), save :: bond_k  ! = 15.0_PREC
   real(PREC), save :: bond_r0  ! = 5.9_PREC

   logical, save :: flg_angl_ReB
   integer, save :: nangl
   integer, allocatable, save :: angl_mp(:,:)
   real(PREC), allocatable, save :: angl_para(:,:)
   real(PREC), save :: angl_k  ! = 10.0_PREC
   real(PREC), save :: angl_t0  ! = 2.618_PREC
   
   integer, save :: ndih
   logical, save :: flg_dih_cos
   logical, save :: flg_dih_exp
   integer, allocatable, save :: dih_mp(:,:)
   real(PREC), allocatable, save :: dih_para(:,:)
   real(PREC), save :: dih_k
   real(PREC), save :: dih_w
   real(PREC), save :: dih_p0

   ! Nearest neighbour parameters
   real(PREC), allocatable :: NN_dG(:)  ! (1:NNT%MAX) if bp_model == 4
   real(PREC), allocatable :: NN_dH(:)  ! (1:NNT%MAX) if bp_model == 4
   real(PREC), allocatable :: NN_dS(:)  ! (1:NNT%MAX) if bp_model == 4

   ! Basepair
   integer,    save :: max_bp_per_nt
   integer,    save :: bp_min_loop  ! = 4 (in the original CAG work), = 3 (for mRNA)
   integer,    save :: bp_model

   integer, allocatable :: bp_map_0(:,:)
   integer, allocatable :: bp_map(:,:)
   real(PREC), allocatable :: bp_map_dG(:,:,:)   ! (nmp, nmp, nrep_proc)
   real(PREC) :: bp_cutoff_energy  ! 0.01 kcal/mol
   real(PREC) :: bp_cutoff_dist
   integer,    save :: bp_seqdep
      ! = 0 (Default): No sequence dependence. Only U0_GC, U0_AU, U0_GU are required.
      ! = 1: Sequence dependent parameters. All possible combinations of trinucleotide-dimer are required.

   real(PREC), save :: coef_dG
   real(PREC), save :: dH0
   real(PREC), save :: dS0

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
   integer, allocatable :: nbp(:)
   integer, save :: nbp_max  ! This defines the size of neighbor list
   integer, allocatable, save :: bp_mp(:,:,:)  ! (3, nbp, nrep_proc) 1: imp1, 2: imp2, 3: bp type

   ! WCA parameters
   real(PREC), save :: wca_sigma  ! = 10.0
   real(PREC), save :: wca_eps    ! = 2.0

   ! WCA list
   integer, allocatable :: nwca(:)
   integer, save :: nwca_max  ! This defines the size of neighbor list
   integer, allocatable, save :: wca_mp(:,:,:)  ! (1:2, nwca, nrep_proc) or (1:2, nwca_max, nrep_proc),
                                                ! 1:imp1, 2:imp2

   ! Electrostatic parameters
   logical, save :: flg_ele
   integer, save :: ele_cutoff_type
   real(PREC), save :: ele_cutoff_inp
   real(PREC), allocatable :: ele_cutoff(:)  ! (nrep_proc)
   real(PREC), allocatable :: ele_coef(:)  ! (nrep_proc)

   ! Electrostatic list
   integer, allocatable :: nele(:)
   integer, save :: nele_max
   integer, allocatable, save :: ele_mp(:, :, :)  ! (2, nele_max, nrep_proc)

   ! Stage potential parameters
   logical, save :: flg_stage
   real(PREC), save :: stage_sigma ! = wca_sigma = 10.0 Angstrom for now
   real(PREC), save :: stage_eps ! = 1.2 kcal/mol


end module var_potential
