module var_state

   use, intrinsic :: ISO_FORTRAN_ENV, only: INT64
   use const
   use const_idx, only : ENE
   use mt_stream, only : mt_state

   implicit none
  
   ! ----------------------------------------------------------------
   logical, save :: restarted
   logical, save :: reset_step

   integer, save :: job

   ! PRNG Me
   !integer(INT64) :: rng_seed
   integer :: rng_seed
   type(mt_state), allocatable :: mts(:)
   type(mt_state) :: mts_rep

   ! Condition
   real(PREC), save :: tempK
   real(PREC), save :: kT

   real(PREC), allocatable :: xyz(:,:,:)     ! (3, nmp, nrep_proc)
   real(PREC), allocatable :: velos(:,:,:)   ! (3, nmp, nrep_proc)
   real(PREC), allocatable :: accels(:,:,:)  ! (3, nmp, nrep_proc)
   !real(PREC), allocatable :: forces(:,:,:)  ! (3, nmp, nrep_proc)

   ! Energies
   real(PREC), allocatable :: energies(:,:)  ! (0:ENE%MAX, nrep_proc)
   real(PREC), allocatable :: Ekinetic(:)    ! (nrep_proc)

   ! Base pairs
   integer, save :: nstep_bp_MC
   logical, save :: flg_bp_MC
   logical, save :: flg_bp_energy
   logical, allocatable, save :: bp_status(:,:)  ! (nbp_max, nrep_proc)
   logical, allocatable, save :: bp_status_MC(:,:)  ! (nbp_max, nrep_proc)
   integer, allocatable, save :: nt_bp_excess(:)
   real(PREC), allocatable, save :: ene_bp(:,:)  ! (nbp_max, nrep_proc)
   real(PREC), allocatable, save :: for_bp(:,:,:)

   ! MD
   integer, save :: integrator
   real(PREC), save :: viscosity_Pas
   real(PREC), save :: dt
   integer(INT64), save :: nstep
   integer(INT64), save :: stop_wall_time_sec
   integer(INT64), save :: nstep_check_stop
   integer, save :: fix_com_origin

   integer(INT64), save :: nstep_save
   integer(INT64), save :: nstep_save_rst
   real(PREC), save :: nl_margin
   integer(INT64) :: istep

   ! Electrostatic
   real(PREC) :: ionic_strength
   real(PREC) :: length_per_charge
   real(PREC), allocatable :: lambdaD(:)           ! (nrep)
   real(PREC), allocatable :: diele(:)             ! (nrep)

   ! Annealing
   integer, save :: opt_anneal
   integer, save :: nanneal
   integer(INT64), allocatable, save :: anneal_step(:)
   real(PREC), allocatable, save :: anneal_tempK(:)
   integer :: ianneal
   integer(INT64) :: istep_anneal_next

   ! Temperature independent potential
   integer :: temp_independent
   real(PREC) :: tempK_ref
   real(PREC), allocatable :: diele_dTcoef(:)

   real(PREC), allocatable :: rg(:)

#ifdef DUMPFORCE
   logical, save :: flg_step_dump_force
#endif

end module var_state
