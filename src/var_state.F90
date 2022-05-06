module var_state

   use,intrinsic :: ISO_FORTRAN_ENV, only: INT64
   use const
   use const_idx, only : ENE

   implicit none
  
   ! ----------------------------------------------------------------
   integer, save :: nthreads

   logical, save :: restarted

   integer, save :: job

   ! Condition
   real(PREC), save :: tempK
   real(PREC), save :: kT
   integer(INT64), save :: rng_seed

   real(PREC), allocatable, save :: xyz(:,:)
   real(PREC), allocatable, save :: velos(:,:)
   real(PREC), allocatable, save :: accels(:,:)
   real(PREC), allocatable, save :: forces(:,:)

   ! Energies
   real(PREC), save :: energies(0:ENE%MAX)
   real(PREC), save :: Ekinetic

   ! Base pairs
   logical, save :: flg_bp_energy
   logical, allocatable, save :: bp_status(:)
   real(PREC), allocatable, save :: ene_bp(:)
   real(PREC), allocatable, save :: for_bp(:,:,:)

   ! MD
   integer, save :: integrator
   real(PREC), save :: viscosity_Pas
   real(PREC), save :: dt
   integer(INT64), save :: nstep
   integer, save :: nstep_save
   integer, save :: nstep_save_rst
   real(PREC), save :: nl_margin
   integer(INT64) :: istep

   ! Electrostatic
   real(PREC), save :: ionic_strength
   real(PREC), save :: length_per_charge
   real(PREC), save :: lambdaD
   real(PREC), save :: diele

   ! Annealing
   integer, save :: opt_anneal
   integer, save :: nanneal
   integer(INT64), allocatable, save :: anneal_step(:)
   real(PREC), allocatable, save :: anneal_tempK(:)
   integer :: ianneal
   integer(INT64) :: istep_anneal_next

   ! variable box
   logical, save :: flg_variable_box
   integer(INT64), save :: variable_box_step
   real(PREC), save :: variable_box_change(3)

end module var_state
