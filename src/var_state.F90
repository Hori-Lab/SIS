module var_state

   use,intrinsic :: ISO_FORTRAN_ENV, only: INT64
   use const
   use const_idx, only : ENE

   implicit none
  
   ! ----------------------------------------------------------------
   integer, save :: nthreads

   integer, save :: job

   ! Condition
   real(PREC), save :: tempK
   real(PREC), save :: kT

   real(PREC), allocatable, save :: xyz(:,:)
   real(PREC), allocatable, save :: velos(:,:)
   real(PREC), allocatable, save :: accels(:,:)
   real(PREC), allocatable, save :: forces(:,:)

   real(PREC), save :: energies(0:ENE%MAX)
   real(PREC), save :: Ekinetic

   ! MD
   integer, save :: integrator
   real(PREC), save :: viscosity_Pas
   real(PREC), save :: dt
   integer(INT64), save :: nstep
   integer, save :: nstep_save
   real(PREC), save :: nl_margin

   ! variable box
   logical, save :: flg_variable_box
   integer(INT64), save :: variable_box_step
   real(PREC), save :: variable_box_change(3)

end module var_state
