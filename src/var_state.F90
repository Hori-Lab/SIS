module var_state

   use const
   use const_idx, only : ENE

   implicit none
  
   ! ----------------------------------------------------------------
   integer, save :: job

   real(PREC), save :: tempK
   real(PREC), save :: viscosity_Pas
   real(PREC), save :: kT

   real(PREC), allocatable, save :: xyz(:,:)
   real(PREC), save :: energies(0:ENE%MAX)
   real(PREC), allocatable, save :: forces(:,:)

   real(PREC), allocatable, save :: velos(:,:)
   real(PREC), allocatable, save :: accels(:,:)
   real(PREC), allocatable, save :: xyz_move(:,:)

   real(PREC), save :: dt
   integer, save :: nstep
   integer, save :: nstep_save

end module var_state
