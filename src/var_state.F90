module var_state

   use const
   use const_idx, only : ENE

   implicit none
  
   ! ----------------------------------------------------------------
   real(PREC), save :: tempK
   real(PREC), save :: kT

   real(PREC), allocatable, save :: xyz(:,:)
   real(PREC), save :: energies(0:ENE%MAX)

end module var_state
