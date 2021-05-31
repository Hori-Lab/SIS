module var_state

   use const
   use const_idx, only : ENE

   implicit none
  
   ! ----------------------------------------------------------------
   integer, save :: nunit
   integer, save :: nmp
   integer, parameter :: nchains = 64
   integer, parameter :: nrepeat = 47
   integer, parameter :: nmp_chain = 3 * nrepeat

   real(PREC), allocatable, save :: xyz(:,:)
   real(PREC), save :: energies(0:ENE%MAX)

end module var_state
