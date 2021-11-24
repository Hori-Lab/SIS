module var_top

   use const
   use const_idx, only : ENE

   implicit none
  
   ! ----------------------------------------------------------------
   integer, save :: nrepeat

   integer, save :: nunit
   integer, save :: nmp
   integer, save :: nchains
   integer, allocatable, save :: nmp_chain(:)
   integer, allocatable, save :: imp_chain(:, :)
   integer, allocatable, save :: ichain_mp(:)

   integer, allocatable, save :: seq(:,:)
   real(PREC), allocatable, save :: mass(:)

end module var_top
