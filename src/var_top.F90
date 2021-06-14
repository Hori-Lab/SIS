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

   logical, save :: flg_pbc
   real(PREC), save :: pbc_box(3)
   real(PREC), save :: pbc_box_half(3)

end module var_top
