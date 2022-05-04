module var_top

   use const
   use const_idx, only : ENE

   implicit none
  
   ! ----------------------------------------------------------------
   integer, save :: nrepeat

   integer, save :: nunit
   integer, save :: nmp
   integer, save :: nchains
   integer, allocatable, save :: nmp_chain(:)  ! Number of mp in each chain
   integer, allocatable, save :: imp_chain(:, :)   ! imp of the x-th particle of the y-th chain = imp_chain(x,y)
   integer, allocatable, save :: ichain_mp(:)      ! ichain of imp
   integer, allocatable, save :: lmp_mp(:)        ! local imp (in each chain) of imp

   integer, allocatable, save :: seq(:,:)
   real(PREC), allocatable, save :: mass(:)

   integer, allocatable, save :: inp_no_charge(:)
   logical, allocatable, save :: has_charge(:)
   real(PREC), allocatable, save :: charge(:)

end module var_top
