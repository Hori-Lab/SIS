module var_top

   use const
   use const_idx, only : ENE

   implicit none
  
   ! ----------------------------------------------------------------
   integer :: nunit
   integer :: nmp
   integer :: nchains
   integer, allocatable :: nmp_chain(:)  ! Number of mp in each chain
   integer, allocatable :: imp_chain(:, :)   ! imp of the x-th particle of the y-th chain = imp_chain(x,y)
   integer, allocatable :: ichain_mp(:)      ! ichain of imp
   integer, allocatable :: lmp_mp(:)        ! local imp (in each chain) of imp

   integer, allocatable :: moltypes(:)   ! chain ID --> MOLT%XXX
   integer, allocatable :: seq(:,:)
   real(PREC), allocatable :: mass(:)

   integer, allocatable :: inp_no_charge(:)
   logical, allocatable :: has_charge(:)
   logical, save :: dummy_has_charge
   !real(PREC), allocatable :: charge(:,:)  ! (nmp, nrep_all)

   logical :: flg_freeze
   integer, allocatable :: frz_ranges(:,:)
   logical, allocatable :: is_frozen(:)

end module var_top
