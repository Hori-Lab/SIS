module var_potential

   use const

   implicit none
  
   integer, save :: nbond
   integer, allocatable, save :: bond_mp(:,:)
   real(PREC), allocatable, save :: bond_para(:,:)

end module var_potential
