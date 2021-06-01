subroutine list_local()

   use const
   use var_top, only : nmp_chain, imp_chain, nchains
   use var_potential, only : nbond, bond_mp, nangl, angl_mp

   implicit none
  
   integer :: ichain
   integer :: n, i, imp
   integer :: ibond, iangl

   do n = 1, 2
  
      ibond = 0
      iangl = 0

       do ichain = 1, nchains
    
          do i = 1, nmp_chain(ichain)-1
    
             imp = imp_chain(i, ichain)
             
             ibond = ibond + 1

             if (i /= nmp_chain(ichain) - 1) then
                iangl = iangl + 1
             endif

             if (n == 2) then

                bond_mp(1, ibond) = imp
                bond_mp(2, ibond) = imp + 1

                if (i /= nmp_chain(ichain) - 1) then
                   angl_mp(1, iangl) = imp
                   angl_mp(2, iangl) = imp + 1
                   angl_mp(3, iangl) = imp + 2
                endif

             endif

          enddo
    
       enddo

       if (n == 1) then
          nbond = ibond
          nangl = iangl
          allocate(bond_mp(2, nbond))
          allocate(angl_mp(3, nangl))
       endif
   enddo

   write(*,*) '#nbond: ', nbond
   write(*,*) '#nangl: ', nangl

end subroutine list_local
