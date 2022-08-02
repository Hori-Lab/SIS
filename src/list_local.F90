subroutine list_local()

   use const
   use var_top, only : nmp_chain, imp_chain, nchains
   use var_potential, only : nbond, bond_mp, nangl, angl_mp, ndihedral, dihedral_mp

   implicit none
  
   integer :: ichain
   integer :: n, i, imp
   integer :: ibond, iangl, idihedral

   print '(a)', 'Making lists of local interactions.'

   do n = 1, 2
  
      ibond = 0
      iangl = 0
      idihedral = 0
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
                !CHECK HERE - -1, -2 or -3
                if (i /= nmp_chain(ichain) - 1) then
                   dihedral_mp(1, idihedral) = imp
                   dihedral_mp(2, idihedral) = imp + 1
                   dihedral_mp(3, idihedral) = imp + 2
                   dihedral_mp(4, idihedral) = imp + 3
                endif

          enddo
    
       enddo

       if (n == 1) then
          nbond = ibond
          nangl = iangl
          ndihedral = idihedral
          allocate(bond_mp(2, nbond))
          allocate(angl_mp(3, nangl))
          allocate(dihedral_mp(4, ndihedral))
       endif
   enddo

   print '(a,i10)', "# Number of bonds: ", nbond
   print '(a,i10)', "# Number of angles: ", nangl
   print '(a,i10)', "# Number of dihedrals: ", ndihedral
   print *

end subroutine list_local
