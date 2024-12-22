subroutine list_local()

   use const
   use const_idx, only : MOLT
   use var_top, only : nmp_chain, imp_chain, nchains, moltypes
   use var_potential, only : nbond, bond_mp, nangl, angl_mp, ndih, dih_mp, flg_dih_cos, flg_dih_exp

   implicit none

   integer :: ichain
   integer :: n, i, imp, imp_ini, imp_end
   integer :: ibond, iangl, idih

   print '(a)', 'Making lists of local interactions.'

   do n = 1, 2

      ibond = 0
      iangl = 0
      idih = 0
      do ichain = 1, nchains

         do i = 1, nmp_chain(ichain)-1

            imp = imp_chain(i, ichain)

            ibond = ibond + 1

            if (i /= nmp_chain(ichain) - 1) then
               iangl = iangl + 1
            endif

            if (i < nmp_chain(ichain) - 2) then
               idih = idih + 1
            endif

            if (n == 2) then

               bond_mp(1, ibond) = imp
               bond_mp(2, ibond) = imp + 1

               if (i /= nmp_chain(ichain) - 1) then
                  angl_mp(1, iangl) = imp
                  angl_mp(2, iangl) = imp + 1
                  angl_mp(3, iangl) = imp + 2
               endif

               if ((flg_dih_cos .or. flg_dih_exp) .and. i < nmp_chain(ichain) - 2) then
                  dih_mp(1, idih) = imp
                  dih_mp(2, idih) = imp + 1
                  dih_mp(3, idih) = imp + 2
                  dih_mp(4, idih) = imp + 3
               endif
            endif

         enddo

         ! Circular
         if (moltypes(ichain) == MOLT%CIRCRNA) then

            if (n == 1) then
               ibond = ibond + 1
               iangl = iangl + 2
               idih = idih + 3

            else
               imp_ini = imp_chain(1, ichain)
               imp_end = imp_chain(nmp_chain(ichain), ichain)

               ibond = ibond + 1
               bond_mp(1, ibond) = imp_end
               bond_mp(2, ibond) = imp_ini

               iangl = iangl + 1
               angl_mp(1, iangl) = imp_end - 1
               angl_mp(2, iangl) = imp_end
               angl_mp(3, iangl) = imp_ini

               iangl = iangl + 1
               angl_mp(1, iangl) = imp_end
               angl_mp(2, iangl) = imp_ini
               angl_mp(3, iangl) = imp_ini + 1

               if (flg_dih_cos .or. flg_dih_exp) then
                  idih = idih + 1
                  dih_mp(1, idih) = imp_end - 2
                  dih_mp(2, idih) = imp_end - 1
                  dih_mp(3, idih) = imp_end
                  dih_mp(4, idih) = imp_ini

                  idih = idih + 1
                  dih_mp(1, idih) = imp_end - 1
                  dih_mp(2, idih) = imp_end
                  dih_mp(3, idih) = imp_ini
                  dih_mp(4, idih) = imp_ini + 1

                  idih = idih + 1
                  dih_mp(1, idih) = imp_end
                  dih_mp(2, idih) = imp_ini
                  dih_mp(3, idih) = imp_ini + 1
                  dih_mp(4, idih) = imp_ini + 2
               endif
            endif
         endif

      enddo

      if (n == 1) then
         nbond = ibond
         nangl = iangl
         allocate(bond_mp(2, nbond))
         allocate(angl_mp(3, nangl))

         if (flg_dih_cos .or. flg_dih_exp) then
            ndih = idih
            allocate(dih_mp(4, ndih))
         endif
      endif
   enddo

   print '(a,i10)', "# Number of bonds: ", nbond
   print '(a,i10)', "# Number of angles: ", nangl
   print '(a,i10)', "# Number of dihedrals: ", ndih
   print *

end subroutine list_local
