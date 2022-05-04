subroutine list_ele()

   use const
   use var_top, only : nmp_chain, imp_chain, nchains, has_charge
   use var_potential, only : nele, ele_mp

   implicit none

   integer :: ichain, jchain
   integer :: n, i, j, j_start, imp, jmp
   integer :: iele

   if (allocated(ele_mp)) then
      deallocate(ele_mp)
   endif

   do n = 1, 2

      iele = 0

      do ichain = 1, nchains 

         do jchain = ichain, nchains 

            do i = 1, nmp_chain(ichain)

               imp = imp_chain(i, ichain)

               if (.not. has_charge(imp)) cycle

               if (ichain == jchain) then
                  j_start = i + 3
               else
                  j_start = 1
               endif

               do j = j_start, nmp_chain(jchain)

                  jmp = imp_chain(j, jchain)

                  if (.not. has_charge(jmp)) cycle

                  iele = iele + 1

                  if (n == 2) then
                     ele_mp(1,iele) = imp
                     ele_mp(2,iele) = jmp
                  endif

               enddo
            enddo

         enddo
      enddo

      if (n == 1) then
         nele = iele
         allocate(ele_mp(2, nele))
      endif
   enddo

   write(*,*) '#nele: ', nele

end subroutine list_ele
