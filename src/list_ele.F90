subroutine list_ele()

   use const
   use var_top, only : nmp_chain, imp_chain, nchains, has_charge
   use var_potential, only : nele, ele_mp, ele_exclude_covalent_bond_pairs, ele_exclude_covalent_angle_pairs
   use var_replica, only : nrep_proc

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
                  j_start = i + 1
               else
                  j_start = 1
               endif

               do j = j_start, nmp_chain(jchain)
                  if (j == i + 1 .and. ele_exclude_covalent_bond_pairs) cycle
                  if (j == i + 2 .and. ele_exclude_covalent_angle_pairs) cycle

                  jmp = imp_chain(j, jchain)

                  if (.not. has_charge(jmp)) cycle

                  iele = iele + 1

                  if (n == 2) then
                     ele_mp(1,iele,1:nrep_proc) = imp
                     ele_mp(2,iele,1:nrep_proc) = jmp
                  endif

               enddo
            enddo

         enddo
      enddo

      if (n == 1) then
         allocate(nele(nrep_proc))
         allocate(ele_mp(2, iele, nrep_proc))
         nele(:) = iele
         ele_mp(:,:,:) = 0
      endif
   enddo

   print '(a,i8)', '# nele: ', nele(1)

end subroutine list_ele
