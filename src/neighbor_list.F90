subroutine neighbor_list()
  
   use const, only : PREC
   use const_idx, only : SEQT
   use pbc, only : flg_pbc, pbc_vec, pbc_wrap
   use var_top, only : nmp_chain, seq, imp_chain, nchains, nmp
   use var_state, only : xyz
   use var_potential, only : wca_nl_cut2, nwca, nwca_max, wca_mp, &
                             bp_nl_cut2, bp_mp, nbp, nbp_max, bp_min_loop

   implicit none

   integer :: ichain, jchain, i, imp, j, jmp, j_start
   integer :: iwca, ibp
   integer :: nhb
   real(PREC) :: d2, v(3)

   if (flg_pbc) then
      call pbc_wrap()
   end if

   if (allocated(wca_mp)) then
      wca_mp(:,:) = 0
   else
      nwca_max = 5 * nmp
      allocate(wca_mp(2, nwca_max))
      wca_mp(:,:) = 0
   endif

   if (allocated(bp_mp)) then
      bp_mp(:,:) = 0
   else
      nbp_max = 5 * nmp
      allocate(bp_mp(3, nbp_max))
      bp_mp(:,:) = 0
   endif

   iwca = 0
   ibp = 0

   do ichain = 1, nchains 

      do jchain = ichain, nchains 

         do i = 1, nmp_chain(ichain)

            imp = imp_chain(i, ichain)
   
            if (ichain == jchain) then
               j_start = i + 1
            else
               j_start = 1
            endif
   
            do j = j_start, nmp_chain(jchain)
               
               jmp = imp_chain(j, jchain)
   
               v(:) = pbc_vec(xyz(:,imp) - xyz(:,jmp))
               d2 = dot_product(v,v)

               if (d2 <= wca_nl_cut2) then
                  if (ichain /= jchain .or. i+2 < j) then
                     iwca = iwca + 1
                     if (iwca > nwca_max) then
                        call reallocate_wca_mp()
                     endif

                     wca_mp(1,iwca) = imp
                     wca_mp(2,iwca) = jmp
                  endif
               endif

               if (d2 <= bp_nl_cut2) then
                  if (ichain /= jchain .or. i+bp_min_loop < j) then
                     nhb = 0
                     if ((seq(i,ichain) == SEQT%G .and. seq(j, jchain) == SEQT%C) .or. &
                         (seq(i,ichain) == SEQT%C .and. seq(j, jchain) == SEQT%G) ) then
                        nhb = 3
                     else if ((seq(i,ichain) == SEQT%A .and. seq(j, jchain) == SEQT%U) .or. &
                              (seq(i,ichain) == SEQT%U .and. seq(j, jchain) == SEQT%A) ) then
                        nhb = 2
                     else if ((seq(i,ichain) == SEQT%G .and. seq(j, jchain) == SEQT%U) .or. &
                              (seq(i,ichain) == SEQT%U .and. seq(j, jchain) == SEQT%G) ) then
                        nhb = 2
                     endif
                    
                     if (nhb > 0) then
                        ibp = ibp + 1
                        if (ibp > nbp_max) then
                           call reallocate_bp_mp()
                        endif
                        bp_mp(1,ibp) = imp
                        bp_mp(2,ibp) = jmp
                        bp_mp(3,ibp) = nhb
                     endif
                  endif
               endif
   
            enddo
         enddo
   
      enddo
   enddo

   nwca = iwca
   nbp = ibp

contains
   
   subroutine reallocate_wca_mp()
      
      integer :: old_max
      integer :: tmp(2, nwca_max)

      old_max = nwca_max   
      tmp(1:2, 1:old_max) = wca_mp(1:2, 1:old_max)

      deallocate(wca_mp)

      nwca_max = int(nwca_max * 1.2)

      allocate(wca_mp(2, nwca_max))

      wca_mp(:, :) = 0
      wca_mp(1:2, 1:old_max) = tmp(1:2, 1:old_max)

   endsubroutine reallocate_wca_mp

   subroutine reallocate_bp_mp()
      
      integer :: old_max
      integer :: tmp(3, nbp_max)

      old_max = nbp_max   
      tmp(1:3, 1:old_max) = bp_mp(1:3, 1:old_max)

      deallocate(bp_mp)

      nbp_max = int(nbp_max * 1.2)

      allocate(bp_mp(3, nbp_max))

      bp_mp(:, :) = 0
      bp_mp(1:3, 1:old_max) = tmp(1:3, 1:old_max)

   endsubroutine reallocate_bp_mp

endsubroutine neighbor_list
