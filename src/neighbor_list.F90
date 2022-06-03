subroutine neighbor_list()
  
   use const, only : PREC
   use const_idx, only : SEQT, BPT
   use pbc, only : flg_pbc, pbc_vec_d, pbc_wrap
   use var_top, only : nmp_chain, seq, imp_chain, nchains, nmp, has_charge
   use var_state, only : xyz, bp_status, ene_bp, for_bp
   use var_potential, only : wca_nl_cut2, nwca, nwca_max, wca_mp, &
                             bp_nl_cut2, bp_mp, nbp, nbp_max, bp_min_loop, &
                             nele, nele_max, ele_mp, ele_nl_cut2, flg_ele

   implicit none

   integer :: ichain, jchain, i, imp, j, jmp, j_start
   integer :: iwca, ibp, iele
   integer :: bptype
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

   if (.not. allocated(bp_mp)) then
      nbp_max = nmp / 2
      allocate(bp_mp(3, nbp_max))
      allocate(bp_status(nbp_max))
      allocate(ene_bp(nbp_max))
      allocate(for_bp(3, 6, nbp_max))
   endif
   bp_mp(:,:) = 0
   bp_status(:) = .False.
   ene_bp(:) = 0.0_PREC
   for_bp(:,:,:) = 0.0_PREC

   if (flg_ele) then
      if (allocated(ele_mp)) then
         ele_mp(:,:) = 0
         !ele_coef(:) = 0.0_PREC
      else
         nele_max = 5 * nmp
         allocate(ele_mp(2, nele_max))
         ele_mp(:,:) = 0
         !allocate(ele_coef(nele_max))
         !ele_coef(:) = 0.0_PREC
      endif
   endif

   iwca = 0
   ibp = 0
   iele = 0

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
   
               v(:) = pbc_vec_d(xyz(:,imp), xyz(:,jmp))
               d2 = dot_product(v,v)

               ! WCA
               if (d2 <= wca_nl_cut2) then
                  if (ichain /= jchain .or. i+2 < j) then
                     iwca = iwca + 1
                     if (iwca > nwca_max) then
                        !write(*,*) 'Error: iwca > nwca_max. iwca =', iwca, 'nwca_max = ', nwca_max
                        call reallocate_wca_mp()
                     endif

                     wca_mp(1,iwca) = imp
                     wca_mp(2,iwca) = jmp
                  endif
               endif

               ! Electrostatic
               if (flg_ele .and. d2 <= ele_nl_cut2) then

                  if (has_charge(imp) .and. has_charge(jmp)) then

                     iele = iele + 1
                     if (iele > nele_max) then
                        call reallocate_ele_mp()
                     endif

                     ele_mp(1,iele) = imp
                     ele_mp(2,iele) = jmp
                     !ele_coef(iele) = charge(i) * charge(j)
                  endif
               endif

               ! Basepairs
               if (d2 <= bp_nl_cut2) then

                  if (i == 1 .or. j == 1 .or. i == nmp_chain(ichain) .or. j == nmp_chain(ichain)) then
                     ! Terminal nucleotides do not form base pairs.
                     continue

                  else if (ichain /= jchain .or. i+bp_min_loop < j) then

                     bptype = BPT%UNDEF
                     !coef = 0.0_PREC

                     if ((seq(i,ichain) == SEQT%G .and. seq(j, jchain) == SEQT%C) .or. &
                         (seq(i,ichain) == SEQT%C .and. seq(j, jchain) == SEQT%G) ) then
                        bptype = BPT%GC

                     else if ((seq(i,ichain) == SEQT%A .and. seq(j, jchain) == SEQT%U) .or. &
                              (seq(i,ichain) == SEQT%U .and. seq(j, jchain) == SEQT%A) ) then
                        bptype = BPT%AU

                     else if ((seq(i,ichain) == SEQT%G .and. seq(j, jchain) == SEQT%U) .or. &
                              (seq(i,ichain) == SEQT%U .and. seq(j, jchain) == SEQT%G) ) then
                        bptype = BPT%GU
                     endif

                     if (bptype /= BPT%UNDEF) then
                        ibp = ibp + 1
                        if (ibp > nbp_max) then
                           !write(*,*) 'Error: ibp > nbp_max. ibp =', ibp, 'nbp_max = ', nbp_max
                           call reallocate_bp_mp()
                        endif
                        bp_mp(1,ibp) = imp
                        bp_mp(2,ibp) = jmp
                        bp_mp(3,ibp) = bptype
                     endif
                  endif
               endif
   
            enddo
         enddo
   
      enddo
   enddo

   nwca = iwca
   nbp = ibp
   nele = iele

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

   subroutine reallocate_ele_mp()

      integer :: old_max
      integer :: tmp(2, nele_max)
      !real(PREC) :: tmp2(nele_max)

      old_max = nele_max   
      tmp(1:2, 1:old_max) = ele_mp(1:2, 1:old_max)
      !tmp2(1:old_max) = ele_coef(1:old_max)

      deallocate(ele_mp)
      !deallocate(ele_coef)

      nele_max = int(nele_max * 1.2)

      allocate(ele_mp(2, nele_max))
      !allocate(ele_coef(nele_max))

      ele_mp(:, :) = 0
      ele_mp(1:2, 1:old_max) = tmp(1:2, 1:old_max)

      !ele_coef(:) = 0.0_PREC
      !ele_coef(1:old_max) = tmp2(1:old_max)

   endsubroutine reallocate_ele_mp

   subroutine reallocate_bp_mp()
      
      integer :: old_max
      integer :: tmp(3, nbp_max)

      old_max = nbp_max   
      tmp(1:3, 1:old_max) = bp_mp(1:3, 1:old_max)

      deallocate(bp_mp)
      deallocate(bp_status)
      deallocate(ene_bp)
      deallocate(for_bp)

      nbp_max = int(nbp_max * 1.2)

      allocate(bp_mp(3, nbp_max))
      allocate(bp_status(nbp_max))
      allocate(ene_bp(nbp_max))
      allocate(for_bp(3, 6, nbp_max))

      bp_mp(:, :) = 0
      bp_mp(1:3, 1:old_max) = tmp(1:3, 1:old_max)

      bp_status(:) = .False.
      ene_bp(:) = 0.0_PREC
      for_bp(:,:,:) = 0.0_PREC

   endsubroutine reallocate_bp_mp

endsubroutine neighbor_list
