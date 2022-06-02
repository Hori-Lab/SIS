subroutine list_bp()

   use const
   use const_idx, only : SEQT, BPT
   use var_top, only : nmp_chain, seq, imp_chain, nchains
   use var_state, only : bp_status, ene_bp, for_bp
   use var_potential, only : nbp, bp_mp, bp_min_loop

   implicit none
  
   integer :: ichain, jchain
   integer :: n, i, j, j_start, imp, jmp
   integer :: ibp

   do n = 1, 2
   
      ibp = 0

      do ichain = 1, nchains 

         do jchain = ichain, nchains 

            do i = 2, nmp_chain(ichain)-1 

               imp = imp_chain(i, ichain)
    
               if (ichain == jchain) then
                  j_start = i + bp_min_loop + 1
               else
                  j_start = 2
               endif
    
               do j = j_start, nmp_chain(jchain)-1
                  
                  jmp = imp_chain(j, jchain)
    
                  if ((seq(i,ichain) == SEQT%G .and. seq(j, jchain) == SEQT%C) .or. &
                      (seq(i,ichain) == SEQT%C .and. seq(j, jchain) == SEQT%G) ) then

                     ibp = ibp + 1

                     if (n == 2) then
                        bp_mp(1,ibp) = imp
                        bp_mp(2,ibp) = jmp
                        bp_mp(3,ibp) = BPT%GC
                        !bp_U0(ibp) = bp_U0_GC
                     endif
    
                  else if ((seq(i,ichain) == SEQT%A .and. seq(j, jchain) == SEQT%U) .or. &
                           (seq(i,ichain) == SEQT%U .and. seq(j, jchain) == SEQT%A) ) then

                     ibp = ibp + 1

                     if (n == 2) then
                        bp_mp(1,ibp) = imp
                        bp_mp(2,ibp) = jmp
                        bp_mp(3,ibp) = BPT%AU
                        !bp_U0(ibp) = bp_U0_AU
                     endif
    
                  else if ((seq(i,ichain) == SEQT%G .and. seq(j, jchain) == SEQT%U) .or. &
                           (seq(i,ichain) == SEQT%U .and. seq(j, jchain) == SEQT%G) ) then

                     ibp = ibp + 1

                     if (n == 2) then
                        bp_mp(1,ibp) = imp
                        bp_mp(2,ibp) = jmp
                        bp_mp(3,ibp) = BPT%GU
                        !bp_U0(ibp) = bp_U0_GU
                     endif

                  endif
    
               enddo
            enddo
    
         enddo
      enddo

      if (n == 1) then
         nbp = ibp
         allocate(bp_mp(3, nbp))
         !allocate(bp_U0(nbp))
         allocate(bp_status(nbp))
         allocate(ene_bp(nbp))
         allocate(for_bp(3, 6, nbp))

         bp_mp(:,:) = 0
         !bp_U0(:) = 0.0_PREC
         bp_status(:) = .False.
         ene_bp(:) = 0.0_PREC
         for_bp(:,:,:) = 0.0_PREC
      endif
   enddo

   write(*,*) '#nbp: ', nbp

end subroutine list_bp
