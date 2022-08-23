subroutine list_bp()

   use const
   use const_idx, only : SEQT, BPT
   use var_top, only : nmp_chain, seq, imp_chain, nchains, nmp
   use var_state, only : bp_status, ene_bp, for_bp, nt_bp_excess
   use var_potential, only : nbp, bp_mp, bp_min_loop
   use var_replica, only : nrep_proc

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
                        bp_mp(1, ibp, 1:nrep_proc) = imp
                        bp_mp(2, ibp, 1:nrep_proc) = jmp
                        bp_mp(3, ibp, 1:nrep_proc) = BPT%GC
                        !bp_U0(ibp) = bp_U0_GC
                     endif
    
                  else if ((seq(i,ichain) == SEQT%A .and. seq(j, jchain) == SEQT%U) .or. &
                           (seq(i,ichain) == SEQT%U .and. seq(j, jchain) == SEQT%A) ) then

                     ibp = ibp + 1

                     if (n == 2) then
                        bp_mp(1, ibp, 1:nrep_proc) = imp
                        bp_mp(2, ibp, 1:nrep_proc) = jmp
                        bp_mp(3, ibp, 1:nrep_proc) = BPT%AU
                        !bp_U0(ibp) = bp_U0_AU
                     endif
    
                  else if ((seq(i,ichain) == SEQT%G .and. seq(j, jchain) == SEQT%U) .or. &
                           (seq(i,ichain) == SEQT%U .and. seq(j, jchain) == SEQT%G) ) then

                     ibp = ibp + 1

                     if (n == 2) then
                        bp_mp(1, ibp, 1:nrep_proc) = imp
                        bp_mp(2, ibp, 1:nrep_proc) = jmp
                        bp_mp(3, ibp, 1:nrep_proc) = BPT%GU
                        !bp_U0(ibp) = bp_U0_GU
                     endif

                  endif
    
               enddo
            enddo
    
         enddo
      enddo

      if (n == 1) then
         allocate(nbp(nrep_proc))
         allocate(bp_mp(3, ibp, nrep_proc))
         !allocate(bp_U0(nbp))
         allocate(bp_status(ibp, nrep_proc))
         allocate(ene_bp(ibp))
         allocate(for_bp(3, 6, ibp))
         allocate(nt_bp_excess(nmp))

         nbp(:) = ibp
         bp_mp(:,:,:) = 0
         !bp_U0(:) = 0.0_PREC
         bp_status(:,:) = .False.
         ene_bp(:) = 0.0_PREC
         for_bp(:,:,:) = 0.0_PREC
         nt_bp_excess(:) = 0
      endif
   enddo

   print '(a,i8)', '# nbp: ', nbp(1)

end subroutine list_bp
