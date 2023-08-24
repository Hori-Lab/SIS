subroutine list_bp2()

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

                  if (.not. is_complement(seq(i-1, ichain), seq(j+1, jchain)) .and. &
                      .not. is_complement(seq(i+1, ichain), seq(j-1, jchain)) ) then
                     cycle
                  endif
    
                  ibp = ibp + 1

                  if (n == 2) then
                     bp_mp(1,ibp) = imp
                     bp_mp(2,ibp) = jmp
                  endif
    
               enddo
            enddo
    
         enddo
      enddo

      if (n == 1) then
         nbp = ibp
         allocate(bp_mp(2, nbp))
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

contains
   logical function is_complement(s1, s2)
   
      integer, intent(in) :: s1, s2

      is_complement = .False.

      if (s1 == SEQT%A) then
         if (s2 == SEQT%U) then
            is_complement = .True.
         endif

      else if (s1 == SEQT%U) then
         if (s2 == SEQT%A .or. s2 == SEQT%G) then
            is_complement = .True.
         endif

      else if (s1 == SEQT%G) then
         if (s2 == SEQT%C .or. s2 == SEQT%U) then
            is_complement = .True.
         endif

      else if (s1 == SEQT%C) then
         if (s2 == SEQT%G) then
            is_complement = .True.
         endif

      endif

   end function is_complement

end subroutine list_bp2
