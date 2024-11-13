subroutine init_bp()

   use, intrinsic :: iso_fortran_env, Only : output_unit
   use const, only : PREC
   use const_idx, only : SEQT, BPT, seqt2nnt, NNENDT, is_complement
   use var_io, only : flg_in_ct, flg_in_bpseq
   use var_top, only : nmp, seq, lmp_mp, ichain_mp, nmp_chain
   use var_potential, only : bp_model, bp_map, bp_min_loop, & !bp_map_dG, bp_map_0, &
                             bp_paras, bp_cutoff_energy, bp_cutoff_dist, bp3_map, &
                             NN_dH, NN_dS, dH0, dS0, bp3_dH, bp3_dS, &
                             flg_NNend, NNend_dH, NNend_dS, dHend0, dSend0
   !use var_replica, only : nrep_proc
   use var_parallel

   implicit none

   integer :: imp, jmp, bptype
   integer :: i, j, h, ichain, jchain
   integer :: w, x, u, z, y, v
   real(PREC) :: cutoff
   real(PREC) :: dH, dS
   logical :: comp_wz, comp_uv
   integer :: bp3_hash(1:4**6)   ! 4**6 = 4096

   allocate(bp_map(nmp, nmp))
   allocate(bp3_map(nmp, nmp))
   bp_map(:,:) = 0
   bp3_map(:,:) = 0

   !if (bp_model == 4 .or. bp_model == 5) then
   !   allocate(bp_map_dG(nmp, nmp, nrep_proc))
   !   bp_map_dG(:,:,:) = 0.0_PREC
   !endif


   if (bp_model == 5) then

      bp3_hash(:) = 0
      bp3_dH(:) = 0.0_PREC
      bp3_dS(:) = 0.0_PREC

      h = 0

      !! seq1: 5'-w x u-3'
      !! seq2: 3'-z y v-5'

      !! Centre base pair
      do x = 0, 2        ! A,         U,      G
         do y = x+1, 3   ! (U, G, C), (G, C), C

            !! Centre has to be paired.
            if (.not. is_complement(x,y)) cycle
               
            !! Upper base pair
            do w = 0, 3   ! A, U, G, C
               do z = 0, 3   ! A, U, G, C

                  comp_wz = is_complement(w, z)
                  if (comp_wz) then
                     dH = 0.5_PREC * (NN_dH(seqt2nnt(w, x, z, y)) - dH0)
                     dS = 0.5_PREC * (NN_dS(seqt2nnt(w, x, z, y)) - dS0)
                  else
                     dH = 0.0_PREC
                     dS = 0.0_PREC
                  endif

                  !! Lower base pair
                  do u = 0, 3   ! A, U, G, C
                     do v = 0, 3   ! A, U, G, C

                        comp_uv = is_complement(u, v)

                        !! At least either side has to be a complementary base pair
                        if ((.not. comp_wz) .and. (.not. comp_uv)) cycle

                        h = h + 1

                        i = bp3_hash_key(w, x, u, z, y, v)
                        bp3_hash(i) = h

                        if (comp_uv) then
                           bp3_dH(h) = dH + 0.5_PREC * (NN_dH(seqt2nnt(x, u, y, v)) - dH0)
                           bp3_dS(h) = (dS + 0.5_PREC * (NN_dS(seqt2nnt(x, u, y, v)) - dS0)) * 1.0e-3_PREC
                        else
                           bp3_dH(h) = dH
                           bp3_dS(h) = dS * 1.0e-3_PREC
                        endif

                        ! End effects
                        if (flg_NNend) then
                           ! u-v is not paired
                           if (comp_wz .and. (.not. comp_uv)) then
                              print *, x
                              print *, y
                              ! AU on AU/CG/GU
                              if ((x == SEQT%A .and. y == SEQT%U) .or. (x == SEQT%U .and. y == SEQT%A)) then
                                 if ((w == SEQT%A .and. z == SEQT%U) .or. (w == SEQT%U .and. z == SEQT%A)) then
                                    bp3_dH(h) = bp3_dH(h) + NNend_dH(NNENDT%AUonAU) - dHend0
                                    bp3_dS(h) = bp3_dS(h) + (NNend_dS(NNENDT%AUonAU) - dSend0) * 1.0e-3_PREC
                                 else if ((w == SEQT%G .and. z == SEQT%C) .or. (w == SEQT%C .and. z == SEQT%G)) then
                                    bp3_dH(h) = bp3_dH(h) + NNend_dH(NNENDT%AUonCG) - dHend0
                                    bp3_dS(h) = bp3_dS(h) + (NNend_dS(NNENDT%AUonCG) - dSend0) * 1.0e-3_PREC
                                 else if ((w == SEQT%G .and. z == SEQT%U) .or. (w == SEQT%U .and. z == SEQT%G)) then
                                    bp3_dH(h) = bp3_dH(h) + NNend_dH(NNENDT%AUonGU) - dHend0
                                    bp3_dS(h) = bp3_dS(h) + (NNend_dS(NNENDT%AUonGU) - dSend0) * 1.0e-3_PREC
                                 else
                                    print '(a)', 'logical error in init_bp'
                                    call sis_abort()
                                 endif

                              ! GU on AU/CG/GU
                              else if ((x == SEQT%G .and. y == SEQT%U) .or. (x == SEQT%U .and. y == SEQT%G)) then
                                 if ((w == SEQT%A .and. z == SEQT%U) .or. (w == SEQT%U .and. z == SEQT%A)) then
                                    bp3_dH(h) = bp3_dH(h) + NNend_dH(NNENDT%GUonAU) - dHend0
                                    bp3_dS(h) = bp3_dS(h) + (NNend_dS(NNENDT%GUonAU) - dSend0) * 1.0e-3_PREC
                                 else if ((w == SEQT%G .and. z == SEQT%C) .or. (w == SEQT%C .and. z == SEQT%G)) then
                                    bp3_dH(h) = bp3_dH(h) + NNend_dH(NNENDT%GUonCG) - dHend0
                                    bp3_dS(h) = bp3_dS(h) + (NNend_dS(NNENDT%GUonCG) - dSend0) * 1.0e-3_PREC
                                 else if ((w == SEQT%G .and. z == SEQT%U) .or. (w == SEQT%U .and. z == SEQT%G)) then
                                    bp3_dH(h) = bp3_dH(h) + NNend_dH(NNENDT%GUonGU) - dHend0
                                    bp3_dS(h) = bp3_dS(h) + (NNend_dS(NNENDT%GUonGU) - dSend0) * 1.0e-3_PREC
                                 else
                                    print '(a)', 'logical error in init_bp'
                                 endif
                              endif

                           ! w-z is not paired
                           else if ((.not. comp_wz) .and. comp_uv) then
                              ! AU on AU/CG/GU
                              if ((x == SEQT%A .and. y == SEQT%U) .or. (x == SEQT%U .and. y == SEQT%A)) then
                                 if ((u == SEQT%A .and. v == SEQT%U) .or. (u == SEQT%U .and. v == SEQT%A)) then
                                    bp3_dH(h) = bp3_dH(h) + NNend_dH(NNENDT%AUonAU) - dHend0
                                    bp3_dS(h) = bp3_dS(h) + (NNend_dS(NNENDT%AUonAU) - dSend0) * 1.0e-3_PREC
                                 else if ((u == SEQT%G .and. v == SEQT%C) .or. (u == SEQT%C .and. v == SEQT%G)) then
                                    bp3_dH(h) = bp3_dH(h) + NNend_dH(NNENDT%AUonCG) - dHend0
                                    bp3_dS(h) = bp3_dS(h) + (NNend_dS(NNENDT%AUonCG) - dSend0) * 1.0e-3_PREC
                                 else if ((u == SEQT%G .and. v == SEQT%U) .or. (u == SEQT%U .and. v == SEQT%G)) then
                                    bp3_dH(h) = bp3_dH(h) + NNend_dH(NNENDT%AUonGU) - dHend0
                                    bp3_dS(h) = bp3_dS(h) + (NNend_dS(NNENDT%AUonGU) - dSend0) * 1.0e-3_PREC
                                 else
                                    print '(a)', 'logical error in init_bp'
                                 endif

                              ! GU on AU/CG/GU
                              else if ((x == SEQT%G .and. y == SEQT%U) .or. (x == SEQT%U .and. y == SEQT%G)) then
                                 if ((u == SEQT%A .and. v == SEQT%U) .or. (u == SEQT%U .and. v == SEQT%A)) then
                                    bp3_dH(h) = bp3_dH(h) + NNend_dH(NNENDT%GUonAU) - dHend0
                                    bp3_dS(h) = bp3_dS(h) + (NNend_dS(NNENDT%GUonAU) - dSend0) * 1.0e-3_PREC
                                 else if ((u == SEQT%G .and. v == SEQT%C) .or. (u == SEQT%C .and. v == SEQT%G)) then
                                    bp3_dH(h) = bp3_dH(h) + NNend_dH(NNENDT%GUonCG) - dHend0
                                    bp3_dS(h) = bp3_dS(h) + (NNend_dS(NNENDT%GUonCG) - dSend0) * 1.0e-3_PREC
                                 else if ((u == SEQT%G .and. v == SEQT%U) .or. (u == SEQT%U .and. v == SEQT%G)) then
                                    bp3_dH(h) = bp3_dH(h) + NNend_dH(NNENDT%GUonGU) - dHend0
                                    bp3_dS(h) = bp3_dS(h) + (NNend_dS(NNENDT%GUonGU) - dSend0) * 1.0e-3_PREC
                                 else
                                    print '(a)', 'logical error in init_bp'
                                 endif
                              endif
                           endif
                        endif ! flg_NNend

                     enddo
                  enddo
               enddo
            enddo
         enddo
      enddo

      !print *, h
      !! h must be 468
   endif


   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   ! Specific pairs given in CT file
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   if (flg_in_ct .or. flg_in_bpseq) then

      if (.not. flg_in_ct .and. .not. flg_in_bpseq) then
         print '(a)', 'Error: either .ct or .bpseq file is required for [Basepair] model = 2.'
         flush(output_unit)
         error stop
      endif

      call read_ss()   ! set bp_map according to ct/bpseq files

   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   ! All pairwise
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   else

      do imp = 1, nmp-1
         i = lmp_mp(imp)
         ichain = ichain_mp(imp)

         ! Either 5' or 3' end
         if (i == 1 .or. i == nmp_chain(ichain)) cycle

         do jmp = imp+1, nmp
            j = lmp_mp(jmp)
            jchain = ichain_mp(jmp)

            ! Either 5' or 3' end
            if (j == 1 .or. j == nmp_chain(jchain)) cycle

            ! Minimum loop length
            if (ichain == jchain .and. i + bp_min_loop >= j) cycle

            if (bp_model == 3 .or. bp_model == 4 .or. bp_model == 5) then
               ! Isolated base pair not allowed
               if (.not. is_complement(seq(i-1, ichain), seq(j+1, jchain)) .and. &
                   .not. is_complement(seq(i+1, ichain), seq(j-1, jchain)) ) then
                  cycle
               endif

               if ((i-1 == 1 .or. j+1 == nmp_chain(jchain)) .and. &
                  (.not. is_complement(seq(i+1, ichain), seq(j-1, jchain)))) then
                  cycle
               endif
               if ((j-1 == 1 .or. i+1 == nmp_chain(ichain)) .and. &
                  (.not. is_complement(seq(i-1, ichain), seq(j+1, jchain)))) then
                  cycle
               endif

            endif

            if (seq(i,ichain) == SEQT%G .and. seq(j, jchain) == SEQT%C) then
               bp_map(imp, jmp) = BPT%GC
               bp_map(jmp, imp) = BPT%GC
            else if (seq(i,ichain) == SEQT%C .and. seq(j, jchain) == SEQT%G) then
               bp_map(imp, jmp) = BPT%CG
               bp_map(jmp, imp) = BPT%CG

            else if (seq(i,ichain) == SEQT%A .and. seq(j, jchain) == SEQT%U) then
               bp_map(imp, jmp) = BPT%AU
               bp_map(jmp, imp) = BPT%AU
            else if (seq(i,ichain) == SEQT%U .and. seq(j, jchain) == SEQT%A) then
               bp_map(imp, jmp) = BPT%UA
               bp_map(jmp, imp) = BPT%UA

            else if (seq(i,ichain) == SEQT%G .and. seq(j, jchain) == SEQT%U) then
               bp_map(imp, jmp) = BPT%GU
               bp_map(jmp, imp) = BPT%GU
            else if (seq(i,ichain) == SEQT%U .and. seq(j, jchain) == SEQT%G) then
               bp_map(imp, jmp) = BPT%UG
               bp_map(jmp, imp) = BPT%UG
            
            else if (seq(i,ichain) == SEQT%A .and. seq(j, jchain) == SEQT%A) then
               bp_map(imp, jmp) = BPT%AA
               bp_map(jmp, imp) = BPT%AA
            else if (seq(i,ichain) == SEQT%A .and. seq(j, jchain) == SEQT%A) then
               bp_map(imp, jmp) = BPT%AA
               bp_map(jmp, imp) = BPT%AA
            endif

         enddo
      enddo
   endif

   if (bp_model == 5) then

      do imp = 1, nmp
         do jmp = imp+1, nmp

            if (bp_map(imp, jmp) == 0) cycle

            i = lmp_mp(imp)
            ichain = ichain_mp(imp)
            j = lmp_mp(jmp)
            jchain = ichain_mp(jmp)

            !if (ichain == jchain .and. (j+1) + bp_min_loop >= (i-1)) then
            !   ! (i-1) and (j+1) should be treated as unpaired (represented as U-U).
            !   bp3_map(imp, jmp) = bp3_hash(bp3_hash_key(SEQT%U, seq(i, ichain), seq(i+1, ichain), &
            !                                             SEQT%U, seq(j, jchain), seq(j-1, jchain)))

            if (i-1 == 1 .or. j+1 == nmp_chain(jchain)) then
               ! (i-1) and (j+1) should be treated as unpaired (represented as U-U).
               bp3_map(imp, jmp) = bp3_hash(bp3_hash_key(SEQT%U, seq(i, ichain), seq(i+1, ichain), &
                                                         SEQT%U, seq(j, jchain), seq(j-1, jchain)))

            else if (j-1 == 1 .or. i+1 == nmp_chain(ichain)) then
               bp3_map(imp, jmp) = bp3_hash(bp3_hash_key(seq(i-1, ichain), seq(i, ichain), SEQT%U, &
                                                         seq(j+1, jchain), seq(j, jchain), SEQT%U))

            else if (ichain == jchain .and. (i+1) + bp_min_loop >= (j-1)) then
               ! (i+1) and (j-1) should be treated as unpaired (represented as U-U).
               bp3_map(imp, jmp) = bp3_hash(bp3_hash_key(seq(i-1, ichain), seq(i, ichain), SEQT%U, &
                                                         seq(j+1, jchain), seq(j, jchain), SEQT%U))

            else
               bp3_map(imp, jmp) = bp3_hash(bp3_hash_key(seq(i-1, ichain), seq(i, ichain), seq(i+1, ichain), &
                                                         seq(j+1, jchain), seq(j, jchain), seq(j-1, jchain)))

            endif

         enddo
      enddo

   else
      bp3_map(:,:) = bp_map(:,:)

   endif

   !call set_bp_map()


   ! Calcuate BP cutoff
   ! If bp_cutoff_energy is not specified in ff, the default value is 0.001 (kcal/mol).
   if (abs(bp_cutoff_energy) <= epsilon(bp_cutoff_energy)) then
      ! When bp_cutoff_energy = 0.0, treat it as in the original way Hung did.
      bp_cutoff_dist = 18.0_PREC

      do bptype = 1, BPT%MAX
         bp_paras(bptype)%cutoff_ddist = bp_cutoff_dist - 13.8_PREC
      enddo

   else
      bp_cutoff_dist = 0.0_PREC

      do bptype = 1, BPT%MAX

         !bp_paras(bptype)%cutoff_ddist = sqrt(log(abs(bp_paras(bptype)%U0 / bp_cutoff_energy)) / bp_paras(bptype)%bond_k)
         ! Now U0 depends on the neighbouring nucleotides and the temperature.
         ! To be on the safe side, let us assume U0 = -10.0 kcal/mol.
         bp_paras(bptype)%cutoff_ddist = sqrt(log(abs(10.0_PREC / bp_cutoff_energy)) / bp_paras(bptype)%bond_k)

         ! To get the maximum (bond_r + cutoff_ddist) as the neighbour list cutoff
         cutoff = bp_paras(bptype)%bond_r + bp_paras(bptype)%cutoff_ddist
         if (cutoff > bp_cutoff_dist) then
            bp_cutoff_dist = cutoff
         endif
      enddo

   endif

   call write_bpcoef()

contains

   function bp3_hash_key(w, x, u, z, y, v) result (i)

      integer, intent(in) :: w, x, u, z, y, v
      integer :: a, b, c, d, e, f
      integer :: i

      if (x < y) then
         a = w
         b = x
         c = u
         d = z
         e = y
         f = v
      else
         f = w
         e = x
         d = u
         c = z
         b = y
         a = v
      endif

      !i = a*(4**5) + b*(4**4) + c*(4**3) + d*(4**2) + e*4 + f + 1
      i = 1 + a + b*4 + c*(4**2) + d*(4**3) + e*(4**4) + f*(4**5)
      !! Add 1 because the index starts with 1 but not 0.

   end function bp3_hash_key

endsubroutine init_bp
