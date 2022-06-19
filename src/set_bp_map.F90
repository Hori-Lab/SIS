subroutine set_bp_map()

   use, intrinsic :: iso_fortran_env, Only : output_unit
   use const
   use const_idx, only : SEQT, BPT, seqt2char
   use var_io, only : flg_in_ct, flg_in_bpseq, cfile_ct_in, cfile_bpseq_in, iopen_hdl
   use var_top, only : nmp, seq, lmp_mp, ichain_mp, nmp_chain
   use var_potential, only : bp_model, bp_map, bp_min_loop

   implicit none

   integer :: imp, jmp
   integer :: i, j, ichain, jchain
   integer :: l, n, idummy
   integer :: istat, hdl
   character(1) :: nt

   allocate(bp_map(nmp, nmp))
   bp_map(:,:) = 0

   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   ! All pairwise
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   if (bp_model == 1) then

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

            if ((seq(i,ichain) == SEQT%G .and. seq(j, jchain) == SEQT%C) .or. &
                (seq(i,ichain) == SEQT%C .and. seq(j, jchain) == SEQT%G) ) then
               bp_map(imp, jmp) = BPT%GC
               bp_map(jmp, imp) = BPT%GC

            else if ((seq(i,ichain) == SEQT%A .and. seq(j, jchain) == SEQT%U) .or. &
                     (seq(i,ichain) == SEQT%U .and. seq(j, jchain) == SEQT%A) ) then
               bp_map(imp, jmp) = BPT%AU
               bp_map(jmp, imp) = BPT%AU

            else if ((seq(i,ichain) == SEQT%G .and. seq(j, jchain) == SEQT%U) .or. &
                     (seq(i,ichain) == SEQT%U .and. seq(j, jchain) == SEQT%G) ) then
               bp_map(imp, jmp) = BPT%GU
               bp_map(jmp, imp) = BPT%GU
            endif
         enddo
      enddo

   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   ! Specific pairs given in CT file
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   else if (bp_model == 2) then

      if (.not. flg_in_ct .and. .not. flg_in_bpseq) then
         print '(a)', 'Error: either .ct or .bpseq file is required for [Basepair] model = 2.'
         flush(output_unit)
         error stop
      endif

      iopen_hdl = iopen_hdl + 1
      hdl = iopen_hdl

      if (flg_in_ct) then
         print '(2a)', "Reading CT file: ", trim(cfile_ct_in)
         open(hdl, file=cfile_ct_in, status='old', action='read', iostat=istat)

         if (istat /= 0) then
            print '(2a)', 'Error: failed to open the CT file. ', trim(cfile_ct_in)
            flush(output_unit)
            error stop
         endif

         ! Header, number of nucleotides
         read(hdl, *) n

         if (n /= nmp) then
            print '(a)', 'Error: CT file format error. The total number of nucleotide does not match.'
            flush(output_unit)
            error stop
         endif

      else ! BPSEQ
         print '(2a)', "Reading BPSEQ file: ", trim(cfile_bpseq_in)
         open(hdl, file=cfile_bpseq_in, status='old', action='read', iostat=istat)

         if (istat /= 0) then
            print '(2a)', 'Error: failed to open the BPSEQ file. ', trim(cfile_bpseq_in)
            flush(output_unit)
            error stop
         endif
      endif

      ! Main
      do l = 1, nmp

         if (flg_in_ct) then
            ! e.g.)   2 G     1     3   218     5
            read(hdl, *, iostat=istat) imp, nt, idummy, idummy, jmp, idummy

            if (istat /= 0) then
               print '(a,i8)', 'Error: CT file format error. Line can not be read for Nucleotide ', l
               flush(output_unit)
               error stop
            end if

            if (imp /= l) then
               print '(a,i8)', 'Error: CT file format error. Nucleotide number does not match for Nucleotide ', l
               flush(output_unit)
               error stop
            endif

         else ! BPSEQ
            ! e.g.)   2 G     218
            read(hdl, *, iostat=istat) imp, nt, jmp

            if (istat /= 0) then
               print '(a,i8)', 'Error: BPSEQ file format error. Line can not be read for Nucleotide ', l
               flush(output_unit)
               error stop
            end if

            if (imp /= l) then
               print '(a,i8)', 'Error: BPSEQ file format error. Nucleotide number does not match for Nucleotide ', l
               flush(output_unit)
               error stop
            endif
         endif

         if (jmp /= 0) then

            if (jmp < imp) then
               idummy = jmp
               jmp = imp
               imp = idummy
            endif

            i = lmp_mp(imp)
            ichain = ichain_mp(imp)
            j = lmp_mp(jmp)
            jchain = ichain_mp(jmp)

            ! Either 5' or 3' end
            if (i == 1 .or. i == nmp_chain(ichain)) then
               print '(a)', 'Warning: The following pair in CT/BPSEQ file will not be considered because (i) is a chain end.'
               print '(a,i5,a,i3,a,a)', '         Nucleotide i ',  i, ' of chain ', ichain, ' - ', seqt2char(seq(i, ichain))
               print '(a,i5,a,i3,a,a)', '         Nucleotide j ',  j, ' of chain ', jchain, ' - ', seqt2char(seq(j, jchain))
               cycle
            endif

            ! Either 5' or 3' end
            if (j == 1 .or. j == nmp_chain(jchain)) then
               print '(a)', 'Warning: The following pair in CT/BPSEQ file will not be considered because (j) is a chain end.'
               print '(a,i5,a,i3,a,a)', '         Nucleotide i ',  i, ' of chain ', ichain, ' - ', seqt2char(seq(i, ichain))
               print '(a,i5,a,i3,a,a)', '         Nucleotide j ',  j, ' of chain ', jchain, ' - ', seqt2char(seq(j, jchain))
               cycle
            endif

            ! Minimum loop length
            if (ichain == jchain .and. i + bp_min_loop >= j) then
               print '(a)', 'Warning: The following pair in CT/BPSEQ file will not be considered due to the minimum loop length required.'
               print '(a,i5,a,i3,a,a)', '         Nucleotide i ',  i, ' of chain ', ichain, ' - ', seqt2char(seq(i, ichain))
               print '(a,i5,a,i3,a,a)', '         Nucleotide j ',  j, ' of chain ', jchain, ' - ', seqt2char(seq(j, jchain))
               cycle
            endif

            if ((seq(i, ichain) == SEQT%G .and. seq(j, jchain) == SEQT%C) .or. &
                (seq(i, ichain) == SEQT%C .and. seq(j, jchain) == SEQT%G)) then
               bp_map(imp, jmp) = BPT%GC
               bp_map(jmp, imp) = BPT%GC

            else if ((seq(i, ichain) == SEQT%A .and. seq(j, jchain) == SEQT%U) .or. &
                     (seq(i, ichain) == SEQT%U .and. seq(j, jchain) == SEQT%A)) then
               bp_map(imp, jmp) = BPT%AU
               bp_map(jmp, imp) = BPT%AU

            else if ((seq(i, ichain) == SEQT%G .and. seq(j, jchain) == SEQT%U) .or. &
                     (seq(i, ichain) == SEQT%U .and. seq(j, jchain) == SEQT%G)) then
               bp_map(imp, jmp) = BPT%GU
               bp_map(jmp, imp) = BPT%GU

            else
               print '(a)', 'Warning: The following pair in CT/BPSEQ file does not form any known types of base pairs.'
               print '(a,i5,a,i3,a,a)', '         Nucleotide i ',  i, ' of chain ', ichain, ' - ', seqt2char(seq(i, ichain))
               print '(a,i5,a,i3,a,a)', '         Nucleotide j ',  j, ' of chain ', jchain, ' - ', seqt2char(seq(j, jchain))

            endif
         endif
      enddo

      close(hdl)
      iopen_hdl = iopen_hdl - 1

      print '(a)', 'Done: reading CT/BPSEQ file'
      print *
      flush(output_unit)
   endif

endsubroutine set_bp_map
