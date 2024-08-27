subroutine read_ss()

   use, intrinsic :: iso_fortran_env, Only : OUTPUT_UNIT

   use const_idx, only : SEQT, BPT, seqt2char, is_complement
   use var_top, only : nmp, nmp_chain, ichain_mp, seq, lmp_mp
   use var_io, only : flg_in_ct, flg_in_bpseq, flg_in_bpl, &
                      cfile_ct_in, cfile_bpseq_in, cfile_bpl_in, iopen_hdl
   use var_potential, only : bp_model, bp_map, bp_min_loop
   use var_parallel

   implicit none

   ! ---------------------------------------------------------------------
   integer :: i, j, n, l, imp, jmp, ichain, jchain, nline
   integer :: istat, hdl
   integer :: idummy

   character(len=1) :: nt

   if (.not. (flg_in_ct .or. flg_in_bpseq .or. flg_in_bpl)) then
      return
   endif

   if (myrank == 0) then

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

         nline = nmp

      else if (flg_in_bpl) then ! base pair list
         print '(2a)', "Reading BPL file: ", trim(cfile_bpl_in)
         open(hdl, file=cfile_bpl_in, status='old', action='read', iostat=istat)

         if (istat /= 0) then
            print '(2a)', 'Error: failed to open the BPL file. ', trim(cfile_bpl_in)
            flush(output_unit)
            error stop
         endif

         n = 0
         do
            read(hdl, *, iostat=istat) imp, jmp
            if (istat < 0) then
               exit
            else if (istat > 0) then
               print '(a)', 'Error: failed to read the BPL file. '//trim(cfile_bpl_in)
               call sis_abort()
            end if
            n = n + 1
         enddo
         rewind(hdl)
         
         nline = n

      else ! BPSEQ
         print '(2a)', "Reading BPSEQ file: ", trim(cfile_bpseq_in)
         open(hdl, file=cfile_bpseq_in, status='old', action='read', iostat=istat)

         if (istat /= 0) then
            print '(2a)', 'Error: failed to open the BPSEQ file. ', trim(cfile_bpseq_in)
            flush(output_unit)
            error stop
         endif

         nline = nmp

      endif

      ! Main
      do l = 1, nline

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

         else if (flg_in_bpl) then
            read(hdl, *, iostat=istat) imp, jmp

            if (istat /= 0) then
               print '(a,i8)', 'Error: BPL file format error. Line can not be read for Nucleotide ', l
               flush(output_unit)
               error stop
            end if

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
               print '(a,i5,a,i3,a,i6,a,a)', '         Nucleotide i ',  i, ' of chain ', ichain, ' (imp= ', imp, ')  - ', seqt2char(seq(i, ichain))
               print '(a,i5,a,i3,a,i6,a,a)', '         Nucleotide j ',  j, ' of chain ', jchain, ' (imp= ', jmp, ')  - ', seqt2char(seq(j, jchain))
               cycle
            endif

            ! Either 5' or 3' end
            if (j == 1 .or. j == nmp_chain(jchain)) then
               print '(a)', 'Warning: The following pair in CT/BPSEQ file will not be considered because (j) is a chain end.'
               print '(a,i5,a,i3,a,i6,a,a)', '         Nucleotide i ',  i, ' of chain ', ichain, ' (imp= ', imp, ')  - ', seqt2char(seq(i, ichain))
               print '(a,i5,a,i3,a,i6,a,a)', '         Nucleotide j ',  j, ' of chain ', jchain, ' (imp= ', jmp, ')  - ', seqt2char(seq(j, jchain))
               cycle
            endif

            ! Minimum loop length
            if (ichain == jchain .and. i + bp_min_loop >= j) then
               print '(a)', 'Warning: The following pair in CT/BPSEQ file will not be considered due to the minimum loop length required.'
               print '(a,i5,a,i3,a,i6,a,a)', '         Nucleotide i ',  i, ' of chain ', ichain, ' (imp= ', imp, ')  - ', seqt2char(seq(i, ichain))
               print '(a,i5,a,i3,a,i6,a,a)', '         Nucleotide j ',  j, ' of chain ', jchain, ' (imp= ', jmp, ')  - ', seqt2char(seq(j, jchain))
               cycle
            endif

            if (bp_model == 3 .or. bp_model == 4 .or. bp_model == 5) then
               ! Isolated base pair not allowed
               if (.not. is_complement(seq(i-1, ichain), seq(j+1, jchain)) .and. &
                   .not. is_complement(seq(i+1, ichain), seq(j-1, jchain)) ) then
                  print '(a)', 'Warning: The following pair in CT/BPSEQ file will not be considered because it is an isolated base pair.'
                  print '(a,i5,a,i3,a,i6,a,a)', '         Nucleotide i ',  i, ' of chain ', ichain, ' (imp= ', imp, ') - ', seqt2char(seq(i, ichain))
                  print '(a,i5,a,i3,a,i6,a,a)', '         Nucleotide j ',  j, ' of chain ', jchain, ' (imp= ', jmp, ') - ', seqt2char(seq(j, jchain))
                  cycle
               endif
            endif

            if (seq(i, ichain) == SEQT%G .and. seq(j, jchain) == SEQT%C) then
               bp_map(imp, jmp) = BPT%GC
               bp_map(jmp, imp) = BPT%GC
            else if (seq(i, ichain) == SEQT%C .and. seq(j, jchain) == SEQT%G) then
               bp_map(imp, jmp) = BPT%CG
               bp_map(jmp, imp) = BPT%CG

            else if (seq(i, ichain) == SEQT%A .and. seq(j, jchain) == SEQT%U) then
               bp_map(imp, jmp) = BPT%AU
               bp_map(jmp, imp) = BPT%AU
            else if (seq(i, ichain) == SEQT%U .and. seq(j, jchain) == SEQT%A) then
               bp_map(imp, jmp) = BPT%UA
               bp_map(jmp, imp) = BPT%UA

            else if (seq(i, ichain) == SEQT%G .and. seq(j, jchain) == SEQT%U) then
               bp_map(imp, jmp) = BPT%GU
               bp_map(jmp, imp) = BPT%GU
            else if (seq(i, ichain) == SEQT%U .and. seq(j, jchain) == SEQT%G) then
               bp_map(imp, jmp) = BPT%UG
               bp_map(jmp, imp) = BPT%UG

            else
               print '(a)', 'Warning: The following pair in CT/BPSEQ file does not form any known types of base pairs.'
               print '(a,i5,a,i3,a,i6,a,a)', '         Nucleotide i ',  i, ' of chain ', ichain, ' (imp= ', imp, ') - ', seqt2char(seq(i, ichain))
               print '(a,i5,a,i3,a,i6,a,a)', '         Nucleotide j ',  j, ' of chain ', jchain, ' (imp= ', jmp, ') - ', seqt2char(seq(j, jchain))

            endif

         endif
      enddo

      close(hdl)
      iopen_hdl = iopen_hdl - 1

      print '(a)', 'Done: reading CT/BPSEQ/BPL file'
      print *
      flush(output_unit)
   endif

#ifdef PAR_MPI

   if (myrank == 0) then
      print '(a)', 'Sending secondary structure (bp_map) data via MPI.'
   else
      print '(a)', 'Receiving secondary structure (bp_map) data via MPI.'
   endif
   flush(OUTPUT_UNIT)

   call MPI_BCAST(bp_map, nmp*nmp, MPI_INTEGER, 0, MPI_COMM_WORLD, istat)

   if (myrank == 0) then
      print '(a)', 'Done: sending secondary structure (bp_map) data via MPI.'
   else
      print '(a)', 'Done: receiving secondary structure (bp_map) data via MPI.'
   endif
   print *
   flush(OUTPUT_UNIT)

#endif
   
end subroutine read_ss
