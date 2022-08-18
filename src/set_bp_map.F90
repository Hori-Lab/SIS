subroutine set_bp_map()

   use, intrinsic :: iso_fortran_env, Only : output_unit
   use const
   use const_idx, only : SEQT, BPT, seqt2char
   use var_io, only : flg_in_ct, flg_in_bpseq, cfile_ct_in, cfile_bpseq_in, iopen_hdl
   use var_top, only : nmp, seq, lmp_mp, ichain_mp, nmp_chain
   use var_potential, only : bp_model, bp_map, bp_min_loop, bp_map_dG

   implicit none

   integer :: imp, jmp
   integer :: i, j, ichain, jchain
   integer :: l, n, idummy
   integer :: istat, hdl
   real(PREC) :: dG
   real(PREC), allocatable :: NN_dG(:,:,:,:)
   character(len=1) :: nt
   character(len=5) :: cNN

   flush(output_unit)

   allocate(bp_map(nmp, nmp))
   bp_map(:,:) = 0

   if (bp_model == 4) then
      allocate(NN_dG(4,4,4,4))
      NN_dG(:,:,:,:) = 0.0_PREC

      iopen_hdl = iopen_hdl + 1
      hdl = iopen_hdl
      open(hdl, file='../NN.txt', status='old', action='read', iostat=istat)
      do
         read (hdl, *, iostat = istat) cNN, dG

         if (istat < 0) then
            exit
         else if (istat > 0) then
            error stop 'Error: cannot read NN file'
         end if

         call store_NN_dG(cNN, dG)
      enddo

      close(hdl)
      iopen_hdl = iopen_hdl - 1

      allocate(bp_map_dG(nmp, nmp))
      bp_map_dG(:,:) = 0.0_PREC
   endif

   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   ! All pairwise
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   if (bp_model == 1 .or. bp_model == 3 .or. bp_model == 4) then

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

            if (bp_model == 3 .or. bp_model == 4) then
               ! Isolated base pair not allowed
               if (.not. is_complement(seq(i-1, ichain), seq(j+1, jchain)) .and. &
                   .not. is_complement(seq(i+1, ichain), seq(j-1, jchain)) ) then
                  cycle
               endif
            endif

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

            if (bp_model == 4) then
               if (bp_map(imp, jmp) > 0) then
                   
                  dG = 0.0_PREC
                  if (i == 2 .or. j+1 == nmp_chain(jchain)) then
                     continue
                  else
                     dG = 0.5 * NN_dG(seq(i-1, ichain), seq(i, ichain), seq(j+1, jchain), seq(j, jchain))
                  endif

                  if (i+1 == nmp_chain(ichain) .or. j == 2) then
                     continue
                  else
                     dG = dG + 0.5 * NN_dG(seq(i, ichain), seq(i+1, ichain), seq(j, jchain), seq(j-1, jchain))
                  endif

                  bp_map_dG(imp, jmp) = dG
                  bp_map_dG(jmp, imp) = dG
               endif
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

   if (bp_model == 4) then
      deallocate(NN_dG)
   endif

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

   subroutine store_NN_dG(cNN, dG)

      character(len=5), intent(in) :: cNN
      real(PREC), intent(in) :: dG
      
      integer :: i
      integer :: s(4)

      do i = 1, 4
         if (cNN(i:i) == 'A') then
            s(i) = SEQT%A
         else if (cNN(i:i) == 'U') then
            s(i) = SEQT%U
         else if (cNN(i:i) == 'G') then
            s(i) = SEQT%G
         else if (cNN(i:i) == 'C') then
            s(i) = SEQT%C
         endif
      enddo

      NN_dG(s(1), s(2), s(3), s(4)) = dG
      NN_dG(s(4), s(3), s(2), s(1)) = dG
   
   endsubroutine store_NN_dG

endsubroutine set_bp_map
