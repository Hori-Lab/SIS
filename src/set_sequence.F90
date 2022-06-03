subroutine set_sequence()

   use const_idx, only : SEQT, seqt2char
   use var_io, only : cfile_fasta_in
   use var_top, only : nmp, nchains, nmp_chain, seq, imp_chain, ichain_mp, nrepeat, lmp_mp

   implicit none

   integer :: i, j, k, imp

   if (allocated(cfile_fasta_in)) then

      call read_fasta()

   else if (nrepeat > 0) then

      allocate(nmp_chain(nchains))
      nmp_chain(:) = 3 * nrepeat
      nmp = sum(nmp_chain)
      write(*, '(a,i5)') '#Nrepeat: ', nrepeat
      allocate(seq(3*nrepeat, nchains))
      allocate(imp_chain(3*nrepeat, nchains))
      allocate(ichain_mp(nmp))
      allocate(lmp_mp(nmp))
      imp = 0
      do i = 1, nchains
         do j = 1, nrepeat
            seq(3*(j-1)+1, i) = SEQT%C
            seq(3*(j-1)+2, i) = SEQT%A
            seq(3*(j-1)+3, i) = SEQT%G
            imp_chain(3*(j-1)+1, i) = imp+1
            imp_chain(3*(j-1)+2, i) = imp+2
            imp_chain(3*(j-1)+3, i) = imp+3
            ichain_mp(imp+1) = i
            ichain_mp(imp+2) = i
            ichain_mp(imp+3) = i
            lmp_mp(imp+1) = 3*(j-1)+1
            lmp_mp(imp+2) = 3*(j-1)+2
            lmp_mp(imp+3) = 3*(j-1)+3
            imp = imp + 3
         enddo
         !write(*,'(a,141(i1))') '# ', seq(:,i)
         !write(*,*) '# ', imp_chain(1,i), imp_chain((nrepeat-1)*3+3, i)
      enddo

   else
      print *,'Error: either FASTA or [repeat] is required.'
      stop (2)

   endif

   print '(a)', '############ System ############'
   print '(a,i8)', 'Nchain: ', nchains
   print *
   do i = 1, nchains
      print '(a, i4)', 'Chain ', i
      print '(a, i10)', 'Nnt: ', nmp_chain(i)
      k = 0
      do j = 1, nmp_chain(i)
         write(6, '(a)', advance='no') seqt2char(seq(j,i))
         k = k + 1
         if (mod(k,100) == 0) then
            write(6, *) ''
            k = 0
         endif
      enddo
      if (k /= 0) then
         write(6, *) ''
      endif
   enddo
   print '(a)', '################################'
   print *

endsubroutine set_sequence
