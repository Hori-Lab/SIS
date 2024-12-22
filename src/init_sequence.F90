subroutine init_sequence()

   use const_idx, only : SEQT, seqt2char, MOLT, molt2char
   use var_io, only : flg_in_fasta
   use var_top, only : nmp, nchains, nmp_chain, seq, imp_chain, ichain_mp, lmp_mp, moltypes
   use var_parallel

   implicit none

   integer :: i, j, k, imp
   integer :: inp_nchains

   ! If [Molecules] is written in the input toml
   if (allocated(moltypes)) then
      inp_nchains = nchains
   endif

   if (flg_in_fasta) then

      call read_fasta()

   else
      print *,'Error: FASTA file is required.'
      call sis_abort()
   endif

   if (allocated(moltypes)) then
      if (inp_nchains /= nchains) then
         print *,'Error: Number of chains in the FASTA file is not consistent with [Molecues] n_chain in the input file.'
         call sis_abort()
      endif
   else
      allocate(moltypes(nchains))
      moltypes(:) = MOLT%RNA
   endif

   print '(a)', '############ System ############'
   print '(a,i8)', 'Nchain: ', nchains
   print *
   do i = 1, nchains
      print '(a, i0)', 'Chain ', i
      print '(a, a)',  'Molecule: ', molt2char(moltypes(i))
      print '(a, i0)', 'Nnt: ', nmp_chain(i)
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
      write(6, *) ''
   enddo
   print '(a)', '################################'
   print *

endsubroutine init_sequence
