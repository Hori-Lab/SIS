subroutine read_fasta()

   use, intrinsic :: iso_fortran_env, Only : iostat_end

   use const
   use const_idx, only : SEQT, char2seqt
   use var_top, only : nmp, nchains, nmp_chain, seq, imp_chain, ichain_mp, lmp_mp
   use var_io, only : cfile_fasta_in, iopen_hdl

   implicit none

   ! ---------------------------------------------------------------------
   integer :: i, n, imp, iseq
   integer :: ichain
   integer :: istat, hdl, nlen
   logical :: flg_reading

   character(len=CHAR_FILE_LINE) :: cline


   iopen_hdl = iopen_hdl + 1
   hdl = iopen_hdl

   write(6,'(a)') "Reading FASTA file: " // trim(cfile_fasta_in)
   open(hdl, file=cfile_fasta_in, status='old', action='read', iostat=istat)

   if (istat /= 0) then
      write(6, '(a)') 'Error: failed to open the FASTA file. '//trim(cfile_fasta_in)
      stop (2)
   endif

   ! Count the number of chains
   nchains = 0
   do
      read (hdl, '(A)', iostat = istat) cline

      if (istat < 0) then
         exit
      else if (istat > 0) then
         error stop 'Error: cannot read fasta file'
      end if

      if (cline(1:1) == '>') then
         nchains = nchains + 1
      endif

   enddo

   allocate(nmp_chain(nchains))

   ! Count the number of nucleotides
   rewind(hdl)
   ichain = 0
   flg_reading = .False.
   do
      read (hdl, '(A)', iostat = istat) cline

      if (istat < 0) then
         exit
      else if (istat > 0) then
         error stop 'Error: cannot read fasta file'
      end if

      if (cline(1:1) == '>') then
         
         if (flg_reading) then
            nmp_chain(ichain) = n
         endif
            
         ichain = ichain + 1
         flg_reading = .True.
         n = 0

      else if (flg_reading) then

         nlen = len(trim(cline))
         do i = 1, nlen
             if (char2seqt(cline(i:i)) /= SEQT%UNDEF) then
                n = n + 1
             endif
         enddo

      endif

   enddo
   
   if (flg_reading) then
      nmp_chain(ichain) = n
   endif

   nmp = sum(nmp_chain)

   allocate(seq(maxval(nmp_chain), nchains))
   allocate(imp_chain(maxval(nmp_chain), nchains))
   allocate(ichain_mp(nmp))
   allocate(lmp_mp(nmp))

   seq(:,:) = SEQT%UNDEF
   imp_chain(:,:) = 0
   ichain_mp(:) = 0
   lmp_mp(:) = 0

   ! Read sequences
   rewind(hdl)
   ichain = 0
   imp = 0
   flg_reading = .False.
   do
      read (hdl, '(A)', iostat = istat) cline

      if (istat < 0) then
         exit
      else if (istat > 0) then
         error stop 'Error: cannot read fasta file'
      end if

      if (cline(1:1) == '>') then
            
         ichain = ichain + 1
         flg_reading = .True.
         n = 0

      else if (flg_reading) then

         nlen = len(trim(cline))
         do i = 1, nlen
             iseq = char2seqt(cline(i:i))
             if (iseq /= SEQT%UNDEF) then
                imp = imp + 1
                n = n + 1
                seq(n, ichain) = iseq
                imp_chain(n, ichain) = imp
                ichain_mp(imp) = ichain
                lmp_mp(imp) = n
             endif
         enddo

      endif

   enddo

   close(hdl)
   iopen_hdl = iopen_hdl - 1

   write(6, '(a)') 'Done: reading FASTA file'
   write(6, *)
   flush(6)
   
end subroutine read_fasta
