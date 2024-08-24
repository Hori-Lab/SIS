subroutine read_restraint()

   use, intrinsic :: ISO_FORTRAN_ENV, only: OUTPUT_UNIT

   use const, only : PREC, CHAR_FILE_LINE
   use var_potential, only : flg_rest_sigb, nrest_sigb, &
                             rest_sigb_id, rest_sigb_rcut, rest_sigb_para
   use var_io, only : cfile_rest_in, iopen_hdl
   use var_parallel

   implicit none

   ! ---------------------------------------------------------------------
   integer :: irest_sigb
   integer :: istat, hdl
   integer :: i, imp, jmp
   real(PREC) :: eps, d, s, rcut
   character(CHAR_FILE_LINE) :: cline

   if (myrank == 0) then
      print '(a)', 'Reading restraint file: '//trim(cfile_rest_in)

      iopen_hdl = iopen_hdl + 1
      hdl = iopen_hdl

      open(hdl, file=cfile_rest_in, status='old', action='read', iostat=istat)

      if (istat /= 0) then
         print *, 'Error: failed to open the restraint file. '//trim(cfile_rest_in)
         call sis_abort()
      endif

      flg_rest_sigb = .False.
      nrest_sigb = 0

      ! Count the number of restraint
      do
         read (hdl, '(a)', iostat=istat) cline

         if (istat < 0) then
            exit

         else if (istat > 0) then
            print *, 'Error: cannot read the restraint file. ' // trim(cfile_rest_in)
            call sis_abort()

         else if (cline(1:1) == '#') then
            cycle

         else if (len(trim(cline)) == 0) then
            cycle
         end if

         if (cline(1:15) == 'Sigmoid-to-bead') then
            nrest_sigb = nrest_sigb + 1

         else
            print *, 'Error: Unknown restraint type in the restraint file'
            print *, trim(cline)
            call sis_abort()

         endif

      enddo

      if (nrest_sigb > 0) then
         flg_rest_sigb = .True.

         allocate(rest_sigb_id(2, nrest_sigb))
         allocate(rest_sigb_rcut(nrest_sigb))
         allocate(rest_sigb_para(3, nrest_sigb))
         rest_sigb_id(:,:) = 0
         rest_sigb_rcut(:) = 0.0_PREC
         rest_sigb_para(:,:) = 0.0_PREC
      endif

      ! Read the steps and temperatures
      rewind(hdl)
      irest_sigb = 0
      do
         read (hdl, '(a)', iostat=istat) cline

         if (istat < 0) then
            exit

         else if (istat > 0) then
            print *, 'Error: cannot read the restraint file. ' // trim(cfile_rest_in)
            call sis_abort()

         else if (cline(1:1) == '#') then
            cycle

         else if (len(trim(cline)) == 0) then
            cycle
         end if

         if (cline(1:15) == 'Sigmoid-to-bead') then
            irest_sigb = irest_sigb + 1

            read(cline(16:), *) imp, jmp, eps, d, s, rcut
            rest_sigb_id(1, irest_sigb) = imp 
            rest_sigb_id(2, irest_sigb) = jmp 
            rest_sigb_para(1, irest_sigb) = eps
            rest_sigb_para(2, irest_sigb) = d
            rest_sigb_para(3, irest_sigb) = s
            rest_sigb_rcut(irest_sigb) = rcut

         else
            print *, 'Error: Unknown restraint type in the restraint file'
            print *, trim(cline)
            call sis_abort()

         endif
      enddo

      close(hdl)
      iopen_hdl = iopen_hdl - 1
   endif  ! myrank = 0

#ifdef PAR_MPI

   if (myrank == 0) then
      print '(a)', 'Sending restraint data via MPI.'
      flush(OUTPUT_UNIT)
   else
      print '(a)', 'Receiving restraint data via MPI.'
      flush(OUTPUT_UNIT)
   endif

   call MPI_BCAST(flg_rest_sigb, 1, MPI_LOGICAL, 0, MPI_COMM_WORLD, istat)
   call MPI_BCAST(nrest_sigb, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, istat)

   if (flg_rest_sigb) then
      if (myrank /= 0) then
         allocate(rest_sigb_id(2, nrest_sigb))
         allocate(rest_sigb_para(3, nrest_sigb))
         allocate(rest_sigb_rcut(nrest_sigb))
         rest_sigb_id(:,:) = 0
         rest_sigb_para(:,:) = 0.0_PREC
         rest_sigb_rcut(:) = 0.0_PREC
      endif

      call MPI_BCAST(rest_sigb_id, 2*nrest_sigb, MPI_INTEGER, 0, MPI_COMM_WORLD, istat)
      call MPI_BCAST(rest_sigb_para, 3*nrest_sigb, PREC_MPI, 0, MPI_COMM_WORLD, istat)
      call MPI_BCAST(rest_sigb_rcut, nrest_sigb, PREC_MPI, 0, MPI_COMM_WORLD, istat)
   endif

#endif

   if (flg_rest_sigb) then
      print '(a)', '# Sigmoid-to-bead:   ID   subject reference  epsilon     d        s     r_cut'
      do irest_sigb = 1, nrest_sigb
         print '(a, i4, 2(1x,i8), 4(1x,f8.2))', '# Sigmoid-to-bead: ', &
                    irest_sigb, (rest_sigb_id(i, irest_sigb), i=1,2), &
                    (rest_sigb_para(i, irest_sigb), i=1,3), &
                    rest_sigb_rcut(irest_sigb)
      enddo
   endif

   print *
   flush(OUTPUT_UNIT)

end subroutine read_restraint
