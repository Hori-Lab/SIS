subroutine read_restraint()

   use, intrinsic :: ISO_FORTRAN_ENV, only: OUTPUT_UNIT

   use const, only : PREC, CHAR_FILE_LINE
   use var_potential, only : nrest_sigb, rest_sigb_id, rest_sigb_para
   use var_io, only : cfile_rest_in, iopen_hdl

   implicit none

   ! ---------------------------------------------------------------------
   integer :: irest_sigb
   integer :: istat, hdl
   integer :: i, imp, jmp
   real(PREC) :: eps, rcut, delta
   character(CHAR_FILE_LINE) :: cline

   print '(a)', 'Reading restraint file: '//trim(cfile_rest_in)

   iopen_hdl = iopen_hdl + 1
   hdl = iopen_hdl

   open(hdl, file=cfile_rest_in, status='old', action='read', iostat=istat)

   if (istat /= 0) then
      print *, 'Error: failed to open the restraint file. '//trim(cfile_rest_in)
      call sis_abort()
   endif

   ! Count the number of restraint
   nrest_sigb = 0
   do
      read (hdl, '(a)', iostat=istat) cline

      if (istat < 0) then
         exit

      else if (istat > 0) then
         print *, 'Error: cannot read the restraint file. ' // trim(cfile_rest_in)
         call sis_abort()

      else if (cline(1:1) == '#') then
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

   allocate(rest_sigb_id(2, nrest_sigb))
   allocate(rest_sigb_para(3, nrest_sigb))

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
      end if

      if (cline(1:15) == 'Sigmoid-to-bead') then
         irest_sigb = irest_sigb + 1
         
         read(cline(16:), *) imp, jmp, eps, rcut, delta
         rest_sigb_id(1, irest_sigb) = imp 
         rest_sigb_id(2, irest_sigb) = jmp 
         rest_sigb_para(1, irest_sigb) = eps
         rest_sigb_para(2, irest_sigb) = rcut
         rest_sigb_para(3, irest_sigb) = delta

      else
         print *, 'Error: Unknown restraint type in the restraint file'
         print *, trim(cline)
         call sis_abort()

      endif
   enddo

   close(hdl)
   iopen_hdl = iopen_hdl - 1

   if (nrest_sigb > 0) then
      print '(a)', '# Sigmoid-to-bead:   ID   subject reference  epsilon  r_cut    delta'
      do irest_sigb = 1, nrest_sigb
         print '(a, i4, 2(1x,i8), 3(1x,f8.2))', '# Sigmoid-to-bead: ', &
               irest_sigb, (rest_sigb_id(i, irest_sigb), i=1,2), (rest_sigb_para(i, irest_sigb), i=1,3)
      enddo
   endif

   print *
   flush(OUTPUT_UNIT)

end subroutine read_restraint
