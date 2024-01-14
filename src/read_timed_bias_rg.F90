subroutine read_timed_bias_rg()

   use, intrinsic :: ISO_FORTRAN_ENV, only: INT64, OUTPUT_UNIT

   use const
   use var_potential, only : ntimed_bias_rg, timed_bias_rg_step, timed_bias_rg_k, timed_bias_rg_0
   use var_io, only : cfile_bias_rg_in, iopen_hdl

   implicit none

   ! ---------------------------------------------------------------------
   integer :: istep
   integer :: istat, hdl
   integer(INT64) :: s
   real(PREC) :: k ,r0

   print '(a)', 'Reading bias-Rg-schedule file: '//trim(cfile_bias_rg_in)

   iopen_hdl = iopen_hdl + 1
   hdl = iopen_hdl

   open(hdl, file=cfile_bias_rg_in, status='old', action='read', iostat=istat)

   if (istat /= 0) then
      print *, 'Error: failed to open the bias-Rg-schedule file. '//trim(cfile_bias_rg_in)
      call sis_abort()
   endif

   ! Count the number of schedules
   ntimed_bias_rg = 0
   do
      read (hdl, *, iostat = istat) s, k, r0

      if (istat < 0) then
         exit

      else if (istat > 0) then
         print *, 'Error: cannot read the bias-Rg-schedule file'
         call sis_abort()

      end if

      if (s >= 0 .and. k >= 0 .and. r0 > 0.0) then
         ntimed_bias_rg = ntimed_bias_rg + 1

      else
         print *, 'Error: cannot read the bias-Rg-schedule file'
         call sis_abort()

      endif

   enddo

   allocate(timed_bias_rg_step(ntimed_bias_rg))
   allocate(timed_bias_rg_k(ntimed_bias_rg))
   allocate(timed_bias_rg_0(ntimed_bias_rg))

   ! Read the steps and temperatures
   rewind(hdl)
   istep = 0
   do
      read (hdl, *, iostat = istat) s, k, r0

      if (istat < 0) then
         exit

      else if (istat > 0) then
         print *, 'Error: cannot read the bias-Rg-schedule file'
         call sis_abort()

      end if

      istep = istep + 1
      timed_bias_rg_step(istep) = s
      timed_bias_rg_k(istep) = k
      timed_bias_rg_0(istep) = r0
   enddo

   close(hdl)
   iopen_hdl = iopen_hdl - 1

   do istep = 1, ntimed_bias_rg
      print '(a,i4,1x,i12,1x,f8.2,1x,f8.2)', '# Timed-bias-Rg: ', &
            istep, timed_bias_rg_step(istep), timed_bias_rg_k(istep), timed_bias_rg_0(istep)
   enddo

   print *
   flush(OUTPUT_UNIT)

end subroutine read_timed_bias_rg
