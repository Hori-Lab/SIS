subroutine read_anneal(stat)

   use, intrinsic :: iso_fortran_env, Only : iostat_end

   use const
   use var_state, only : nanneal, anneal_step, anneal_tempK
   use var_io, only : cfile_anneal_in, iopen_hdl

   implicit none

   logical, intent(out) :: stat

   ! ---------------------------------------------------------------------
   integer :: istep
   integer :: istat, hdl
   integer :: s
   real(PREC) :: t

   stat = .False.

   write(*,*) 'Reading annealing-schedule file: '//trim(cfile_anneal_in)

   iopen_hdl = iopen_hdl + 1
   hdl = iopen_hdl

   open(hdl, file=cfile_anneal_in, status='old', action='read', iostat=istat)

   if (istat /= 0) then
      write(*,*) 'Error: failed to open the annealing-schedule file. '//trim(cfile_anneal_in)
      stop (2)
   endif

   ! Count the number of annealing steps
   nanneal = 0
   do
      read (hdl, *, iostat = istat) s, t

      if (istat < 0) then
         exit

      else if (istat > 0) then
         write(*,*) 'Error: cannot read the annealing-schedule file'
         return

      end if

      if (s >= 0 .and. t > 0.0) then
         nanneal = nanneal + 1

      else
         write(*,*) 'Error: cannot read the annealing-schedule file'
         return

      endif

   enddo

   allocate(anneal_tempK(nanneal))
   allocate(anneal_step(nanneal))

   ! Read the steps and temperatures
   rewind(hdl)
   istep = 0
   do
      read (hdl, *, iostat = istat) s, t

      if (istat < 0) then
         exit

      else if (istat > 0) then
         write(*,*) 'Error: cannot read the annealing-schedule file'
         return

      end if

      istep = istep + 1
      anneal_step(istep) = s
      anneal_tempK(istep) = t
   enddo

   close(hdl)
   iopen_hdl = iopen_hdl - 1
   
end subroutine read_anneal
