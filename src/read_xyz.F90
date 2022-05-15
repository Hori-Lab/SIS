subroutine read_xyz(cfilepath, nmp, xyz)

   use const
   use var_io, only : iopen_hdl

   implicit none

   character(len=*), intent(in) :: cfilepath
   integer, intent(in) :: nmp
   real(PREC), intent(out) :: xyz(3, nmp)

   ! ---------------------------------------------------------------------
   integer :: imp, n
   real(PREC) :: x, y, z

   integer :: istat, hdl
   character(10) :: c10

   imp = 0

   iopen_hdl = iopen_hdl + 1
   hdl = iopen_hdl

   open(hdl, file=cfilepath, status='old', action='read', iostat=istat)

   if (istat /= 0) then
      print '(2a)', 'Error: failed to open the XYZ file. ', trim(cfilepath)
      error stop
   endif


   read (hdl, *, iostat = istat) n     ! Number of coordinates
   if (istat /= 0) then
      print '(2a)', 'Error: in reading XYZ file (N) ', trim(cfilepath)
      error stop
   end if

   if (n /= nmp) then
      print '(2a)', "The number of coordinates in the XYZ file is not consitent (N != nmp). ", trim(cfilepath)
      error stop
   endif

   read (hdl, *, iostat = istat)   ! Ignore a comment line
   if (istat /= 0) then
      print '(2a)', 'Error: in reading XYZ file (comment) ', trim(cfilepath)
      error stop
   end if

   do imp = 1, n
      read (hdl, *, iostat = istat) c10, x, y, z

      if (istat /= 0) then
         print '(2a)', 'Error: in reading XYZ file (coordinates) ', trim(cfilepath)
         error stop
      end if

      xyz(1, imp) = x
      xyz(2, imp) = y
      xyz(3, imp) = z
   end do

   close(hdl)
   iopen_hdl = iopen_hdl - 1

end subroutine read_xyz
