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
      error stop 'Error: failed to open the XYZ file. '//trim(cfilepath)
   endif


   read (hdl, *, iostat = istat) n     ! Number of coordinates
   write(*,*) n
   if (istat /= 0) then
      error stop 'Error: in reading XYZ file (N) ' // trim(cfilepath)
   end if

   if (n /= nmp) then
      error stop "The number of coordinates in the XYZ file is not consitent (N != nmp). " // trim(cfilepath)
   endif

   read (hdl, *, iostat = istat)   ! Ignore a comment line
   if (istat /= 0) then
      error stop 'Error: in reading XYZ file (comment) ' // trim(cfilepath)
   end if

   do imp = 1, n
      read (hdl, *, iostat = istat) c10, x, y, z
         write(*,*) c10, x, y, z

      if (istat /= 0) then
         error stop 'Error: in reading XYZ file (coordinates) ' // trim(cfilepath)
      end if

      xyz(1, imp) = x
      xyz(2, imp) = y
      xyz(3, imp) = z
   end do

   close(hdl)
   iopen_hdl = iopen_hdl - 1

end subroutine read_xyz
