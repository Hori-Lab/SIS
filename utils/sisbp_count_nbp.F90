program sisbp_count_nbp

   use, intrinsic :: iso_fortran_env, Only : iostat_end
   implicit none

   integer, parameter :: PREC = 4
   integer, parameter :: hdl_bp = 10
   integer(4) :: kind_int, kind_real
   character(len=1000) :: cfile_bp
   integer :: i, j
   real(PREC) :: e
   integer :: nbp
   integer :: istat
   logical :: flg_frame

   if (command_argument_count() == 1) then
      
      call get_command_argument(1, cfile_bp)  

   else
      write(*,*) 'Usage: PROGRAM (input bp file)'
      stop
   endif

   open(hdl_bp, file=cfile_bp, status='old', action='read', form='unformatted',access='stream')

   read(hdl_bp) kind_int
   read(hdl_bp) kind_real

   loop_read: do while (.True.)

      flg_frame = .True.
      nbp = 0

      loop_frame: do while (flg_frame) 

         call read_bp(hdl_bp, i, j, e, istat)
         if (istat == iostat_end) then
            exit loop_read
         endif

         if (i == 0 .and. j == 0) then
            flg_frame = .False.
         else
            nbp = nbp + 1
         endif

      enddo loop_frame
      write(*,*) nbp
   enddo loop_read


contains
   subroutine read_bp(hdl, i, j, e, istat)

      integer, intent(in) :: hdl
      integer, intent(out) :: i, j
      real, intent(out) :: e
      integer(2) :: i2, j2
      real(4) :: e4
      integer, intent(out) :: istat

      if (kind_int == 2 .and. kind_real == 4) then
         read(hdl, iostat=istat) i2, j2, e4
         i = int(i2)
         j = int(j2)
         e = real(e4, kind=PREC)
      endif

   endsubroutine read_bp

endprogram sisbp_count_nbp
