program sisbp_count_nbp_nt

   use, intrinsic :: iso_fortran_env, Only : iostat_end
   implicit none

   integer, parameter :: PREC = 4
   integer, parameter :: hdl_bp = 10
   integer(4) :: kind_int, kind_real
   character(len=1000) :: cfile_bp, c_in
   integer :: nnt
   integer :: i, j
   real(PREC) :: e
   integer :: nbp
   integer :: istat
   logical :: flg_frame
   integer, allocatable :: nbp_nt(:)
   integer :: i_nt, n
   integer :: nnt_count_bp(20)  ! Maximum 20 basepairs considered for each nt (should be much more than enough)
   integer :: max_bp

   if (command_argument_count() == 2) then
      
      call get_command_argument(1, cfile_bp)  
      call get_command_argument(2, c_in)
      read(c_in,*) nnt

   else
      write(*,*) 'Usage: PROGRAM (input bp file) (#nt)'
      stop
   endif

   allocate(nbp_nt(nnt))

   open(hdl_bp, file=cfile_bp, status='old', action='read', form='unformatted',access='stream')

   read(hdl_bp) kind_int
   read(hdl_bp) kind_real

   loop_read: do while (.True.)

      flg_frame = .True.
      nbp = 0
      nbp_nt(:) = 0
      nnt_count_bp(:) = 0

      loop_frame: do while (flg_frame) 

         call read_bp(hdl_bp, i, j, e, istat)
         if (istat == iostat_end) then
            exit loop_read
         endif

         if (i == 0 .and. j == 0) then
            flg_frame = .False.
         else
            nbp = nbp + 1
            nbp_nt(i) = nbp_nt(i) + 1
            nbp_nt(j) = nbp_nt(j) + 1
         endif

      enddo loop_frame

      !! Statistics
      max_bp = 0
      do i_nt = 1, nnt

         n = nbp_nt(i_nt)

         if (n > 20) then
            write(*,*) 'Error: nbp_nt(i_nt) > 20. i_nt = ', i_nt, ' nbp_nt(i_nt) = ', n
            stop
         else if (n > max_bp) then
            max_bp = n
         endif

         nnt_count_bp(n) = nnt_count_bp(n) + 1
      enddo

      write(*,*) max_bp, (nnt_count_bp(i), i=1,max_bp)

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

endprogram sisbp_count_nbp_nt
