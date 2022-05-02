subroutine write_rst()

   use const_idx, only : RSTBLK
   use var_io,    only : cfile_prefix, iopen_hdl
   use var_top, only : nmp
   use var_state, only : istep, &
                         opt_anneal, ianneal, &
                         xyz, velos, accels

   implicit none

   integer :: i, imp
   integer :: istep_sim
   integer :: grep
   integer :: lunout
   integer :: istat
   integer :: nblock_size

   ! Open
   lunout = iopen_hdl
   iopen_hdl = iopen_hdl + 1
   open(lunout, file=trim(cfile_prefix)//'.rst', status='replace', action='write', &
        iostat=istat, form = 'unformatted', access = 'stream')

   if(istat > 0) then
      error stop 'Error: cannot open the file: ' // cfile_prefix // '.rst'
   end if

   grep = 1   ! Reserved for replica later
   istep_sim = 1   ! Unused but reserved

   ! Step 
   write(lunout) RSTBLK%STEP 
   nblock_size = calc_size(1, 1, 0, 0)
   write(lunout) nblock_size
   write(lunout) istep_sim    ! M_INT
   write(lunout) istep        ! L_INT

   ! Annealing
   if (opt_anneal > 0) then
      write(lunout) RSTBLK%ANNEAL
      nblock_size = calc_size(1, 0, 0, 0)
      write(lunout) nblock_size
      write(lunout) ianneal      ! M_INT
   endif

   ! Coordinate
   write(lunout) RSTBLK%XYZ
   nblock_size = calc_size(2, 0, nmp*3, 0)
   write(lunout) nblock_size
   write(lunout) grep     ! M_INT
   write(lunout) nmp      ! M_INT
   do imp = 1, nmp
      write(lunout) (xyz(i,imp),i=1,3)   ! PREC
   enddo

   ! Velocity
   write(lunout) RSTBLK%VELO
   nblock_size = calc_size(2, 0, nmp*3, 0)
   write(lunout) nblock_size
   write(lunout) grep    ! M_INT
   write(lunout) nmp     ! M_INT
   do imp = 1, nmp
      write(lunout) (velos(i,imp),i=1,3)  ! PREC
   enddo

   ! Acceleration
   write(lunout) RSTBLK%ACCEL
   nblock_size = calc_size(2, 0, nmp*3, 0)
   write(lunout) nblock_size
   write(lunout) grep    ! M_INT
   write(lunout) nmp     ! M_INT
   do imp = 1, nmp
      write(lunout) (accels(i,imp),i=1,3) ! PREC
   enddo

   close(lunout)
   iopen_hdl = iopen_hdl - 1

contains
   integer function calc_size(mi, li, dr, lo)

      use const, only : PREC, M_INT, L_INT, LOGIC

      integer, intent(in) :: mi  !< # of medium-size integer (default integer)
      integer, intent(in) :: li  !< # of long integer
      integer, intent(in) :: dr  !< # of double-precision real
      integer, intent(in) :: lo  !< # of logical

      calc_size = M_INT * mi + L_INT * li + PREC * dr + LOGIC * lo
   endfunction calc_size

endsubroutine write_rst
