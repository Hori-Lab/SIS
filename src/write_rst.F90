subroutine write_rst()

   use mt_stream, only : save !, print (for debug)
   use const, only : MTS_SIZE
   use const_idx, only : RSTBLK
   use var_io,    only : hdl_rst, cfile_rst
   use var_top, only : nmp
   use var_state, only : istep, mts_rep, mts, &
                         opt_anneal, ianneal, &
                         xyz, velos, accels
   use var_potential, only : flg_twz, ntwz_FR, twz_FR_init
   use var_replica, only : nrep_proc, nrep_all, rep2lab, irep2grep
   use var_parallel, only : myrank

   implicit none

   integer :: i, imp, irep, ipair
   integer :: istep_sim
   integer :: grep
   integer :: lunout
   integer :: istat
   integer :: nblock_size

   do irep = 1, nrep_proc

      grep = irep2grep(irep)

      ! Open
      lunout = hdl_rst(irep)
      open(lunout, file=trim(cfile_rst(irep)), status='replace', action='write', &
           iostat=istat, form = 'unformatted', access = 'stream')

      if(istat > 0) then
         print '(2a)', 'Error: cannot open the file: ', trim(cfile_rst(irep))
         error stop
      end if

      istep_sim = 1   ! Unused but reserved

      ! replica
      if (myrank == 0 .and. irep == 1) then
         if (nrep_all > 1) then
            write(lunout) RSTBLK%REPLICA
            nblock_size = calc_size(1+2*nrep_all, 0, 0, 0)
            write(lunout) nblock_size
            write(lunout) nrep_all  ! M_INT
            do i = 1, nrep_all
               write(lunout) i, rep2lab(i)  ! M_INT
            enddo
         endif
      endif

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
         write(lunout) (xyz(i,imp,irep),i=1,3)   ! PREC
      enddo

      ! Velocity
      write(lunout) RSTBLK%VELO
      nblock_size = calc_size(2, 0, nmp*3, 0)
      write(lunout) nblock_size
      write(lunout) grep    ! M_INT
      write(lunout) nmp     ! M_INT
      do imp = 1, nmp
         write(lunout) (velos(i,imp,irep),i=1,3)  ! PREC
      enddo

      ! Acceleration
      write(lunout) RSTBLK%ACCEL
      nblock_size = calc_size(2, 0, nmp*3, 0)
      write(lunout) nblock_size
      write(lunout) grep    ! M_INT
      write(lunout) nmp     ! M_INT
      do imp = 1, nmp
         write(lunout) (accels(i,imp,irep),i=1,3) ! PREC
      enddo

      ! Tweezers
      if (flg_twz) then
         if (ntwz_FR > 0) then
            write(lunout) RSTBLK%TWZ
            nblock_size = calc_size(2, 0, ntwz_FR*2*3, 0)
            write(lunout) nblock_size
            write(lunout) grep    ! M_INT
            write(lunout) ntwz_FR     ! M_INT
            do ipair = 1, ntwz_FR
               write(lunout) (twz_FR_init(i, 1, ipair), i=1,3)  ! PREC
               write(lunout) (twz_FR_init(i, 2, ipair), i=1,3)  ! PREC
            enddo
         endif
      endif

      ! PRNGREP, mts_rep for replica exchange
      if (myrank == 0 .and. irep == 1) then
         if (nrep_all > 1) then
            write(lunout) RSTBLK%PRNGREP
            nblock_size = MTS_SIZE
            write(lunout) nblock_size
            call save(mts_rep, lunout)
            !call print(mts_rep)
         endif
      endif

      ! PRNG, mts for the process
      write(lunout) RSTBLK%PRNG
      nblock_size = calc_size(1, 0, 0, 0) + MTS_SIZE
      write(lunout) nblock_size
      write(lunout) grep    ! M_INT
      call save(mts(irep), lunout)
      !call print(mts(irep))

      close(lunout)

   enddo

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
