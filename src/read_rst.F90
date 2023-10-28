subroutine read_rst(itype_wanted)

   use mt_stream, only : read, mt_state  !, print (for debug)
   use const
   use const_idx
   use var_io,  only : hdl_in_rst
   use var_top, only : nmp
   use var_state, only : xyz, velos, accels, istep, ianneal, mts, mts_rep
   use var_replica, only : nrep_all, rep2lab, lab2rep, grep2rank, grep2irep
#ifdef PAR_MPI
   use var_replica, only : nrep_proc
#endif
   use var_parallel

   implicit none

   integer, intent(in) :: itype_wanted

   integer :: i, imp, itype
   integer :: n, m
   integer :: irep, grep, rank
   integer :: idummy
   integer :: istat
   integer :: nblock_size
   logical :: flg_done(1:nrep_all)
   real(PREC), allocatable :: temp_array(:,:)
#ifdef PAR_MPI
   integer, parameter :: TAG = 1
#endif

   if (myrank == 0) then

      rewind(hdl_in_rst)

      ! For checking completion
      flg_done(:) = .false.

      ! Do-loop for reading restart file
      do
         ! Read block-identifier
         read (hdl_in_rst, iostat=istat) itype
         if (istat < 0) then
            exit
         else if (istat > 0) then
            error stop 'Error: cannot read the restart file.'
         endif 

         read (hdl_in_rst) nblock_size

         if (itype /= itype_wanted) then
            call sub_skip_block(hdl_in_rst, nblock_size)
            cycle
         endif

         !#############################################################################
         select case (itype)

         !----------------------------
         ! step
         !----------------------------
         case(RSTBLK%STEP)
            read (hdl_in_rst) idummy
            read (hdl_in_rst) istep
#ifdef PAR_MPI
            call MPI_BCAST(istep, L_INT, MPI_BYTE, 0, MPI_COMM_WORLD, istat)
#endif
            flg_done(:) = .true.

            write(6, '(a)') '## RESTART: step has been loaded.'
            flush(6)
            exit

         !----------------------------
         ! Annealing
         !----------------------------
         case(RSTBLK%ANNEAL)
            read (hdl_in_rst) ianneal
            flg_done(:) = .true.

            write(6, '(a)') '## RESTART: annealing information has been loaded.'
            flush(6)
            exit

         !----------------------------
         ! xyz
         !----------------------------
         case(RSTBLK%XYZ)
            read (hdl_in_rst) grep
            read (hdl_in_rst) n
            if (n /= nmp) then
               print '(a,i10,a,i10)', 'Error: nmp is not consistent in the restart file. n=',n,' nmp=',nmp
               error stop
            endif

            allocate(temp_array(3, nmp))
            do imp = 1, nmp
               read (hdl_in_rst) (temp_array(i,imp), i=1,3)
            enddo

            rank = grep2rank(grep)
            irep = grep2irep(grep)
            if (rank == 0) then
               xyz(1:3, 1:nmp, irep) = temp_array(1:3, 1:nmp)
            endif

#ifdef PAR_MPI
            if (rank == 0) then
               ! Add local communication when necessary.
               continue
            else
               call MPI_SEND(irep, 1, MPI_INTEGER, rank, TAG, MPI_COMM_WORLD, istat)
               call MPI_SEND(temp_array, 3*nmp, PREC_MPI, rank, TAG, MPI_COMM_WORLD, istat)
            endif
#endif
            deallocate(temp_array)

            flg_done(grep) = .true.
            if (all(flg_done)) then
               write(6, '(a)') '## RESTART: coordinates data has been loaded.'
               flush(6)
               exit
            endif

         !----------------------------
         ! velos
         !----------------------------
         case(RSTBLK%VELO)
            read (hdl_in_rst) grep
            read (hdl_in_rst) n
            if (n /= nmp) then
               print '(a,i10,a,i10)', 'Error: nmp is not consistent. n=',n,' nmp=',nmp
               error stop
            endif

            allocate(temp_array(3, nmp))
            do imp = 1, nmp
               read (hdl_in_rst) (temp_array(i,imp), i=1,3)
            enddo

            rank = grep2rank(grep)
            irep = grep2irep(grep)
            if (rank == 0) then
               velos(1:3, 1:nmp, irep) = temp_array(1:3, 1:nmp)
            endif

#ifdef PAR_MPI
            if (rank == 0) then
               ! Add local communication when necessary.
               continue
            else
               call MPI_SEND(irep, 1, MPI_INTEGER, rank, TAG, MPI_COMM_WORLD, istat)
               call MPI_SEND(temp_array, 3*nmp, PREC_MPI, rank, TAG, MPI_COMM_WORLD, istat)
            endif
#endif
            deallocate(temp_array)

            flg_done(grep) = .true.
            if (all(flg_done)) then
               write(6, '(a)') '## RESTART: velocity data has been loaded.'
               flush(6)
               exit
            endif

         !----------------------------
         ! accels
         !----------------------------
         case(RSTBLK%ACCEL)
            read (hdl_in_rst) grep
            read (hdl_in_rst) n
            if (n /= nmp) then
               print '(a,i10,a,i10)', 'Error: nmp is not consistent in the restart file. n=',n,' nmp=',nmp
               error stop
            endif

            allocate(temp_array(3, nmp))
            do imp = 1, nmp
               read (hdl_in_rst) (temp_array(i,imp), i=1,3)
            enddo

            rank = grep2rank(grep)
            irep = grep2irep(grep)
            if (rank == 0) then
               accels(1:3, 1:nmp, irep) = temp_array(1:3, 1:nmp)
            endif

#ifdef PAR_MPI
            if (rank == 0) then
               ! Add local communication when necessary.
               continue
            else
               call MPI_SEND(irep, 1, MPI_INTEGER, rank, TAG, MPI_COMM_WORLD, istat)
               call MPI_SEND(temp_array, 3*nmp, PREC_MPI, rank, TAG, MPI_COMM_WORLD, istat)
            endif
#endif
            deallocate(temp_array)

            flg_done(grep) = .true.
            if (all(flg_done)) then
               write(6, '(a)') '## RESTART: acceleration data has been loaded.'
               flush(6)
               exit
            endif

         !----------------------------
         ! replica
         !----------------------------
         case(RSTBLK%REPLICA)
            read (hdl_in_rst) n
            if (n /= nrep_all) then
               print '(a,i10,a,i10)', 'Error: nrep_all is not consistent in the restart file. n=',n,' nrep_all=', nrep_all
               error stop
            endif

            do i=1, nrep_all
               read (hdl_in_rst) n, m
               rep2lab(n) = m
            enddo

#ifdef PAR_MPI
            call MPI_BCAST(rep2lab, nrep_all, MPI_INTEGER, 0, MPI_COMM_WORLD, istat)
#endif
            do irep = 1, nrep_all
               i = rep2lab(irep)
               lab2rep(i) = irep
            enddo

            flg_done(:) = .true.
            write(6, '(a)') '## RESTART: replica data has been loaded.'
            flush(6)
            exit

         !----------------------------
         ! PRNGREP (mts_rep for replica exchange)
         !----------------------------
         case(RSTBLK%PRNGREP)
            ! mts_rep for replica exchange
            call read(mts_rep, hdl_in_rst)
            !print *, 'mts_rep'
            !call print(mts_rep)

#ifdef PAR_MPI
            call MPI_BCAST(mts_rep, MTS_SIZE, MPI_BYTE, 0, MPI_COMM_WORLD, istat)
#endif
            flg_done(:) = .true.
            write(6, '(a)') '## RESTART: mts_rep has been loaded.'
            flush(6)
            exit

         !----------------------------
         ! PRNG (mts for each process)
         !----------------------------
         case(RSTBLK%PRNG)
            read (hdl_in_rst) grep

            rank = grep2rank(grep)
            irep = grep2irep(grep)

            if (rank == 0) then
               call read(mts(irep), hdl_in_rst)
               !print *, 'mts(', irep, ')'
               !call print(mts(irep))
            else
               call read(mts(0), hdl_in_rst)
               !print *, 'mts(', 0, ')'
               !call print(mts(0))
            endif
#ifdef PAR_MPI
            if (rank == 0) then
               ! Add local communication when necessary.
               continue
            else
               call MPI_SEND(irep, 1, MPI_INTEGER, rank, TAG, MPI_COMM_WORLD, istat)
               call MPI_SEND(mts(0), MTS_SIZE, MPI_BYTE, rank, TAG, MPI_COMM_WORLD, istat)
            endif
#endif
            flg_done(grep) = .true.
            if (all(flg_done)) then
               write(6, '(a)') '## RESTART: mts data has been loaded.'
               flush(6)
               exit
            endif

         case default
            print '(a,i3)', 'Error: Unknown block-identifier in the restart file. itype=',itype
            error stop

         endselect
         !#############################################################################

      enddo

      ! check completion
      if (.not. all(flg_done)) then
         call sub_not_found(itype_wanted)
      endif

#ifdef PAR_MPI
   else ! myrank /= 0

      !#############################################################################
      select case (itype_wanted)

      !----------------------------
      ! step
      !----------------------------
      case(RSTBLK%STEP)
         call MPI_BCAST(istep, L_INT, MPI_BYTE, 0, MPI_COMM_WORLD, istat)
         write(6, '(a)') '## RESTART: annealing information has been received.'
         flush(6)

      !----------------------------
      ! xyz
      !----------------------------
      case(RSTBLK%XYZ)

         do i = 1, nrep_proc
            call MPI_RECV(irep, 1, MPI_INTEGER, 0, TAG, MPI_COMM_WORLD, istats_mpi, istat)
            call MPI_RECV(xyz(:,:,irep), 3*nmp, PREC_MPI, 0, TAG, MPI_COMM_WORLD, istats_mpi, istat)
         enddo

         ! Add local communication when necessary.

         write(6, '(a)') '## RESTART: coordinates data has been received.'
         flush(6)

      !----------------------------
      ! velos
      !----------------------------
      case(RSTBLK%VELO)

         do i = 1, nrep_proc
            call MPI_RECV(irep, 1, MPI_INTEGER, 0, TAG, MPI_COMM_WORLD, istats_mpi, istat)
            call MPI_RECV(velos(:,:,irep), 3*nmp, PREC_MPI, 0, TAG, MPI_COMM_WORLD, istats_mpi, istat)
         enddo

         ! Add local communication when necessary.

         write(6, '(a)') '## RESTART: velocity data has been received.'
         flush(6)

      !----------------------------
      ! accels
      !----------------------------
      case(RSTBLK%ACCEL)

         do i = 1, nrep_proc
            call MPI_RECV(irep, 1, MPI_INTEGER, 0, TAG, MPI_COMM_WORLD, istats_mpi, istat)
            call MPI_RECV(accels(:,:,irep), 3*nmp, PREC_MPI, 0, TAG, MPI_COMM_WORLD, istats_mpi, istat)
         enddo

         ! Add local communication when necessary.

         write(6, '(a)') '## RESTART: acceleration data has been received.'
         flush(6)

      !----------------------------
      ! replica
      !----------------------------
      case(RSTBLK%REPLICA)

         call MPI_BCAST(rep2lab, nrep_all, MPI_INTEGER, 0, MPI_COMM_WORLD, istat)

         do irep = 1, nrep_all
            i = rep2lab(irep)
            lab2rep(i) = irep
         enddo

         write(6, '(a)') '## RESTART: replica data has been received.'
         flush(6)

      !----------------------------
      ! PRNGREP
      !----------------------------
      case(RSTBLK%PRNGREP)

         call MPI_BCAST(mts_rep, MTS_SIZE, MPI_BYTE, 0, MPI_COMM_WORLD, istat)
         !print *, 'mts_rep'
         !call print(mts_rep)

         write(6, '(a)') '## RESTART: mts_rep data has been received.'
         flush(6)

      !----------------------------
      ! PRNG
      !----------------------------
      case(RSTBLK%PRNG)

         call MPI_RECV(irep, 1, MPI_INTEGER, 0, TAG, MPI_COMM_WORLD, istats_mpi, istat)
         call MPI_RECV(mts(irep), MTS_SIZE, MPI_BYTE, 0, TAG, MPI_COMM_WORLD, istats_mpi, istat)
         !print *, 'mts(', irep, ')'
         !call print(mts(irep))

         write(6, '(a)') '## RESTART: mts data has been received.'
         flush(6)

      case default
         print '(a,i3)', 'Error: Logical error undefined itype_wanted used.' ,itype_wanted
         error stop

      endselect
      !#############################################################################

#endif
   endif

#ifdef PAR_MPI
   call MPI_BARRIER(MPI_COMM_WORLD, istat)
#endif

contains

   subroutine sub_not_found(itype_wanted)

      use const_idx, only : RSTBLK
      integer, intent(in) :: itype_wanted

      select case (itype_wanted)
      case(RSTBLK%STEP)
         print '(a)', 'Absence or invalid format of "STEP" block in restart file.'
      case(RSTBLK%ANNEAL)
         print '(a)', 'Absence or invalid format of "ANNEAL" block in restart file.'
      case(RSTBLK%XYZ)
         print '(a)', 'Absence or invalid format of "XYZ" block in restart file.'
      case(RSTBLK%VELO)
         print '(a)', 'Absence or invalid format of "VELO" block in restart file.'
      case(RSTBLK%ACCEL)
         print '(a)', 'Absence or invalid format of "ACCEL" block in restart file.'
      case default
         continue
      endselect

      error stop

   endsubroutine sub_not_found

   subroutine sub_skip_block(hdl_in_rst, nblock_size)
      integer, intent(in) :: hdl_in_rst
      integer, intent(in) :: nblock_size
      character(nblock_size) :: a
      read (hdl_in_rst) a
   endsubroutine sub_skip_block

endsubroutine read_rst
