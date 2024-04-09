subroutine read_rst(itype_wanted, rst_status)

   use mt_stream, only : read, mt_state  !, print !(for debug)
   use const, only : PREC
   use const_idx, only : RSTBLK
   use var_io,  only : hdl_in_rst
   use var_top, only : nmp
   use var_state, only : xyz, velos, accels, istep, ianneal, mts, mts_rep
   use var_replica, only : nrep_all, rep2lab, lab2rep, grep2rank, grep2irep
   use var_potential, only : ntwz_FR, twz_FR_init
#ifdef PAR_MPI
   use const, only : L_INT
   use var_replica, only : nrep_proc
#endif
   use var_parallel

   implicit none

   integer, intent(in) :: itype_wanted
   integer, intent(out) :: rst_status

   integer :: i, imp, itype, ipair
   integer :: n, m
   integer :: irep, grep, rank
   integer :: idummy
   integer :: istat
   integer :: nblock_size
   logical :: flg_done(1:nrep_all)
   real(PREC), allocatable :: temp_array(:,:)
   real(PREC), allocatable :: temp_array3(:,:,:)
#ifdef PAR_MPI
   integer, parameter :: TAG = 1
#endif

   rst_status = 0

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
            ! This is to signal that PRNGREP was successfully loaded
            call MPI_BCAST(0, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, istat)

            !call MPI_BCAST(mts_rep, MTS_SIZE, MPI_BYTE, 0, MPI_COMM_WORLD, istat)
            ! This does not work because mts%state is an array of pointers

            call MPI_BCAST(mts_rep%i, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, istat)
            call MPI_BCAST(mts_rep%stream_id, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, istat)
            call MPI_BCAST(mts_rep%istatus, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, istat)
            call MPI_BCAST(mts_rep%nn, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, istat)
            call MPI_BCAST(mts_rep%mm, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, istat)
            call MPI_BCAST(mts_rep%rr, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, istat)
            call MPI_BCAST(mts_rep%ww, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, istat)
            call MPI_BCAST(mts_rep%aaa, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, istat)
            call MPI_BCAST(mts_rep%wmask, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, istat)
            call MPI_BCAST(mts_rep%umask, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, istat)
            call MPI_BCAST(mts_rep%lmask, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, istat)
            call MPI_BCAST(mts_rep%shift0, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, istat)
            call MPI_BCAST(mts_rep%shift1, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, istat)
            call MPI_BCAST(mts_rep%maskB, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, istat)
            call MPI_BCAST(mts_rep%maskC, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, istat)
            call MPI_BCAST(mts_rep%shiftB, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, istat)
            call MPI_BCAST(mts_rep%shiftC, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, istat)
            call MPI_BCAST(mts_rep%mag(0:1), 2, MPI_INTEGER, 0, MPI_COMM_WORLD, istat)
            call MPI_BCAST(mts_rep%state(0:623), 624, MPI_INTEGER, 0, MPI_COMM_WORLD, istat)
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

               !call MPI_SEND(mts(0), MTS_SIZE, MPI_BYTE, rank, TAG, MPI_COMM_WORLD, istat)
               ! This does not work because mts%state is an array of pointers

               call MPI_SEND(mts(0)%i, 1, MPI_INTEGER, rank, TAG, MPI_COMM_WORLD, istat)
               call MPI_SEND(mts(0)%stream_id, 1, MPI_INTEGER, rank, TAG, MPI_COMM_WORLD, istat)
               call MPI_SEND(mts(0)%istatus, 1, MPI_INTEGER, rank, TAG, MPI_COMM_WORLD, istat)
               call MPI_SEND(mts(0)%nn, 1, MPI_INTEGER, rank, TAG, MPI_COMM_WORLD, istat)
               call MPI_SEND(mts(0)%mm, 1, MPI_INTEGER, rank, TAG, MPI_COMM_WORLD, istat)
               call MPI_SEND(mts(0)%rr, 1, MPI_INTEGER, rank, TAG, MPI_COMM_WORLD, istat)
               call MPI_SEND(mts(0)%ww, 1, MPI_INTEGER, rank, TAG, MPI_COMM_WORLD, istat)
               call MPI_SEND(mts(0)%aaa, 1, MPI_INTEGER, rank, TAG, MPI_COMM_WORLD, istat)
               call MPI_SEND(mts(0)%wmask, 1, MPI_INTEGER, rank, TAG, MPI_COMM_WORLD, istat)
               call MPI_SEND(mts(0)%umask, 1, MPI_INTEGER, rank, TAG, MPI_COMM_WORLD, istat)
               call MPI_SEND(mts(0)%lmask, 1, MPI_INTEGER, rank, TAG, MPI_COMM_WORLD, istat)
               call MPI_SEND(mts(0)%shift0, 1, MPI_INTEGER, rank, TAG, MPI_COMM_WORLD, istat)
               call MPI_SEND(mts(0)%shift1, 1, MPI_INTEGER, rank, TAG, MPI_COMM_WORLD, istat)
               call MPI_SEND(mts(0)%maskB, 1, MPI_INTEGER, rank, TAG, MPI_COMM_WORLD, istat)
               call MPI_SEND(mts(0)%maskC, 1, MPI_INTEGER, rank, TAG, MPI_COMM_WORLD, istat)
               call MPI_SEND(mts(0)%shiftB, 1, MPI_INTEGER, rank, TAG, MPI_COMM_WORLD, istat)
               call MPI_SEND(mts(0)%shiftC, 1, MPI_INTEGER, rank, TAG, MPI_COMM_WORLD, istat)
               call MPI_SEND(mts(0)%mag(0:1), 2, MPI_INTEGER, rank, TAG, MPI_COMM_WORLD, istat)
               call MPI_SEND(mts(0)%state(0:623), 624, MPI_INTEGER, rank, TAG, MPI_COMM_WORLD, istat)
            endif
#endif
            flg_done(grep) = .true.
            if (all(flg_done)) then
               write(6, '(a)') '## RESTART: mts data has been loaded.'
               flush(6)
               exit
            endif

         !----------------------------
         ! Tweezers
         !----------------------------
         case(RSTBLK%TWZ)
            read (hdl_in_rst) grep

            if (ntwz_FR > 0) then
               read (hdl_in_rst) n

               if (n /= ntwz_FR) then
                  print '(a,i10,a,i10)', 'Error: ntwz_FR is not consistent in the restart file. n=',n,' ntwz_FR=',ntwz_FR
                  error stop
               endif

               allocate(temp_array3(3, 2, ntwz_FR))
               do ipair = 1, ntwz_FR
                  read (hdl_in_rst) (temp_array3(i, 1, ipair), i=1,3)
                  read (hdl_in_rst) (temp_array3(i, 2, ipair), i=1,3)
               enddo

               rank = grep2rank(grep)
               irep = grep2irep(grep)
               if (rank == 0) then
                  twz_FR_init(1:3, 1:2, 1:ntwz_FR) = temp_array3(1:3, 1:2, 1:ntwz_FR)
               endif

#ifdef PAR_MPI
               if (rank == 0) then
                  ! Add local communication when necessary.
                  continue
               else
                  call MPI_SEND(irep, 1, MPI_INTEGER, rank, TAG, MPI_COMM_WORLD, istat)
                  call MPI_SEND(temp_array3, 3*2*ntwz_FR, PREC_MPI, rank, TAG, MPI_COMM_WORLD, istat)
               endif
#endif
               deallocate(temp_array3)
            endif

            flg_done(grep) = .true.
            if (all(flg_done)) then
               write(6, '(a)') '## RESTART: tweezers data has been loaded.'
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

         ! This is to signal that PRNGREP was successfully loaded
         call MPI_BCAST(i, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, istat)
         if (i == -1) then
            rst_status = -1
            return
         endif

         !call MPI_BCAST(mts_rep, MTS_SIZE, MPI_BYTE, 0, MPI_COMM_WORLD, istat)
         ! This does not work because mts%state is an array of pointers

         call MPI_BCAST(mts_rep%i, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, istat)
         call MPI_BCAST(mts_rep%stream_id, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, istat)
         call MPI_BCAST(mts_rep%istatus, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, istat)
         call MPI_BCAST(mts_rep%nn, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, istat)
         call MPI_BCAST(mts_rep%mm, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, istat)
         call MPI_BCAST(mts_rep%rr, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, istat)
         call MPI_BCAST(mts_rep%ww, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, istat)
         call MPI_BCAST(mts_rep%aaa, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, istat)
         call MPI_BCAST(mts_rep%wmask, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, istat)
         call MPI_BCAST(mts_rep%umask, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, istat)
         call MPI_BCAST(mts_rep%lmask, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, istat)
         call MPI_BCAST(mts_rep%shift0, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, istat)
         call MPI_BCAST(mts_rep%shift1, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, istat)
         call MPI_BCAST(mts_rep%maskB, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, istat)
         call MPI_BCAST(mts_rep%maskC, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, istat)
         call MPI_BCAST(mts_rep%shiftB, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, istat)
         call MPI_BCAST(mts_rep%shiftC, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, istat)
         call MPI_BCAST(mts_rep%mag(0:1), 2, MPI_INTEGER, 0, MPI_COMM_WORLD, istat)
         call MPI_BCAST(mts_rep%state(0:623), 624, MPI_INTEGER, 0, MPI_COMM_WORLD, istat)
         !print *, 'mts_rep'
         !call print(mts_rep)

         write(6, '(a)') '## RESTART: mts_rep data has been received.'
         flush(6)

      !----------------------------
      ! PRNG
      !----------------------------
      case(RSTBLK%PRNG)

         do i = 1, nrep_proc
            call MPI_RECV(irep, 1, MPI_INTEGER, 0, TAG, MPI_COMM_WORLD, istats_mpi, istat)
            if (irep == -1) then
               rst_status = -1
               return
            endif

            !call MPI_RECV(mts(irep), MTS_SIZE, MPI_BYTE, 0, TAG, MPI_COMM_WORLD, istats_mpi, istat)
            ! This does not work because mts%state is an array of pointers

            call MPI_RECV(mts(irep)%i, 1, MPI_INTEGER, 0, TAG, MPI_COMM_WORLD, istats_mpi, istat)
            call MPI_RECV(mts(irep)%stream_id, 1, MPI_INTEGER, 0, TAG, MPI_COMM_WORLD, istats_mpi, istat)
            call MPI_RECV(mts(irep)%istatus, 1, MPI_INTEGER, 0, TAG, MPI_COMM_WORLD, istats_mpi, istat)
            call MPI_RECV(mts(irep)%nn, 1, MPI_INTEGER, 0, TAG, MPI_COMM_WORLD, istats_mpi, istat)
            call MPI_RECV(mts(irep)%mm, 1, MPI_INTEGER, 0, TAG, MPI_COMM_WORLD, istats_mpi, istat)
            call MPI_RECV(mts(irep)%rr, 1, MPI_INTEGER, 0, TAG, MPI_COMM_WORLD, istats_mpi, istat)
            call MPI_RECV(mts(irep)%ww, 1, MPI_INTEGER, 0, TAG, MPI_COMM_WORLD, istats_mpi, istat)
            call MPI_RECV(mts(irep)%aaa, 1, MPI_INTEGER, 0, TAG, MPI_COMM_WORLD, istats_mpi, istat)
            call MPI_RECV(mts(irep)%wmask, 1, MPI_INTEGER, 0, TAG, MPI_COMM_WORLD, istats_mpi, istat)
            call MPI_RECV(mts(irep)%umask, 1, MPI_INTEGER, 0, TAG, MPI_COMM_WORLD, istats_mpi, istat)
            call MPI_RECV(mts(irep)%lmask, 1, MPI_INTEGER, 0, TAG, MPI_COMM_WORLD, istats_mpi, istat)
            call MPI_RECV(mts(irep)%shift0, 1, MPI_INTEGER, 0, TAG, MPI_COMM_WORLD, istats_mpi, istat)
            call MPI_RECV(mts(irep)%shift1, 1, MPI_INTEGER, 0, TAG, MPI_COMM_WORLD, istats_mpi, istat)
            call MPI_RECV(mts(irep)%maskB, 1, MPI_INTEGER, 0, TAG, MPI_COMM_WORLD, istats_mpi, istat)
            call MPI_RECV(mts(irep)%maskC, 1, MPI_INTEGER, 0, TAG, MPI_COMM_WORLD, istats_mpi, istat)
            call MPI_RECV(mts(irep)%shiftB, 1, MPI_INTEGER, 0, TAG, MPI_COMM_WORLD, istats_mpi, istat)
            call MPI_RECV(mts(irep)%shiftC, 1, MPI_INTEGER, 0, TAG, MPI_COMM_WORLD, istats_mpi, istat)
            call MPI_RECV(mts(irep)%mag(0:1), 2, MPI_INTEGER, 0, TAG, MPI_COMM_WORLD, istats_mpi, istat)
            call MPI_RECV(mts(irep)%state(0:623), 624, MPI_INTEGER, 0, TAG, MPI_COMM_WORLD, istats_mpi, istat)
            !print *, 'mts(', irep, ')'
            !call print(mts(irep))

            write(6, '(a,i5,a)') '## RESTART: mts data has been received for irep = ',irep,'.'
            flush(6)
         enddo

      !----------------------------
      ! Tweezers
      !----------------------------
      case(RSTBLK%TWZ)

         if (ntwz_FR > 0) then
            do i = 1, nrep_proc
               call MPI_RECV(irep, 1, MPI_INTEGER, 0, TAG, MPI_COMM_WORLD, istats_mpi, istat)
               call MPI_RECV(twz_FR_init(:,:,:), 3*2*ntwz_FR, PREC_MPI, 0, TAG, MPI_COMM_WORLD, istats_mpi, istat)
            enddo

            ! Add local communication when necessary.

            write(6, '(a)') '## RESTART: tweezers data has been received.'
            flush(6)
         endif

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

#ifdef PAR_MPI
      integer :: rank
#endif

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
      case(RSTBLK%REPLICA)
         print '(a)', 'Absence or invalid format of "REPLICA" block in restart file.'
      case(RSTBLK%PRNG)
         print '(a)', 'Absence or invalid format of "PRNG" block in restart file.'
      case(RSTBLK%PRNGREP)
         print '(a)', 'Absence or invalid format of "PRNGREP" block in restart file.'
      case default
         print '(a)', 'Absence or invalid format of unknown block in restart file.'
      endselect
      flush(6)

      if (itype_wanted == RSTBLK%PRNGREP) then
#ifdef PAR_MPI
         call MPI_BCAST(-1, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, istat)
#endif
         rst_status = -1
         return

      else if (itype_wanted == RSTBLK%PRNG) then
#ifdef PAR_MPI
         do rank = 1, nprocs-1
            call MPI_SEND(-1, 1, MPI_INTEGER, rank, TAG, MPI_COMM_WORLD, istat)
         enddo
#endif
         rst_status = -1
         return

      else
         call sis_abort()
      endif

   endsubroutine sub_not_found

   subroutine sub_skip_block(hdl_in_rst, nblock_size)
      integer, intent(in) :: hdl_in_rst
      integer, intent(in) :: nblock_size
      character(nblock_size) :: a
      read (hdl_in_rst) a
   endsubroutine sub_skip_block

endsubroutine read_rst
