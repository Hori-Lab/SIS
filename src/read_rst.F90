subroutine read_rst(itype_wanted)

   use const
   use const_idx
   use var_io,  only : hdl_rst
   use var_top, only : nmp
   use var_state, only : xyz, velos, accels, istep, ianneal

   implicit none

   integer, intent(in) :: itype_wanted

   integer :: i, imp, itype
   integer :: n
   integer :: grep
   integer :: idummy
   integer :: istat
   integer :: nblock_size
   logical, allocatable :: flg_done(:)

   !if (myrank == 0) then

      rewind(hdl_rst)

      ! For checking completion
      allocate(flg_done(1:1))
      flg_done(:) = .false.

      ! Do-loop for reading restart file
      do 
         ! Read block-identifier
         read (hdl_rst, iostat=istat) itype
         if (istat < 0) then
            exit
         else if (istat > 0) then
            error stop 'Error: cannot read the restart file.'
         endif 

         read (hdl_rst) nblock_size

         if (itype /= itype_wanted) then
            call sub_skip_block(hdl_rst, nblock_size)
            cycle
         endif

         !#############################################################################
         select case (itype)

         !----------------------------
         ! step
         !----------------------------
         case(RSTBLK%STEP)
            read (hdl_rst) idummy
            read (hdl_rst) istep
            flg_done(:) = .true.

            write(6, *) '## RESTART: annealing information has been loaded.'
            flush(6)
            exit

         !----------------------------
         ! Annealing
         !----------------------------
         case(RSTBLK%ANNEAL)
            read (hdl_rst) ianneal
            flg_done(:) = .true.

            write(6, *) '## RESTART: step has been loaded.'
            flush(6)
            exit

         !----------------------------
         ! xyz
         !----------------------------
         case(RSTBLK%XYZ)
            read (hdl_rst) grep
            read (hdl_rst) n
            if (n /= nmp) then
               print '(a,i,a,i)', 'Error: nmp is not consistent. n=',n,' nmp=',nmp
               error stop
            endif

            do imp = 1, nmp
               read (hdl_rst) (xyz(i,imp), i=1,3)
            enddo

            flg_done(:) = .true.
            write(6,*) '## RESTART: coordinates data has been loaded.'
            flush(6)
            exit

         !----------------------------
         ! velos
         !----------------------------
         case(RSTBLK%VELO)
            read (hdl_rst) grep
            read (hdl_rst) n
            if (n /= nmp) then
               print '(a,i,a,i)', 'Error: nmp is not consistent. n=',n,' nmp=',nmp
               error stop
            endif

            do imp = 1, nmp
               read (hdl_rst) (velos(i,imp), i=1,3)
            enddo

            flg_done(:) = .true.
            write(6,*) '## RESTART: velocity data has been loaded.'
            flush(6)
            exit

         !----------------------------
         ! accels
         !----------------------------
         case(RSTBLK%ACCEL)
            read (hdl_rst) grep
            read (hdl_rst) n
            if (n /= nmp) then
               print '(a,i,a,i)', 'Error: nmp is not consistent. n=',n,' nmp=',nmp
               error stop
            endif

            do imp = 1, nmp
               read (hdl_rst) (accels(i,imp), i=1,3)
            enddo

            flg_done(:) = .true.
            write(6,*) '## RESTART: acceleration data has been loaded.'
            flush(6)
            exit

         case default
            print '(a,i)', 'Error: Unknown block-identifier in restart file. itype=',itype
            error stop

         endselect
         !#############################################################################

      enddo

      ! check completion
      if (.not. all(flg_done)) then
         deallocate(flg_done)
         call sub_not_found(itype_wanted)
      endif
      deallocate(flg_done)

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

   subroutine sub_skip_block(hdl_rst, nblock_size)
      integer, intent(in) :: hdl_rst
      integer, intent(in) :: nblock_size
      character(nblock_size) :: a
      read (hdl_rst) a
   endsubroutine sub_skip_block

endsubroutine read_rst
