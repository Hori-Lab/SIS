program show_rst

   use mt_stream, only : read, print, set_mt19937, new
   use const
   use const_idx, only : RSTBLK
   use mt_stream, only : mt_state

   implicit none

   integer :: iarg
   integer :: iargc
   character(1000) :: cfile_input

   integer :: i, imp, itype
   integer :: nblock_size
   integer :: istat
   integer :: nmp
   integer(L_INT) :: i_ll
   real(PREC) :: r1,r2,r3
   integer, parameter :: luninp = 10
   type(mt_state) :: mts

   iarg = iargc()

   if (iarg /= 1) then
      write(6, *) 'Usage: % show_rst [restart file]'
      stop
   end if

   call getarg(1, cfile_input) 

   ! Open
   open (luninp, file=trim(cfile_input), status='OLD', action='READ', iostat = istat, &
             form = 'unformatted', access = 'stream')
   ! exception
   if(istat > 0) then
      print '(2a)', 'Error: cannot open the file in input ', trim(cfile_input)
      error stop
   end if

   do 
      read (luninp, iostat=istat) itype
      if (istat < 0) then
         exit
      else if (istat > 0) then
         print '(2a)', 'Error: cannot open the file in input ', trim(cfile_input)
         error stop
      endif

      read (luninp) nblock_size

      select case (itype)
      case(RSTBLK%REPLICA)
         write(*,*) '# Replica'
         read (luninp) n
         write(*,*) 'nrep_all:', n
         do i = 1, n
            read (luninp) i2, i3
            write(*,*) i2, i3
         enddo

      case(RSTBLK%STEP)
         write(*,*) '# Step'
         read (luninp) i
         write(*,*) 'istep_sim:', i
         read (luninp) i_ll
         write(*,*) 'istep:', i_ll

      case(RSTBLK%ANNEAL)
         write(*,*) '# Annealing'
         read (luninp) i
         write(*,*) 'ianneal:', i

      case(RSTBLK%XYZ)
         write(*,*) '# Coordinates'
         read (luninp) i
         write(*,*) 'replica:', i
         read (luninp) nmp
         write(*,*) 'nmp:',nmp
         do imp = 1, nmp
            read (luninp) r1,r2,r3
            write(*,*) r1,r2,r3
         enddo

      case(RSTBLK%VELO)
         write(*,*) '# Velocities'
         read (luninp) i
         write(*,*) 'replica:', i
         read (luninp) nmp
         write(*,*) 'nmp:',nmp
         do imp = 1, nmp
            read (luninp) r1,r2,r3
            write(*,*) r1,r2,r3
         enddo

      case(RSTBLK%ACCEL)
         write(*,*) '# Accelerations'
         read (luninp) i
         write(*,*) 'replica:', i
         read (luninp) nmp
         write(*,*) 'nmp:',nmp
         do imp = 1, nmp
            read (luninp) r1,r2,r3
            write(*,*) r1,r2,r3
         enddo

      case(RSTBLK%PRNGREP)
         write(*,*) '# PRNGREP, mts_rep for replica exchange'
         call set_mt19937()
         call new(mts)
         call read(mts, luninp)
         call print(mts)

      case(RSTBLK%PRNG)
         write(*,*) '# PRNG, mts for each process'
         read (luninp) i
         write(*,*) 'replica:', i
         call set_mt19937()
         call new(mts)
         call read(mts, luninp)
         call print(mts)

      case default
         write(*,*) '#######################'
         write(*,*) 'unknown type of block'
         write(*,*) 'block identifier:', itype
         write(*,*) 'exit'
         exit

      endselect
   enddo

      
   close(luninp)

endprogram show_rst
