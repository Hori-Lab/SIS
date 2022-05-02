program show_rst

   use const
   use const_idx, only : RSTBLK

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
      error stop 'Error: cannot open the file in input ' // trim(cfile_input)
   end if

   do 
      read (luninp, iostat=istat) itype
      if (istat < 0) then
         exit
      else if (istat > 0) then
         error stop 'Error: cannot read the file in input ' // trim(cfile_input)
      endif

      read (luninp) nblock_size

      select case (itype)
      case(RSTBLK%STEP)
         write(*,*) '# step'
         read (luninp) i
         write(*,*) 'istep_sim:', i
         read (luninp) i_ll
         write(*,*) 'istep:', i_ll

      case(RSTBLK%ANNEAL)
         write(*,*) '# annealing'
         read (luninp) i
         write(*,*) 'ianneal:', i

      case(RSTBLK%XYZ)
         write(*,*) '# xyz'
         read (luninp) i
         write(*,*) 'replica:', i
         read (luninp) nmp
         write(*,*) 'nmp:',nmp
         do imp = 1, nmp
            read (luninp) r1,r2,r3
            write(*,*) r1,r2,r3
         enddo

      case(RSTBLK%VELO)
         write(*,*) '# velocities'
         read (luninp) i
         write(*,*) 'replica:', i
         read (luninp) nmp
         write(*,*) 'nmp:',nmp
         do imp = 1, nmp
            read (luninp) r1,r2,r3
            write(*,*) r1,r2,r3
         enddo

      case(RSTBLK%ACCEL)
         write(*,*) '# accelerations'
         read (luninp) i
         write(*,*) 'replica:', i
         read (luninp) nmp
         write(*,*) 'nmp:',nmp
         do imp = 1, nmp
            read (luninp) r1,r2,r3
            write(*,*) r1,r2,r3
         enddo

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
