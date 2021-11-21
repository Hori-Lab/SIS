module dcd

   use, intrinsic :: iso_fortran_env, Only : iostat_end
   use const, only : PREC, CHAR_FILE_PATH

   implicit none

   type type_open_mode
      integer :: READ
      integer :: WRITE
   endtype type_open_mode
   type(type_open_mode), parameter :: DCD_OPEN_MODE = type_open_mode(1,2)

   type file_dcd
      character(len=CHAR_FILE_PATH) :: filepath

      integer :: hdl  ! file handle
      logical :: flg_opened
      logical :: flg_header_read
      logical :: flg_reach_end

      logical :: flg_unitcell
      real(PREC) :: box(3)

      integer :: nunit
      integer :: nmp_all

      contains
      procedure :: open_to_read => open_to_read
      procedure :: open_to_write => open_to_write
      procedure :: read_header => read_header
      procedure :: skip_header => skip_header
      procedure :: copy_header => copy_header
      procedure :: write_header => write_header
      procedure :: read_nmp => read_nmp
      procedure :: read_onestep => read_onestep
      procedure :: write_onestep => write_onestep
      procedure :: close => close
      procedure :: rewind => rewind

   end type file_dcd

   interface file_dcd
      module procedure init_file_dcd
   end interface file_dcd

   contains

   type(file_dcd) function init_file_dcd(handle, filepath, mode)
      integer, intent(in) :: handle
      character(len=*), intent(in) :: filepath
      integer, intent(in) :: mode

      init_file_dcd%flg_opened = .False.
      init_file_dcd%flg_header_read = .False.
      init_file_dcd%flg_reach_end = .False.
      init_file_dcd%flg_unitcell = .False.
      init_file_dcd%box(:) = -1.0
      init_file_dcd%nmp_all = -1
      init_file_dcd%nunit = -1

      if (mode == DCD_OPEN_MODE%READ) then
         call init_file_dcd%open_to_read(handle, filepath)

      else if (mode == DCD_OPEN_MODE%WRITE) then
         call init_file_dcd%open_to_write(handle, filepath)

      endif

   end function init_file_dcd

   subroutine open_to_read(self, handle, filepath)

      class(file_dcd), intent(inout) :: self
      integer, intent(in) :: handle
      character(len=*), intent(in) :: filepath

      integer :: iopen_status

      if (self%flg_opened) then
         write(*,*) 'Error: dcd file is already opened: ' // trim(filepath)
         stop
      end if

      self%hdl = handle
      self%filepath = filepath

      open(self%hdl, file = self%filepath, status = 'old', action = 'read', iostat=iopen_status, &
          form = 'unformatted', access = 'stream')
      if(iopen_status > 0) then 
         write(*,*) 'Error: cannot open the file: ' // trim(filepath)
         stop
      end if

      self%flg_opened = .True.

   end subroutine open_to_read

   subroutine open_to_write(self, handle, filepath)

      class(file_dcd), intent(inout) :: self
      integer, intent(in) :: handle
      character(len=*), intent(in) :: filepath

      integer :: iopen_status

      if (self%flg_opened) then
         write(*,*) 'Error: dcd file is already opened: ' // trim(filepath)
         stop
      end if

      self%hdl = handle
      self%filepath = filepath

      open(self%hdl, file = self%filepath, status = 'replace', action = 'write', iostat=iopen_status, &
           form = 'unformatted', access = 'stream')
      if(iopen_status > 0) then 
         write(*,*) 'Error: cannot open the file: ' // trim(filepath)
         stop
      end if

      self%flg_opened = .True.

   end subroutine open_to_write


   subroutine rewind(self)

      class(file_dcd), intent(inout) :: self

      rewind(self%hdl)

   end subroutine rewind

  
   subroutine close(self)

      class(file_dcd), intent(inout) :: self

      self%flg_opened = .False.
      close(self%hdl)

   end subroutine close


   subroutine read_header(self)
      implicit none

      class(file_dcd), intent(inout) :: self
    
      integer :: i
      integer :: idummy
      integer :: nblock_size
            
      ! ---------------------------------------------------------------------
      rewind(self%hdl)
    
      ! Check if there is unitcell
      self%flg_unitcell = .false.
      read (self%hdl) nblock_size
      do i = 1, nblock_size, 4
        read(self%hdl) idummy
        if (i == 45 .and. idummy == 1) then
           self%flg_unitcell = .true.
        endif
      enddo

      rewind(self%hdl)

      read (self%hdl) nblock_size
      do i = 1, nblock_size, 4
        read(self%hdl) idummy
      enddo
      read (self%hdl) nblock_size
      
      read (self%hdl) nblock_size
      do i = 1, nblock_size, 4
        read(self%hdl) idummy
      enddo
      read (self%hdl) nblock_size
    
      read (self%hdl) nblock_size
      do i = 1, nblock_size, 4
        read(self%hdl) idummy
      enddo
      read (self%hdl) nblock_size
    
      self%flg_header_read = .True.
    
   end subroutine read_header
  
 
   subroutine skip_header(self)
      implicit none

      class(file_dcd), intent(inout) :: self
    
      integer :: i
      integer :: idummy
      integer :: nblock_size
            
      ! ---------------------------------------------------------------------
      rewind(self%hdl)
    
      read (self%hdl) nblock_size
      do i = 1, nblock_size, 4
        read(self%hdl) idummy
      enddo
      read (self%hdl) nblock_size
      
      read (self%hdl) nblock_size
      do i = 1, nblock_size, 4
        read(self%hdl) idummy
      enddo
      read (self%hdl) nblock_size
    
      read (self%hdl) nblock_size
      do i = 1, nblock_size, 4
        read(self%hdl) idummy
      enddo
      read (self%hdl) nblock_size
    
   end subroutine skip_header


   subroutine copy_header(self, other)
      implicit none

      class(file_dcd), intent(inout) :: self
      class(file_dcd), intent(inout) :: other
    
      integer :: i
      integer :: idummy
      integer :: nblock_size
            
      ! ---------------------------------------------------------------------
      rewind(self%hdl)
    
      ! Check if there is unitcell
      self%flg_unitcell = .false.
      read (self%hdl) nblock_size
      do i = 1, nblock_size, 4
        read(self%hdl) idummy
        if (i == 45 .and. idummy == 1) then
           self%flg_unitcell = .true.
        endif
      enddo

      ! copy
      other%flg_unitcell = self%flg_unitcell
    
      rewind(self%hdl)
      rewind(other%hdl)
    
      read(self%hdl) nblock_size
      write(other%hdl) nblock_size
      do i = 1, nblock_size, 4
        read(self%hdl) idummy
        write(other%hdl) idummy
      enddo
      read(self%hdl) nblock_size
      write(other%hdl) nblock_size
      
      read(self%hdl) nblock_size
      write(other%hdl) nblock_size
      do i = 1, nblock_size, 4
        read(self%hdl) idummy
        write(other%hdl) idummy
      enddo
      read(self%hdl) nblock_size
      write(other%hdl) nblock_size
    
      read(self%hdl) nblock_size
      write(other%hdl) nblock_size
      do i = 1, nblock_size, 4
        read(self%hdl) idummy
        write(other%hdl) idummy
      enddo
      read(self%hdl) nblock_size
      write(other%hdl) nblock_size
    
      self%flg_header_read = .True.
      flush(other%hdl)
    
   end subroutine copy_header

   subroutine write_header(self, natom)
      implicit none

      class(file_dcd), intent(inout) :: self
      integer, intent(in) :: natom

      integer :: i
      integer :: ntitle, nblock_size
      integer, parameter :: idummy = 0
      character(4), parameter :: ctype = 'CORD'
      character(80) :: title

      ! ---------------------------------------------------------------------
      rewind(self%hdl)

      ! copy
      nblock_size = 84
      write(self%hdl) nblock_size
      write(self%hdl) ctype
      do i = 1, 20
         write(self%hdl) idummy
      enddo
      write(self%hdl) nblock_size

      ntitle = 1
      nblock_size = 4 + 80 * ntitle
      write(self%hdl) nblock_size
      write(self%hdl) ntitle
      title(1:80) = ' '
      do i = 1, ntitle
        write(self%hdl) title
      enddo
      write(self%hdl) nblock_size

      nblock_size = 4
      write(self%hdl) nblock_size
      write(self%hdl) natom
      write(self%hdl) nblock_size

      self%flg_header_read = .True.
      flush(self%hdl)

   end subroutine write_header

   subroutine read_nmp(self, nmp, istat)
    
      implicit none
    
      class(file_dcd), intent(inout) :: self
      integer, intent(out) :: nmp
      integer, intent(out) :: istat

      integer :: i
      integer :: idummy
      real(8) :: rdummy
      integer :: nblock_size
      ! ---------------------------------------------------------------------

      rewind(self%hdl)
      
      call self%skip_header()

      if (self%flg_unitcell) then
         read (self%hdl, iostat=istat) nblock_size
   
         select case (istat)
         case (0)
            continue
         case (iostat_end)
            self%flg_reach_end = .True.
            return
         case default
            write(*,*) 'unknown error in read_onestep of dcd'
            stop
         end select
         
         if (nblock_size == 48) then
            read(self%hdl) self%box(1)
            read(self%hdl) rdummy
            read(self%hdl) self%box(2)
            read(self%hdl) rdummy
            read(self%hdl) rdummy
            read(self%hdl) self%box(3)
         else
            do i = 1, nblock_size, 4
               read(self%hdl) idummy
            enddo
         endif
         read (self%hdl) nblock_size
   
         read (self%hdl) nmp
         nmp = nmp / 4

         if (self%nmp_all < 0) then
            self%nmp_all = nmp
         endif
   
      else
         read (self%hdl, iostat=istat) nmp
         nmp = nmp / 4

         if (self%nmp_all < 0) then
            self%nmp_all = nmp
         endif
   
         select case (istat)
         case (0)
            continue
         case (iostat_end)
            self%flg_reach_end = .True.
            return
         case default
            write(*,*) 'unknown error in read_onestep of dcd'
            stop
         end select
            
      endif

      istat = 0

      rewind(self%hdl)

      call self%skip_header()
    
   end subroutine read_nmp


   subroutine read_onestep(self, nmp, xyz, istat)
    
      implicit none
    
      class(file_dcd), intent(inout) :: self
      integer, intent(in) :: nmp
      real(PREC), intent(out) :: xyz(3,nmp)
      integer, intent(out) :: istat

      integer :: i
      integer :: imp
      integer :: idummy
      real(8) :: rdummy
      integer :: num, nblock_size
      real(4), allocatable :: xyz4(:,:)
      ! ---------------------------------------------------------------------

      if (self%flg_reach_end) then
         istat = iostat_end
         return
      endif

      if (self%flg_unitcell) then
         read (self%hdl, iostat=istat) nblock_size
   
         select case (istat)
         case (0)
            continue
         case (iostat_end)
            self%flg_reach_end = .True.
            return
         case default
            write(*,*) 'unknown error in read_onestep of dcd'
            stop
         end select
         
         if (nblock_size == 48) then
            read(self%hdl) self%box(1)
            read(self%hdl) rdummy
            read(self%hdl) self%box(2)
            read(self%hdl) rdummy
            read(self%hdl) rdummy
            read(self%hdl) self%box(3)
         else
            do i = 1, nblock_size, 4
               read(self%hdl) idummy
            enddo
         endif
         read (self%hdl) nblock_size
   
         read (self%hdl) num

         if (self%nmp_all < 0) then
            self%nmp_all = num / 4
         endif
   
      else
         read (self%hdl, iostat=istat) num

         if (self%nmp_all < 0) then
            self%nmp_all = num / 4
         endif
   
         select case (istat)
         case (0)
            continue
         case (iostat_end)
            self%flg_reach_end = .True.
            return
         case default
            write(*,*) 'unknown error in read_onestep of dcd'
            stop
         end select
            
      endif

      !write(*,*) 'nmp_all', self%nmp_all
      allocate(xyz4(3, self%nmp_all))

      !num = nmp_all*4
      !read (self%hdl) num
      read (self%hdl) (xyz4(1,imp),imp=1,self%nmp_all)
      read (self%hdl) num
      read (self%hdl) num
      read (self%hdl) (xyz4(2,imp),imp=1,self%nmp_all)
      read (self%hdl) num
      read (self%hdl) num
      read (self%hdl) (xyz4(3,imp),imp=1,self%nmp_all)
      read (self%hdl) num
    
      xyz(1:3, 1:self%nmp_all) = real(xyz4(1:3, 1:self%nmp_all), kind=PREC)
    
      deallocate(xyz4)

      istat = 0
    
   end subroutine read_onestep

   subroutine write_onestep(self, nmp, xyz)
    
      implicit none
    
      class(file_dcd), intent(inout) :: self
      integer, intent(in) :: nmp
      real(PREC), intent(in) :: xyz(3,nmp)

      integer :: imp
      integer :: num, nblock_size
      ! ---------------------------------------------------------------------

      if (self%flg_unitcell) then

         if (self%box(1) < 0.0) then
            write(*,*) 'Error: box information is not set in write_onestep, dcd.F90'
            stop
         endif

         nblock_size = 48 
         write(self%hdl) nblock_size
   
         write(self%hdl) self%box(1)
         write(self%hdl) real(0.0, kind=PREC)
         write(self%hdl) self%box(2)
         write(self%hdl) real(0.0, kind=PREC)
         write(self%hdl) real(0.0, kind=PREC)
         write(self%hdl) self%box(3)

         write(self%hdl) nblock_size

      endif
   
      num = nmp*4
      write(self%hdl) num
      write(self%hdl) (real(xyz(1,imp), kind=4),imp=1, nmp)
      write(self%hdl) num
      write(self%hdl) num
      write(self%hdl) (real(xyz(2,imp), kind=4),imp=1, nmp)
      write(self%hdl) num
      write(self%hdl) num
      write(self%hdl) (real(xyz(3,imp), kind=4),imp=1, nmp)
      write(self%hdl) num
    
      flush(self%hdl)

   end subroutine write_onestep

end module dcd
