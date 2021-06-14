subroutine read_input(cfilepath, stat)
      
   use tomlf

   use const
   use const_phys
   use const_idx
   use var_io, only : iopen_hdl, cfile_ff, cfile_prefix, cfile_dcd_in
   use var_state, only : job
   use var_top, only : nrepeat, nchains
  
   implicit none

   character(len=CHAR_FILE_PATH), intent(in) :: cfilepath
   logical, intent(out) :: stat

   !======= TOML
    type(toml_table), allocatable :: table
    type(toml_table), pointer :: group, node
    type(toml_key), allocatable :: list(:)
   !======= 

   integer :: i

   integer :: istat
   integer :: hdl
   !character(CHAR_FILE_LINE) :: cline
   character(len=:), allocatable :: cline

   iopen_hdl = iopen_hdl + 1
   hdl = iopen_hdl

   open(hdl, file=cfilepath, status='old', action='read', iostat=istat)

   call toml_parse(table, hdl)

   close(hdl)
   iopen_hdl = iopen_hdl - 1

   if (.not. allocated(table)) then
      stat = .False.
      return
   endif

   call get_value(table, "title", cline)
   write(*,*) cline

   !################# job #################
   call get_value(table, "job", group)
   call get_value(group, "type", cline)

   !cline = trim(cline)
   if (cline == 'DCD') then
      job = JOBT%DCD
   else
      write(*,*) 'Unknown job type: '//trim(cline)
   endif

   !################# files #################
   call get_value(table, "files", group)

   if (.not. associated(group)) then
      write(*,*) 'no files gorup in input'
      stat = .False.
      return
   endif

   call get_value(group, "in", node)
   if (associated(node)) then
      call get_value(node, "ff", cfile_ff)
      write(*,*) 'cfile_ff', cfile_ff

      !call node%get_keys(list)
      !do i = 1, size(list)
      !   !call get_value(node, list(i)%key, cline, stat=stat)
      !   call get_value(node, list(i)%key, cline)
      !   write(*,*) list(i)%key, cline
      !enddo

      if (job == JOBT%DCD) then
         call get_value(node, "dcd", cfile_dcd_in)
      endif

   else
      write(*,*) 'no files.in'
      stop
   endif

   call get_value(group, "out", node)
   if (associated(node)) then
      call get_value(node, "prefix", cfile_prefix)

      !call node%get_keys(list)
      !do i = 1, size(list)
      !   !call get_value(node, list(i)%key, cline, stat=stat)
      !   call get_value(node, list(i)%key, cline)
      !   write(*,*) list(i)%key, cline
      !enddo

   else
      write(*,*) 'no files.out'
      stop
   endif

   call get_value(table, "repeat", group)
   call get_value(group, "n_repeat", nrepeat)
   call get_value(group, "n_chain", nchains)


   call table%destroy

   stat = .True.

end subroutine read_input
