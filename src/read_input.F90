subroutine read_input(cfilepath, stat)
      
   use tomlf

   use const
   use const_phys
   use const_idx
   use var_io, only : iopen_hdl, &
                      flg_out_bp, flg_out_bpe, flg_out_bpall, &
                      cfile_ff, cfile_dcd_in, &
                      cfile_prefix, cfile_pdb_ini, cfile_fasta_in
   use var_state, only : job, tempK
   use var_top, only : nrepeat, nchains
  
   implicit none

   character(len=CHAR_FILE_PATH), intent(in) :: cfilepath
   logical, intent(out) :: stat


   !======= TOML
   type(toml_table), allocatable :: table
   type(toml_table), pointer :: group, node
   type(toml_array), pointer :: array
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

   else if (cline == 'CHECK_FORCE') then
      job = JOBT%CHECK_FORCE

   else if (cline == 'MD') then
      job = JOBT%MD

   else
      write(*,*) 'Unknown job type: '//trim(cline)
      stat = .False.
      return
   endif

   !################# input files #################
   call get_value(table, "files", group)

   if (.not. associated(group)) then
      write(*,*) 'no files gorup in input'
      stat = .False.
      return
   endif

   call get_value(group, "in", node)
   if (associated(node)) then
      call get_value(node, "ff", cfile_ff)

      if (job == JOBT%DCD) then
         call get_value(node, "dcd", cfile_dcd_in)
      endif

      call get_value(node, "pdb_ini", cfile_pdb_ini)

      call get_value(node, "fasta", cfile_fasta_in)

   else
      write(*,*) 'no files.in'
      stop
   endif

   !################# output files #################
   call get_value(group, "out", node)
   if (associated(node)) then
      call get_value(node, "prefix", cfile_prefix)

      !call node%get_keys(list)
      !do i = 1, size(list)
      !   !call get_value(node, list(i)%key, cline, stat=stat)
      !   call get_value(node, list(i)%key, cline)
      !   write(*,*) list(i)%key, cline
      !enddo

      call get_value(node, "types", array)

      do i = 1, len(array)
         call get_value(array, i, cline)

         if (cline == "bp") then
            flg_out_bp = .True.
         else if (cline == "bpall") then
            flg_out_bpall = .True.
         else if (cline == "bpe") then
            flg_out_bpe = .True.
         else
            write(*,*) 'Unknown output type: '//trim(cline)
            stat = .False.
            return
         endif
      enddo

   else
      write(*,*) 'no files.out'
      stop
   endif

   !################# Condition #################
   call get_value(table, "condition", group)
   call get_value(group, "tempK", tempK)


   !################# Repeat sequence #################
   if (.not. allocated(cfile_fasta_in)) then
      call get_value(table, "repeat", group)
      if (associated(group)) then
         call get_value(group, "n_repeat", nrepeat)
         call get_value(group, "n_chain", nchains)
      endif
   else
      nrepeat = 0
   endif

   call table%destroy

   stat = .True.

end subroutine read_input
