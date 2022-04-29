subroutine read_input(cfilepath, stat)
      
   use tomlf

   use const
   use const_phys
   use const_idx, only : JOBT, INTGRT
   use pbc, only : flg_pbc, pbc_box, set_pbc_size
   use var_io, only : iopen_hdl, &
                      flg_progress, step_progress, &
                      flg_out_bp, flg_out_bpe, flg_out_bpall, &
                      cfile_ff, cfile_dcd_in, &
                      cfile_prefix, cfile_pdb_ini, cfile_fasta_in, cfile_anneal_in
   use var_state, only : job, tempK, kT, viscosity_Pas, opt_anneal, &
                         nstep, dt, nstep_save, integrator, nl_margin, &
                         flg_variable_box, variable_box_step, variable_box_change, &
                         rng_seed
   use var_top, only : nrepeat, nchains
  
   implicit none

   character(len=*), intent(in) :: cfilepath
   logical, intent(out) :: stat


   !======= TOML
   type(toml_table), allocatable :: table
   type(toml_table), pointer :: group, node
   type(toml_array), pointer :: array
   !======= 

   integer :: i
   integer :: istat
   integer :: hdl
   real(PREC) :: v(3)
   character(len=:), allocatable :: cline

   stat = .False.

   iopen_hdl = iopen_hdl + 1
   hdl = iopen_hdl

   write(*,*) "Reading input file: ", trim(cfilePath)
   open(hdl, file=cfilepath, status='old', action='read', iostat=istat)

   if (istat /= 0) then
      write(*,*) 'Error: failed to open the input file. '//trim(cfilepath)
      stop (2)
   endif

   call toml_parse(table, hdl)

   close(hdl)
   iopen_hdl = iopen_hdl - 1

   if (.not. allocated(table)) then
      return
   endif

   call get_value(table, "title", cline)
   write(*,*) '# title: ', trim(cline)

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
      write(*,*) 'Error: Unknown job type, '//trim(cline)
      return
   endif
   write(*,*) '# job type: ', trim(cline), ' (job=', job,')'

   !################# input files #################
   call get_value(table, "files", group)

   if (.not. associated(group)) then
      write(*,*) 'Error in input file: no files gorup in input.'
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
      call get_value(node, "anneal", cfile_anneal_in)

   else
      write(*,*) 'Error in input file: no files.in.'
      return
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
            write(*,*) 'Error in input file: Unknown output type, '//trim(cline)
            return
         endif
      enddo

   else
      write(*,*) 'Error in input file: no files.out.'
      return
   endif

   !################# Condition #################
   call get_value(table, "condition", group)

   call get_value(group, "rng_seed", rng_seed)

   opt_anneal = 0
   call get_value(group, "opt_anneal", opt_anneal)

   if (opt_anneal > 0 .and. .not. allocated(cfile_anneal_in)) then
      write(*,*) 'Error: opt_anneal requires anneal in [files.in].'
      return
   endif

   tempK = -1.0
   call get_value(group, "tempK", tempK)

   if (opt_anneal == 0 .and. tempK < 0.0) then
      write(*,*) 'Error: tempK is invalid or undefined in [condition].'
      return
   endif

   kT = BOLTZ_KCAL_MOL * tempK
   write(*,*) '# tempK: ', tempK
   write(*,*) '# rng_seed: ', rng_seed

   !################# Repeat sequence #################
   if (.not. allocated(cfile_fasta_in)) then
      call get_value(table, "repeat", group)
      if (associated(group)) then
         call get_value(group, "n_repeat", nrepeat)
         call get_value(group, "n_chain", nchains)
         write(*,*) '# repeat n_repeat: ', nrepeat
         write(*,*) '# repeat n_chains: ', nchains
      endif
   else
      nrepeat = 0
   endif

   !################# MD #################
   if (job == JOBT%MD) then

      call get_value(table, "MD", group, requested=.false.)

      if (.not. associated(group)) then 
         write(*,*) 'Error: [MD] field required.'
         return
      endif

      !###### integrator #######
      call get_value(group, "integrator", cline)

      if (.not. allocated(cline)) then
         write(*,*) 'Error: integrator is required in [MD].'
         return

      else if (cline == 'GJF-2GJ') then
         integrator = INTGRT%LD_GJF2GJ

      else
         write(*,*) 'Error: Unknown integrator type, '//trim(cline)
         return
      endif
      write(*,*) '# MD integrator: ', trim(cline)

      !###### dt #######
      dt = -1.0
      call get_value(group, "dt", dt, stat=istat)
      if (istat /= 0 .or. dt < 0.0) then
         write(*,*) 'Error: invalid value for dt in [MD].'
         return
      endif
      write(*,*) '# MD dt: ', dt

      !###### nstep #######
      nstep = -1
      call get_value(group, "nstep", nstep, stat=istat)
      if (istat /= 0 .or. nstep < 0) then
         write(*,*) 'Error: invalid value for nstep in [MD].'
         return
      endif
      write(*,*) '# MD nstep: ', nstep

      !###### nstep_save #######
      nstep_save = -1
      call get_value(group, "nstep_save", nstep_save, stat=istat)
      if (istat /= 0 .or. nstep_save < 0) then
         write(*,*) 'Error: invalid value for nstep_save in [MD].'
         return
      endif
      write(*,*) '# MD nstep_save: ', nstep_save

      !###### neighbor_list_margin ######
      nl_margin = -1.0
      call get_value(group, "neighbor_list_margin", nl_margin)
      if (istat /= 0) then
         write(*,*) 'Error: invalid value for neighbor_list_margin in [MD].'
         return
      else if (nl_margin < 0.0) then
         nl_margin = 10.0_PREC
         write(*,*) 'Warning: neighbor_list_margin is not specified in [MD] field. The default value will be used.'
      endif
      write(*,*) '# MD neighbor_list_margin: ', nl_margin

      !###### viscosity_Pas ######
      viscosity_Pas = -1.0
      call get_value(group, "viscosity_Pas", viscosity_Pas)
      if (istat /= 0) then
         write(*,*) 'Error: invalid value for viscosity_Pas in [MD].'
         return
      else if (viscosity_Pas < 0.0) then
         viscosity_Pas = 0.00001_PREC
         write(*,*) 'Warning: viscosity_Pas is not specified in [MD] field. The default value will be used.'
      endif
      write(*,*) '# MD viscosity_Pas: ', viscosity_Pas

   endif

   !################# box #################
   flg_pbc = .False.
   call get_value(table, "PBC_box", group, requested=.False.)
   if (associated(group)) then 
      flg_pbc = .True.
      call get_value(group, "x", v(1))
      call get_value(group, "y", v(2))
      call get_value(group, "z", v(3))
      call set_pbc_size(v)
      write(*,*) '# pbc_box x: ', pbc_box(1)
      write(*,*) '# pbc_box y: ', pbc_box(2)
      write(*,*) '# pbc_box z: ', pbc_box(3)
   endif

   !################# variable box #################
   flg_variable_box = .False.
   call get_value(table, "variable_box", group, requested=.False.)
   if (associated(group)) then 
      flg_variable_box = .True.
      call get_value(group, "step", variable_box_step)
      call get_value(group, "change_x", variable_box_change(1))
      call get_value(group, "change_y", variable_box_change(2))
      call get_value(group, "change_z", variable_box_change(3))
      write(*,*) '# variable_box change_x: ', variable_box_change(1)
      write(*,*) '# variable_box change_y: ', variable_box_change(2)
      write(*,*) '# variable_box change_z: ', variable_box_change(3)
   endif

   !################# Progress #################
   flg_progress = .False.
   call get_value(table, "Progress", group, requested=.False.)
   if (associated(group)) then
      flg_progress = .True.
      call get_value(group, "step", step_progress)
      write(*,*) '# Progress step: ', step_progress
   endif

   call table%destroy

   write(6,*) 'Done: reading input file'
   write(6,*) ''
   flush(6)
   stat = .True.

end subroutine read_input
