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
                      cfile_prefix, cfile_pdb_ini, cfile_xyz_ini, cfile_fasta_in, cfile_anneal_in
   use var_state, only : job, tempK, kT, viscosity_Pas, opt_anneal, &
                         nstep, dt, nstep_save, nstep_save_rst, integrator, nl_margin, &
                         flg_variable_box, variable_box_step, variable_box_change, &
                         rng_seed, &
                         ionic_strength, length_per_charge
   use var_potential, only : flg_ele, ele_cutoff_type, ele_cutoff_inp, &
                             bp_min_loop, max_bp_per_nt
   use var_top, only : nrepeat, nchains, inp_no_charge
  
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

   write(6, '(a)') "Reading input file: " // trim(cfilePath)
   open(hdl, file=cfilepath, status='old', action='read', iostat=istat)

   if (istat /= 0) then
      print '(2a)', 'Error: failed to open the input file. ', trim(cfilepath)
      error stop
   endif

   call toml_parse(table, hdl)

   close(hdl)
   iopen_hdl = iopen_hdl - 1

   if (.not. allocated(table)) then
      return
   endif

   call get_value(table, "title", cline)
   write(6, '(a)') '# title: ' // trim(cline)

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
      call get_value(node, "xyz_ini", cfile_xyz_ini)
      call get_value(node, "fasta", cfile_fasta_in)
      call get_value(node, "anneal", cfile_anneal_in)

   else
      write(*,*) 'Error in input file: no files.in.'
      return
   endif

   if (job == JOBT%MD .or. job == JOBT%DCD) then
      if (len(cfile_pdb_ini) < 1 .and. len(cfile_xyz_ini) < 1) then
         error stop 'Initial structure is not specified. Either XYZ or PDB is required.'

      else if (len(cfile_pdb_ini) > 0 .and. len(cfile_xyz_ini) > 0) then
         error stop 'Both XYZ and PDB are specified for the initial structure. Please use only one of them.'

      endif
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

      !###### nstep_save_rst #######
      nstep_save_rst = nstep_save
      call get_value(group, "nstep_save_rst", nstep_save_rst, stat=istat)
      if (istat /= 0) then
         write(*,*) 'Error: invalid value for nstep_save_rst in [MD].'
         return
      endif
      write(*,*) '# MD nstep_save_rst: ', nstep_save_rst

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
      write(*,*)

   endif

   !################# Basepair #################
   call get_value(table, "Basepair", group, requested=.False.)

   if (associated(group)) then

      write(6, '(a)') '# Basepair'

      ! max_bp_per_nt
      max_bp_per_nt = INVALID_INT_VALUE
      call get_value(group, "max_bp_per_nt", max_bp_per_nt)
      if (max_bp_per_nt > INVALID_INT_JUDGE) then
         write(6, '(a)') '#### max_bp_per_nt is not specified in the input file. Default value applies.'
         max_bp_per_nt = -1   ! default
      endif

      ! min_loop
      bp_min_loop = -1
      call get_value(group, "min_loop", bp_min_loop)
      if (bp_min_loop < 0) then
         write(6, '(a)') '#### min_loop is not specified in the input file. Default value applies.'
         bp_min_loop = 3    ! default
      endif

   else
      write(6, '(a)') '#### [Basepair] section does not exist in the input file. Default values will be used.'
      max_bp_per_nt = 1  ! default
      bp_min_loop = 3    ! default

   endif

   write(6, '(a,i6)') '# Basepair, max_bp_per_nt: ', max_bp_per_nt
   write(6, '(a,i6)') '# Basepair, min_loop: ', bp_min_loop
   write(6,*)


   !################# Electrostatic #################
   flg_ele = .False.
   call get_value(table, "Electrostatic", group, requested=.False.)

   if (associated(group)) then

      write(6, '(a)') '# Electrostatic: On'
      flg_ele = .True.

      ! ionic_strength
      ionic_strength = INVALID_VALUE
      call get_value(group, "ionic_strength", ionic_strength)
      if (ionic_strength > INVALID_JUDGE) then
         write(6, '(a)') "Error: Invalid value for ionic_strength in [Electrostatic]."
         return
      else
         write(6, '(a,g15.8)') '# Electrostatic, ionic strength: ', ionic_strength
      endif

      ! cutoff type   1: distance-based (default) 
      !               2: x Debye length
      ele_cutoff_type = 1
      call get_value(group, "cutoff_type", ele_cutoff_type)
      if (ele_cutoff_type == 1) then
         write(6, '(a)') "# Electrostatic, cutoff type: 1 (distance-based)"
      else if (ele_cutoff_type == 2) then
         write(6, '(a)') "# Electrostatic, cutoff type: 2 cutoff will be multiplied by the Debye length"
      else
         write(6, '(a)') "Error: Invalid values for cutoff_type in [Electrostatic]."
         return
      endif

      ! cutoff
      ele_cutoff_inp = INVALID_VALUE
      call get_value(group, "cutoff", ele_cutoff_inp)
      if (ele_cutoff_inp > INVALID_JUDGE) then
         write(6, '(a)') "Error: Invalid values or absence of cutoff in [Electrostatic]."
         return
      else
         write(6, '(a,g15.8)') '# Electrostatic, cutoff: ', ele_cutoff_inp
      endif

      ! length_per_charge
      length_per_charge = INVALID_VALUE
      call get_value(group, "length_per_charge", length_per_charge)
      write(6, '(a,g15.8)') '# Electrostatic, length_per_charge: ', length_per_charge
      if (length_per_charge > INVALID_JUDGE) then
         write(6, '(a)') "Error: Invalid value for length_per_charge in [Electrostatic]."
         return
      else
         write(6, '(a,g15.8)') '# Electrostatic, length per charge: ', length_per_charge
      endif

      ! (optional) No charge particles
      call get_value(group, "no_charge", array)

      if (len(array) > 0) then
         write(6, '(a)') "# No charges on the following particles:"
         allocate(inp_no_charge(len(array)))
         do i = 1, len(array)
            call get_value(array, i, inp_no_charge(i))
            write(6, '(i8, 1x, i8)') i, inp_no_charge(i)
         enddo
      endif

      write(6,*)
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
