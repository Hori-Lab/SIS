subroutine read_input(cfilepath, stat)
      
   use tomlf

   use const
   use const_phys
   use const_idx, only : JOBT, INTGRT, REPT
   use pbc, only : flg_pbc, pbc_box, set_pbc_size
   use var_io, only : iopen_hdl, &
                      flg_progress, step_progress, &
                      flg_out_bp, flg_out_bpe, flg_out_bpall, &
                      flg_in_ct, flg_in_bpseq, &
                      cfile_ff, cfile_dcd_in, &
                      cfile_prefix, cfile_pdb_ini, cfile_xyz_ini, cfile_fasta_in, cfile_anneal_in, &
                      cfile_ct_in, cfile_bpseq_in
   use var_state, only : job, tempK, kT, viscosity_Pas, opt_anneal, temp_independent,  temp_ref, &
                         nstep, dt, nstep_save, nstep_save_rst, integrator, nl_margin, &
                         flg_variable_box, variable_box_step, variable_box_change, &
                         rng_seed, stop_wall_time_sec, fix_com_origin, &
                         ionic_strength, length_per_charge
   use var_potential, only : flg_ele, ele_cutoff_type, ele_cutoff_inp, &
                             bp_min_loop, max_bp_per_nt, bp_model
   use var_top, only : nrepeat, nchains, inp_no_charge
   use var_replica, only : n_replica_temp, nstep_rep_exchange, nstep_rep_save, flg_exchange, &
                           replica_values, flg_replica
  
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
   integer(INT64) :: idummy
   real(PREC) :: rdummy, v(3)
   character(len=:), allocatable :: cline
   character(len=5) :: cquery

   stat = .False.

   iopen_hdl = iopen_hdl + 1
   hdl = iopen_hdl

   print '(2a)', "Reading input file: ", trim(cfilePath)
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
   print '(2a)', '# title: ', trim(cline)

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
      print '(2a)', 'Error: Unknown job type, ', trim(cline)
      return
   endif
   print '(3a,i3,a)', '# job type: ', trim(cline), ' (job=', job,')'

   !################# input files #################
   call get_value(table, "files", group)

   if (.not. associated(group)) then
      print '(a)', 'Error in input file: no files gorup in input.'
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
      call get_value(node, "ct", cfile_ct_in)
      call get_value(node, "bpseq", cfile_bpseq_in)
      call get_value(node, "anneal", cfile_anneal_in)

   else
      print '(a)', 'Error in input file: no files.in.'
      return
   endif

   if (allocated(cfile_ct_in)) flg_in_ct = .True.
   if (allocated(cfile_bpseq_in)) flg_in_bpseq = .True.

   if (flg_in_ct .and. flg_in_bpseq) then
      print '(a)', 'Error: only one of ct and bpseq files can be specified.'
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
      !   print *, list(i)%key, cline
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
            print '(a)', 'Error in input file: Unknown output type, '//trim(cline)
            return
         endif
      enddo

   else
      print '(a)', 'Error in input file: no files.out.'
      return
   endif

   !################# Replica #################
   flg_replica = .False.
   call get_value(table, "replica", group, requested=.False.)

   if (associated(group)) then

      flg_replica = .True.

      n_replica_temp = 0
      call get_value(group, "n_replica_temp", n_replica_temp)
      call get_value(group, "nstep_exchange", nstep_rep_exchange)
      call get_value(group, "nstep_save", nstep_rep_save)

      flg_exchange = .True.
      call get_value(group, "exchange", flg_exchange)

      print '(a,i16)', '# Replica, n_replica: ', n_replica_temp
      print '(a,i16)', '# Replica, nstep_exchange: ', nstep_rep_exchange
      print '(a,i16)', '# Replica, nstep_save: ', nstep_rep_save
      if (flg_exchange) then
         print '(a)', '# Replica, exchange: True'
      else
         print '(a)', '# Replica, exchange: False'
      endif
      print '(a)', '#'

      if (n_replica_temp > 0) then
         call get_value(group, "temperature", node, requested=.False.)
         if (associated(node)) then
            do i = 1, n_replica_temp
               write(cquery, '(i0)') i
               call get_value(node, cquery, replica_values(i, REPT%TEMP))
            enddo
         else
            print '(a)', 'Error in input file: [replcia.temperature] is needed.'
            return
         endif
      endif

   else
      n_replica_temp = 1
      flg_exchange = .False.
   endif

   !################# Condition #################
   call get_value(table, "condition", group)

   call get_value(group, "rng_seed", rng_seed)

   opt_anneal = 0
   call get_value(group, "opt_anneal", opt_anneal)

   if (opt_anneal > 0 .and. .not. allocated(cfile_anneal_in)) then
      print '(a)', 'Error: opt_anneal requires anneal in [files.in].'
      return
   endif

   tempK = -1.0
   call get_value(group, "tempK", tempK)

   if (opt_anneal == 0 .and. tempK < 0.0) then
      print '(a)', 'Error: tempK is invalid or undefined in [condition].'
      return
   endif

   kT = BOLTZ_KCAL_MOL * tempK
   print '(a,g15.8)', '# Condition, tempK: ', tempK
   print '(a,i16)', '# Condition, rng_seed: ', rng_seed
   print '(a)', '#'

   temp_independent = 0
   call get_value(group, "temp_independent", temp_independent)

   if (temp_independent > 0) then
      call get_value(group, "temp_ref", temp_ref)
   endif

   !################# Repeat sequence #################
   if (.not. allocated(cfile_fasta_in)) then
      call get_value(table, "repeat", group)
      if (associated(group)) then
         call get_value(group, "n_repeat", nrepeat)
         call get_value(group, "n_chain", nchains)
         print '(a,i10)', '# Repeat n_repeat: ', nrepeat
         print '(a,i10)', '# Repeat n_chains: ', nchains
         print '(a)', '#'
      endif
   else
      nrepeat = 0
   endif

   !################# MD #################
   if (job == JOBT%MD) then

      call get_value(table, "MD", group, requested=.false.)

      if (.not. associated(group)) then 
         print '(a)', 'Error: [MD] field required.'
         return
      endif

      !###### integrator #######
      call get_value(group, "integrator", cline)

      if (.not. allocated(cline)) then
         print '(a)', 'Error: integrator is required in [MD].'
         return

      else if (cline == 'GJF-2GJ') then
         integrator = INTGRT%LD_GJF2GJ

      else
         print '(2a)', 'Error: Unknown integrator type, ', trim(cline)
         return
      endif
      print '(2a)', '# MD integrator: ', trim(cline)

      !###### dt #######
      dt = -1.0
      call get_value(group, "dt", dt, stat=istat)
      if (istat /= 0 .or. dt < 0.0) then
         print '(a)', 'Error: invalid value for dt in [MD].'
         return
      endif
      print '(a,g15.8)', '# MD dt: ', dt

      !###### nstep #######
      nstep = -1
      call get_value(group, "nstep", nstep, stat=istat)
      if (istat /= 0 .or. nstep < 0) then
         print '(a)', 'Error: invalid value for nstep in [MD].'
         return
      endif
      print '(a,i16)', '# MD nstep: ', nstep

      !###### nstep_save #######
      nstep_save = -1
      call get_value(group, "nstep_save", nstep_save, stat=istat)
      if (istat /= 0 .or. nstep_save < 0) then
         print '(a)', 'Error: invalid value for nstep_save in [MD].'
         return
      endif
      print '(a,i16)', '# MD nstep_save: ', nstep_save

      !###### nstep_save_rst #######
      nstep_save_rst = nstep_save
      call get_value(group, "nstep_save_rst", nstep_save_rst, stat=istat)
      if (istat /= 0) then
         print '(a)', 'Error: invalid value for nstep_save_rst in [MD].'
         return
      endif
      print '(a,i16)', '# MD nstep_save_rst: ', nstep_save_rst

      !###### neighbor_list_margin ######
      nl_margin = -1.0
      call get_value(group, "neighbor_list_margin", nl_margin, stat=istat)
      if (istat /= 0) then
         print '(a)', 'Error: invalid value for neighbor_list_margin in [MD].'
         return
      else if (nl_margin < 0.0) then
         nl_margin = 10.0_PREC
         print '(a)', 'Warning: neighbor_list_margin is not specified in [MD] field. The default value will be used.'
      endif
      print '(a,g15.8)', '# MD neighbor_list_margin: ', nl_margin

      !###### viscosity_Pas ######
      viscosity_Pas = -1.0
      call get_value(group, "viscosity_Pas", viscosity_Pas, stat=istat)
      if (istat /= 0) then
         print '(a)', 'Error: invalid value for viscosity_Pas in [MD].'
         return
      else if (viscosity_Pas < 0.0) then
         viscosity_Pas = 0.00001_PREC
         print '(a)', 'Warning: viscosity_Pas is not specified in [MD] field. The default value will be used.'
      endif
      print '(a,f12.8)', '# MD viscosity_Pas: ', viscosity_Pas

      !###### stop_wall_time_hour ######
      rdummy = -1.0
      call get_value(group, "stop_wall_time_hour", rdummy, stat=istat)

      if (istat /= 0) then
         ! Try integer
         call get_value(group, "stop_wall_time_hour", idummy, stat=istat)

         if (istat /= 0) then
            print '(a)', 'Error: invalid value for stop_wall_time_hour in [MD].'
            return
         endif

         rdummy = real(idummy, kind=PREC)
      endif

      if (rdummy < 0.0) then
         stop_wall_time_sec = -1
         print '(a,g15.8,a)', '# MD stop_wall_time_hour: ', rdummy, ' (wall time limit not set)'

      else
         stop_wall_time_sec = int(rdummy * 3600.0_PREC, kind=INT64)
         print '(a,g15.8)', '# MD stop_wall_time_hour: ', rdummy
      endif

      !###### fix_com_origin ######
      fix_com_origin = 0
      call get_value(group, "fix_com_origin", fix_com_origin, stat=istat)

      if (istat /= 0) then
         print '(a)', 'Error: invalid value for fix_com_origin in [MD].'
         return
      endif
      print '(a,i5)', '# MD fix_com_origin: ', fix_com_origin

      print '(a)', '#'

   endif

   !################# Basepair #################
   call get_value(table, "Basepair", group, requested=.False.)

   if (associated(group)) then

      ! bp_model
      bp_model = INVALID_INT_VALUE
      call get_value(group, "model", bp_model)
      if (bp_model > INVALID_INT_JUDGE) then
         print '(a)', '# [Basepair] model is not specified in the input file. Default value applies.'
         bp_model = 1   ! default
      endif

      ! max_bp_per_nt
      max_bp_per_nt = INVALID_INT_VALUE
      call get_value(group, "max_bp_per_nt", max_bp_per_nt)
      if (max_bp_per_nt > INVALID_INT_JUDGE) then
         print '(a)', '# [Basepair] max_bp_per_nt is not specified in the input file. Default value applies.'
         max_bp_per_nt = -1   ! default
      endif

      ! min_loop
      bp_min_loop = -1
      call get_value(group, "min_loop", bp_min_loop)
      if (bp_min_loop < 0) then
         print '(a)', '# [Basepair] min_loop is not specified in the input file. Default values are used.'
         bp_min_loop = 3    ! default
      endif

   else
      print '(a)', '# [Basepair] section does not exist in the input file. Default values are used.'
      max_bp_per_nt = 1  ! default
      bp_min_loop = 3    ! default

   endif

   print '(a,i6)', '# Basepair, model: ', bp_model
   print '(a,i6)', '# Basepair, max_bp_per_nt: ', max_bp_per_nt
   print '(a,i6)', '# Basepair, min_loop: ', bp_min_loop
   print '(a)', '#'


   !################# Electrostatic #################
   flg_ele = .False.
   call get_value(table, "Electrostatic", group, requested=.False.)

   if (associated(group)) then

      print '(a)', '# Electrostatic: On'
      flg_ele = .True.

      ! ionic_strength
      ionic_strength = INVALID_VALUE
      call get_value(group, "ionic_strength", ionic_strength)
      if (ionic_strength > INVALID_JUDGE) then
         print '(a)', "Error: Invalid value for ionic_strength in [Electrostatic]."
         return
      else
         print '(a,g15.8)', '# Electrostatic, ionic strength: ', ionic_strength
      endif

      ! cutoff type   1: distance-based (default) 
      !               2: x Debye length
      ele_cutoff_type = 1
      call get_value(group, "cutoff_type", ele_cutoff_type)
      if (ele_cutoff_type == 1) then
         print '(a)', "# Electrostatic, cutoff type: 1 (distance-based)"
      else if (ele_cutoff_type == 2) then
         print '(a)', "# Electrostatic, cutoff type: 2 cutoff will be multiplied by the Debye length"
      else
         print '(a)', "Error: Invalid values for cutoff_type in [Electrostatic]."
         return
      endif

      ! cutoff
      ele_cutoff_inp = INVALID_VALUE
      call get_value(group, "cutoff", ele_cutoff_inp)
      if (ele_cutoff_inp > INVALID_JUDGE) then
         print '(a)', "Error: Invalid values or absence of cutoff in [Electrostatic]."
         return
      else
         print '(a,g15.8)', '# Electrostatic, cutoff: ', ele_cutoff_inp
      endif

      ! length_per_charge
      length_per_charge = INVALID_VALUE
      call get_value(group, "length_per_charge", length_per_charge)
      if (length_per_charge > INVALID_JUDGE) then
         print '(a)', "Error: Invalid value for length_per_charge in [Electrostatic]."
         return
      else
         print '(a,g15.8)', '# Electrostatic, length per charge: ', length_per_charge
      endif

      ! (optional) No charge particles
      call get_value(group, "no_charge", array)

      if (len(array) > 0) then
         print '(a)', "# Electrostatic, no charges on the following particles:"
         allocate(inp_no_charge(len(array)))
         do i = 1, len(array)
            call get_value(array, i, inp_no_charge(i))
            print '(a, i8, 1x, i8)', '#         ', i, inp_no_charge(i)
         enddo
      endif

      print '(a)', '#'
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
      print '(a,g15.8)', '# pbc_box x: ', pbc_box(1)
      print '(a,g15.8)', '# pbc_box y: ', pbc_box(2)
      print '(a,g15.8)', '# pbc_box z: ', pbc_box(3)
      print '(a)', '#'
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
      print '(a,i16)', '# variable_box step: ', variable_box_step
      print '(a,g15.8)', '# variable_box change_x: ', variable_box_change(1)
      print '(a,g15.8)', '# variable_box change_y: ', variable_box_change(2)
      print '(a,g15.8)', '# variable_box change_z: ', variable_box_change(3)
      print '(a)', '#'
   endif

   !################# Progress #################
   flg_progress = .False.
   call get_value(table, "Progress", group, requested=.False.)
   if (associated(group)) then
      flg_progress = .True.
      call get_value(group, "step", step_progress)
      print '(a,i16)', '# Progress step: ', step_progress
      print '(a)', '#'
   endif

   call table%destroy

   print '(a)', 'Done: reading input file'
   print *
   flush(6)

   stat = .True.

end subroutine read_input
