subroutine read_input(cfilepath, stat)
      
   use tomlf

   use const
   use const_phys
   use const_idx, only : JOBT, INTGRT
   use pbc, only : flg_pbc, pbc_box, set_pbc_size
   use var_io, only : iopen_hdl, &
                      flg_progress, step_progress, &
                      flg_out_bp, flg_out_bpe, flg_out_bpall, &
                      flg_in_ct, flg_in_bpseq, &
                      cfile_ff, cfile_dcd_in, &
                      cfile_prefix, cfile_pdb_ini, cfile_xyz_ini, cfile_fasta_in, cfile_anneal_in, &
                      cfile_ct_in, cfile_bpseq_in
   use var_state, only : job, tempK, kT, viscosity_Pas, opt_anneal, temp_independent,  tempK_ref, &
                         nstep, dt, nstep_save, nstep_save_rst, integrator, nl_margin, &
                         flg_variable_box, variable_box_step, variable_box_change, &
                         rng_seed, stop_wall_time_sec, fix_com_origin, &
                         ionic_strength, length_per_charge
   use var_potential, only : flg_ele, ele_cutoff_type, ele_cutoff_inp, &
                             bp_min_loop, max_bp_per_nt, bp_model, &
                             flg_stage, stage_sigma, stage_eps
   use var_top, only : nrepeat, nchains, inp_no_charge
  
   implicit none

   character(len=*), intent(in) :: cfilepath
   logical, intent(out) :: stat


   !======= TOML
   type(toml_table), allocatable :: table
   type(toml_table), pointer :: group, node
   type(toml_array), pointer :: array
   type(toml_context) :: context
   type(toml_error), allocatable :: tm_err
   !======= 

   integer :: i
   integer :: istat, origin
   integer :: hdl
   real(PREC) :: rdummy, v(3)
   character(len=:), allocatable :: cline

   stat = .False.

   iopen_hdl = iopen_hdl + 1
   hdl = iopen_hdl

   print '(2a)', "Reading input file: ", trim(cfilePath)
   open(hdl, file=cfilepath, status='old', action='read', iostat=istat)

   if (istat /= 0) then
      print '(2a)', 'Error: failed to open the input file. ', trim(cfilepath)
      return
   endif

   call toml_load(table, hdl, context=context, error=tm_err)

   close(hdl)
   iopen_hdl = iopen_hdl - 1

   if (allocated(tm_err)) then
      print '(a)', tm_err%message
      flush(6)
      return
   endif


   !################# title #################
   ! (optional)
   call get_value(table, "title", cline)
   if (allocated(cline)) then
      print '(2a)', '# title: ', trim(cline)
      flush(6)
   endif

   !################# [Job] #################
   call get_value(table, "Job", group, requested=.False.)
   if (.not. associated(group)) then
      call get_value(table, "job", group, requested=.False.)
      if (.not. associated(group)) then
         print '(a)', context%report("[Job] section is required.", 0)
         return
      else
         print '(a)', context%report("[job] is deprecated. Please use [Job].", 0, level=toml_level%warning)
      endif
   endif

   !----------------- type -----------------
   call get_value(group, "type", cline, stat=istat, origin=origin)

   if (istat /= 0) then
      print '(a)', context%report('invalid type in [Job].', origin, "expected either DCD, CHECK_FORCE, or MD.")
      return
   endif

   if (cline == 'DCD') then
      job = JOBT%DCD

   else if (cline == 'CHECK_FORCE') then
      job = JOBT%CHECK_FORCE

   else if (cline == 'MD') then
      job = JOBT%MD

   else
      print '(a)', context%report('invalid type in [Job].', origin, "expected either DCD, CHECK_FORCE, or MD.")
      return
   endif
   print '(3a,i3,a)', '# Job type: ', trim(cline), ' (type = ', job,')'
   flush(6)

   !################# [Files] #################
   call get_value(table, "Files", group, requested=.False.)
   if (.not. associated(group)) then
      call get_value(table, "files", group, origin=origin, requested=.False.)

      if (.not. associated(group)) then
         print '(a)', context%report("[Files] section required.", 0)
         return
      else
         print '(a)', context%report("[files] is deprecated. Please use [Files].", origin, level=toml_level%warning)
      endif
   endif

   !================= [Files.In] =================
   call get_value(group, "In", node, requested=.False.)
   if (.not. associated(node)) then
      call get_value(group, "in", node, origin=origin, requested=.False.)

      if (.not. associated(group)) then
         print '(a)', context%report("[Files.In] section required.", 0)
         return
      else
         print '(a)', context%report("[files.in] is deprecated. Please use [Files.In].", origin, level=toml_level%warning)
      endif
   endif

   if (associated(node)) then
      call get_value(node, "ff", cfile_ff, stat=istat, origin=origin)
      if (istat /= 0 .or. .not. allocated(cfile_ff)) then
         print '(a)', context%report("invalid ff file name in [Files.In].", origin)
         return
      endif

      if (job == JOBT%DCD) then
         call get_value(node, "dcd", cfile_dcd_in, stat=istat, origin=origin)

         if (istat /= 0 .or. .not. allocated(cfile_dcd_in)) then
            print '(a)', context%report("invalid dcd file name in [Files.In]. Required when Job type = DCD.", origin)
            return
         endif
      endif

      call get_value(node, "pdb_ini", cfile_pdb_ini)
      call get_value(node, "xyz_ini", cfile_xyz_ini)
      call get_value(node, "fasta", cfile_fasta_in)
      call get_value(node, "ct", cfile_ct_in)
      call get_value(node, "bpseq", cfile_bpseq_in)
      call get_value(node, "anneal", cfile_anneal_in)

   else
      print '(a)', context%report("Files.In section required.", 0)
      return
   endif

   if (allocated(cfile_ct_in)) flg_in_ct = .True.
   if (allocated(cfile_bpseq_in)) flg_in_bpseq = .True.

   if (flg_in_ct .and. flg_in_bpseq) then
      print '(a)', 'Error: only one of ct and bpseq files can be specified.'
      return
   endif

   if (job == JOBT%MD) then
      if (len(cfile_pdb_ini) < 1 .and. len(cfile_xyz_ini) < 1) then
         print '(a)', 'Error: Initial structure is not specified. Either XYZ or PDB is required.'
         return

      else if (len(cfile_pdb_ini) > 0 .and. len(cfile_xyz_ini) > 0) then
         print '(a)', 'Error: Both XYZ and PDB are specified for the initial structure. Please use only one of them.'
         return

      endif
   endif


   !================= [Files.Out] =================
   call get_value(group, "Out", node, requested=.False.)
   if (.not. associated(node)) then
      call get_value(group, "out", node, origin=origin, requested=.False.)

      if (.not. associated(group)) then
         print '(a)', context%report("[Files.Out] section required.", 0)
         return
      else
         print '(a)', context%report("[files.out] is deprecated. Please use [Files.Out].", origin, level=toml_level%warning)
      endif
   endif

   if (associated(node)) then

      !----------------- prefix -----------------
      call get_value(node, "prefix", cfile_prefix, stat=istat, origin=origin)
      if (istat /= 0 .or. .not. allocated(cfile_prefix)) then
         print '(a)', context%report("invalid prefix in [Files.Out].", origin)
         return
      endif

      !call node%get_keys(list)
      !do i = 1, size(list)
      !   !call get_value(node, list(i)%key, cline, stat=stat)
      !   call get_value(node, list(i)%key, cline)
      !   print *, list(i)%key, cline
      !enddo

      !----------------- types -----------------
      call get_value(node, "types", array)

      do i = 1, len(array)
         call get_value(array, i, cline, stat=istat, origin=origin)

         if (cline == "bp") then
            flg_out_bp = .True.
         else if (cline == "bpall") then
            flg_out_bpall = .True.
         else if (cline == "bpe") then
            flg_out_bpe = .True.
         else
            print '(a)', context%report("invalid output type.", origin, "expected either bp, bpall, or bpe.")
            return
         endif
      enddo

   else
      print '(a)', context%report("Files.Out section required.", 0)
      return
   endif

   !################# [Condition] #################
   call get_value(table, "Condition", group, requested=.False.)
   if (.not. associated(group)) then
      call get_value(table, "condition", group, requested=.False.)

      if (.not. associated(group)) then
         print '(a)', context%report("[Condition] section is required.", 0)
         return
      else
         print '(a)', context%report("[condition] is deprecated. Please use [Condition].", 0, level=toml_level%warning)
      endif
   endif

   !----------------- rng_seed -----------------
   ! (optional, highly recommended)
   ! Seed value for the pseudorandom number generator.
   ! Set a 64-bit integer. If omitted, it will be set based on SYSTEM_CLOCK.
   call get_value(group, "rng_seed", rng_seed, stat=istat, origin=origin)
   if (istat /= 0) then
      print '(a)', context%report("rng_seed (random-number seed value) is not set in [condition].", 0, "a value is set by SYSTEM_CLOCK.", level=toml_level%warning)
      print '(a)', context%report("a value is set by SYSTEM_CLOCK.", origin=0, level=toml_level%warning)
      call SYSTEM_CLOCK(i)
      rng_seed = int(i, kind=INT64)
   endif
   print '(a,i16)', '# Condition, rng_seed: ', rng_seed

   !----------------- opt_anneal -----------------
   ! (optional) flag for simulated annealing.
   ! 0: No annealing (default)
   ! 1: Annealing ("anneal" is required in [files.in])
   call get_value(group, "opt_anneal", opt_anneal, 0, stat=istat, origin=origin)
   if (istat /= 0 .or. all(opt_anneal /= (/0, 1/))) then
      print '(a)', context%report("invalid opt_anneal in [condition].", origin, "expected either 0 or 1")
      return
   endif

   if (opt_anneal == 1 .and. .not. allocated(cfile_anneal_in)) then
      print '(a)', 'Error: opt_anneal requires anneal in [files.in].'
      return
   endif

   !----------------- tempK -----------------
   ! Specify the simulation temperature in Kelvin.
   ! This line will be ignored if opt_anneal = 1.
   tempK = -1.0
   call get_value(group, "tempK", tempK, stat=istat, origin=origin)
   if (istat /= 0 .or. (opt_anneal == 0 .and. tempK < 0.0)) then
      print '(a)', context%report("Error: invalid tempK in [condition]", origin, "expected real positive value")
      return
   endif

   kT = BOLTZ_KCAL_MOL * tempK
   print '(a,g15.8)', '# Condition, tempK: ', tempK

   !----------------- temp_independent -----------------
   ! (optional) temperature independent potential
   ! 0: Use the original (temperature-dependent) potential. (default)
   ! 1: Use temperature-independent potential. temp_ref is required.
   call get_value(group, "temp_independent", temp_independent, 0, stat=istat, origin=origin)

   if (istat /= 0 .or. all(temp_independent /= (/0, 1/))) then
      print '(a)', context%report("Error: invalid temp_independent in [condition]", origin, "expected either 0 or 1.")
      return
   endif

   if (temp_independent > 0) then
      print '(a,g15.8)', '# Condition, temp_independent: ', tempK

      call get_value(group, "tempK_ref", tempK_ref, stat=istat, origin=origin)
      if (istat /= 0 .or. tempK_ref <= 0.0_PREC) then
         print '(a)', context%report("Error: invalid tempK_ref in [condition]", origin, "expected positive real value.")
         return
      endif

   endif
   print '(a)', '#'

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

      !----------------- integrator -----------------
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

      !----------------- dt -----------------
      call get_value(group, "dt", dt, stat=istat, origin=origin)
      if (istat /= 0 .or. dt < 0.0) then
         print '(a)', context%report("invalid dt value in [MD].", origin, "expected a positive real value.")
         return
      endif
      print '(a,g15.8)', '# MD dt: ', dt

      !----------------- nstep -----------------
      call get_value(group, "nstep", nstep, stat=istat, origin=origin)
      if (istat /= 0 .or. nstep < 0) then
         print '(a)', context%report("invalid nstep value in [MD].", origin, "expected a positive integer value.")
         return
      endif
      print '(a,i16)', '# MD nstep: ', nstep

      !----------------- nstep_save -----------------
      call get_value(group, "nstep_save", nstep_save, stat=istat, origin=origin)
      if (istat /= 0 .or. nstep_save < 0) then
         print '(a)', context%report("invalid nstep_save value in [MD].", origin, "expected a positive integer value.")
         return
      endif
      print '(a,i16)', '# MD nstep_save: ', nstep_save

      !----------------- nstep_save_rst -----------------
      call get_value(group, "nstep_save_rst", nstep_save_rst, nstep_save, stat=istat, origin=origin)
      if (istat /= 0 .or. nstep_save_rst < 0) then
         print '(a)', context%report("invalid nstep_save_rst value in [MD].", origin, "expected a positive integer value.")
         return
      endif
      print '(a,i16)', '# MD nstep_save_rst: ', nstep_save_rst

      !----------------- neighbor_list_margin -----------------
      call get_value(group, "neighbor_list_margin", nl_margin, -1.0_PREC, stat=istat, origin=origin)
      if (istat /= 0) then
         print '(a)', context%report("invalid neighbor_list_margin value in [MD].", origin, "expected a positive real value.")
         return
      else if (nl_margin < 0.0) then
         nl_margin = 10.0_PREC
         print '(a)', context%report("neighbor_list_margin is not specified in [MD].", origin, "Default value 10.0 is used.", level=toml_level%warning)
      endif
      print '(a,g15.8)', '# MD neighbor_list_margin: ', nl_margin

      !----------------- viscosity_Pas -----------------
      call get_value(group, "viscosity_Pas", viscosity_Pas, -1.0_PREC, stat=istat, origin=origin)
      if (istat /= 0) then
         print '(a)', context%report("invalid viscosity_Pas value in [MD].", origin, "expected a positive real value.")
         return
      else if (viscosity_Pas < 0.0) then
         viscosity_Pas = 0.00001_PREC
         print '(a)', context%report("viscosity_Pas is not specified in [MD].", origin, "Default value 0.00001 is used.", level=toml_level%warning)
      endif
      print '(a,f12.8)', '# MD viscosity_Pas: ', viscosity_Pas

      !----------------- stop_wall_time_hour -----------------
      call get_value(group, "stop_wall_time_hour", rdummy, -1.0_PREC, stat=istat, origin=origin)

      if (istat /= 0) then
         print '(a)', context%report("Error: invalid stop_wall_time_hour in [MD]", origin, "expected real value.")
         return
      endif
      if (rdummy < 0.0) then
         stop_wall_time_sec = -1
         print '(a,g15.8,a)', '# MD stop_wall_time_hour: ', rdummy, ' (wall time limit not set)'

      else
         stop_wall_time_sec = int(rdummy * 3600.0_PREC, kind=INT64)
         print '(a,g15.8)', '# MD stop_wall_time_hour: ', rdummy
      endif

      !----------------- fix_com_origin -----------------
      call get_value(group, "fix_com_origin", fix_com_origin, 0, stat=istat, origin=origin)

      if (istat /= 0) then
         print '(a)', 'Error: invalid value for fix_com_origin in [MD].'
         return
      endif
      print '(a,i5)', '# MD fix_com_origin: ', fix_com_origin

      print '(a)', '#'

   endif

   !################# [Basepair] #################
   call get_value(table, "Basepair", group, requested=.False.)

   if (associated(group)) then

      ! bp_model
      call get_value(group, "model", bp_model, stat=istat, origin=origin)
      if (istat /= 0 .or. all(bp_model /= (/1, 3, 4, 5/))) then
         print '(a)', context%report("model is not specified in [Basepair].", origin, "expected either 1, 3, 4, or 5.")
         return
      else
         print '(a,i6)', '# Basepair, model: ', bp_model
      endif

      !----------------- max_bp_per_nt -----------------
      max_bp_per_nt = INVALID_INT_VALUE
      call get_value(group, "max_bp_per_nt", max_bp_per_nt, stat=istat, origin=origin)
      if (istat /= 0) then
         print '(a)', context%report('max_bp_per_nt is not specified in [Basepair]. Default value applies.', origin=0, level=toml_level%warning)
         max_bp_per_nt = -1   ! default
      endif

      !----------------- min_loop -----------------
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

   print '(a,i6)', '# Basepair, max_bp_per_nt: ', max_bp_per_nt
   print '(a,i6)', '# Basepair, min_loop: ', bp_min_loop
   print '(a)', '#'


   !################# [Electrostatic] #################
   ! (optional)
   ! If this section is given, electrostatic interaction is enabled.
   flg_ele = .False.
   call get_value(table, "Electrostatic", group, requested=.False.)

   if (associated(group)) then

      print '(a)', '# Electrostatic: On'
      flg_ele = .True.

      !----------------- ionic_strength -----------------
      ! Ionic strength of the monovalent-ions in molar units.
      call get_value(group, "ionic_strength", ionic_strength, stat=istat, origin=origin)
      if (istat /= 0 .or. ionic_strength <= 0.0_PREC) then
         print '(a)', context%report("invalid ionic_strength in [Electrostatic].", origin, "expected a positive real value.")
         return
      endif
      print '(a,g15.8)', '# Electrostatic, ionic strength: ', ionic_strength

      !----------------- cutoff_type -----------------
      ! How to specify the cutoff distance for electrostatic interactions.
      !     = 1: Distance-based.
      !          The cutoff will be specified as distance in Angstrom. (default)
      !     = 2: Multiple of the Debye length. 
      !          The cutoff will be a factor to be multiplied by Debye length.
      ele_cutoff_type = 1
      call get_value(group, "cutoff_type", ele_cutoff_type, 1, stat=istat, origin=origin)
      if (istat /= 0 .or. all(ele_cutoff_type /= (/1, 2/))) then
         print '(a)', context%report("invalid ele_cutoff_type in [Electrostatic].", origin, "expected either 1 or 2.")
         return
      endif

      if (ele_cutoff_type == 1) then
         print '(a)', "# Electrostatic, cutoff type: 1 (distance-based)"
      else if (ele_cutoff_type == 2) then
         print '(a)', "# Electrostatic, cutoff type: 2 cutoff will be multiplied by the Debye length"
      endif

      !----------------- cutoff -----------------
      ! Either diestance (cutoff_type = 1) or multiple (cutoff_type = 2),
      ! depending on the choice of cutoff_type.
      call get_value(group, "cutoff", ele_cutoff_inp, stat=istat, origin=origin)
      if (istat /= 0 .or. ele_cutoff_inp < 0.0_PREC) then
         print '(a)', context%report("invalid cutoff value in [Electrostatic].", origin, "expected a positive real value.")
         return
      endif
      print '(a,g15.8)', '# Electrostatic, cutoff: ', ele_cutoff_inp

      !----------------- length_per_charge -----------------
      ! Paremeter in ion-condensation theory in Angstrom.
      call get_value(group, "length_per_charge", length_per_charge, stat=istat, origin=origin)
      if (istat /= 0 .or. length_per_charge < 0.0_PREC) then
         print '(a)', context%report("invalid length_per_charge value in [Electrostatic].", origin, "expected a positive real value.")
         return
      endif
      print '(a,g15.8)', '# Electrostatic, length per charge: ', length_per_charge

      !----------------- no_charge -----------------
      ! (optional) array of positive integeres.
      ! Particles having no charges.
      call get_value(group, "no_charge", array)

      if (len(array) > 0) then
         print '(a)', "# Electrostatic, no charges on the following particles:"
         allocate(inp_no_charge(len(array)))
         do i = 1, len(array)
            call get_value(array, i, inp_no_charge(i), stat=istat, origin=origin)
            if (istat /= 0 .or. inp_no_charge(i) < 1) then
               print '(a)', context%report("invalid particle ID in no_charge, [Electrostatic].", origin, "expected a positive integer.")
               return
            endif
            print '(a, i8, 1x, i8)', '#         ', i, inp_no_charge(i)
         enddo
      endif

      print '(a)', '#'
   endif

   !################# [PBC_box] #################
   ! (optional)
   ! If this section is given, periodic boundary condition is enabled.
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

   !################# [variable_box] #################
   ! (optional)
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

   !################# [Stage] #################
   ! (optional)
   flg_stage = .False.
   call get_value(table, "Stage", group, requested=.False.)
   if (.not. associated(group)) then 
      call get_value(table, "stage", group, requested=.False.)
      if (associated(group)) then 
         print '(a)', context%report("[stage] is deprecated. Please use [Stage].", 0, level=toml_level%warning)
      endif
   endif

   if (associated(group)) then 
      flg_stage = .True.
      call get_value(group, "sigma", stage_sigma, stat=istat, origin=origin)
      if (istat /= 0) then
         print '(a)', context%report("invalid sigma value in [Stage].", origin, "expected a real value.")
         return
      endif

      call get_value(group, "epsilon", stage_eps, stat=istat, origin=origin)
      if (istat /= 0) then
         print '(a)', context%report("invalid epsilon value in [Stage].", origin, "expected a real value.")
         return
      endif

      print '(a,g15.8)', '# Stage sigma: ', stage_sigma
      print '(a,g15.8)', '# Stage epsilon: ', stage_eps
      print '(a)', '#'
   endif

   !################# [Progress] #################
   ! (optional)
   ! If this section is given, progress information will be output to STDOUT. 
   flg_progress = .False.
   call get_value(table, "Progress", group, requested=.False.)
   if (associated(group)) then
      flg_progress = .True.
      call get_value(group, "step", step_progress, stat=istat, origin=origin)
      if (istat /= 0 .or. step_progress < 1) then
         print '(a)', context%report("invalid step value in [Progress].", origin, "expected an integer value.")
         return
      endif
      print '(a,i16)', '# Progress step: ', step_progress
      print '(a)', '#'
   endif

   call table%destroy

   print '(a)', 'Done: reading input file'
   print *
   flush(6)

   stat = .True.

end subroutine read_input
