subroutine read_input(cfilepath)

   use,intrinsic :: ISO_FORTRAN_ENV, only: OUTPUT_UNIT, INT64
   use tomlf

   use const, only : PREC, L_INT, CHAR_FILE_PATH, MAX_REPLICA, MAX_REP_PER_DIM, MAX_REP_DIM
   use const_phys, only : BOLTZ_KCAL_MOL, JOUL2KCAL_MOL, &
                          INVALID_JUDGE, INVALID_VALUE, INVALID_INT_JUDGE, INVALID_INT_VALUE
   use const_idx, only : JOBT, INTGRT, REPT
   use pbc, only : flg_pbc, set_pbc_size
   use var_io, only : iopen_hdl, & !flg_gen_init_struct, &
                      flg_progress, step_progress, &
                      flg_out_bpcoef, flg_out_bp, flg_out_bpe, flg_out_bpall, &
                      flg_in_ct, flg_in_bpseq, flg_in_fasta, flg_in_pdb, flg_in_xyz, &
                      cfile_ff, cfile_dcd_in, &
                      cfile_prefix, cfile_pdb_ini, cfile_xyz_ini, cfile_fasta_in, cfile_anneal_in, &
                      cfile_ct_in, cfile_bpseq_in
   use var_state, only : job, tempK, kT, viscosity_Pas, opt_anneal, temp_independent,  tempK_ref, &
                         nstep, dt, nstep_save, nstep_save_rst, integrator, nl_margin, &
                         flg_variable_box, variable_box_step, variable_box_change, &
                         rng_seed, stop_wall_time_sec, nstep_check_stop, fix_com_origin, &
                         ionic_strength, length_per_charge, nstep_bp_MC
   use var_potential, only : flg_ele, ele_cutoff_type, ele_cutoff_inp, ele_exclude_covalent_bond_pairs, &
                             bp_min_loop, max_bp_per_nt, bp_model, &
                             flg_stage, stage_sigma, stage_eps, &
                             flg_twz, ntwz_DCF, twz_DCF_pairs, twz_DCF_direction, &
                             flg_bias_ss, bias_ss_force
   use var_top, only : nrepeat, nchains, inp_no_charge, &
                       flg_freeze, frz_ranges
   use var_replica, only : nrep, nstep_rep_exchange, nstep_rep_save, flg_exchange, &
                           replica_values, flg_replica, flg_repvar
   use var_parallel

   implicit none

   character(len=*), intent(in) :: cfilepath


   !======= TOML
   type(toml_table), allocatable :: table
   type(toml_table), pointer :: group, node
   type(toml_array), pointer :: array, nested_array
   type(toml_context) :: context
   type(toml_error), allocatable :: tm_err
   !======= 

   integer :: i
   integer :: nrepdim
   integer :: istat, origin
   integer :: hdl
   integer :: len_array
#ifdef PAR_MPI
   integer :: strlen
#endif
   real(PREC) :: rdummy, boxsize(3)
   character(len=:), allocatable :: cline
   character(len=5) :: cquery

   flush(OUTPUT_UNIT)

   if (myrank == 0) then

      iopen_hdl = iopen_hdl + 1
      hdl = iopen_hdl

      print '(2a)', "Reading input file: ", trim(cfilePath)
      open(hdl, file=cfilepath, status='old', action='read', iostat=istat)

      if (istat /= 0) then
         print '(2a)', 'Error: failed to open the input file. ', trim(cfilepath)
         call sis_abort()
      endif

      call toml_load(table, hdl, context=context, error=tm_err)

      close(hdl)
      iopen_hdl = iopen_hdl - 1

      if (allocated(tm_err)) then
         print '(a)', tm_err%message
         call sis_abort()
      endif

      !################# title #################
      ! (optional)
      call get_value(table, "title", cline)
      if (allocated(cline)) then
         print '(2a)', '# Title: ', trim(cline)
         flush(6)
      endif

      !################# [Job] #################
      call get_value(table, "Job", group, requested=.False.)
      if (.not. associated(group)) then
         print '(a)', context%report("[Job] section is required.", 0)
         call sis_abort()
      endif

      !----------------- type -----------------
      call get_value(group, "type", cline, stat=istat, origin=origin)

      if (istat /= 0) then
         print '(a)', context%report('[Job] type value is invalid.', origin, "expected either DCD, CHECK_FORCE, or MD.")
         call sis_abort()
      endif

      if (cline == 'DCD') then
         job = JOBT%DCD

      else if (cline == 'CHECK_FORCE') then
         job = JOBT%CHECK_FORCE

      else if (cline == 'MD') then
         job = JOBT%MD

      else
         print '(a)', context%report('[Job] type value is invalid.', origin, "expected either DCD, CHECK_FORCE, or MD.")
         call sis_abort()
      endif

      !################# [Files] #################
      call get_value(table, "Files", group, requested=.False.)
      if (.not. associated(group)) then
         print '(a)', context%report("[Files] section required.", 0)
         call sis_abort()
      endif

      !================= [Files.In] =================
      call get_value(group, "In", node, requested=.False.)
      if (.not. associated(node)) then
         print '(a)', context%report("[Files.In] section required.", 0)
         call sis_abort()
      endif

      if (associated(node)) then
         call get_value(node, "ff", cfile_ff, stat=istat, origin=origin)
         if (istat /= 0 .or. .not. allocated(cfile_ff)) then
            print '(a)', context%report("invalid ff file name in [Files.In].", origin)
            call sis_abort()
         endif

         if (job == JOBT%DCD) then
            call get_value(node, "dcd", cfile_dcd_in, stat=istat, origin=origin)

            if (istat /= 0 .or. .not. allocated(cfile_dcd_in)) then
               print '(a)', context%report("invalid dcd file name in [Files.In]. Required when Job type = DCD.", origin)
               call sis_abort()
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
         call sis_abort()
      endif

      if (allocated(cfile_pdb_ini)) flg_in_pdb = .True.
      if (allocated(cfile_xyz_ini)) flg_in_xyz = .True.
      if (allocated(cfile_fasta_in)) flg_in_fasta = .True.
      if (allocated(cfile_ct_in)) flg_in_ct = .True.
      if (allocated(cfile_bpseq_in)) flg_in_bpseq = .True.

      if (flg_in_ct .and. flg_in_bpseq) then
         print '(a)', 'Error: only one of ct and bpseq files can be specified.'
         call sis_abort()
      endif

      if (job == JOBT%MD) then
         if (len(cfile_pdb_ini) < 1 .and. len(cfile_xyz_ini) < 1) then
            print '(a)', 'Error: Initial structure is not specified. Either XYZ or PDB is required.'
            call sis_abort()

         else if (len(cfile_pdb_ini) > 0 .and. len(cfile_xyz_ini) > 0) then
            print '(a)', 'Error: Both XYZ and PDB are specified for the initial structure. Please use only one of them.'
            call sis_abort()

         endif
      endif


      !================= [Files.Out] =================
      call get_value(group, "Out", node, requested=.False.)
      if (.not. associated(node)) then
         print '(a)', context%report("[Files.Out] section required.", 0)
         call sis_abort()
      endif

      if (associated(node)) then

         !----------------- prefix -----------------
         call get_value(node, "prefix", cfile_prefix, stat=istat, origin=origin)
         if (istat /= 0 .or. .not. allocated(cfile_prefix)) then
            print '(a)', context%report("invalid prefix in [Files.Out].", origin)
            call sis_abort()
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

            if (cline == "bpcoef") then
               flg_out_bpcoef = .True.
            else if (cline == "bp") then
               flg_out_bp = .True.
            else if (cline == "bpall") then
               flg_out_bpall = .True.
            else if (cline == "bpe") then
               flg_out_bpe = .True.
            else
               print '(a)', context%report("invalid output type.", origin, "expected either bpcoef, bp, bpall, or bpe.")
               call sis_abort()
            endif
         enddo

      else
         print '(a)', context%report("Files.Out section required.", 0)
         call sis_abort()
      endif

      !################# [Condition] #################
      call get_value(table, "Condition", group, requested=.False.)
      if (.not. associated(group)) then
         print '(a)', context%report("[Condition] section is required.", 0)
         call sis_abort()
      endif

      !----------------- rng_seed -----------------
      ! (optional, highly recommended)
      ! Seed value for the pseudorandom number generator.
      ! Set a 64-bit integer. If omitted, it will be set based on SYSTEM_CLOCK.
      call get_value(group, "rng_seed", rng_seed, stat=istat, origin=origin)
      if (istat /= 0) then
         print '(a)', context%report("[Condition] rng_seed (random-number seed value) is not set.", 0, "a value is set by SYSTEM_CLOCK.", level=toml_level%warning)
         print '(a)', context%report("[Condition] rng_seed value is set by SYSTEM_CLOCK.", origin=0, level=toml_level%info)
         call SYSTEM_CLOCK(rng_seed)
      endif

      !----------------- opt_anneal -----------------
      ! (optional) flag for simulated annealing.
      ! 0: No annealing (default)
      ! 1: Annealing ("anneal" is required in [files.in])
      call get_value(group, "opt_anneal", opt_anneal, 0, stat=istat, origin=origin)
      if (istat /= 0 .or. all(opt_anneal /= (/0, 1/))) then
         print '(a)', context%report("invalid opt_anneal in [Condition].", origin, "expected either 0 or 1")
         call sis_abort()
      endif

      if (opt_anneal == 1 .and. .not. allocated(cfile_anneal_in)) then
         print '(a)', 'Error: opt_anneal requires anneal in [files.in].'
         call sis_abort()
      endif

      !----------------- tempK -----------------
      ! Specify the simulation temperature in Kelvin.
      ! This line will be ignored if opt_anneal = 1.
      tempK = -1.0
      call get_value(group, "tempK", tempK, stat=istat, origin=origin)
      if (istat /= 0 .or. (opt_anneal == 0 .and. tempK < 0.0)) then
         print '(a)', context%report("Error: invalid tempK in [Condition]", origin, "expected real positive value")
         call sis_abort()
      endif

      !----------------- temp_independent -----------------
      ! (optional) temperature independent potential
      ! 0: Use the original (temperature-dependent) potential. (default)
      ! 1: Use temperature-independent potential. temp_ref is required.
      call get_value(group, "temp_independent", temp_independent, 0, stat=istat, origin=origin)

      if (istat /= 0 .or. all(temp_independent /= (/0, 1/))) then
         print '(a)', context%report("Error: invalid temp_independent in [Condition]", origin, "expected either 0 or 1.")
         call sis_abort()
      endif

      if (temp_independent > 0) then
         call get_value(group, "tempK_ref", tempK_ref, stat=istat, origin=origin)
         if (istat /= 0 .or. tempK_ref <= 0.0_PREC) then
            print '(a)', context%report("Error: invalid tempK_ref in [Condition]", origin, "expected positive real value.")
            call sis_abort()
         endif
      endif

      !################# Repeat sequence #################
      if (.not. allocated(cfile_fasta_in)) then
         call get_value(table, "Repeat", group)
         if (associated(group)) then
            call get_value(group, "n_repeat", nrepeat)
            call get_value(group, "n_chain", nchains)
         endif
      else
         nrepeat = 0
      endif

      !################# MD #################
      if (job == JOBT%MD) then

         call get_value(table, "MD", group, requested=.false.)

         if (.not. associated(group)) then 
            print '(a)', 'Error: [MD] field required.'
            call sis_abort()
         endif

         !----------------- integrator -----------------
         call get_value(group, "integrator", cline)

         !----------------- dt -----------------
         call get_value(group, "dt", dt, stat=istat, origin=origin)
         if (istat /= 0 .or. dt < 0.0) then
            print '(a)', context%report("[MD] dt value is invalid.", origin, "expected a positive real value.")
            call sis_abort()
         endif

         !----------------- nstep -----------------
         call get_value(group, "nstep", nstep, stat=istat, origin=origin)
         if (istat /= 0 .or. nstep < 0) then
            print '(a)', context%report("[MD] nstep value is invalid.", origin, "expected a positive integer value.")
            call sis_abort()
         endif

         !----------------- nstep_save -----------------
         call get_value(group, "nstep_save", nstep_save, stat=istat, origin=origin)
         if (istat /= 0 .or. nstep_save < 0) then
            print '(a)', context%report("[MD] nstep_save value is invalid.", origin, "expected a positive integer value.")
            call sis_abort()
         endif

         !----------------- nstep_save_rst -----------------
         call get_value(group, "nstep_save_rst", nstep_save_rst, nstep_save, stat=istat, origin=origin)
         if (istat /= 0 .or. nstep_save_rst < 0) then
            print '(a)', context%report("[MD] nstep_save_rst value is invalid.", origin, "expected a positive integer value.")
            call sis_abort()
         endif

         !----------------- neighbor_list_margin -----------------
         call get_value(group, "neighbor_list_margin", nl_margin, -1.0_PREC, stat=istat, origin=origin)
         if (istat /= 0) then
            print '(a)', context%report("[MD] neighbor_list_margin value is invalid.", origin, "expected a positive real value.")
            call sis_abort()
         else if (nl_margin < 0.0) then
            nl_margin = 10.0_PREC
            print '(a)', context%report("[MD] neighbor_list_margin is not specified.", origin, "Default value 10.0 is used.", level=toml_level%warning)
         endif

         !----------------- viscosity_Pas -----------------
         call get_value(group, "viscosity_Pas", viscosity_Pas, -1.0_PREC, stat=istat, origin=origin)
         if (istat /= 0) then
            print '(a)', context%report("[MD] invalid viscosity_Pas value.", origin, "expected a positive real value.")
            call sis_abort()
         else if (viscosity_Pas < 0.0) then
            viscosity_Pas = 0.00001_PREC
            print '(a)', context%report("[MD] viscosity_Pas is not specified.", origin, "Default value 0.00001 is used.", level=toml_level%warning)
         endif

         !----------------- stop_wall_time_hour -----------------
         call get_value(group, "stop_wall_time_hour", rdummy, -1.0_PREC, stat=istat, origin=origin)

         if (istat /= 0) then
            print '(a)', context%report("Error: [MD] steop_wall_time_hour value is invalid", origin, "expected real value.")
            call sis_abort()
         endif
         if (rdummy < 0.0) then
            stop_wall_time_sec = -1
         else
            stop_wall_time_sec = int(rdummy * 3600.0_PREC, kind=INT64)
         endif

         nstep_check_stop = 100
         if (stop_wall_time_sec > 0) then
            call get_value(group, "nstep_check_stop", nstep_check_stop, int(100, kind=INT64), stat=istat, origin=origin)
         endif

         !----------------- fix_com_origin -----------------
         call get_value(group, "fix_com_origin", fix_com_origin, 0, stat=istat, origin=origin)

         if (istat /= 0) then
            print '(a)', 'Error: [MD] fix_com_origin value is invalid.'
            call sis_abort()
         endif

      endif

      !################# Replica #################
      flg_replica = .False.
      replica_values(:,:) = INVALID_VALUE
      nrep(:) = -1

      call get_value(table, "Replica", group, requested=.False.)

      if (associated(group)) then

         flg_replica = .True.

         nrep(REPT%TEMP) = 0
         call get_value(group, "nrep_temp", nrep(REPT%TEMP))

         if (nrep(REPT%TEMP) > MAX_REP_PER_DIM) then
            print '(a)', 'Error in input file: Number of replicas (nrep_temp) exceeds MAX_REP_PER_DIM in [Replica].'
            print '(a)', 'Please reduce the number of replicas or increase MAX_REP_PER_DIM in const.F90.'
            call sis_abort()
         endif

         nrep(REPT%TWZDCF) = 0
         call get_value(group, "nrep_force", nrep(REPT%TWZDCF))

         if (nrep(REPT%TWZDCF) > MAX_REP_PER_DIM) then
            print '(a)', 'Error in input file: Number of replicas (nrep_force) exceeds MAX_REP_PER_DIM in [Replica].'
            print '(a)', 'Please reduce the number of replicas or increase MAX_REP_PER_DIM in const.F90.'
            call sis_abort()
         endif

         nrep(REPT%ION) = 0
         call get_value(group, "nrep_ion", nrep(REPT%ION))

         if (nrep(REPT%ION) > MAX_REP_PER_DIM) then
            print '(a)', 'Error in input file: Number of replicas (nrep_ion) exceeds MAX_REP_PER_DIM in [Replica].'
            print '(a)', 'Please reduce the number of replicas or increase MAX_REP_PER_DIM in const.F90.'
            call sis_abort()
         endif

         if (nrep(REPT%TEMP) < 2 .and. nrep(REPT%TWZDCF) <2 .and. nrep(REPT%ION) < 2) then
            print '(a)', 'Error in input file: There should be more than one replica specified in any of nrep_temp, nrep_force, or nrep_ion'
            call sis_abort()
         endif

         call get_value(group, "nstep_exchange", nstep_rep_exchange)
         call get_value(group, "nstep_save", nstep_rep_save)

         flg_exchange = .True.
         call get_value(group, "exchange", flg_exchange)

         flg_repvar(:)= .False.

         nrepdim = 0
         !----------------- Replica.Temperature -----------------
         if (nrep(REPT%TEMP) > 1) then

            flg_repvar(REPT%TEMP) = .True.
            nrepdim = nrepdim + 1

            call get_value(group, "Temperature", node, requested=.False.)
            if (associated(node)) then
               do i = 1, nrep(REPT%TEMP)
                  write(cquery, '(i0)') i
                  call get_value(node, cquery, replica_values(i, REPT%TEMP))

                  if (replica_values(i, REPT%TEMP) > INVALID_JUDGE) then
                     print '(a,i4,a)', 'Error: Invalid value for replica(', i, ') in [Replica.Temperature].'
                     call sis_abort()
                  endif
               enddo
            else
               print '(a)', 'Error in input file: [Replica.Temperature] is required.'
               call sis_abort()
            endif

         else
            ! Ignore if nrep_temp = 0 or 1
            nrep(REPT%TEMP) = 1
         endif

         !----------------- Replica.Force -----------------
         if (nrep(REPT%TWZDCF) > 1) then

            flg_repvar(REPT%TWZDCF) = .True.
            nrepdim = nrepdim + 1

            call get_value(group, "Force", node, requested=.False.)
            if (associated(node)) then
               do i = 1, nrep(REPT%TWZDCF)
                  write(cquery, '(i0)') i
                  rdummy = INVALID_VALUE
                  call get_value(node, cquery, rdummy)

                  if (rdummy > INVALID_JUDGE) then
                     print '(a,i4,a)', 'Error: Invalid value for replica(', i, ') in [Replica.Force].'
                     call sis_abort()
                  endif

                  ! Convert the unit from pN to kcal/mol/A
                  replica_values(i, REPT%TWZDCF) = rdummy * (JOUL2KCAL_MOL * 1.0e-22)
               enddo
            else
               print '(a)', 'Error in input file: [Replica.Force] is required.'
               call sis_abort()
            endif

         else
            ! Ignore if nrep_force = 0 or 1
            nrep(REPT%TWZDCF) = 1
         endif

         !----------------- Replica.Ionic_strength -----------------
         if (nrep(REPT%ION) > 1) then

            flg_repvar(REPT%ION) = .True.
            nrepdim = nrepdim + 1

            call get_value(group, "Ionic_strength", node, requested=.False.)
            if (associated(node)) then
               do i = 1, nrep(REPT%ION)
                  write(cquery, '(i0)') i
                  call get_value(node, cquery, replica_values(i, REPT%ION))

                  if (replica_values(i, REPT%ION) > INVALID_JUDGE) then
                     print '(a,i4,a)', 'Error: Invalid value for replica(', i, ') in [Replica.Ionic_strength].'
                     call sis_abort()
                  endif
               enddo
            else
               print '(a)', 'Error in input file: [Replica.Ionic_strength] is required.'
               call sis_abort()
            endif

         else
            ! Ignore if nrep_temp = 0 or 1
            nrep(REPT%ION) = 1
         endif


         if (nrepdim > MAX_REP_DIM) then
            print '(a)', 'Error in input file: Number of replica dimensions exceeds MAX_REP_DIM in [Replica].'
            print '(a)', 'Please increase MAX_REP_DIM in const.F90.'
            call sis_abort()
         endif
      else
         nrep(:) = 1
         flg_exchange = .False.
      endif

      !################# [Basepair] #################
      call get_value(table, "Basepair", group, requested=.False.)

      if (associated(group)) then

         ! bp_model
         call get_value(group, "model", bp_model, stat=istat, origin=origin)
         if (istat /= 0 .or. all(bp_model /= (/1, 3, 4, 5/))) then
            print '(a)', context%report("[Basepair] model value is not specified.", origin, "expected either 1, 3, 4, or 5.")
            call sis_abort()
         endif

         !----------------- nstep_MC -----------------
         call get_value(group, "nstep_MC", nstep_bp_MC, stat=istat, origin=origin)
         if (istat /= 0) then
            print '(a)', context%report('[Basepair] nstep_MC is not specified. Default value applies.', origin=0, level=toml_level%warning)
            nstep_bp_MC = 1  ! Default
         else if (nstep_bp_MC < 0) then
            print '(a)', context%report("[Basepair] nstep_MC value is invalid.", origin, "nstep_MC has to be 0 or positive value.")
            call sis_abort()
         endif

         !----------------- max_bp_per_nt -----------------
         max_bp_per_nt = 0   ! Set 0 if no MC
         if (nstep_bp_MC > 0) then
            call get_value(group, "max_bp_per_nt", max_bp_per_nt, stat=istat, origin=origin)
            if (istat /= 0) then
               print '(a)', context%report('[Basepair] max_bp_per_nt is not specified. Default value applies.', origin=0, level=toml_level%warning)
               max_bp_per_nt = 1  ! Default
            else if (max_bp_per_nt <= 0) then
               print '(a)', context%report('[Basepair] max_bp_per_nt value is invalid.', origin=0, level=toml_level%error)
               call sis_abort()
            endif
         endif

         !----------------- min_loop -----------------
         call get_value(group, "min_loop", bp_min_loop, stat=istat, origin=origin)
         if (istat /= 0) then
            print '(a)', context%report('[Basepair] min_loop is not specified. Default value applies.', origin=0, level=toml_level%warning)
            bp_min_loop = 3  ! Default
         else if (bp_min_loop < 0) then
            print '(a)', context%report('[Basepair] min_loop value is invalid.', origin=0, level=toml_level%warning)
            call sis_abort()
         endif

      else
         print '(a)', '# [Basepair] section does not exist in the input file. Default values are used.'
         bp_model = 5
         nstep_bp_MC = 1
         max_bp_per_nt = 1  ! default
         bp_min_loop = 3    ! default
      endif

      !################# [Electrostatic] #################
      ! (optional)
      ! If this section is given, electrostatic interaction is enabled.
      flg_ele = .False.
      call get_value(table, "Electrostatic", group, requested=.False.)

      if (associated(group)) then
         flg_ele = .True.

         !----------------- ionic_strength -----------------
         ! Ionic strength of the monovalent-ions in molar units.
         call get_value(group, "ionic_strength", ionic_strength, stat=istat, origin=origin)
         if (istat /= 0 .or. ionic_strength <= 0.0_PREC) then
            print '(a)', context%report("[Electrostatic] ionic_strength value is invalid.", origin, "expected a positive real value.")
            call sis_abort()
         endif

         !----------------- cutoff_type -----------------
         ! How to specify the cutoff distance for electrostatic interactions.
         !     = 1: Distance-based.
         !          The cutoff will be specified as distance in Angstrom. (default)
         !     = 2: Multiple of the Debye length. 
         !          The cutoff will be a factor to be multiplied by Debye length.
         call get_value(group, "cutoff_type", ele_cutoff_type, 1, stat=istat, origin=origin)
         if (istat /= 0 .or. all(ele_cutoff_type /= (/1, 2/))) then
            print '(a)', context%report("[Electrostatic] ele_cutoff_type value is invalid.", origin, "expected either 1 or 2.")
            call sis_abort()
         endif

         !----------------- cutoff -----------------
         ! Either diestance (cutoff_type = 1) or multiple (cutoff_type = 2),
         ! depending on the choice of cutoff_type.
         call get_value(group, "cutoff", ele_cutoff_inp, stat=istat, origin=origin)
         if (istat /= 0 .or. ele_cutoff_inp < 0.0_PREC) then
            print '(a)', context%report("[Electrostatic] cutoff value is invalid.", origin, "expected a positive real value.")
            call sis_abort()
         endif

         !----------------- length_per_charge -----------------
         ! Paremeter in ion-condensation theory in Angstrom.
         call get_value(group, "length_per_charge", length_per_charge, stat=istat, origin=origin)
         if (istat /= 0 .or. length_per_charge < 0.0_PREC) then
            print '(a)', context%report("[Electrostatic] length_per_charge value is invalid.", origin, "expected a positive real value.")
            call sis_abort()
         endif

         !----------------- no_charge -----------------
         ! (optional) array of positive integeres.
         ! Particles having no charges.
         call get_value(group, "no_charge", array)

         if (len(array) > 0) then
            allocate(inp_no_charge(len(array)))
            do i = 1, len(array)
               call get_value(array, i, inp_no_charge(i), stat=istat, origin=origin)
               if (istat /= 0 .or. inp_no_charge(i) < 1) then
                  print '(a)', context%report("[Electrostatic] no_charge has invalid particle ID.", origin, "expected a positive integer.")
                  call sis_abort()
               endif
            enddo
         endif

         !----------------- exclude_covalent_bond_pairs -----------------
         ! (optional) true / false
         call get_value(group, "exclude_covalent_bond_pairs", ele_exclude_covalent_bond_pairs, .True., stat=istat, origin=origin)

         if (istat /= 0) then
            print '(a)', context%report("[Electrostatic] ele_exclude_covalent_bond_pairs value is invalid.", origin, "expected either true or false.")
            call sis_abort()
         endif

      endif

      !################# [PBC_box] #################
      ! (optional)
      ! If this section is given, periodic boundary condition is enabled.
      flg_pbc = .False.
      call get_value(table, "PBC_box", group, requested=.False.)
      if (associated(group)) then 
         flg_pbc = .True.
         call get_value(group, "x", boxsize(1))
         call get_value(group, "y", boxsize(2))
         call get_value(group, "z", boxsize(3))
      endif

      !################# [variable_box] #################
      ! (optional)
      flg_variable_box = .False.
      call get_value(table, "Variable_box", group, requested=.False.)
      if (associated(group)) then 
         flg_variable_box = .True.
         call get_value(group, "step", variable_box_step)
         call get_value(group, "change_x", variable_box_change(1))
         call get_value(group, "change_y", variable_box_change(2))
         call get_value(group, "change_z", variable_box_change(3))
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
            print '(a)', context%report("[Progress] step value is invalid.", origin, "expected an integer value.")
            call sis_abort()
         endif
      endif

      !################# [Stage] #################
      ! (optional)
      flg_stage = .False.
      call get_value(table, "Stage", group, requested=.False.)

      if (associated(group)) then 
         flg_stage = .True.
         call get_value(group, "sigma", stage_sigma, stat=istat, origin=origin)
         if (istat /= 0) then
            print '(a)', context%report("[Stage] sigma value is invalid.", origin, "expected a real value.")
            call sis_abort()
         endif

         call get_value(group, "epsilon", stage_eps, stat=istat, origin=origin)
         if (istat /= 0) then
            print '(a)', context%report("[Stage] epsilon value is invalid.", origin, "expected a real value.")
            call sis_abort()
         endif
      endif

      !################# [Freeze] #################
      ! (optional)
      flg_freeze = .False.
      call get_value(table, "Freeze", group, requested=.False.)

      if (associated(group)) then 
         flg_freeze = .True.

         call get_value(group, "id_ranges", array, stat=istat, origin=origin)
         if (istat /= 0) then
            print '(a)', context%report("invalid id_ranges value in [Freeze].", origin, "expected an array of integer pairs.")
            call sis_abort()
         endif
         if (associated(array)) then
            len_array = len(array)

            if (len_array > 0) then
               allocate(frz_ranges(2,len_array))
               do i = 1, len_array
                  call get_value(array, i, nested_array)
                  if (associated(nested_array)) then
                     if (len(nested_array) /= 2) then
                        print '(a)', context%report("invalid id_ranges value in [Freeze].", &
                                     origin, "array(s) have to have two elements each.")
                        call sis_abort()
                     endif
                     call get_value(nested_array, 1, frz_ranges(1,i))
                     call get_value(nested_array, 2, frz_ranges(2,i))
                     if ((frz_ranges(1,i) < 1) .or. (frz_ranges(2,i) < 1) .or. &
                         (frz_ranges(2,i) <= frz_ranges(1,i))) then
                        print '(a)', context%report("invalid id_ranges value in [Freeze].", &
                                     origin, "A range has to be in the form of [i,j] where 0 < i < j.")
                        call sis_abort()
                     endif
                  else
                     print '(a)', context%report("invalid id_ranges value in [Freeze].", &
                                  origin, "expected an array of integer pairs.")
                     call sis_abort()
                  end if
               end do
            else
               flg_freeze = .False.
            endif
         else
            flg_freeze = .False.
         end if
      endif

      !################# [Tweezers] #################
      ! (optional)
      ! If nrep_force > 1, this is necessary.
      flg_twz = .False.
      ntwz_DCF = 0
      call get_value(table, "Tweezers", group, requested=.False.)

      if (associated(group)) then

         !================= [Tweezers.Dual_Constnat_Force] =================
         call get_value(group, "Dual_Constant_Force", node, requested=.False.)
         if (associated(node)) then

            flg_twz = .True.

            call get_value(node, "id_pairs", array, stat=istat, origin=origin)
            if (istat /= 0) then
               print '(a)', context%report("invalid id_pairs value in [Tweezers.Dual_Constant_Force].", &
                            origin, "expected an array of integer pairs [imp1, imp2].")
               call sis_abort()
            endif
            if (associated(array)) then
               ntwz_DCF = len(array)

               if (ntwz_DCF > 0) then
                  allocate(twz_DCF_pairs(2, ntwz_DCF))
                  do i = 1, ntwz_DCF
                     call get_value(array, i, nested_array)
                     if (associated(nested_array)) then
                        if (len(nested_array) /= 2) then
                           print '(a)', context%report("invalid id_pairs value in [Tweezers.Dual_Constant_Force].", &
                                        origin, "array(s) have to have two elements each.")
                           call sis_abort()
                        endif
                        call get_value(nested_array, 1, twz_DCF_pairs(1, i))
                        call get_value(nested_array, 2, twz_DCF_pairs(2, i))
                        if ((twz_DCF_pairs(1, i) < 1) .or. (twz_DCF_pairs(2, i) < 1)) then
                           print '(a)', context%report("invalid id_pairs value in [Tweezers.Dual_Constant_Force].", &
                              origin, "expected an array of integer pairs.")
                           call sis_abort()
                        endif
                     else
                        print '(a)', context%report("invalid id_pairs value in [Tweezers.Dual_Constant_Force].", &
                                     origin, "expected an array of integer pairs.")
                        call sis_abort()
                     endif
                  enddo
               else
                  flg_twz = .False.
               endif
            else
               flg_twz = .False.
            endif

            ! Forces
            if (flg_twz) then
               call get_value(node, "forces_pN", array, stat=istat, origin=origin)
               if (istat /= 0) then
                  print '(a)', context%report("invalid forces_pN in [Tweezers.Dual_Constant_Force].", &
                               origin, "expected an array of vectors.")
                  call sis_abort()
               endif

               if (associated(array)) then
                  if (len(array) /= ntwz_DCF) then
                     print '(a)', context%report("invalid forces_pN in [Tweezers.Dual_Constant_Force].", &
                                  origin, "the array has to have the same numbe of vectors (fx, fy, yz) as id_pairs.")
                     call sis_abort()
                  endif

                  allocate(twz_DCF_direction(3, ntwz_DCF))

                  do i = 1, ntwz_DCF
                     call get_value(array, i, nested_array)
                     if (associated(nested_array)) then
                        if (len(nested_array) /= 3) then
                           print '(a)', context%report("invalid forces_pN vector in [Tweezers.Dual_Constant_Force].", &
                                        origin, "array(s) have to be a vector (fx, fy, fz).")
                           call sis_abort()
                        endif
                        call get_value(nested_array, 1, twz_DCF_direction(1,i))
                        call get_value(nested_array, 2, twz_DCF_direction(2,i))
                        call get_value(nested_array, 3, twz_DCF_direction(3,i))

                     else
                        print '(a)', context%report("invalid forces_pN vector in [Tweezers.Dual_Constant_Force].", &
                                     origin, "array(s) have to be a vector (fx, fy, fz).")
                        call sis_abort()
                     endif
                  enddo

                  ! Convert the unit from pN to kcal/mol/A
                  twz_DCF_direction(:,:) = twz_DCF_direction(:,:) * (JOUL2KCAL_MOL * 1.0e-22)
               endif
            endif
         endif

      endif

      ! If f-REMD
      if (flg_repvar(REPT%TWZDCF) .and. ntwz_DCF /= 1) then
         print '(a)', 'Error in input file: A force pair should be specified in [Tweezers.Dual_Constant_Force].'
         call sis_abort()
      endif

      !################# [Bias_SS] #################
      ! (optional)
      flg_bias_ss = .False.
      call get_value(table, "Bias_SS", group, requested=.False.)

      if (associated(group)) then
         flg_bias_ss = .True.

         if (.not. (flg_in_ct .or. flg_in_bpseq)) then
            print '(a)', "Incompatible setup in [Bias_SS]. To use Bias_SS, either a ct or bpseq file is required in [Files.In]."
            call sis_abort()
         endif

         !----------------- force_pN -----------------
         call get_value(group, "force_pN", bias_ss_force, stat=istat, origin=origin)
         if (istat /= 0 .or. bias_ss_force < 0.0_PREC) then
            print '(a)', context%report("invalid bias_ss_force value in [Bias_SS].", origin, "expected a positive real value.")
            call sis_abort()
         endif

         ! Convert the unit from pN to kcal/mol/A
         bias_ss_force = bias_ss_force * (JOUL2KCAL_MOL * 1.0e-22)
      endif

      call table%destroy

   endif ! myrank == 0


#ifdef PAR_MPI

   call MPI_BCAST(job, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, istat)

   !! File names (cfile_ff, cfile_dcd_in, cfile_pdb_ini, cfile_xyz_ini, cfile_fast_in, cfile_ct_in,
   !! cfile_bpseq_in, cfile_anneal_in) do not need to be sent becuase it will be read by myrank = 0

   call MPI_BCAST(flg_in_fasta, 1, MPI_LOGICAL, 0, MPI_COMM_WORLD, istat)
   !call MPI_BCAST(flg_gen_init_struct, 1, MPI_LOGICAL, 0, MPI_COMM_WORLD, istat)
   call MPI_BCAST(flg_in_pdb, 1, MPI_LOGICAL, 0, MPI_COMM_WORLD, istat)
   call MPI_BCAST(flg_in_xyz, 1, MPI_LOGICAL, 0, MPI_COMM_WORLD, istat)
   call MPI_BCAST(flg_in_ct, 1, MPI_LOGICAL, 0, MPI_COMM_WORLD, istat)
   call MPI_BCAST(flg_in_bpseq, 1, MPI_LOGICAL, 0, MPI_COMM_WORLD, istat)

   ! Because the length of cfile_prefix is unknown...
   if (myrank == 0) strlen = len(cfile_prefix)
   call MPI_BCAST(strlen, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, istat)
   if (myrank /= 0) allocate(character(len=strlen)::cfile_prefix)

   call MPI_BCAST(cfile_prefix, strlen, MPI_CHARACTER, 0, MPI_COMM_WORLD, istat)
   call MPI_BCAST(flg_out_bpcoef, 1, MPI_LOGICAL, 0, MPI_COMM_WORLD, istat)
   call MPI_BCAST(flg_out_bp, 1, MPI_LOGICAL, 0, MPI_COMM_WORLD, istat)
   call MPI_BCAST(flg_out_bpe, 1, MPI_LOGICAL, 0, MPI_COMM_WORLD, istat)
   call MPI_BCAST(flg_out_bpall, 1, MPI_LOGICAL, 0, MPI_COMM_WORLD, istat)

   call MPI_BCAST(flg_replica, 1, MPI_LOGICAL, 0, MPI_COMM_WORLD, istat)
   call MPI_BCAST(flg_repvar, REPT%MAX, MPI_LOGICAL, 0, MPI_COMM_WORLD, istat)
   call MPI_BCAST(nrep, REPT%MAX, MPI_INTEGER, 0, MPI_COMM_WORLD, istat)
   call MPI_BCAST(nstep_rep_exchange, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, istat)
   call MPI_BCAST(nstep_rep_save, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, istat)
   call MPI_BCAST(flg_exchange, 1, MPI_LOGICAL, 0, MPI_COMM_WORLD, istat)

   call MPI_BCAST(replica_values, MAX_REP_PER_DIM*REPT%MAX, PREC_MPI, 0, MPI_COMM_WORLD, istat)

   call MPI_BCAST(rng_seed, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, istat)
   call MPI_BCAST(opt_anneal, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, istat)
   call MPI_BCAST(tempK, 1, PREC_MPI, 0, MPI_COMM_WORLD, istat)
   call MPI_BCAST(temp_independent, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, istat)
   call MPI_BCAST(tempK_ref, 1, PREC_MPI, 0, MPI_COMM_WORLD, istat)

   call MPI_BCAST(nrepeat, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, istat)
   call MPI_BCAST(nchains, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, istat)

   call MPI_BCAST(integrator, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, istat)
   call MPI_BCAST(dt, 1, PREC_MPI, 0, MPI_COMM_WORLD, istat)
   call MPI_BCAST(nstep, L_INT, MPI_BYTE, 0, MPI_COMM_WORLD, istat)
   call MPI_BCAST(nstep_save, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, istat)
   call MPI_BCAST(nstep_save_rst, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, istat)
   call MPI_BCAST(stop_wall_time_sec, L_INT, MPI_BYTE, 0, MPI_COMM_WORLD, istat)
   call MPI_BCAST(nstep_check_stop, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, istat)

   call MPI_BCAST(nl_margin, 1, PREC_MPI, 0, MPI_COMM_WORLD, istat)
   call MPI_BCAST(viscosity_Pas, 1, PREC_MPI, 0, MPI_COMM_WORLD, istat)

   call MPI_BCAST(fix_com_origin, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, istat)

   call MPI_BCAST(bp_model, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, istat)
   call MPI_BCAST(nstep_bp_MC, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, istat)
   call MPI_BCAST(max_bp_per_nt, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, istat)
   call MPI_BCAST(bp_min_loop, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, istat)

   call MPI_BCAST(flg_ele, 1, MPI_LOGICAL, 0, MPI_COMM_WORLD, istat)
   call MPI_BCAST(ionic_strength, 1, PREC_MPI, 0, MPI_COMM_WORLD, istat)
   call MPI_BCAST(ele_cutoff_type, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, istat)
   call MPI_BCAST(ele_cutoff_inp, 1, PREC_MPI, 0, MPI_COMM_WORLD, istat)
   call MPI_BCAST(length_per_charge, 1, PREC_MPI, 0, MPI_COMM_WORLD, istat)
   call MPI_BCAST(ele_exclude_covalent_bond_pairs, 1, MPI_LOGICAL, 0, MPI_COMM_WORLD, istat)

   if (myrank == 0) then
      i = 0
      if (allocated(inp_no_charge)) then
         i = size(inp_no_charge)
      endif
   endif
   call MPI_BCAST(i, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, istat)

   if (i > 0) then
      if (myrank /= 0) then
         allocate(inp_no_charge(i))
      endif
      call MPI_BCAST(inp_no_charge, i, MPI_INTEGER, 0, MPI_COMM_WORLD, istat)
   endif

   call MPI_BCAST(flg_pbc, 1, MPI_LOGICAL, 0, MPI_COMM_WORLD, istat)
   call MPI_BCAST(boxsize, 3, PREC_MPI, 0, MPI_COMM_WORLD, istat)

   call MPI_BCAST(flg_variable_box, 1, MPI_LOGICAL, 0, MPI_COMM_WORLD, istat)
   call MPI_BCAST(variable_box_step, L_INT, MPI_BYTE, 0, MPI_COMM_WORLD, istat)
   call MPI_BCAST(variable_box_change, 3, PREC_MPI, 0, MPI_COMM_WORLD, istat)

   call MPI_BCAST(flg_progress, 1, MPI_LOGICAL, 0, MPI_COMM_WORLD, istat)
   call MPI_BCAST(step_progress, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, istat)

   call MPI_BCAST(flg_stage, 1, MPI_LOGICAL, 0, MPI_COMM_WORLD, istat)
   call MPI_BCAST(stage_sigma, 1, PREC_MPI, 0, MPI_COMM_WORLD, istat)
   call MPI_BCAST(stage_eps, 1, PREC_MPI, 0, MPI_COMM_WORLD, istat)

   call MPI_BCAST(flg_freeze, 1, MPI_LOGICAL, 0, MPI_COMM_WORLD, istat)
   if (flg_freeze) then
      call MPI_BCAST(frz_ranges, size(frz_ranges), MPI_INTEGER, 0, MPI_COMM_WORLD, istat)
   endif

   call MPI_BCAST(flg_twz, 1, MPI_LOGICAL, 0, MPI_COMM_WORLD, istat)
   if (flg_twz) then
      call MPI_BCAST(ntwz_DCF, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, istat)
      if (ntwz_DCF > 0) then
         if (myrank /= 0) then
            allocate(twz_DCF_pairs(2, ntwz_DCF))
            allocate(twz_DCF_direction(3, ntwz_DCF))
         endif
         call MPI_BCAST(twz_DCF_pairs, 2*ntwz_DCF, MPI_INTEGER, 0, MPI_COMM_WORLD, istat)
         call MPI_BCAST(twz_DCF_direction, 3*ntwz_DCF, PREC_MPI, 0, MPI_COMM_WORLD, istat)
      endif
   endif

   call MPI_BCAST(flg_bias_ss, 1, MPI_LOGICAL, 0, MPI_COMM_WORLD, istat)
   if (flg_bias_ss) then
      call MPI_BCAST(bias_ss_force, 1, PREC_MPI, 0, MPI_COMM_WORLD, istat)
   endif

#endif

   kT = BOLTZ_KCAL_MOL * tempK

   if (flg_pbc) call set_pbc_size(boxsize)

   if (myrank /= 0) then
      print '(a)', "Input data received via MPI."
   endif

   if (job == JOBT%DCD) then
      print '(a,i3,a)', '# Job type: DCD (job = ', job,')'
   else if (job == JOBT%CHECK_FORCE) then
      print '(a,i3,a)', '# Job type: CHECK_FORCE (job = ', job,')'
   else if (job == JOBT%MD) then
      print '(a,i3,a)', '# Job type: MD (job = ', job,')'
   endif
   print '(a)', '#'

   if (flg_replica) then
      print '(a,i16)', '# Replica, nrep_temp: ', nrep(REPT%TEMP)
      print '(a,i16)', '# Replica, nrep_force: ', nrep(REPT%TWZDCF)
      print '(a,i16)', '# Replica, nrep_ion: ', nrep(REPT%ION)
      print '(a,i16)', '# Replica, nstep_exchange: ', nstep_rep_exchange
      print '(a,i16)', '# Replica, nstep_save: ', nstep_rep_save
      print '(a,L1)', '# Replica, exchange: ', flg_exchange
      print '(a)', '#'
   endif

   print '(a,a)', '# Files.In, ff: ', trim(cfile_ff)
   if (job == JOBT%DCD) then
      print '(a,a)', '# Files.In, dcd: ', trim(cfile_dcd_in)
   endif
   if (flg_in_pdb) print '(a,a)', '# Files.In, pdb: ', trim(cfile_pdb_ini)
   if (flg_in_xyz) print '(a,a)', '# Files.In, xyz: ', trim(cfile_xyz_ini)
   if (flg_in_fasta) print '(a,a)', '# Files.In, fasta: ', trim(cfile_fasta_in)
   if (flg_in_ct) print '(a,a)', '# Files.In, ct: ', trim(cfile_ct_in)
   if (flg_in_bpseq) print '(a,a)', '# Files.In, bpseq: ', trim(cfile_bpseq_in)
   if (myrank == 0) then
      if (opt_anneal > 0) print '(a,a)', '# Files.In, anneal: ', trim(cfile_anneal_in)
   endif
   print '(a)', '#'

   print '(a,a)', '# Files.Out, prefix: ', trim(cfile_prefix)
   print '(a,L)', '# Files.Out, bpcoef: ', flg_out_bpcoef
   print '(a,L)', '# Files.Out, bp: ', flg_out_bp
   print '(a,L)', '# Files.Out, bpall: ', flg_out_bpall
   print '(a,L)', '# Files.Out, bpe: ', flg_out_bpe
   print '(a)', '#'

   print '(a,g15.8)', '# Condition, tempK: ', tempK
   print '(a,i16)', '# Condition, rng_seed: ', rng_seed
   print '(a,i16)', '# Condition, temp_independent: ', temp_independent
   if (temp_independent /= 0) then
      print '(a,g15.8)', '# Condition, tempK_ref: ', tempK_ref
   endif
   print '(a)', '#'

   if (nrepeat > 0) then
      print '(a,i10)', '# Repeat n_repeat: ', nrepeat
      print '(a,i10)', '# Repeat n_chains: ', nchains
      print '(a)', '#'
   endif

   if (integrator == INTGRT%LD_GJF2GJ) then
      print '(a)', '# MD integrator: GJF-2GJ'
   endif
   print '(a,g15.8)', '# MD dt: ', dt
   print '(a,i16)', '# MD nstep: ', nstep
   print '(a,i16)', '# MD nstep_save: ', nstep_save
   print '(a,i16)', '# MD nstep_save_rst: ', nstep_save_rst
   print '(a,g15.8)', '# MD neighbor_list_margin: ', nl_margin
   print '(a,f12.8)', '# MD viscosity_Pas: ', viscosity_Pas

   if (stop_wall_time_sec < 0) then
      print '(a,g15.8,a)', '# MD stop_wall_time_hour: -1 (wall time limit not set)'
   else
      print '(a,g15.8)', '# MD stop_wall_time_hour: ', real(stop_wall_time_sec, kind=PREC) / 3600.0_PREC
      print '(a,g15.8)', '# MD nstep_check_stop: ', nstep_check_stop
   endif

   print '(a,i5)', '# MD fix_com_origin: ', fix_com_origin
   print '(a)', '#'

   print '(a,i6)', '# Basepair, model: ', bp_model
   print '(a,i6)', '# Basepair, nstep_MC: ', nstep_bp_MC
   print '(a,i6)', '# Basepair, max_bp_per_nt: ', max_bp_per_nt
   print '(a,i6)', '# Basepair, min_loop: ', bp_min_loop
   print '(a)', '#'

   if (flg_ele) then
      print '(a,g15.8)', '# Electrostatic, ionic strength: ', ionic_strength

      if (ele_cutoff_type == 1) then
         print '(a)', "# Electrostatic, cutoff type: 1 (distance-based)"
      else if (ele_cutoff_type == 2) then
         print '(a)', "# Electrostatic, cutoff type: 2 cutoff will be multiplied by the Debye length"
      endif

      print '(a,g15.8)', '# Electrostatic, cutoff: ', ele_cutoff_inp
      print '(a,g15.8)', '# Electrostatic, length per charge: ', length_per_charge
      if (allocated(inp_no_charge)) then
         print '(a)', "# Electrostatic, no charges on the following particles:"
         do i = 1, size(inp_no_charge)
            print '(a, i8, 1x, i8)', '#         ', i, inp_no_charge(i)
         enddo
      endif
      print '(a,l1)', '# Electrostatic, exclude_covalent_bond_pairs: ', ele_exclude_covalent_bond_pairs
      print '(a)', '#'
   endif

   if (flg_pbc) then
      print '(a,g15.8)', '# pbc_box x: ', boxsize(1)
      print '(a,g15.8)', '# pbc_box y: ', boxsize(2)
      print '(a,g15.8)', '# pbc_box z: ', boxsize(3)
      print '(a)', '#'
   endif

   if (flg_variable_box) then
      print '(a,i16)', '# variable_box step: ', variable_box_step
      print '(a,g15.8)', '# variable_box change_x: ', variable_box_change(1)
      print '(a,g15.8)', '# variable_box change_y: ', variable_box_change(2)
      print '(a,g15.8)', '# variable_box change_z: ', variable_box_change(3)
      print '(a)', '#'
   endif

   if (flg_progress) then
      print '(a,i16)', '# Progress step: ', step_progress
      print '(a)', '#'
   endif

   if (flg_stage) then
      print '(a,g15.8)', '# Stage sigma: ', stage_sigma
      print '(a,g15.8)', '# Stage epsilon: ', stage_eps
      print '(a)', '#'
   endif

   if (flg_freeze) then
      print '(a)', '# Freeze: On'
      print '(a,i5)', '# Freeze, number of specified ranges:', size(frz_ranges, 2)
      do i = 1, size(frz_ranges,2)
         print '(a,i5,a,i5,x,i5)', '# Freeze, range(',i,'): ', frz_ranges(1,i), frz_ranges(2,i)
      enddo
      print '(a)', '#'
   endif

   if (flg_twz) then
      print '(a)', '# Tweezers: On'
      print '(a,i5)', '# Tweezers, number of pairs:', ntwz_DCF
      print '(a)', '# Tweezers, mode pair  imp1  imp2   fx      fy      fz'
      do i = 1, ntwz_DCF
         print '(a,i3,x,i5,x,i5,3(x,f7.2))', '# Tweezers,  DCF ',i, twz_DCF_pairs(1,i), twz_DCF_pairs(2,i),&
               twz_DCF_direction(1,i) / (JOUL2KCAL_MOL*1.0e-22), twz_DCF_direction(2,i) / (JOUL2KCAL_MOL*1.0e-22), &
               twz_DCF_direction(3,i) / (JOUL2KCAL_MOL*1.0e-22)  ! output in pN
      enddo
      print '(a)', '#'
   endif

   if (flg_bias_ss) then
      print '(a)', '# Bias_SS: On'
      print '(a,x,f7.2)', '# Bias_SS, force', bias_ss_force / (JOUL2KCAL_MOL*1.0e-22)  ! in pN
      print '(a)', '#'
   endif

   if (myrank == 0) then
      print '(a)', 'Done: reading input file'
   else
      print '(a)', 'Done: receiving input data'
   endif
   print *
   flush(OUTPUT_UNIT)

end subroutine read_input
