subroutine read_input(cfilepath)

   use,intrinsic :: ISO_FORTRAN_ENV, only: OUTPUT_UNIT
   use tomlf

   use const, only : PREC, L_INT, CHAR_FILE_PATH
   use const_phys, only : BOLTZ_KCAL_MOL, INVALID_JUDGE, INVALID_VALUE, INVALID_INT_JUDGE, INVALID_INT_VALUE
   use const_idx, only : JOBT, INTGRT, REPT
   use pbc, only : flg_pbc, set_pbc_size
   use var_io, only : iopen_hdl, &
                      flg_progress, step_progress, &
                      flg_out_bp, flg_out_bpe, flg_out_bpall, &
                      flg_in_ct, flg_in_bpseq, flg_in_fasta, flg_in_pdb, flg_in_xyz, &
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
   use var_parallel

   implicit none

   character(len=*), intent(in) :: cfilepath


   !======= TOML
   type(toml_table), allocatable :: table
   type(toml_table), pointer :: group, node
   type(toml_array), pointer :: array
   !======= 

   integer :: i
   integer :: istat
   integer :: hdl
   integer(L_INT) :: idummy
   real(PREC) :: rdummy, boxsize(3)
   character(len=:), allocatable :: cline
   character(len=5) :: cquery

   if (myrank == 0) then

      iopen_hdl = iopen_hdl + 1
      hdl = iopen_hdl

      print '(2a)', "Reading input file: ", trim(cfilePath)
      open(hdl, file=cfilepath, status='old', action='read', iostat=istat)

      if (istat /= 0) then
         print '(2a)', 'Error: failed to open the input file. ', trim(cfilepath)
         call sis_abort()
      endif

      call toml_parse(table, hdl)

      close(hdl)
      iopen_hdl = iopen_hdl - 1

      if (.not. allocated(table)) then
         print '(2a)', 'Error: could not obtain data from the toml input file.'
         call sis_abort()
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
         call sis_abort()
      endif

      !################# input files #################
      call get_value(table, "files", group)

      if (.not. associated(group)) then
         print '(a)', 'Error in input file: no files gorup in input.'
         call sis_abort()
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
         call sis_abort()
      endif

      if (allocated(cfile_fasta_in)) flg_in_fasta = .True.
      if (allocated(cfile_xyz_ini)) flg_in_xyz = .True.
      if (allocated(cfile_pdb_ini)) flg_in_pdb = .True.
      if (allocated(cfile_ct_in)) flg_in_ct = .True.
      if (allocated(cfile_bpseq_in)) flg_in_bpseq = .True.

      if (flg_in_ct .and. flg_in_bpseq) then
         print '(a)', 'Error: only one of ct and bpseq files can be specified.'
         call sis_abort()
      endif

      if (job == JOBT%MD .or. job == JOBT%DCD) then
         if (.not. flg_in_pdb .and. .not. flg_in_xyz) then
            print '(a)', 'Error: Initial structure is not specified. Either XYZ or PDB is required.'
            call sis_abort()

         else if (flg_in_pdb .and. flg_in_xyz) then
            print '(a)', 'Error: Both XYZ and PDB are specified for the initial structure. Please use only one of them.'
            call sis_abort()

         endif
      endif


      !################# output files #################
      call get_value(group, "out", node)
      if (associated(node)) then
         call get_value(node, "prefix", cline)
         cfile_prefix = trim(cline)

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
               call sis_abort()
            endif
         enddo

      else
         print '(a)', 'Error in input file: no files.out.'
         call sis_abort()
      endif

      !################# Replica #################
      flg_replica = .False.
      replica_values(:,:) = INVALID_VALUE
      n_replica_temp = -1

      call get_value(table, "replica", group, requested=.False.)

      if (associated(group)) then

         flg_replica = .True.

         n_replica_temp = 0
         call get_value(group, "n_replica_temp", n_replica_temp)
         call get_value(group, "nstep_exchange", nstep_rep_exchange)
         call get_value(group, "nstep_save", nstep_rep_save)

         flg_exchange = .True.
         call get_value(group, "exchange", flg_exchange)

         if (n_replica_temp > 0) then
            call get_value(group, "temperature", node, requested=.False.)
            if (associated(node)) then
               do i = 1, n_replica_temp
                  write(cquery, '(i0)') i
                  call get_value(node, cquery, replica_values(i, REPT%TEMP))
               enddo
            else
               print '(a)', 'Error in input file: [replica.temperature] is needed.'
               call sis_abort()
            endif

         else
            print '(a)', 'Error in input file: n_replica_temp has to be more than zero in [replica].'
            call sis_abort()
         endif

         do i = 1, n_replica_temp
            if (replica_values(i, REPT%TEMP) > INVALID_JUDGE) then
               print '(a,i4,a)', 'Error: Invalid value for replica(', i, ') in [replica.temperature].'
               call sis_abort()
            endif
         enddo

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
         call sis_abort()
      endif

      tempK = -1.0
      call get_value(group, "tempK", tempK)

      if (opt_anneal == 0 .and. tempK < 0.0) then
         print '(a)', 'Error: tempK is invalid or undefined in [condition].'
         call sis_abort()
      endif

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

         !###### integrator #######
         call get_value(group, "integrator", cline)

         if (.not. allocated(cline)) then
            print '(a)', 'Error: integrator is required in [MD].'
            call sis_abort()

         else if (cline == 'GJF-2GJ') then
            integrator = INTGRT%LD_GJF2GJ

         else
            print '(2a)', 'Error: Unknown integrator type, ', trim(cline)
            call sis_abort()
         endif

         !###### dt #######
         dt = -1.0
         call get_value(group, "dt", dt, stat=istat)
         if (istat /= 0 .or. dt < 0.0) then
            print '(a)', 'Error: invalid value for dt in [MD].'
            call sis_abort()
         endif

         !###### nstep #######
         nstep = -1
         call get_value(group, "nstep", nstep, stat=istat)
         if (istat /= 0 .or. nstep < 0) then
            print '(a)', 'Error: invalid value for nstep in [MD].'
            call sis_abort()
         endif

         !###### nstep_save #######
         nstep_save = -1
         call get_value(group, "nstep_save", nstep_save, stat=istat)
         if (istat /= 0 .or. nstep_save < 0) then
            print '(a)', 'Error: invalid value for nstep_save in [MD].'
            call sis_abort()
         endif

         !###### nstep_save_rst #######
         nstep_save_rst = nstep_save
         call get_value(group, "nstep_save_rst", nstep_save_rst, stat=istat)
         if (istat /= 0) then
            print '(a)', 'Error: invalid value for nstep_save_rst in [MD].'
            call sis_abort()
         endif

         !###### neighbor_list_margin ######
         nl_margin = -1.0
         call get_value(group, "neighbor_list_margin", nl_margin, stat=istat)
         if (istat /= 0) then
            print '(a)', 'Error: invalid value for neighbor_list_margin in [MD].'
            call sis_abort()
         else if (nl_margin < 0.0) then
            nl_margin = 10.0_PREC
            print '(a)', 'Warning: neighbor_list_margin is not specified in [MD] field. The default value will be used.'
         endif

         !###### viscosity_Pas ######
         viscosity_Pas = -1.0
         call get_value(group, "viscosity_Pas", viscosity_Pas, stat=istat)
         if (istat /= 0) then
            print '(a)', 'Error: invalid value for viscosity_Pas in [MD].'
            call sis_abort()
         else if (viscosity_Pas < 0.0) then
            viscosity_Pas = 0.00001_PREC
            print '(a)', 'Warning: viscosity_Pas is not specified in [MD] field. The default value will be used.'
         endif

         !###### stop_wall_time_hour ######
         rdummy = -1.0
         call get_value(group, "stop_wall_time_hour", rdummy, stat=istat)

         if (istat /= 0) then
            ! Try integer
            call get_value(group, "stop_wall_time_hour", idummy, stat=istat)

            if (istat /= 0) then
               print '(a)', 'Error: invalid value for stop_wall_time_hour in [MD].'
               call sis_abort()
            endif

            rdummy = real(idummy, kind=PREC)
         endif

         if (rdummy < 0.0) then
            stop_wall_time_sec = -1
         else
            stop_wall_time_sec = int(rdummy * 3600.0_PREC, kind=L_INT)
         endif

         !###### fix_com_origin ######
         fix_com_origin = 0
         call get_value(group, "fix_com_origin", fix_com_origin, stat=istat)

         if (istat /= 0) then
            print '(a)', 'Error: invalid value for fix_com_origin in [MD].'
            call sis_abort()
         endif

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
            call sis_abort()
         endif

         ! cutoff type   1: distance-based (default) 
         !               2: x Debye length
         ele_cutoff_type = 1
         call get_value(group, "cutoff_type", ele_cutoff_type)
         if (ele_cutoff_type /= 1 .and. ele_cutoff_type /= 2) then
            print '(a)', "Error: Invalid values for cutoff_type in [Electrostatic]."
            call sis_abort()
         endif

         ! cutoff
         ele_cutoff_inp = INVALID_VALUE
         call get_value(group, "cutoff", ele_cutoff_inp)
         if (ele_cutoff_inp > INVALID_JUDGE) then
            print '(a)', "Error: Invalid values or absence of cutoff in [Electrostatic]."
            call sis_abort()
         endif

         ! length_per_charge
         length_per_charge = INVALID_VALUE
         call get_value(group, "length_per_charge", length_per_charge)
         if (length_per_charge > INVALID_JUDGE) then
            print '(a)', "Error: Invalid value for length_per_charge in [Electrostatic]."
            call sis_abort()
         endif

         ! (optional) No charge particles
         call get_value(group, "no_charge", array)

         if (len(array) > 0) then
            allocate(inp_no_charge(len(array)))
            do i = 1, len(array)
               call get_value(array, i, inp_no_charge(i))
            enddo
         endif

      endif

      !################# box #################
      flg_pbc = .False.
      call get_value(table, "PBC_box", group, requested=.False.)
      if (associated(group)) then 
         flg_pbc = .True.
         call get_value(group, "x", boxsize(1))
         call get_value(group, "y", boxsize(2))
         call get_value(group, "z", boxsize(3))
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
      endif

      !################# Progress #################
      flg_progress = .False.
      call get_value(table, "Progress", group, requested=.False.)
      if (associated(group)) then
         flg_progress = .True.
         call get_value(group, "step", step_progress)
      endif

      call table%destroy

   endif ! myrank == 0

#ifdef PAR_MPI

   call MPI_BCAST(job, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, istat)

   !! File names (cfile_ff, cfile_dcd_in, cfile_pdb_ini, cfile_xyz_ini, cfile_fast_in, cfile_ct_in,
   !! cfile_bpseq_in, cfile_anneal_in) do not need to be sent becuase it will be read by myrank = 0

   call MPI_BCAST(flg_in_fasta, 1, MPI_LOGICAL, 0, MPI_COMM_WORLD, istat)
   call MPI_BCAST(flg_in_pdb, 1, MPI_LOGICAL, 0, MPI_COMM_WORLD, istat)
   call MPI_BCAST(flg_in_xyz, 1, MPI_LOGICAL, 0, MPI_COMM_WORLD, istat)
   call MPI_BCAST(flg_in_ct, 1, MPI_LOGICAL, 0, MPI_COMM_WORLD, istat)
   call MPI_BCAST(flg_in_bpseq, 1, MPI_LOGICAL, 0, MPI_COMM_WORLD, istat)

   call MPI_BCAST(cfile_prefix, CHAR_FILE_PATH, MPI_CHARACTER, 0, MPI_COMM_WORLD, istat)
   call MPI_BCAST(flg_out_bp, 1, MPI_LOGICAL, 0, MPI_COMM_WORLD, istat)
   call MPI_BCAST(flg_out_bpe, 1, MPI_LOGICAL, 0, MPI_COMM_WORLD, istat)
   call MPI_BCAST(flg_out_bpall, 1, MPI_LOGICAL, 0, MPI_COMM_WORLD, istat)

   call MPI_BCAST(flg_replica, 1, MPI_LOGICAL, 0, MPI_COMM_WORLD, istat)
   call MPI_BCAST(n_replica_temp, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, istat)
   call MPI_BCAST(nstep_rep_exchange, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, istat)
   call MPI_BCAST(nstep_rep_save, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, istat)
   call MPI_BCAST(flg_exchange, 1, MPI_LOGICAL, 0, MPI_COMM_WORLD, istat)

   call MPI_BCAST(replica_values, n_replica_temp, PREC_MPI, 0, MPI_COMM_WORLD, istat)

   call MPI_BCAST(rng_seed, L_INT, MPI_BYTE, 0, MPI_COMM_WORLD, istat)
   call MPI_BCAST(opt_anneal, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, istat)
   call MPI_BCAST(tempK, 1, PREC_MPI, 0, MPI_COMM_WORLD, istat)
   call MPI_BCAST(temp_independent, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, istat)
   call MPI_BCAST(temp_ref, 1, PREC_MPI, 0, MPI_COMM_WORLD, istat)

   call MPI_BCAST(nrepeat, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, istat)
   call MPI_BCAST(nchains, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, istat)

   call MPI_BCAST(integrator, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, istat)
   call MPI_BCAST(dt, 1, PREC_MPI, 0, MPI_COMM_WORLD, istat)
   call MPI_BCAST(nstep, L_INT, MPI_BYTE, 0, MPI_COMM_WORLD, istat)
   call MPI_BCAST(nstep_save, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, istat)
   call MPI_BCAST(nstep_save_rst, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, istat)
   call MPI_BCAST(stop_wall_time_sec, L_INT, MPI_BYTE, 0, MPI_COMM_WORLD, istat)

   call MPI_BCAST(nl_margin, 1, PREC_MPI, 0, MPI_COMM_WORLD, istat)
   call MPI_BCAST(viscosity_Pas, 1, PREC_MPI, 0, MPI_COMM_WORLD, istat)

   call MPI_BCAST(fix_com_origin, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, istat)

   call MPI_BCAST(bp_model, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, istat)
   call MPI_BCAST(max_bp_per_nt, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, istat)
   call MPI_BCAST(bp_min_loop, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, istat)

   call MPI_BCAST(flg_ele, 1, MPI_LOGICAL, 0, MPI_COMM_WORLD, istat)
   call MPI_BCAST(ionic_strength, 1, PREC_MPI, 0, MPI_COMM_WORLD, istat)
   call MPI_BCAST(ele_cutoff_type, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, istat)
   call MPI_BCAST(ele_cutoff_inp, 1, PREC_MPI, 0, MPI_COMM_WORLD, istat)
   call MPI_BCAST(length_per_charge, 1, PREC_MPI, 0, MPI_COMM_WORLD, istat)

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

#endif

   kT = BOLTZ_KCAL_MOL * tempK


   if (flg_pbc) call set_pbc_size(boxsize)

   if (myrank /= 0) then
      print '(a)', "Input data received via MPI."
   endif

   if (job == JOBT%DCD) then
      print '(a,i3,a)', '# job type: DCD (job = ', job,')'
   else if (job == JOBT%CHECK_FORCE) then
      print '(a,i3,a)', '# job type: CHECK_FORCE (job = ', job,')'
   else if (job == JOBT%MD) then
      print '(a,i3,a)', '# job type: MD (job = ', job,')'
   endif

   if (flg_replica) then
      print '(a,i16)', '# Replica, n_replica: ', n_replica_temp
      print '(a,i16)', '# Replica, nstep_exchange: ', nstep_rep_exchange
      print '(a,i16)', '# Replica, nstep_save: ', nstep_rep_save
      if (flg_exchange) then
         print '(a)', '# Replica, exchange: True'
      else
         print '(a)', '# Replica, exchange: False'
      endif
      print '(a)', '#'
   endif

   print '(a,g15.8)', '# Condition, tempK: ', tempK
   print '(a,i16)', '# Condition, rng_seed: ', rng_seed
   print '(a,i16)', '# Condition, temp_independent: ', temp_independent
   if (temp_independent /= 0) then
      print '(a,g15.8)', '# Condition, temp_ref: ', temp_ref
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
   endif

   print '(a,i5)', '# MD fix_com_origin: ', fix_com_origin
   print '(a)', '#'

   print '(a,i6)', '# Basepair, model: ', bp_model
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


   if (myrank == 0) then
      print '(a)', 'Done: reading input file'
   else
      print '(a)', 'Done: receiving input data'
   endif
   print *
   flush(OUTPUT_UNIT)

end subroutine read_input
