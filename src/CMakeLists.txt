set(SIS_EXE sis)

option(BP_HALT_IEEE_EXCEPTIONS "Use ieee halting mode functions in write_bp" OFF)
if (BP_HALT_IEEE_EXCEPTIONS)
    set(CMAKE_Fortran_FLAGS_DEBUG   "${CMAKE_Fortran_FLAGS_DEBUG} -DBP_HALT_IEEE_EXCEPTIONS")
    set(CMAKE_Fortran_FLAGS_RELEASE   "${CMAKE_Fortran_FLAGS_RELEASE} -DBP_HALT_IEEE_EXCEPTIONS")
    set(CMAKE_Fortran_FLAGS_RELWITHDEBINFO   "${CMAKE_Fortran_FLAGS_RELWITHDEBINFO} -DBP_HALT_IEEE_EXCEPTIONS")
endif (BP_HALT_IEEE_EXCEPTIONS)

option(ALLOW_ARGUMENT_MISMATCH "Use -fallow-argument-mismatch option" OFF)
if (ALLOW_ARGUMENT_MISMATCH)
    set(CMAKE_Fortran_FLAGS_DEBUG   "${CMAKE_Fortran_FLAGS_DEBUG} -fallow-argument-mismatch")
    set(CMAKE_Fortran_FLAGS_RELEASE   "${CMAKE_Fortran_FLAGS_RELEASE} -fallow-argument-mismatch")
    set(CMAKE_Fortran_FLAGS_RELWITHDEBINFO   "${CMAKE_Fortran_FLAGS_RELWITHDEBINFO} -fallow-argument-mismatch")
endif (ALLOW_ARGUMENT_MISMATCH)

add_executable(${SIS_EXE}
    const.F90 const_idx.F90 const_phys.F90
    rnd_boxmuller.F90
    dcd.F90
    pbc.F90 progress.F90
    var_parallel.F90 var_replica.F90
    var_top.F90 var_state.F90 var_potential.F90 var_io.F90
    read_force_field.F90 read_input.F90 read_pdb.F90 read_xyz.F90 read_fasta.F90
    read_anneal.F90 read_ss.F90 read_timed_bias_rg.F90 read_restraint.F90
    write_out.F90
    write_rep_table.F90
    write_rst.F90 read_rst.F90
    write_bp.F90
    write_bpcoef.F90
    init_mts.F90 init_out_files.F90 init_replica.F90 init_ele.F90 init_freeze.F90
    init_tweezers.F90 init_bias_rg.F90
    init_sequence.F90 init_structure.F90 gen_random_coil.F90
    list_local.F90 list_bp.F90
    list_exv.F90 list_ele.F90
    init_bp.F90 set_ele.F90
    replica_exchange.F90
    energy_replica.F90 energy_sumup.F90
    energy_kinetic.F90
    energy_bond.F90 energy_angl.F90 energy_angl_ReB.F90
    energy_bp.F90 energy_bp_limit.F90 energy_bp_triplet.F90 energy_bp_limit_triplet.F90
    energy_dih_cos.F90 energy_dih_exp.F90
    energy_wca.F90 energy_ele_DH.F90 energy_stage.F90 energy_tweezers.F90
    energy_rg.F90 energy_restraint.F90
    force.F90
    force_bond.F90 force_angl.F90 force_angl_ReB.F90
    force_bp.F90 force_bp_limit.F90 force_bp_triplet.F90 force_bp_limit_triplet.F90
    force_dih_cos.F90 force_dih_exp.F90
    force_wca.F90 force_ele_DH.F90 force_stage.F90 force_tweezers.F90
    force_rg.F90 force_restraint.F90
    neighbor_list.F90
    job_md.F90 job_dcd.F90 job_check_force.F90
    sis_abort.F90
    warning_stage.F90
    main.F90
)

add_dependencies(${SIS_EXE} toml-f-lib)
add_dependencies(${SIS_EXE} mt_stream-lib)

target_link_libraries(${SIS_EXE} toml-f-lib) 
target_link_libraries(${SIS_EXE} githash)
target_link_libraries(${SIS_EXE} mt_stream-lib)

add_executable(show_rst
    const.F90 const_idx.F90
    show_rst.F90
)
add_dependencies(show_rst mt_stream-lib)
target_link_libraries(show_rst mt_stream-lib)
