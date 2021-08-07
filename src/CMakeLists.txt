set(SIS_EXE sis)

add_executable(${SIS_EXE}
    const.F90 const_idx.F90 const_phys.F90
    var_top.F90 var_state.F90 var_potential.F90 var_io.F90
    dcd.F90 read_force_field.F90 read_input.F90 read_pdb.F90
    list_local.F90 list_bp.F90
    energy.F90
    energy_bond.F90 energy_angl.F90 energy_bp.F90 read_sisinfo.F90
    force.F90
    force_bond.F90 force_angl.F90
    job_dcd.F90 job_check_force.F90
    main.F90
)

add_dependencies(${SIS_EXE} toml-f-lib)
target_link_libraries(${SIS_EXE} toml-f-lib) 
