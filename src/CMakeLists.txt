set(SIS_EXE sis)

add_executable(${SIS_EXE}
    const.F90 const_idx.F90 const_phys.F90
    mt19937-64.F90 rnd_boxmuller.F90
    dcd.F90
    pbc.F90 progress.F90
    var_top.F90 var_state.F90 var_potential.F90 var_io.F90
    read_force_field.F90 read_input.F90 read_pdb.F90 read_xyz.F90 read_fasta.F90
    read_anneal.F90
    write_rst.F90 read_rst.F90
    list_local.F90 list_bp.F90 list_bp2.F90
    list_exv.F90 list_ele.F90
    set_sequence.F90 set_bp_map.F90
    set_ele.F90
    energy.F90
    energy_kinetic.F90
    energy_bond.F90 energy_angl.F90 energy_angl_ReB.F90
    energy_bp.F90 energy_bp_limit.F90 energy_bp_limit_triplet.F90
    energy_dih_cos.F90 energy_dih_exp.F90
    energy_wca.F90 energy_ele_DH.F90
    force.F90
    force_bond.F90 force_angl.F90 force_angl_ReB.F90
    force_bp.F90 force_bp_limit.F90 force_bp_limit_triplet.F90
    force_dih_cos.F90 force_dih_exp.F90
    force_wca.F90 force_ele_DH.F90
    neighbor_list.F90
    job_md.F90 job_dcd.F90 job_check_force.F90
    main.F90
)

add_dependencies(${SIS_EXE} toml-f-lib)
target_link_libraries(${SIS_EXE} toml-f-lib) 
target_link_libraries(${SIS_EXE} githash)

add_executable(show_rst
    const.F90 const_idx.F90
    show_rst.F90
)
