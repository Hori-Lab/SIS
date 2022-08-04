subroutine energy()

   use const
   use const_idx, only : ENE
   use var_state, only : energies
   use var_potential, only : flg_ele, max_bp_per_nt, flg_dihedral

   implicit none

   energies(:) = 0.0e0_PREC

   call energy_bond(energies(ENE%BOND))

   call energy_angl(energies(ENE%ANGL))
   
   if (flg_dihedral) call energy_dihedral(energies(ENE%DIHE))

   if (max_bp_per_nt < 1) then
      call energy_bp(energies(ENE%BP))
   else
      call energy_bp_limit(energies(ENE%BP))
   endif

   call energy_wca(energies(ENE%EXV))

   if (flg_ele) call energy_ele_DH(energies(ENE%ELE))

   energies(ENE%TOTAL) = sum(energies(1:ENE%MAX))

end subroutine energy
