subroutine energy()

   use const
   use const_idx, only : ENE
   use var_state, only : energies
   use var_potential, only : flg_ele

   implicit none

   energies(:) = 0.0e0_PREC

   call energy_bond(energies(ENE%BOND))
   call energy_angl(energies(ENE%ANGL))
   call energy_bp(energies(ENE%BP))
   call energy_wca(energies(ENE%EXV))

   if (flg_ele) call energy_ele_DH(energies(ENE%ELE))

   energies(ENE%TOTAL) = sum(energies(1:ENE%MAX))

end subroutine energy
