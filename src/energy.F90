subroutine energy()

   use const
   use const_idx, only : ENE
   use var_state, only : energies

   implicit none

   energies(:) = 0.0

   call energy_bond(energies(ENE%BOND))
   call energy_angl(energies(ENE%ANGL))

   energies(ENE%TOTAL) = sum(energies(1:ENE%MAX))

end subroutine energy
