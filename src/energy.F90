subroutine energy()

   use const
   use const_idx, only : ENE
   use var_state, only : energies
   use var_potential, only : flg_angl_ReB, flg_ele, max_bp_per_nt, flg_dih_cos, flg_dih_exp, bp_model

   implicit none

   energies(:) = 0.0e0_PREC

   call energy_bond(energies(ENE%BOND))

   if (flg_angl_ReB) then
      call energy_angl_ReB(energies(ENE%ANGL))
   else
      call energy_angl(energies(ENE%ANGL))
   endif
   
   if (flg_dih_cos) call energy_dih_cos(energies(ENE%DIHE))

   if (flg_dih_exp) call energy_dih_exp(energies(ENE%DIHE))

   if (max_bp_per_nt < 1) then
      call energy_bp(energies(ENE%BP))
   else
      if (bp_model == 4) then
         call energy_bp_limit_triplet(energies(ENE%BP))
      else
         call energy_bp_limit(energies(ENE%BP))
      endif
   endif

   call energy_wca(energies(ENE%EXV))

   if (flg_ele) call energy_ele_DH(energies(ENE%ELE))

   energies(ENE%TOTAL) = sum(energies(1:ENE%MAX))

end subroutine energy
