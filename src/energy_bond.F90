subroutine energy_bond(irep, Ebd)
      
   use const
   use pbc, only : pbc_vec_d
   use var_state, only : xyz
   use var_potential, only : nbond, bond_mp, bond_k, bond_r0
  
   implicit none
  
   integer, intent(in) :: irep
   real(PREC), intent(inout) :: Ebd

   integer :: ibd, imp1, imp2
   !real(PREC) :: k, r0
   real(PREC) :: d

   do ibd = 1, nbond
      imp1 = bond_mp(1, ibd)
      imp2 = bond_mp(2, ibd)
      !k = bond_para(1, ibd)
      !r0 = bond_para(2, ibd)

      d = norm2(pbc_vec_d(xyz(:, imp1, irep), xyz(:, imp2, irep)))

      !Ebd = Ebd + 0.5 * k * (d - r0) ** 2
      Ebd = Ebd + 0.5 * bond_k * (d - bond_r0) ** 2

   enddo
  
end subroutine energy_bond
