subroutine energy_dihedral(Edihedral)
   
   use const
   use pbc, only : pbc_vec_d
   use var_state, only : xyz
   use var_potential, only : ndihedral, dihedral_mp, angl_kphi, angl_phi0

   implicit none

   real(PREC), intent(inout) :: Edihedral
  
   integer :: idih, imp1, imp2, imp3, imp4
   real(PREC) :: dih
   real(PREC) :: c1(3), c2(3)
   real(PREC) :: v12(3), v23(3), v34(3)

   do idih = 1, ndihedral
      imp1 = dihedral_mp(1, idih)
      imp2 = dihedral_mp(2, idih)
      imp3 = dihedral_mp(3, idih)
      imp4 = dihedral_mp(4, idih)

      v12(:)= pbc_vec_d(xyz(:, imp1), xyz(:, imp2))
      v23(:)= pbc_vec_d(xyz(:, imp2), xyz(:, imp3))
      v34(:)= pbc_vec_d(xyz(:, imp3), xyz(:, imp4))

      c1(1) = v12(2)*v23(3) - v12(3)*v23(2)
      c1(2) = v12(3)*v23(1) - v12(1)*v23(3)
      c1(3) = v12(1)*v23(2) - v12(2)*v23(1)

      c2(1) = v23(2)*v34(3) - v23(3)*v34(2)
      c2(2) = v23(3)*v34(1) - v23(1)*v34(3)
      c2(3) = v23(1)*v34(2) - v23(2)*v34(1)

      dih = atan2(norm2(v23) * dot_product(v12, c2), dot_product(c1, c2))

      Edihedral = Edihedral + angl_kphi * (1.0_PREC + cos(dih + angl_phi0))
    end do
end subroutine energy_dihedral
