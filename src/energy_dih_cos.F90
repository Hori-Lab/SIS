subroutine energy_dih_cos(Edih)
   
   use const
   use pbc, only : pbc_vec_d
   use var_state, only : xyz
   use var_potential, only : ndih, dih_mp, dih_k, dih_p0

   implicit none

   real(PREC), intent(inout) :: Edih
  
   integer :: idih, imp1, imp2, imp3, imp4
   real(PREC) :: dih
   real(PREC) :: c1(3), c2(3)
   real(PREC) :: v12(3), v23(3), v34(3)

   do idih = 1, ndih
      imp1 = dih_mp(1, idih)
      imp2 = dih_mp(2, idih)
      imp3 = dih_mp(3, idih)
      imp4 = dih_mp(4, idih)

      v12(:)= pbc_vec_d(xyz(:, imp2), xyz(:, imp1))
      v23(:)= pbc_vec_d(xyz(:, imp3), xyz(:, imp2))
      v34(:)= pbc_vec_d(xyz(:, imp4), xyz(:, imp3))

      c1(1) = v12(2)*v23(3) - v12(3)*v23(2)
      c1(2) = v12(3)*v23(1) - v12(1)*v23(3)
      c1(3) = v12(1)*v23(2) - v12(2)*v23(1)

      c2(1) = v23(2)*v34(3) - v23(3)*v34(2)
      c2(2) = v23(3)*v34(1) - v23(1)*v34(3)
      c2(3) = v23(1)*v34(2) - v23(2)*v34(1)

      dih = atan2(norm2(v23) * dot_product(v12, c2), dot_product(c1, c2))

      Edih = Edih + dih_k * (1.0_PREC + cos(dih + dih_p0))
    end do
end subroutine energy_dih_cos
