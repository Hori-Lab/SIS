subroutine energy_dih_exp(irep, Edih)

   use const
   use pbc, only : pbc_vec_d
   use var_state, only : xyz
   use var_potential, only : ndih, dih_mp, dih_k, dih_w, dih_p0

   implicit none

   integer, intent(in) :: irep
   real(PREC), intent(inout) :: Edih

   integer :: idih, imp1, imp2, imp3, imp4
   real(PREC) :: dih, coef
   real(PREC) :: c1(3), c2(3)
   real(PREC) :: v12(3), v23(3), v34(3)

   coef = - 0.5_PREC * dih_w

   do idih = 1, ndih
      imp1 = dih_mp(1, idih)
      imp2 = dih_mp(2, idih)
      imp3 = dih_mp(3, idih)
      imp4 = dih_mp(4, idih)

      v12(:)= pbc_vec_d(xyz(:, imp2, irep), xyz(:, imp1, irep))
      v23(:)= pbc_vec_d(xyz(:, imp3, irep), xyz(:, imp2, irep))
      v34(:)= pbc_vec_d(xyz(:, imp4, irep), xyz(:, imp3, irep))

      c1(1) = v12(2)*v23(3) - v12(3)*v23(2)
      c1(2) = v12(3)*v23(1) - v12(1)*v23(3)
      c1(3) = v12(1)*v23(2) - v12(2)*v23(1)

      c2(1) = v23(2)*v34(3) - v23(3)*v34(2)
      c2(2) = v23(3)*v34(1) - v23(1)*v34(3)
      c2(3) = v23(1)*v34(2) - v23(2)*v34(1)

      dih = atan2(norm2(v23) * dot_product(v12, c2), dot_product(c1, c2))

      Edih = Edih - dih_k * exp(coef * (dih - dih_p0)**2)
    end do

end subroutine energy_dih_exp
