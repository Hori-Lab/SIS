subroutine energy_dihedral(Edihedral)
   
   use const
   use pbc, only : pbc_vec_d
   use var_state, only : xyz
   use var_top, only : nmp
   use var_potential, only : ndihedral, dihedral_mp, angl_kphi, angl_phi0

   implicit none

   real(PREC), intent(inout) :: Edihedral
  
   integer :: ibd, imp1, imp2, imp3, imp4
   real(PREC) :: f, delta, cosine, anglsign, d1
   real(PREC) :: c1(3), c2(3), c3(3)
   real(PREC) :: v12(3), v23(3), v43(3), fi(3), fj(3), fk(3), fl(3)

   do ibd = 1, ndihedral
      imp1 = dihedral_mp(1, ibd)
      imp2 = dihedral_mp(2, ibd)
      imp3 = dihedral_mp(3, ibd)
      imp4 = dihedral_mp(4, ibd)

      v12(:)= pbc_vec_d(xyz(:, imp1), xyz(:, imp2))
      v23(:)= pbc_vec_d(xyz(:, imp2), xyz(:, imp3))
      v43(:)= pbc_vec_d(xyz(:, imp4), xyz(:, imp3))

      c1(1) = v12(2)*v23(3) - v12(3)*v23(2)
      c1(2) = v12(3)*v23(1) - v12(1)*v23(3)
      c1(3) = v12(1)*v23(2) - v12(2)*v23(1)

      c2(1) = v43(2)*v23(3) - v43(3)*v23(2)
      c2(2) = v43(3)*v23(1) - v43(1)*v23(3)
      c2(3) = v43(1)*v23(2) - v43(2)*v23(1)

      d1 = dot_product(c1, c2)
      cosine = - d1 /((NORM2(c1)*NORM2(c2)))

      if (abs(cosine) > 1.0e0_PREC) then
            cosine = sign(1.0e0_PREC, cosine)
      endif

      c3(1) = c2(2)*c1(3) - c2(3)*c1(2)
      c3(2) = c2(3)*c1(1) - c2(1)*c1(3)
      c3(3) = c2(1)*c1(2) - c2(2)*c1(1)

      anglsign = dot_product(c3, v23)
      if (anglsign < 0) then
         cosine = - cosine
      end if
      if (anglsign > 0) then
         cosine = abs(cosine)
      end if

      delta = (acos(cosine) - angl_phi0)
      Edihedral = Edihedral + angl_kphi *(1+cos(delta))
    end do
end subroutine energy_dihedral
