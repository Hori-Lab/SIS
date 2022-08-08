subroutine force_dihedral(forces)
   
   use const
   use pbc, only : pbc_vec_d
   use var_state, only : xyz
   use var_top, only : nmp
   use var_potential, only : ndihedral, dihedral_mp, angl_kphi, angl_phi0

   implicit none

   real(PREC), intent(inout) :: forces(3, nmp)
  
   integer :: idih, imp1, imp2, imp3, imp4
   real(PREC) :: dih, akj, akj2, pre, dvijvkj_akj2, dvklvkj_akj2
   real(PREC) :: m(3), n(3)
   real(PREC) :: vij(3), vkj(3), vkl(3), fi(3), fj(3), fk(3), fl(3)

   !$omp   do private(imp1, imp2, imp3, imp4, vij, vkj, vkl, m, n, &
   !$omp&             akj2, akj, dih, pre, fi, fj, fk, fl, dvijvkj_akj2, dvklvkj_akj2)
   do idih = 1, ndihedral
      imp1 = dihedral_mp(1, idih)
      imp2 = dihedral_mp(2, idih)
      imp3 = dihedral_mp(3, idih)
      imp4 = dihedral_mp(4, idih)

      vij(:) = pbc_vec_d(xyz(:, imp1), xyz(:, imp2))
      vkj(:) = pbc_vec_d(xyz(:, imp3), xyz(:, imp2))
      vkl(:) = pbc_vec_d(xyz(:, imp3), xyz(:, imp4))

      ! m = rij x rkj
      m(1) = vij(2)*vkj(3) - vij(3)*vkj(2)
      m(2) = vij(3)*vkj(1) - vij(1)*vkj(3)
      m(3) = vij(1)*vkj(2) - vij(2)*vkj(1)

      ! n = rkj x rkl
      n(1) = vkj(2)*vkl(3) - vkj(3)*vkl(2)
      n(2) = vkj(3)*vkl(1) - vkj(1)*vkl(3)
      n(3) = vkj(1)*vkl(2) - vkj(2)*vkl(1)

      akj2 = dot_product(vkj, vkj)
      akj = sqrt(akj2)

      dih = atan2(-akj * dot_product(vij, n), dot_product(m, n))

      pre = - angl_kphi * sin(dih + angl_phi0) * akj
      fi(:) =  pre / dot_product(m,m) * m(:)
      fl(:) = -pre / dot_product(n,n) * n(:)
      
      dvijvkj_akj2 = dot_product(vij, vkj) / akj2
      dvklvkj_akj2 = dot_product(vkl, vkj) / akj2
      fj(:) = (-1.0_PREC + dvijvkj_akj2) * fi(:) + (          - dvklvkj_akj2) * fl(:)
      fk(:) = (          - dvijvkj_akj2) * fi(:) + (-1.0_PREC + dvklvkj_akj2) * fl(:)

      forces(:, imp1) = forces(:, imp1) + fi(:)
      forces(:, imp2) = forces(:, imp2) + fj(:)
      forces(:, imp3) = forces(:, imp3) + fk(:) 
      forces(:, imp4) = forces(:, imp4) + fl(:)
   end do 
   !$omp end do nowait

end subroutine force_dihedral
