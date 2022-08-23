subroutine force_angl_ReB(irep, forces)

   use const
   use pbc, only : pbc_vec_d
   use var_state, only : xyz
   use var_top, only : nmp
   use var_potential, only : nangl, angl_mp, angl_k, angl_t0

   implicit none

   integer, intent(in) :: irep
   real(PREC), intent(inout) :: forces(3, nmp)

   integer :: ibd, imp1, imp2, imp3
   real(PREC) :: f, cosine, co0
   real(PREC) :: d12, d32, a1232
   real(PREC) :: v12(3), v32(3), f1(3), f3(3)

   co0 = cos(angl_t0)

   !$omp do private(imp1,imp2,imp3,v12,v32,d12,d32,a1232,cosine,f,f1,f3)
   do ibd = 1, nangl
      imp1 = angl_mp(1, ibd)
      imp2 = angl_mp(2, ibd)
      imp3 = angl_mp(3, ibd)

      v12(:) = pbc_vec_d(xyz(:, imp1, irep), xyz(:, imp2, irep))
      v32(:) = pbc_vec_d(xyz(:, imp3, irep), xyz(:, imp2, irep))

      d12 = dot_product(v12, v12)
      d32 = dot_product(v32, v32)
      a1232 = sqrt(d12*d32)

      cosine = dot_product(v12, v32) / a1232

      f = - angl_k * (cosine - co0) * (1.0_PREC - cosine * co0) / (1.0_PREC - cosine**2)**2
      f1(:) = f * (v32(:) / a1232 - cosine * v12(:) / d12)
      f3(:) = f * (v12(:) / a1232 - cosine * v32(:) / d32)

      forces(:, imp1) = forces(:, imp1) + f1(:)
      forces(:, imp2) = forces(:, imp2) - (f1(:) + f3(:))
      forces(:, imp3) = forces(:, imp3) + f3(:)
   enddo
   !$omp end do nowait

end subroutine force_angl_ReB
