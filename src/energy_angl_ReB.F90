subroutine energy_angl_ReB(Eangl)

   use const
   use pbc, only : pbc_vec_d
   use var_state, only : xyz
   use var_potential, only : nangl, angl_mp, angl_k, angl_t0 !angl_para

   implicit none

   real(PREC), intent(inout) :: Eangl

   integer :: ibd, imp1, imp2, imp3
   real(PREC) :: co, co0, v12(3), v32(3)

   co0 = cos(angl_t0)

   do ibd = 1, nangl
      imp1 = angl_mp(1, ibd)
      imp2 = angl_mp(2, ibd)
      imp3 = angl_mp(3, ibd)

      v12(:) = pbc_vec_d(xyz(:, imp1), xyz(:, imp2))
      v32(:) = pbc_vec_d(xyz(:, imp3), xyz(:, imp2))

      co = dot_product(v32,v12) / sqrt(dot_product(v12,v12) * dot_product(v32,v32))

      Eangl = Eangl + 0.5_PREC * angl_k * (co - co0)**2 / (1.0_PREC - co**2)
   enddo

end subroutine energy_angl_ReB
