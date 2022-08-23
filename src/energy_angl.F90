subroutine energy_angl(irep, Eangl)

   use const
   use pbc, only : pbc_vec_d
   use var_state, only : xyz
   use var_potential, only : nangl, angl_mp, angl_k, angl_t0

   implicit none
  
   integer, intent(in) :: irep
   real(PREC), intent(inout) :: Eangl

   integer :: ibd, imp1, imp2, imp3
   !real(PREC) :: k, t0
   real(PREC) :: t
   real(PREC) :: v12(3), v32(3)

   do ibd = 1, nangl
      imp1 = angl_mp(1, ibd)
      imp2 = angl_mp(2, ibd)
      imp3 = angl_mp(3, ibd)

      v12(:) = pbc_vec_d(xyz(:, imp1, irep), xyz(:, imp2, irep))
      v32(:) = pbc_vec_d(xyz(:, imp3, irep), xyz(:, imp2, irep))

      t = dot_product(v32, v12) / sqrt(dot_product(v12,v12) * dot_product(v32,v32))

      if(t > 1.0e0_PREC) then
         t = 1.0e0_PREC
      else if(t < -1.0e0_PREC) then
         t = -1.0e0_PREC
      end if

      t = acos(t)
      
      !k = angl_para(1, ibd)
      !t0 = angl_para(2, ibd)

      !Eangl = Eangl + 0.5 * k * (t - t0) ** 2
      Eangl = Eangl + 0.5 * angl_k * (t - angl_t0) ** 2
   enddo
  
end subroutine energy_angl
