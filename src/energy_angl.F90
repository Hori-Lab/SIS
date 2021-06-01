subroutine energy_angl(Eangl)

   use const
   use var_state, only : xyz
   use var_potential, only : nangl, angl_mp, angl_k, angl_t0 !angl_para

   implicit none
  
   real(PREC), intent(inout) :: Eangl

   integer :: ibd, imp1, imp2, imp3
   !real(PREC) :: k, t0
   real(PREC) :: t

   do ibd = 1, nangl
      imp1 = angl_mp(1, ibd)
      imp2 = angl_mp(2, ibd)
      imp3 = angl_mp(3, ibd)

      t = xyz_angle(xyz(:,imp1), xyz(:,imp2), xyz(:,imp3))
      
      !k = angl_para(1, ibd)
      !t0 = angl_para(2, ibd)

      !Eangl = Eangl + 0.5 * k * (t - t0) ** 2
      Eangl = Eangl + 0.5 * angl_k * (t - angl_t0) ** 2
   enddo
  
contains

   real(PREC) function xyz_angle(x1, x2, x3) result (angle)

      real(PREC), intent(in) :: x1(3), x2(3), x3(3)
      real(PREC) :: v12(3), v32(3)
      real(PREC) :: co

      v12(:) = x1(:) - x2(:)
      v32(:) = x3(:) - x2(:)

      co = dot_product(v32, v12) / sqrt(dot_product(v12,v12) * dot_product(v32,v32))

      if(co > 1.0e0_PREC) then
         co = 1.0e0_PREC
      else if(co < -1.0e0_PREC) then
         co = -1.0e0_PREC
      end if

      angle = acos(co)

   endfunction xyz_angle

end subroutine energy_angl
