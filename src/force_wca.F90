subroutine force_wca()

   use const, only : PREC
   use var_state, only : xyz, forces
   use var_potential, only : wca_sigma, wca_eps, wca_mp, nwca

   implicit none
  
   integer :: iwca, imp1, imp2
   real(PREC) :: v12(3), d2, for(3)
   real(PREC) :: wca_sigma_2, coef

   wca_sigma_2 = wca_sigma**2
   coef = 12.0e0_PREC * wca_eps / wca_sigma_2

   !$omp parallel do private(imp1,imp2,v12,d2,for)
   do iwca = 1, nwca

      imp1 = wca_mp(1, iwca)
      imp2 = wca_mp(2, iwca)

      v12(:) = pbc_vec(xyz(:,imp1) - xyz(:,imp2))
      
      d2 = dot_product(v12, v12)

      if (d2 >= wca_sigma_2) cycle

      d2 = wca_sigma_2 / d2

      for(:) = coef * (d2**7 - d2**4) * v12(:)

      forces(:, imp1) = forces(:, imp1) + for(:)
      forces(:, imp2) = forces(:, imp2) - for(:)
   enddo
   !$omp end parallel do

contains

   function pbc_vec(v) result (new_vec)
      
      use const, only : PREC
      use var_top, only : flg_pbc, pbc_box, pbc_box_half

      real(PREC) :: new_vec(3)
      real(PREC), intent(in) :: v(3)

      integer :: i

      if (.not. flg_pbc) then
         new_vec(:) = v(:)
         return
      endif

      do i = 1, 3
         if(v(i) > pbc_box_half(i)) then
            new_vec(i) = v(i) - pbc_box(i)
   
         else if(v(i) < -pbc_box_half(i)) then
            new_vec(i) = v(i) + pbc_box(i)
   
         else
            new_vec(i) = v(i)
   
         end if
      end do

   end function pbc_vec

end subroutine force_wca
