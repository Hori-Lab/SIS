subroutine energy_wca(Ewca)

   use const
   use var_state, only : xyz
   use var_potential, only : wca_sigma, wca_eps, wca_mp, nwca

   implicit none
  
   real(PREC), intent(inout) :: Ewca

   integer :: iwca, imp, jmp
   real(PREC) :: d
   real(PREC) :: e_wca

   e_wca = 0.0e0_PREC

   !$omp parallel do private(imp,jmp,d) reduction(+:e_wca)
   do iwca = 1, nwca

      imp = wca_mp(1, iwca)
      jmp = wca_mp(2, iwca)
      
      d = mp_distance(imp, jmp)

      if (d >= wca_sigma) cycle

      d = wca_sigma / d
      e_wca = e_wca + wca_eps * (d**12 - 2 * d**6 + 1.0)

   enddo
   !$omp end parallel do

   Ewca = e_wca

contains

   function pbc_vec(v) result (new_vec)
      
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

   function mp_distance(imp1, imp2) result(d)

      real(PREC) :: d
      integer, intent(in) :: imp1, imp2

      d = norm2( pbc_vec(xyz(:,imp1) - xyz(:, imp2)) )

   endfunction mp_distance

end subroutine energy_wca
