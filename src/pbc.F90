module pbc

   use const, only : PREC

   logical, save :: flg_pbc
   real(PREC), save :: pbc_box(3)
   real(PREC), save :: pbc_box_half(3)

contains

   function pbc_vec(v) result (new_vec)
      
      use const, only : PREC

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

   subroutine pbc_wrap()
      
      use const
      use var_top, only : nmp
      use var_state, only : xyz

      implicit none
  
      integer :: imp, i
      real(PREC) :: box_max(3), box_min(3)
      real(PREC) :: x
     
      box_max(1:3) = 0.5*pbc_box(1:3)
      box_min(1:3) = -0.5*pbc_box(1:3)
  
      do imp = 1, nmp
         do i = 1, 3
            x = xyz(i, imp)
            if(x > box_max(i)) then
               xyz(i, imp) = x - pbc_box(i) * (int((x - box_max(i))/pbc_box(i)) + 1)
            else if(x < box_min(i)) then
               xyz(i, imp) = x + pbc_box(i) * (int((box_min(i) - x)/pbc_box(i)) + 1)
            endif
         enddo
      enddo

   endsubroutine pbc_wrap

endmodule pbc
