module pbc

   use const, only : PREC

   logical, save :: flg_pbc
   real(PREC), save :: pbc_box(3)
   real(PREC), save :: pbc_box_half(3)

contains

   subroutine set_pbc_size(boxsize)

      real(PREC), intent(in) :: boxsize(3)

      pbc_box(:) = boxsize(:)
      pbc_box_half(:) = 0.5_PREC * boxsize(:)

   endsubroutine set_pbc_size

   pure function pbc_vec_d(v1, v2) result (new_vec)

      use const, only : PREC

      real(PREC) :: new_vec(3)
      real(PREC), intent(in) :: v1(3), v2(3)

      integer :: i

      new_vec(:) = v1(:) - v2(:)

      if (.not. flg_pbc) then
         return
      endif

      do i = 1, 3
         if(new_vec(i) > pbc_box_half(i)) then
            new_vec(i) = new_vec(i) - pbc_box(i)

         else if(new_vec(i) < -pbc_box_half(i)) then
            new_vec(i) = new_vec(i) + pbc_box(i)

         end if
      end do

   end function pbc_vec_d

   pure function pbc_vec(v) result (new_vec)
      
      use const, only : PREC

      real(PREC) :: new_vec(3)
      real(PREC), intent(in) :: v(3)

      integer :: i

      if (.not. flg_pbc) then
         new_vec(1:3) = v(1:3)
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
     
      !box_max(1:3) =  0.5*pbc_box(1:3)
      !box_min(1:3) = -0.5*pbc_box(1:3)
      box_max(1:3) = pbc_box(1:3)
      box_min(1:3) = 0.0_PREC
  
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
