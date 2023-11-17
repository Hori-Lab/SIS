subroutine init_freeze()
   
   use, intrinsic :: iso_fortran_env, only : output_unit
   use var_top, only : nmp, is_frozen, frz_ranges

   implicit none

   integer :: i, imp_begin, imp_end

   !! Allocation
   allocate(is_frozen(nmp))
   is_frozen(:) = .False.

   do i = 1, size(frz_ranges, 2)
      imp_begin = frz_ranges(1, i)
      imp_end = frz_ranges(2, i)
      is_frozen(imp_begin:imp_end) = .True.
   enddo

endsubroutine init_freeze
