module const_idx

  use,intrinsic :: ISO_FORTRAN_ENV, only: REAL64, INT64
  use :: const

  implicit none

  type energy_types
     integer :: TOTAL  ! 0
     integer :: BOND   ! 1
     integer :: ANGL   ! 2
     integer :: HB     ! 3
     integer :: ELE    ! 4
     integer :: MAX
  endtype energy_types
  type(energy_types), parameter :: ENE = energy_types(0,1,2,3,4,4)

end module const_idx
