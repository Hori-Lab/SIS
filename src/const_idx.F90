module const_idx

  use,intrinsic :: ISO_FORTRAN_ENV, only: REAL64, INT64
  use :: const

  implicit none

  type energy_types
     integer :: TOTAL  ! 0
     integer :: BOND   ! 1
     integer :: ANGL   ! 2
     integer :: BP     ! 3
     integer :: ELE    ! 4
     integer :: MAX
  endtype energy_types
  type(energy_types), parameter :: ENE = energy_types(0,1,2,3,4,4)

  type seq_types
     integer :: A   ! 1
     integer :: U   ! 2
     integer :: G   ! 3
     integer :: C   ! 4
     integer :: MAX
  endtype seq_types
  type(seq_types), parameter :: SEQT = seq_types(1,2,3,4,4)

end module const_idx
