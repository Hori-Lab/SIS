module const_idx

   use,intrinsic :: ISO_FORTRAN_ENV, only: REAL64, INT64
   use :: const

   implicit none

   type energy_types
      integer :: TOTAL  ! 0
      integer :: BOND   ! 1
      integer :: ANGL   ! 2
      integer :: BP     ! 3
      integer :: EXV    ! 4
      integer :: MAX
   endtype energy_types
   type(energy_types), parameter :: ENE = energy_types(0,1,2,3,4,4)

   type seq_types
      integer :: UNDEF  ! 0
      integer :: A   ! 1
      integer :: U   ! 2
      integer :: G   ! 3
      integer :: C   ! 4
      integer :: MAX
   endtype seq_types
   type(seq_types), parameter :: SEQT = seq_types(0,1,2,3,4,4)

   type job_types
      integer :: DEBUG       ! 0
      integer :: CHECK_FORCE ! 1
      integer :: MD          ! 2
      integer :: DCD         ! 3
      integer :: MAX
   endtype job_types
   type(job_types), parameter :: JOBT = job_types(0,1,2,3,3)

   type integrator_types
      integer :: UNDEF       ! 0
      integer :: LD_GJF2GJ   ! 1
      integer :: MAX
   endtype integrator_types
   type(integrator_types), parameter :: INTGRT = integrator_types(0,1,1)

   type bp_types
      integer :: UNDEF    ! 0
      integer :: GC_WCF   ! 1
      integer :: AU_WCF   ! 2
      integer :: GU_WBL   ! 3
      integer :: MAX
   endtype bp_types
   type(bp_types), parameter :: BPT = bp_types(0,1,2,3,3)

contains
   function seqt2char(i) result(c)
      integer, intent(in) :: i
      character(1) :: c

      select case (i)
      case (SEQT%A)
         c = 'A'
      case (SEQT%U)
         c = 'U'
      case (SEQT%G)
         c = 'G'
      case (SEQT%C)
         c = 'C'
      case (SEQT%UNDEF)
         c = '?'
      case default
         c = '?'
      endselect

   endfunction seqt2char

   function char2seqt(c) result (i)
      character(1), intent(in) :: c
      integer :: i

      if (c == 'A' .or. c == 'a') then
         i = SEQT%A

      else if (c == 'U' .or. c == 'u') then
         i = SEQT%U

      else if (c == 'G' .or. c == 'g') then
         i = SEQT%G

      else if (c == 'C' .or. c == 'c') then
         i = SEQT%C

      else
         i = SEQT%UNDEF
      endif

   endfunction char2seqt

end module const_idx
