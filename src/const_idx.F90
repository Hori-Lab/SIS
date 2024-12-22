module const_idx

   use,intrinsic :: ISO_FORTRAN_ENV, only: REAL64, INT64

   implicit none

   type energy_types
      integer :: TOTAL  ! 0
      integer :: BOND   ! 1
      integer :: ANGL   ! 2
      integer :: DIHE   ! 3
      integer :: BP     ! 4
      integer :: EXV    ! 5
      integer :: ELE    ! 6
      integer :: STAGE  ! 7
      integer :: TWZ    ! 8
      integer :: RG     ! 9
      integer :: REST   !10
      integer :: MAX
   endtype energy_types
   type(energy_types), parameter :: ENE = energy_types(0,1,2,3,4,5,6,7,8,9,10,10)

   type mol_types
      integer :: UNDEF    ! 0
      integer :: RNA      ! 1
      integer :: CIRCRNA  ! 2
      integer :: MAX      ! 2
   endtype mol_types
   type(mol_types), parameter :: MOLT = mol_types(0,1,2,2)

   type seq_types
      integer :: UNDEF  ! -1
      integer :: A   ! 0
      integer :: U   ! 1
      integer :: G   ! 2
      integer :: C   ! 3
      integer :: D   ! 4
      integer :: MAX
   endtype seq_types
   type(seq_types), parameter :: SEQT = seq_types(-1,0,1,2,3,4,4)

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
      integer :: BD_EM       ! 2
      integer :: MAX
   endtype integrator_types
   type(integrator_types), parameter :: INTGRT = integrator_types(0,1,2,2)

   type bp_types
      integer :: UNDEF   ! 0
      integer :: GC      ! 1
      integer :: AU      ! 2
      integer :: GU      ! 3
      integer :: CG      ! 4
      integer :: UA      ! 5
      integer :: UG      ! 6
      integer :: MAX
   endtype bp_types
   type(bp_types), parameter :: BPT = bp_types(0,1,2,3,4,5,6,6)

   character(len=*), dimension(6), parameter ::  BPTYPE_CHAR = (/'GC', 'AU', 'GU', 'CG', 'UA', 'UG'/)

   type rst_block
      integer :: STEP    ! 1
      integer :: ANNEAL  ! 2
      integer :: XYZ     ! 3
      integer :: VELO    ! 4
      integer :: ACCEL   ! 5
      integer :: REPLICA ! 6
      integer :: PRNG    ! 7
      integer :: PRNGREP ! 8
      integer :: TWZ     ! 9
      integer :: PBC     ! 10
   endtype rst_block
   type(rst_block), parameter :: RSTBLK = rst_block(1,2,3,4,5,6,7,8,9,10)

   type replica_type
      integer :: TEMP
      integer :: TWZDCF
      integer :: ION
      integer :: RG
      integer :: MAX
   endtype replica_type
   type(replica_type), parameter :: REPT = replica_type(1,2,3,4,4)

   type nn_types
      integer :: GC_CG
      integer :: CC_GG
      integer :: GA_CU
      integer :: CG_GC
      integer :: AC_UG
      integer :: CA_GU
      integer :: AG_UC
      integer :: UA_AU
      integer :: AU_UA
      integer :: AA_UU
      integer :: GC_UG
      integer :: CU_GG
      integer :: GG_CU
      integer :: CG_GU
      integer :: AU_UG
      integer :: GA_UU
      integer :: UG_GU
      integer :: UA_GU
      integer :: GG_UU
      integer :: GU_UG
      integer :: AG_UU
      integer :: MAX
   endtype nn_types
   type(nn_types), parameter :: NNT = nn_types(1,2,3,4,5,6,7,8,9,10,11,12,13,14,15, &
                                               16,17,18,19,20,21,21)

   type nnend_types
      integer :: AUonAU
      integer :: AUonCG
      integer :: AUonGU
      integer :: GUonCG
      integer :: GUonAU
      integer :: GUonGU
   endtype nnend_types
   type(nnend_types), parameter :: NNENDT = nnend_types(1,2,3,4,5,6)

   type potential_types
      integer :: HARMONIC
      integer :: FLATBOTTOMED
   endtype potential_types
   type(potential_types), parameter :: POTT = potential_types(1,2)

   type tomlfstatus_types
      integer :: SUCCESS          !  0
      integer :: FATAL            ! -1
      integer :: DUPLICATE_KEY    ! -2
      integer :: TYPE_MISMATCH    ! -3
      integer :: CONVERSION_ERROR ! -4
      integer :: MISSING_KEY      ! -5
   endtype tomlfstatus_types
   type(tomlfstatus_types), parameter :: TOMLFSTAT = tomlfstatus_types(0,-1,-2,-3,-4,-5)

contains
   function molt2char(i) result(c)
      integer, intent(in) :: i
      character(:), allocatable :: c

      select case (i)
      case (MOLT%UNDEF)
         c = '?'
      case (MOLT%RNA)
         c = 'linearRNA'
      case (MOLT%CIRCRNA)
         c = 'circRNA'
      case default
         c = '?'
      endselect

   endfunction molt2char

   function char2molt(c) result (i)
      character(:), allocatable, intent(in) :: c
      integer :: i

      if (c == 'linearRNA') then
         i = MOLT%RNA

      else if (c == 'circRNA') then
         i = MOLT%CIRCRNA

      else
         i = MOLT%UNDEF
      endif

   endfunction char2molt

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
      case (SEQT%D)
         c = 'D'
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

      else if (c == 'D' .or. c == 'd') then
         i = SEQT%D

      else
         i = SEQT%UNDEF
      endif

   endfunction char2seqt

   function char2nnt(c) result (i)
      character(len=5), intent(in) :: c
      integer :: i

      if (c == 'GC_CG') then
         i = NNT%GC_CG
      else if (c == 'CC_GG' .or. c == 'GG_CC') then
         i = NNT%CC_GG
      else if (c == 'GA_CU' .or. c == 'UC_AG') then
         i = NNT%GA_CU
      else if (c == 'CG_GC') then
         i = NNT%CG_GC
      else if (c == 'AC_UG' .or. c == 'GU_CA') then
         i = NNT%AC_UG
      else if (c == 'CA_GU' .or. c == 'UG_AC') then
         i = NNT%CA_GU
      else if (c == 'AG_UC' .or. c == 'CU_GA') then
         i = NNT%AG_UC
      else if (c == 'UA_AU') then
         i = NNT%UA_AU
      else if (c == 'AU_UA') then
         i = NNT%AU_UA
      else if (c == 'AA_UU' .or. c == 'UU_AA') then
         i = NNT%AA_UU
      else if (c == 'GC_UG' .or. c == 'GU_CG') then
         i = NNT%GC_UG
      else if (c == 'CU_GG' .or. c == 'GG_UC') then
         i = NNT%CU_GG
      else if (c == 'GG_CU' .or. c == 'UC_GG') then
         i = NNT%GG_CU
      else if (c == 'CG_GU' .or. c == 'UG_GC') then
         i = NNT%CG_GU
      else if (c == 'AU_UG' .or. c == 'GU_UA') then
         i = NNT%AU_UG
      else if (c == 'GA_UU' .or. c == 'UU_AG') then
         i = NNT%GA_UU
      else if (c == 'UG_GU') then
         i = NNT%UG_GU
      else if (c == 'UA_GU' .or. c == 'UG_AU') then
         i = NNT%UA_GU
      else if (c == 'GG_UU' .or. c == 'UU_GG') then
         i = NNT%GG_UU
      else if (c == 'GU_UG') then
         i = NNT%GU_UG
      else if (c == 'AG_UU' .or. c == 'UU_GA') then
         i = NNT%AG_UU
      endif

   endfunction char2nnt

   function nnt2char(i) result (c)
      integer, intent(in) :: i
      character(len=5) :: c

      if (i == NNT%GC_CG) then
         c = 'GC_CG'
      else if (i == NNT%CC_GG) then
         c = 'CC_GG'
      else if (i == NNT%GA_CU) then
         c = 'GA_CU'
      else if (i == NNT%CG_GC) then
         c = 'CG_GC'
      else if (i == NNT%AC_UG) then
         c = 'AC_UG'
      else if (i == NNT%CA_GU) then
         c = 'CA_GU'
      else if (i == NNT%AG_UC) then
         c = 'AG_UC'
      else if (i == NNT%UA_AU) then
         c = 'UA_AU'
      else if (i == NNT%AU_UA) then
         c = 'AU_UA'
      else if (i == NNT%AA_UU) then
         c = 'AA_UU'
      else if (i == NNT%GC_UG) then
         c = 'GC_UG'
      else if (i == NNT%CU_GG) then
         c = 'CU_GG'
      else if (i == NNT%GG_CU) then
         c = 'GG_CU'
      else if (i == NNT%CG_GU) then
         c = 'CG_GU'
      else if (i == NNT%AU_UG) then
         c = 'AU_UG'
      else if (i == NNT%GA_UU) then
         c = 'GA_UU'
      else if (i == NNT%UG_GU) then
         c = 'UG_GU'
      else if (i == NNT%UA_GU) then
         c = 'UA_GU'
      else if (i == NNT%GG_UU) then
         c = 'GG_UU'
      else if (i == NNT%GU_UG) then
         c = 'GU_UG'
      else if (i == NNT%AG_UU) then
         c = 'AG_UU'
      endif
   endfunction nnt2char

   function seqt2nnt(i,j,k,l) result (nnt)
      integer, intent(in) :: i, j, k, l
      integer :: nnt
      character(1) :: ci, cj, ck, cl

      ci = seqt2char(i)
      cj = seqt2char(j)
      ck = seqt2char(k)
      cl = seqt2char(l)

      nnt = char2nnt(ci//cj//'_'//ck//cl)
   endfunction seqt2nnt

   logical function is_complement(s1, s2)
      integer, intent(in) :: s1, s2

      is_complement = .False.

      if (s1 == SEQT%A) then
         if (s2 == SEQT%U) then
            is_complement = .True.
         endif

      else if (s1 == SEQT%U) then
         if (s2 == SEQT%A .or. s2 == SEQT%G) then
            is_complement = .True.
         endif

      else if (s1 == SEQT%G) then
         if (s2 == SEQT%C .or. s2 == SEQT%U) then
            is_complement = .True.
         endif

      else if (s1 == SEQT%C) then
         if (s2 == SEQT%G) then
            is_complement = .True.
         endif

      endif
   end function is_complement

end module const_idx
