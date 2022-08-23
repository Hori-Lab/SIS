module const

   use,intrinsic :: ISO_FORTRAN_ENV, only: REAL64, INT64

   implicit none

   integer, save      :: LOGIC    !< byte size of default logical
   integer, save      :: M_INT    !< byte size of default integer   (medium-size integer)
   integer, parameter :: L_INT = INT64 !< byte size of larger integer ("long long int" in C)

#ifdef PREC_INTEL16
   integer, parameter :: PREC = 16 !< precision of floating-point
#else
   integer, parameter :: PREC = REAL64  !< precision of floating-point (default, 4 or 8)
#endif

   integer, parameter :: CHAR_FILE_PATH = 1000
   integer, parameter :: CHAR_FILE_LINE = 500
   integer, parameter :: FILENAME_DIGIT_REPLICA = 4

   integer, parameter :: MAX_REPLICA = 8
contains

   subroutine init_const()

      logical, parameter :: ldummy = .True.
      integer, parameter :: idummy = 1

      LOGIC = sizeof(ldummy)
      M_INT = sizeof(idummy)

   end subroutine init_const

end module const
