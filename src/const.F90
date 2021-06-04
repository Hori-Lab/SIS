module const

  use,intrinsic :: ISO_FORTRAN_ENV, only: REAL64, INT64

  implicit none

#ifdef PREC_INTEL16
  integer, parameter :: PREC = 16 !< precision of floating-point
#else
  integer, parameter :: PREC = REAL64  !< precision of floating-point (default, 4 or 8)
#endif

  integer, save      :: S_REAL    !< byte size of default real   (single precision)
  integer, save      :: M_INT     !< byte size of default integer   (medium-size integer)
  integer, save      :: LOGIC     !< byte size of default logical
  integer, parameter :: L_INT = INT64 !< byte size of larger integer ("long long int" in C)

  integer, parameter :: CHAR_FILE_PATH = 1000
  integer, parameter :: CHAR_FILE_LINE = 500

end module const
