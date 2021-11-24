module const

   use,intrinsic :: ISO_FORTRAN_ENV, only: REAL64

   implicit none

#ifdef PREC_INTEL16
   integer, parameter :: PREC = 16 !< precision of floating-point
#else
   integer, parameter :: PREC = REAL64  !< precision of floating-point (default, 4 or 8)
#endif

   integer, parameter :: CHAR_FILE_PATH = 1000
   integer, parameter :: CHAR_FILE_LINE = 500

end module const
