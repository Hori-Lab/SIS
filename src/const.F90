module const

   use,intrinsic :: ISO_FORTRAN_ENV, only: REAL64, INT64, INT32

   implicit none

#ifdef PREC_INTEL16
   integer, parameter :: PREC = 16 !< precision of floating-point
#else
   integer, parameter :: PREC = REAL64  !< precision of floating-point (default, 4 or 8)
#endif

   integer, save      :: LOGIC    !< byte size of default logical
   integer, save      :: M_INT    !< byte size of default integer   (medium-size integer)
   integer, parameter :: L_INT = INT64 !< byte size of larger integer ("long long int" in C)
   integer, parameter :: MTS_SIZE =INT32 * 19 + INT32 * 624
                         ! in mt_stream.F90
                         !integer(INT32) :: i = -1         ! state vector index
                         !integer(INT32) :: stream_id = -1 ! stream ID
                         !integer(INT32) :: istatus = -1   ! initialization status
                         !integer(INT32) :: nn = -1        ! MT parameter N
                         !integer(INT32) :: mm = -1        ! MT parameter M
                         !integer(INT32) :: rr = -1        ! MT parameter R
                         !integer(INT32) :: ww = -1        ! MT parameter W (width =32)
                         !integer(INT32) :: aaa = 0        ! Companion matrix parameter
                         !integer(INT32) :: wmask = 0      ! 32-bit mask
                         !integer(INT32) :: umask = 0      ! Twist mask x(1)
                         !integer(INT32) :: lmask = 0      ! Twist mask x(0)
                         !integer(INT32) :: shift0 = 0     ! Temparing parameters ...
                         !integer(INT32) :: shift1 = 0
                         !integer(INT32) :: maskB  = 0
                         !integer(INT32) :: maskC  = 0
                         !integer(INT32) :: shiftB = 0
                         !integer(INT32) :: shiftC = 0
                         !integer(INT32) :: mag(0:1) = 0   ! mag(0) = 0, mag(1) = aaa
                         !integer(INT32), pointer :: state(:) => NULL()  ! state vector (0:this%nn-1)
                         !    The length of the last array should be 624 because,
                         !        integer(INT32), parameter :: MT19937_N = 624
                         !        g_mt_master%nn    = MT19937_N

   integer, parameter :: CHAR_FILE_PATH = 1000
   integer, parameter :: CHAR_FILE_LINE = 500
   integer, parameter :: FILENAME_DIGIT_REPLICA = 4

   integer, parameter :: MAX_REP_DIM = 1
   integer, parameter :: MAX_REP_PER_DIM = 256
   integer, parameter :: MAX_REPLICA = MAX_REP_PER_DIM * MAX_REP_DIM
contains

   subroutine init_const()

      logical, parameter :: ldummy = .True.
      integer, parameter :: idummy = 1

      LOGIC = sizeof(ldummy)
      M_INT = sizeof(idummy)

   end subroutine init_const

end module const
