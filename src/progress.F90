module progress

   use,intrinsic :: ISO_FORTRAN_ENV, only: INT64, REAL64

   implicit none

   private
   real(REAL64), parameter :: sec_in_day = 60*60*24.0
   real(REAL64), parameter :: fn_million = 1000000.0

   integer(INT64), save :: step_pre, step_ini, t_rate
   integer(INT64), save :: clock_pre, clock_ini
   real(REAL64), save :: rt_rate

   integer(INT64), save :: wall_clock_0

   logical, save :: flg_init = .False.

   public :: progress_init
   public :: progress_update
   public :: wall_time_sec

contains

   subroutine progress_init(istep)

      use,intrinsic :: ISO_FORTRAN_ENV, only: OUTPUT_UNIT

      implicit none

      integer(INT64), intent(in) :: istep

      integer(INT64) :: clock

      write(OUTPUT_UNIT,'(a)') '# PROGRESS'
      write(OUTPUT_UNIT,'(a)') '#1: step'
      write(OUTPUT_UNIT,'(a)') '#2: million steps / day (latest cycle)'
      write(OUTPUT_UNIT,'(a)') '#3: million steps / day (averaged)'
      write(OUTPUT_UNIT,'(a)') '#4: remaining time / hours'

      call system_clock(clock,t_rate)

      clock_pre = clock
      clock_ini = clock
      step_ini = istep
      step_pre = istep
      rt_rate = 1.0_REAL64 / real(t_rate, kind=REAL64)

      wall_clock_0 = -clock

      flg_init = .True.
   endsubroutine progress_init

   subroutine progress_update(istep, nstep)

      use,intrinsic :: ISO_FORTRAN_ENV, only: OUTPUT_UNIT, INT64, REAL64

      implicit none
   
      integer(INT64), intent(in) :: istep
      integer(INT64), intent(in) :: nstep

      integer(INT64) :: clock
      real(REAL64) :: mstep_per_day_pre, mstep_per_day_ini
      real(REAL64) :: days_from_pre, days_from_ini, remain_hours

      if (.not. flg_init) then
         call progress_init(istep)
         return
      endif

      call system_clock(clock)

      days_from_ini = rt_rate * (clock - clock_ini) / sec_in_day
      days_from_pre = rt_rate * (clock - clock_pre) / sec_in_day

      mstep_per_day_ini = real(istep - step_ini, kind=REAL64) / fn_million / days_from_ini
      mstep_per_day_pre = real(istep - step_pre, kind=REAL64) / fn_million / days_from_pre

      remain_hours = real(nstep - istep, kind=REAL64) / fn_million / mstep_per_day_ini * 24.0

      write(OUTPUT_UNIT,'(i12,1x,g13.6,1x,g13.6,1x,f6.1)') istep, mstep_per_day_pre, mstep_per_day_ini, remain_hours
      flush(OUTPUT_UNIT)

      clock_pre = clock
      step_pre = istep

   endsubroutine progress_update

   integer(INT64) function wall_time_sec()
      integer(INT64) :: clock

      call system_clock(clock)

      wall_time_sec = int((wall_clock_0 + clock) * rt_rate, kind=INT64)

   end function wall_time_sec

endmodule progress
