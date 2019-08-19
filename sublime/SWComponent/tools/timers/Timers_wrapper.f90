!
! Fortran wrapper of Timers
!

subroutine starttime(pilename)
    implicit none
    character(20) :: pilename
    integer(kind=4) :: length

    length = len_trim(pilename)

    call starttimer_f(pilename,length)

end subroutine starttime

subroutine endtime(pilename)
    implicit none
    character(20) :: pilename
    integer(kind=4) :: length

    length = len_trim(pilename)

    call endtimer_f(pilename,length)

end subroutine endtime

subroutine printtime
    implicit none

    call printtimer_f

end subroutine printtime

subroutine getmaxtime
    implicit none

    call maxtimesum_f

end subroutine getmaxtime
