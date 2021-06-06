!*************************************************************************
!>
!  DIVA constants.

    module diva_constants

    use iso_fortran_env, only: wp => real64

    implicit none

    public

    !++ substitute for kdim, maxord, maxstf below
    integer, parameter :: kdim = 20
    integer, parameter :: maxord = 2
    integer, parameter :: maxstf = 1

    real(wp),dimension(5),parameter :: d1mach = &
        [tiny(1.0_wp), &
         huge(1.0_wp), &
         real(radix(1.0_wp),wp)**(-digits(1.0_wp)), &
         epsilon(1.0_wp), &
         log10(real(radix(1.0_wp),wp))]

    end module diva_constants
!*************************************************************************