!*************************************************************************
!>
!  DIVA constants.

    module diva_constants

    use iso_fortran_env, only: wp => real64

    implicit none

    public

    !++S Default KDIM = 16 [for single precision]
    !++  Default KDIM = 20
    !++  Default MAXORD = 2, MAXSTF = 1
    !++  Default INTEGO, VAREQ, OUTPUT, DUMP, GSTOP, EXTRAP
    !++  Default STIFF=.F., ARGM=.F., ERRSTO=.F.

    !++ substitute for kdim, maxord, maxstf below
    integer, parameter :: kdim = 20
    integer, parameter :: maxord = 2
    integer, parameter :: maxstf = 1

    real(wp),dimension(5),parameter :: d1mach = &
        [tiny(1.0_wp), &
         huge(1.0_wp), &
         real(radix(1.0_wp),wp)**(-digits(1.0_wp)), &
         epsilon(1.0_wp), &
         log10(real(radix(1.0_wp),wp))] !! machine constants

    end module diva_constants
!*************************************************************************