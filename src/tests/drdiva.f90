!*************************************************************************
!>
!  Sample driver for DIVA --  Set up to solve two second order equations.
!
!### History
!  * 2010-06-09 DRDIVA  Krogh Used parameters for all dimenssions.
!  * 2001-05-25 DRDIVA  Krogh Minor change for making .f90 version.
!  * 1996-06-14 DRDIVA  Krogh  Small change in output format
!  * 1994-11-02 DRDIVA  Krogh  Changes to use M77CON
!  * 1994-07-18 DRDIVA  Krogh   Last change.

    program DRDIVA

    use diva_module

    implicit none

    integer, parameter :: ineq  = 2
    integer, parameter :: ifdim = 16*ineq+1
    integer, parameter :: ikdim = 6
    integer, parameter :: itdim = 4
    integer, parameter :: iydim = 4*ineq
    integer, parameter :: ndig = 10
    double precision, parameter :: tol = 10.0d0 **(-ndig)

    integer :: neq, kord(ikdim), iopt(6)
    double precision :: tspecs(itdim), y(iydim), t, h, delt, tfinal
    double precision :: f(ifdim)

    neq         = 2
    t           = 0.d0
    h           = 1.d0
    delt        = 6.283185307179586477d0
    tfinal      = 2.d1
    y(1:4)      = [1.d0, 0.d0, 0.d0, 1.d0 ]
    tspecs(1:4) = [t, h, delt, tfinal]

    ! set option for error control, local absolute error < tol.
    iopt(1) = 16
    iopt(2) = 6;    kord(6) = 2
    iopt(3) = 3;    f(3) = tol

    ! group the system to be treated as a single unit, set tolerance value
    ! set option for second order equations
    iopt(4) = 17
    iopt(5) = 2

    ! flag end of options
    iopt(6) = 0

    ! do the integration
    kord(1) = 0
    do
        call diva(tspecs, y, f, kord, neq, divaf, divao, &
                itdim, iydim, ifdim, ikdim, iopt)
        if (kord(1) == 1) exit
    end do

    contains
!*************************************************************************

    !*************************************************************************
        subroutine divaf(t, y, f, kord)

        !! Sample derivative subroutine for use with DIVA
        !! This evaluates derivatives for a simple two body problem.

        implicit none

        integer :: kord
        double precision :: t
        double precision :: y(4)
        double precision :: f(2)

        double precision :: tp

        ! Evaluate the derivatives
        tp = y(1)*y(1) + y(3)*y(3)
        tp = 1.d0 / (tp * sqrt(tp))
        f(1) = -y(1) * tp
        f(2) = -y(3) * tp

        end subroutine divaf
    !*************************************************************************

    !*************************************************************************
        subroutine divao(tspecs, y, f, kord)

        !! Sample output subroutine for use with DIVA.
        !! This subroutine gives output for a simple two body problem.

        implicit none

        integer :: kord
        double precision :: tspecs(4)
        double precision :: y(4)
        double precision :: f(2)

        ! Do the output
        if (kord == 1) then
           write (*, '(12X,A//,8X,A,13X,A,11X,A,9X,A)') &
                     'RESULTS FOR A SIMPLE 2-BODY PROBLEM',&
                     'T','U/V','UP/VP','UPP/VPP'
        end if
        write (*, '(1P,SP,4E15.6/ 15X,3E15.6/ A)') &
                  tspecs(1), y(1), y(2), f(1), y(3), y(4), f(2), ' '

        end subroutine divao
    !*************************************************************************

    end program drdiva
!*************************************************************************


