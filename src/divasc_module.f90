!*************************************************************************
!>
!  The secondary common block for the package.  This contains
!  variables that are required for doing interpolation and is
!  separate to simplify saving the variables that are required
!  when the solution is being dumped (saved).

    module divasc_module

    use diva_constants

    implicit none

    public

    double precision :: tn       !! The value of TSPECS(1) at the conclusion of the last step.
    double precision :: xi(kdim) !! XI(K) = TSPECS(1) - value of TSPECS(1) K steps previous.
    integer          :: iopst    !! Intended for possible use in stiff equations.
    integer          :: kordi    !! Order of differential equation being integrated.  If
                                 !! all orders are the same, this set once at the beginning.
    integer          :: kqmaxd   !! Maximum integration order used for stiff equations.
    integer          :: kqmaxi   !! Maximum integration order used for nonstiff equations.
    integer          :: ldt      !! Used to keep track of state of difference table.
                                 !!
                                 !!  * -5  Used only on first step to indicate that an extra iteration
                                 !!        is desired to get a firm estimate on the error.
                                 !!  * -4  Set on initialization before there is any difference table.
                                 !!  * -3  Set just after predicting, interpolation is not allowed when
                                 !!        this value is set.
                                 !!  * -2  Set when difference table is to be updated to the end of the
                                 !!        current step, but no interpolation is to be done.  (For
                                 !!        dumping the solution.)
                                 !!  *  0  Calculations for current step are complete, it is not o.k. to
                                 !!        update the difference table.
                                 !!  *  1  Difference table has been updated to the end of the current
                                 !!        step, e.g. by doing an interpolation.
                                 !!  *  2  Set when doing a special interpolation during computation of
                                 !!        derivatives.  (For delay equations.)
    integer          :: maxdif   !! Maximum differentiations required for stiff equations.
    integer          :: maxint   !! Maximum integrations required.  (= max. order of
                                 !! differential equations if equations are not stiff.)
    integer          :: nkdko    !! If this is nonzero (option 17), it gives the location
                                 !! in KORD() where a vector defining the order of each equation is
                                 !! specified.
    integer          :: nte      !! Total number of equations being integrated = NEQ.
    integer          :: nyny     !! Location in Y() where the base value for Y() is saved.
    integer          :: ndtf     !! Location in F() where difference table starts.
    integer          :: numdt    !! Maximum allowed number of differences available for doing an integration.

    integer :: ivc1(12) !! Array used for output of variables IOPST to NUMDT in common block DIVASC.
    double precision :: tneq(1) !! Array of dimension 1 equivalenced to TN so that an
                                !! array can be passed to *MESS.

    common / divasc / tn, xi, iopst, kordi, kqmaxd, kqmaxi, ldt, &
                      maxdif, maxint, nkdko, nte, nyny, ndtf, numdt

    equivalence (ivc1(1), iopst)
    equivalence (tneq, tn)

    end module divasc_module
!*************************************************************************
