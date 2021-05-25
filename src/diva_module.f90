!*************************************************************************
!>
!  Modernized version of the DIVA Variable Order Adams Method
!  Ordinary Differential Equation Solver From the MATH77 library.
!
!### Original Copyright
!
! Copyright (c) 1996 California Institute of Technology, Pasadena, CA.
! ALL RIGHTS RESERVED.
! Based on Government Sponsored Research NAS7-03001.


    module diva_module
!*************************************************************************

    use iso_fortran_env, only: wp => real64

    implicit none

    real(wp),dimension(5),parameter :: d1mach = &
        [tiny(1.0_wp), &
         huge(1.0_wp), &
         real(radix(1.0_wp),wp)**(-digits(1.0_wp)), &
         epsilon(1.0_wp), &
         log10(real(radix(1.0_wp),wp))]

    contains
!*************************************************************************

!*************************************************************************
!>
!  Main routine.
!
!--D replaces "?": ?IVA,?IVAA,?IVABU,?IVACO,?IVACR,?IVAEV,?IVAF,?IVAHC,
!-- & ?IVAG,?IVAIN,?IVAMC,?IVAO,?IVAOP,?IVAPR,?IVASC,?IVACE,?IVAIE,
!-- & ?IVAPE,?MESS
!
! Note a "*" at the start of a name is used to indicate "D" for the
! double precision version and "S" for the single precision version.
!
! When converting between precisions, don't forget to change the value
! of KDIM set in parameter statements in a variety of routines, and to
! adjust comments for the data statements associated with EIBND in
! *IVACR, and B in *IVAHC.
!
! Entries
!  *IVA    Main entry for starting the package.
!  *IVAA   Main program inside the package, calls the other routines,
!          and does checks for output, and noise.  Called by the user
!          if reverse communication is used.
!  *IVABU  Back ups the solution to the current base time, if a step
!          that has been started must be taken over for some reason.
!  *IVACO  Called by user to get certain information from the common
!          blocks.
!  *IVACR  Corrects the solution, estimates errors, and selects order.
!  *IVADB  Subroutine to assist in debugging codes.  Called by user to
!          get a formatted list of all the variables used in the
!          integration.  Not required in usual case.
!  *IVADE  Needed only for delay differential equations.  This is called
!          by the user from the derivative subprogram.
!  *IVAG   Required only if the user has G-Stops, i.e. places to call
!          his output subroutine where certain functions have zeroes.
!  *IVAHC  Compute coefficients that depend on the step size history.
!  *IVAIN  Used to interpolate to arbitrary points.
!  *IVAOP  Used to process user option requests.
!  *IVAPR  Used to update the differences and to predict the solution
!          at the end of the current step.
!
! External Routines
!  *1MACH  Not used in the Fortran 95 version.  ("*" is "D" for double
!          and "R" for single precision.) This returns constants that
!          depend on the floating point arithmetic.  Input arguments of
!          1 to 4 give respectively:  underflow limit, overflow limit,
!          smallest relative difference between two floating point
!          numbers, and the largest relative difference between two
!          floating point numbers.
! DERIVS (formal) Name of subroutine to be called for computing

!  OPTCHK  Used in checking storage allocation.
!  *MESS   Used to output error messages and diaganostic messages.
!          (Just MESS if no floating point is output.)
!  *ZERO   Called only if *IVAG is used.  Iterates to find zeros of
!          arbitrary (continuous) functions.
!
! Common blocks -- As a left over from the distant past, some variables
!   are in common so that they would be saved.
!  *IVAEV  Holds variables that depend on the environment.
!  *IVAMC  The main common block for the package.
!  *IVASC  The secondary common block for the package.  This contains
!          variables that are required for doing interpolation and is
!          separate to simplify saving the variables that are required
!          when the solution is being dumped (saved).
!
! Common variables and local variables
! ALPHA  (*IVAMC) Array with I-th entry = (current step size) / XI(I).
!   Used in computing integration coefficients.
! B      (*IVAHC) Array used to get started on computing integration
!   coefficients.  B(K) = 1. / (K*(K+1))
! BAKMIN (*IVADE) The largest delay at the initial point.
! BETA   (*IVAMC) Array with I-th entry = product (K=1,I-1) of
!   (current (XI(K)) / XI(K) from previous step),  BETA(1)=1.  Used in
!    updating the difference tables.
! C      (*IVAIN) Array used to hold integration/interpolation coeffs.
! C0     Parameter = 0. (in *IVAA,DE,CR,A,G,HC,IN,OP,PR)
! C1     Parameter = 1. (in *IVA,A,CR,DA,HC,IN,OP)
! C10    Parameter = 10. (in *IVAA,CR,OP)
! C1000  Parameter = 1000. (in *IVACR)
! C16    Parameter = 16. (in *IVAA,OP)
! C1M3   Parameter = .001 (in *IVAA)
! C1M5   Parameter = .00001 (in *IVAA)
! C1P125 Parameter = 1.125 (in *IVAA,HC,OP)
! C1P3   Parameter = 1.3 (in *IVAA)
! C1P4   Parameter = 1.4 (in *IVACR)
! C2     Parameter = 2. (in *IVAA,DE,BU,CR,IN,OP)
! C20    Parameter = 20. (in *IVACR)
! C2P5M3 Parameter = .0025 (in *IVAA)
! C4     Parameter = 4. (in *IVACR,OP)
! C40    Parameter = 40. (in *IVACR)
! C4096  Parameter = 4096. (in *IVAA)
! C6     Parameter = 6. (in *IVAA)
! C8M3   Parameter = .008 (in *IVAA)
! CM2    Parameter = -2. (in *IVACR)
! CM8    Parameter = -8. (in *IVACR)
! CMP5   Parameter = -.5 (in *IVACR)
! CMP75  Parameter = -.75 (in *IVAOP)
! CP0625 Parameter = .0625 (in *IVAA)
! CP1    Parameter = .1 (in *IVAA,CR,DA,HC)
! CP125  Parameter = .125 (in *IVACR)
! CP25   Parameter = .25 (in *IVAA,CR,DE,OP)
! CP3    Parameter = .3 (in *IVAA,OP)
! CP4    Parameter = .4 (in *IVAA)
! CP5    Parameter = .5 (in *IVAA,CR,DA,DE,HC,OP)
! CP5625 Parameter = .5625 (in *IVAHC)
! CP625  Parameter = .625 (in *IVAOP)
! CP75   Parameter = .75 (in *IVACR,OP)
! CP8    Parameter = .8 (in *IVACR)
! CP875  Parameter = .875 (in *IVAA, OP)
! CP9    Parameter = .9 (in *IVAOP)
! CP9375 Parameter = .9375 (in *IVACR)
! CQ3125 Parameter = .03125 (in *IVACR)
! CRBQI  Parameter = .421875 (in *IVAHC)  Initial val for computing RBQ.
! CSUM   (*IVAIN) Array used to contain partial sums of the integration
!   coefficients.  This is used to corrrect for a difference table that
!   has not yet been updated.
! D      (*IVAMC) Array to be used later to store coefficients for
!   integrating stiff equations.
!   derivatives.  Not used if option 13 is set.
! DISADJ (*IVAA) Value of stepsize when discontinuity is indicated.
! DNOISE (*IVAMC) Used in determining if noise is limiting the
!   precision.  It is usually |highest difference used in correcting|
!   of the equation with the largest error estimate.
! DS     (*IVAMC) Array to be used later to store coefficients for
!   estimating errors when integrating stiff equations.
! DVC2   (*IVADB) Array used for output of variables HC to TOUT in
!   common block *IVAMC.
! E      (*IVACR) (Estimated error) / (Requested accuracy)
! EAVE   (*IVAMC) This is a weighted average of past values of EIMAX.
!   It is adjusted to account for expected changes due to step changes.
! EEPS10 (*IVAEV) = 10. * (machine epsilon).
! EEPS16 (*IVAEV) = 16. * (machine epsilon).
! EEPS2  (*IVAEV) =  2. * (machine epsilon).
! EEPT75 (*IVAEV) = (machine epsilon) ** (.75)
! EI     (*IVACR) Estimate for what E would be if step size increased.
! EIBND  (*IVACR) Array containing limits on the estimated error with
!   the stepsize increased.  This array tends to make the code a little
!   more conservative on step size increases at low order.
! EIMAX  (*IVAMC) Estimate of (error estimate / error requested) if the
!   step size should be increased.
! EIMIN  (*IVAMC) An error estimate is small enough to allow a step
!   increase if the estimate of ((error with the step size increased) /
!   (error requested)) is less than EIMIN.
! EIMINO (*IVAA) Set to C8M3 and never changed.  When step size is being
!   reduced if EIMIN <= EIMINO then the reduction factor is set to
!   CP875.  This variable could be a parameter.
! EMAX   (*IVAMC) Largest value computed for (error estimate) / (error
!   requested).
! EOVEP2 (*IVAEV) = EEPS2 * (largest floating point number).
! EPS    (*IVACR) Current absolute error tolerance.  Also used for
!   temporary storage when computing the desired value of EPS.
! ERCOEF (*IVACR) (Error coefficient from formula) / EPS
! EREP   (*IVAMC) If EMAX > EREP, a step is repeated.  Ordinarily
!   this has the value .3.  This is set < 0 if the error tolerance is
!   specified improperly, and is set to a large value if the user
!   requests complete control over the step size.  EREP is also set
!   < 0 after a user specified discontinuity.
! EROV10 (*IVAEV) = 10. / (largest floating point number).
! ETA    (*IVAIN) Array used in computing integration/interp. coeffs.
! EVC    (*IVADB) Array used for output of variables EEPS2 to EROV10 in
!   common block *IVAEV.
! EXR    (*IVAA) Set to CP1 and never changed.  If it is estimated the
!   the (error estimate) / (error requested) on the next step will be
!   >= EXR then the step size is reduced.  Could be a parameter.
! F      (formal) Array used to store derivative values, the difference
!   tables, error tolerance requests, and values used by some other
!   options. (in *IVA,A,BU,CR,DA,DB,G,IN,PR)
! FDAT  (*IVAMC) Used to store data for error messages.  (Local array in
!   *IVAIN.)
! FOPT  (formal) in *IVAOP.  Passed as place to save floating point data
!   for options.  This package passes F in for FOPT when calling *IVAOP.
! G      (*IVAMC) Integration coefficients used for predicting solution.
!   G(I, J) gives the I-th coefficient for integrating a J-th order
!   differential equation.  G(1, 1) is equal to the step size.
! GAMMA  (*IVAIN) Array used in computing integration/interp. coeffs.
! GG     (*IVAHC) Array of length = max. differential equation order
!   allowed by code - 1.  GG(K) = (HH**(K+1)) / K!
! GNEW   (formal) in *IVAG.  Current value for vector function g, whose
!   zeroes are to be found.
! GOINT  (*IVACR) Used for assigned go to used in computing integration
!   coefficients.
! GOLD   (*IVAG) Previous value for element of G whose zero search is
!   active.
! GS     (*IVAMC) Integration coefficients used in estimating errors.
! GT     (formal) in *IVAG.  Previous value of GNEW.
! HC     (*IVAMC) Ratio of (new step size) / (old step size)
! HDEC   (*IVAMC) Default value to use for HC when reducing the step
!   size.  (Values closer to 1 may be used some of the time.)
! HH     Equivalenced to G(1,1) = current step size in *IVAA,CR,DA,G,HC.
! HI     (*IVAIN) Step length from the base value of the independent
!   variable for the interpolation.
! HINC   (*IVAMC) Default value to use for HC when increasing the step
!   size.  (Values closer to 1 may be used some of the time.)
! HINCC  (*IVAMC) Actual value used for default value of HC when
!   increasing the step size.  Set to HINC after start is considered
!   complete.  During the start HINCC is set to 1.125.
! HMAX   (*IVAMC) Largest value allowed for abs(step size).  Default
!   value is a very large number.
! HMAXP9 (*IVAMC) .9 * HMAX.
! HMIN   (*IVAMC) Smallest value allowed for abs(step size).  Default
!   value is 0.
! HNEW   (*IVADE) Value of step size when iterating at initial point
!   for delay differential equations.
! I      Used for temporary storage. (*IVAA,BU,CR,DA,DE,G,IN,OP,PR)
! IA     (*IVAOP) absolute value of first integer stored for an option.
! ICF    (*IVAMC) Final index for current loop in *IVACR.  Required by
!   option 18.
! ICI    (*IVAIN) Temporary index, = 0 for interpolation, 1 or 0 for
!   differentiation, and d-1, d-2, ... 0 for integration, where d is the
!   order of the differential equation.  Index of first location
!   in C() used is ICI + an offset.
! ICS    (*IVAMC) Starting index for current loop in *IVACR.
! ID     (formal) Array use to contain integer data from common.  Values
!   are returned in locations 1 to 5 as follows.
!   1    KEMAX = Index of equation with largest error estimate
!   2    KSTEP = Current step number
!   3    NUMDT = Number of differences used for each equation
!   4            Reserved for future use
!   5            Reserved for future use
! IDAT   (*IVAMC) Used to store integer for error messages.  (Also used
!   in *IVAA for temporary storage of KORD(2).  (Local array in *IVAIN.)
! IDE    (*IVADE - formal) Array used to contain past information so
!   that delays can stretch back indefinitely.  If the first location is
!   0, then any interpolations requested must be in the range of the
!   current difference tables.  At present, only the value 0 is allowed
!   in IDE(1).  This array is intended for the support of saving long
!   past histories.  IDE(2) must contain the declared dimension of WDE.
! IDEF   (*IVADE -- formal) Flag giving indicaion of what is going on.
!   = 0  User should compute derivatives and return to the main
!        integrator.
!   = 1  Code is computing additional values in order to get past data
!        necessary for starting.  User should compute derivatives and
!        call *IVADE.
!   < 0  Indicates an error condition.  If *IVADE is called without
!        changing the value of IDEF, the integration is stopped and an
!        error message printed.  Possible error flags are:
!    -1  Difference tables do not span back far enough to compute the
!        past values of Y needed.
!    -2  There is not enough space in WDE to get the required starting
!        values.
! IDIMF  (formal) Declared dimension of F().
! IDIMK  (formal) Declared dimension of KORD().
! IDIMT  (formal) Declared dimension of TSPECS().
! IDIMY  (formal) Declared dimension of Y().
! IDT    (*IVAIN) Used as a base index into the difference table.
! IFLAG  (formal in *IVAG) Used for communication with user.
!   = 1  Continue as if *IVAG was not called.
!   = 2  Check KORD(1) as one would do at start of OUTPUT if no G-Stops
!        were present. (Exit if in DERIVS.)
!   = 3  Return to the integrator.
!   = 4  Compute G and return to *IVAG.
!   = 5  A G-Stop has been found, and NSTOP gives its index.  (If NSTOP
!        < 0, the stop was an extrapolating stop.)
!   = 6  Same as 5, but requested accuracy was not met.
!   = 7  Same as 5, but there is a probable error in computing G.
!   = 8  Fatal error of some type.  (An error message has been printed.)
! IG     (*IVAG)  IG = KORD(2) on the initial entry (0 for extrapolating
!   G-Stops, and 1 for interpolating).
! IGFLG  (*IVAMC) Used primarily in *ivag, but also used in *iva to keep
!   track of the state of GSTOP calculations.
!   = -2 Extrapolatory G's initialized, but not the interpolatory.
!   = -1 Interpolatory G's initialized, but not the extrapolatory.
!   =  0 Set when integration is started or restarted, or option setting
!        GSTOP is set.
!   =  1 Iterating to find a GSTOP.
!   =  2 User told that a GSTOP was found.
!   =  3 Checking G's at point where a GSTOP was located.
!   =  4 Checking G's at a T output point.
!   =  5 Usual case, no sign change detected.
! IGSTOP (*IVAMC) IGSTOP(k) is set in *ivag to the index of the last G
!   with a 0, where k is one for an interpolatory G-Stop, and k is two
!   for an extrapolatory G-Stop.
! IGTYPE (*IVAMC) Array with two elements as for IGSTOP, but this saves
!   a flag giving the nature of convergence to the stop.
!   = 0  All known G-stops completely processed.
!   = 4  Need to compute next value while iterating.
!   = 5  Got good convergence.
!   = 6  Got convergence, but not to desired accuracy.
!   = 7  Problem in getting convergence.
!   = 8  A fatal error of some type.
! IHI    (*IVA) Last location used by the current option.
! ILGREP (*IVAMC) Used when correction to keep track of equations that
!   are to use a certain error tolerance.
! ILGROR (*IVACR) Index of last equation in the current group of
!   equations grouped for selecting integration order.
! ILOW   (*IVA) First location used by the current option.
! INCOM  (*IVADE) Array equivalenced to LDT in the common block *IVASC.
!   Used to simplify saving information in the common block.
! INCOP  (*IVAOP) Array containing data giving the amount of space in
!   IOPT used for each of the options.
! INGS   Current index for G-stop being examined in DIVAG.
! INICAS (*IVADE) Used to track the initialization for a delay equation.
!   = 1  Very beginning.
!   = 2  Getting derivative at the very beginning.
!   = 3  Getting derivatives at points prior to the initial point.
!   = 4  Getting derivative at initial point after iteration is started.
! INTCHK (*IVA) Array passed to OPTCHK containing information on storage
!   allocation.  See comments in OPTCHK for details.
! INTEG  (*IVAIN) Number of integrations being done. (<0 for
!   differentiations and =0 for interpolation.)  Also used as counter
!   when computing integration coefficients.
!        (*IVAPR) Number of integrations being done.
! INTEGS (*IVAPR) = -1 for equations that are not stiff, 0 for those
!   that are stiff.
! INTEGZ (*IVAIN) min(INTEG, 0)
! INTERP (*IVAIN) added to the usual integration order to get the order
!   to be used when interpolating: 3-KQMAXI, if HI=0; 1, if
!   |HI| > |XI(1)| and HI * XI(1) < 0; 0, otherwise -- the usual case.
! IOP10  (*IVAMC) Number of times diagnostic output is to be given when
!   leaving *ivacr (the corrector).
! IOP11  (*IVAMC) Gives current step number of the method.  Tells how
!   many of certain coefficients must be computed. (Has nothing to do
!   with options.) = min(max integ order + 1, KDIM).  Also set when
!   starting to flag that certain memory locations must be set to 0.
! IOP12  (*IVAMC) Points to location in F() where user supplied values
!   of HINC, HDEC, HMIN, and HMAX.  (0 if option 12 not used.)
! IOP13  (*IVAMC) If not zero, reverse communication will be used for
!   getting the values of derivatives.  Associated with option 13.
! IOP14  (*IVAMC) If not zero, reverse communication will be used in
!   place of calls to the output routine.  Associated with option 14.
! IOP15  (*IVAMC) If not zero, a return will be made to the user after
!   the initialization.  Associated with option 15.  This might be used
!   to overlay *iva, some of the user's code, and perhaps *ivaop.
! IOP16  (*IVAMC) Points to location in KORD() where information for
!   specifying the error tolerance is specified.  See option 16.
! IOP17  (*IVAMC) Used in initialization for option 17, afterwards this
!   cell is used by KEXIT which is equivalenced to IOP17.
! IOP18  (*IVAMC) Points to location in KORD() where information for
!   specifying a grouping of equations for derivative evaluation is
!   stored.  See option 18.
! IOP19  (*IVAMC) Points to location in KORD() where information for
!   specifying a grouping of equations for integration order control
!   is stored.  See option 19.
! IOP20  (*IVAMC) Used for option 20, gives first location in F where
!   estimated errors are to be stored.  Expected to be useful in a
!   program for solving boundary value problems using multiple shooting.
! IOP21  (*IVAMC) Was used for stiff equations option (never completely
!   coded).  The optional code still uses this (don't activate it!).
!   Now used to flag the location if F where the user has stored the
!    tolerance to use in finding G-Stops.
! IOP21S (*IVAMC) Was used for stiff equations see above.
! IOP22  (*IVAMC) Set aside for possible option for stiff equations.
! IOP3   (*IVAMC) Value set by option 3.
!   =  0 Interpolate to final point. (The default)
!   =  1 Integrate to final point.
!   = -1 Extrapolate to final point.
! IOP4   (*IVAMC) Value set by option 4.  The output routine is called
!   with KORD(1) = 4, every IOP4 steps.  (Default value for IOP4 is a
!   very large number.
! IOP5   (*IVAMC) Value provided by option 5, used to specify extra
!   output points.
! IOP6   (*IVAMC) Value provided by option 6.  If nonzero, the output
!   routine is called at the end of every step.  If > 0, there are
!   IOP6 interpolating G-Stops.
! IOP7   (*IVAMC) Value provided by option 7.  If > 0, there are K7
!   extrapolating G-Stops.
! IOP8   (*IVAMC) Value provided by option 8.  If nonzero, the output
!   routine is called with KORD(1)=8 whenever the step size is changed.
! IOP9   (*IVAMC) Value provided by option 9.  Used to specify that the
!   user wishes to save the solution.
! IOPIVA (*IVA) Used to save length of IOPT vector for error messages.
! IOPST  (*IVASC) Intended for possible use in stiff equations.
! IOPT   (formal *IVA and IVAOP) Used to specify options.
! IOPTC  (*IVAOP) In *IVAOP equivalenced so that IOPTC(3) is equivalent
!   to IOP3.
! IOPTS  (*IVAOP) Array containing the current default values to be
!   stored into IOPTC.
! IORD   (*IVACR) Index of first equation in the current group of
!   equations grouped for selecting integration order.
! IOUTKO (*IVADC) Used in *IVADI to point to KORD to keep track of
!   equation grouping for diagnostic output.
! ISVCOM (*IVADE) Used to save info. in the common block *IVASC.
! ITERS  (*IVADE) Counts iterations in starting delay differential
!   equations.  Max. value for this is arbitrarily 100.
! ITOLEP (*IVAMC) Used for temporary storage, and for the index of a
!   tolerance relative to the start of tolerances.
! IVC1   (*IVADB) Array used for output of variables IOPST to NUMDT in
!   common block *IVASC.
! IVC2   (*IVADB) Array used for output of variables ICF to NY in
!   common block *IVAMC.
! IWB    (*IVADE) Current base index for saving F values in WDE when
!   starting delay differential equations.
! IY     (*IVAMC) Used for the current index to the Y() array.  (Local
!   variable in *IVAIN used in computing IYI.)  Equivalenced to
!   IZFLAG in *IVAG.
! IYI    (*IVAIN) Y(IYI) is currently being computed.
! IYN    (*IVAIN) Y(IYN) is base Y() corresponding to Y(IYI).
! IYNI   (*IVAIN) Used as base index for computing IYN as IY is for INI.
! IYO    (*IVADE) Points to first base value of Y for current
!   interpolation when getting values for a delay differential equation.
! IZFLAG (*IVAG)  Equivalenced to IY.  Set to 0 initially, and later
!   set to the value returned by *ZERO.
!    = 0  Value set on entry at start of search.
!    = 1  Compute next g again.
!    = 2  Normal terminiation.
!    = 3  Normal termination -- error criterion not satisfied.
!    = 4  Apparent discontinuity -- no zero found.
!    = 5  Couldn't find a sign change.
!    = 6  *ZERO was called with a bad value in IZFLAG.
! J      For temporary storage. (In *IVA,A,BU,CR,DA,DB,DE,HC,IN,OP,PR)
! J1     (*IVAA & DA) Used for temporary storage.
! J2     (*IVAA) Used for temporary storage.
! JL     (*IVA) Used for checking storage.
! JLGREP (*IVACR) Contents of first location of KORD (called LGROUP in
!   *IVACR) for the current error tolerance rule.
! JLGROR (*IVACR) Contents of first location of KORD for the current
!   integration order control.
! JLIM   (*IVA) Used for checking second item in KORD list for options
!   16 and 19.
! K      For temporary storage.  (In *IVA,A,BU,CR,DA,DB,DE,HC,IN,OP,PR)
! KDIM   Parameter giving the largest number of differences supported.
!        Used in all the routines.
! KEMAX  (*IVAMC) Index associated with equation giving the largest
!   value for (estimated error) / (requested error).
! KEXIT  (*IVAMC) Equivalenced to IOP17 which is not used after
!   initialization.  Defines actions when KORD2I = -7.  (Referenced in
!   (*IVAA,DA,G).)
!   =  1  Take the step over with reduced H.
!   =  2  Take the step over.
!   =  3  Do the end of step call to OUTPUT.
!   =  4  Reset TMARK, then do same as for KEXIT = 2.
!   =  5  Reset TMARK, then do same as for KEXIT = 3.
!   =  6  Give the fatal error diagnostic.
! KFERR  (*IVA)  Temporary storage in checking for option 16.
! KGO    (*IVA)  Used to tell from whence a check is being done or an
!   error message is being written.
!   = 1 Checking an equation group for variational equations.
!   = 2 Checking an equation group for diagnostic print.
!   = 3 Checking an equation group for integration order control.
!   = 4 Checking an equation group for error control.
!   = 5 Checking an equation group for specifying ODE orders.
!   = 6 Found a problem with output type for printing.
!   = 7 Found a problem with an output group for printing.
!   = 8 Found a problem with input NEQ.
!   = 9 Order specified for the ODE's in the system is out of range.
!   =10 Option 16 was not used (an error).
!   =11 Error tolerance of 0 specified without proper flags.
! KIS    (*IVAMC) Used to check if it is time to dump the solution.
!   The check involves incrementing KIS at the end of the step, and
!   dumping the solution if KIS is 0.
!   = -1  Set in *ivacr when it is time to dump solution
!   =  0  When starting
!   =  2  After being dumped.
!   This is set to 1000 just after a user specified discontinuity, and
!   counted up from that point.
! KMARK  (*IVAMC) Identifies the type of output associated with the next
!   output point specified by TSPECS.
! KONV   (*IVADE) Counts iterations.  Test for convergence if KONV > 1.
! KORD   (formal in *IVA,A,BU,CR,DA,DB,DE,G,IN,PR) KORD(1) is used to
!   return flags for the user to test, and KORD(2) tells what routine
!   the flag is associated with.  See KORD1I and KORD2I below and the
!   write up for the program.  KORD(3) is used for communicating extra
!   information to the user in some cases.  KORD(4) to KORD(NTE+3) are
!   used for integration order for the equations, and the rest of KORD()
!   is available for user options.
! KORD1I (*IVAMC) Helps in defining the state of the integrator.
!   Frequently has the same value as KORD(1).  Meaning depends on the
!   value of KORD(2), or the value about to be assigned to KORD(2).
!   <  0  Happens when preparing to give output with extrapolation.
!   =  0  Happens when checking F at points for noise test.
!   =  1  (KORD(2)=-1)  End of integration has been reached.
!   =  1  (KORD(2)= 0)  Computing first predicted derivative.
!   =  1  (KORD(2)= 1)  Output for initial point.
!   =  2  (KORD(2)=-1)  Giving diagnostic for noise limiting precision.
!   =  2  (KORD(2)= 0)  Computing corrected derivative.
!   =  2  (KORD(2)= 1)  Output for TSPECS(3).
!   =  3  (KORD(2)=-1)  Diagnostic for step size reduction too fast.
!   =  3  (KORD(2)= 0)  Computing variational derivative.
!   =  3  (KORD(2)= 1)  Output for TSPECS(4).
!   =  4  (KORD(2)=-1)  Error, discontinuity.
!   =  4  (KORD(2)= 1)  Output for certain number of steps.
!   =  5  (KORD(2)= 0)  Get initial derivatives for stiff equations.
!   =  5  (KORD(2)= 1)  Extra output from TSPECS.
!   =  6  (KORD(2)= 1)  End of step output.
!   =  7  (KORD(2)= 0)  Evaluate G before extrapolated output point.
!   =  7  (KORD(2)= 1)  Evaluate G before extrapolated output point.
!                       (Also used when checking for other G's after
!                        finding one.)
!   =  8  (KORD(2)= 1)  Tell user step size has changed.
!   =  9  (KORD(2)= 1)  Request for user to save solution.
!   = 11  (KORD(2)=-1)  Error, step size too small at end of start.
!   = 12  (KORD(2)=-1)  Error, step size is too small.
!   = 13  (KORD(2)=-1)  Error, output points specified badly.
!   = 21  (KORD(2)=-1)  H too small to give reasonable change when added
!                       to T.
!   = 22  (KORD(2)=-1)  Error, bad tolerance.
!   = 23  (KORD(2)=-1)  Set after message for a fatal error.
!   = 24  Set on error message in *iva, along with KORD2I = -4.
!   Also used as an index into MLOC in *IVAA when an error is being
!   processsed, see MLOC below.
! KORD2I (*IVAMC) Helps in defining the state of the integrator.
!   Frequently has the same value as KORD(2).
!   = -3  Set in *ivag, to get a derivative evaluation.
!   = -2  Set in *ivag, to get another entry to OUTPUT.
!   = -1  Return to calling program, done, interrupt, or got an error.
!   =  1  Calling OUTPUT or returning to user for OUTPUT type action.
!   =  0  Calling DERIVS or returning to user for DERIVS type action.
!   = -4  Error message in *iva and in *ivaop, along with KORD1I = 24.
!   = -5  Starting
!   = -6  Starting, getting the initial derivative value or derivatives
!         for the noise test.
!   = -7  Done some extrapolation, KEXIT defines the action to take.
!         Set in *ivag to activate KEXIT action in *iva.
!   = -8  Set when user has requested adjustment of the difference
!         tables for a discontinutiy.
! KORDI  (*IVASC) Order of differential equation being integrated.  If
!   all orders are the same, this set once at the beginning.
! KOUNT   (*IVADE) Count of number of points back from the initial point
!   when solving a delay differential equation.
! KOUNTM  (*IVADE) Largest value currrently allowed for KOUNT.
! KOUNTX  (*IVADE) Largest value allowed for KOUNTM.
! KOUTKO  Used in DIVACR to track where output is wanted.
! KPRED  (*IVAMC) Value assigned to KORD1I when getting a predicted
!   derivative.  (1 used now, 5 planned for use with stiff equations.)
! KQD    (*IVACR) = max(2, integration order)
! KQDCON (*IVAMC) Number of coefficients computed with constant step
!   size for stiff equations.
! KQICON (*IVAMC) Number of coefficients computed with constant step
!   size for nonstiff equations.
! KQL    (*IVACR) Integration order at start of (*IVACR)
! KQLORD (*IVACR) Saved value of KQL when equations are grouped for
!   controlling the integration order.
! KQMAXD (*IVASC) Maximum integration order used for stiff equations.
! KQMAXI (*IVASC) Maximum integration order used for nonstiff equations.
! KQMAXS (*IVAMC) Maximum integration order for equations that have
!   some limit on the error that can be committed.
! KQMXDS (*IVAMC) Used to save KQMAXD in case step is repeated and the
!   solution must be dumped.
! KQMXI  (*IVAIN) Maximum integration order used for integration or
!   interpolation, = KQMAXI+INTERP-1.
! KQMXS  (*IVAIN) Maximum step number, = max(KQMXI, KQMAXD).
! KQMXIL (*IVAMC) Value of KQMAXI the last time integration coefficients
!   were computed.
! KQMXIP (*IVAMC) = KQMAXI + MAXINT, for computing integration coeffs.
! KQMXIS (*IVAMC) Used to save KQMAXI in case step is repeated and the
!   solution must be dumped.
! KQN    (*IVACR) Value of integration order at end of *IVACR.
! KQQ    Used for the integration order for current equation.  (Values
!   < 0 are intended for stiff equations.)  (In *IVA,BU,DA,IN,PR)
! KSC    (*IVAMC) Number of steps that have been taken with a constant
!   step size.
! KSOUT  (*IVAMC) When KSTEP reaches this value, the output routine is
!   called with KORD(1) = 4.  The default value is a very large number.
! KSSTRT (*IVAMC) Set when ending one derivative per step to KSTEP + 2.
!   Checked later in *IVAHC to decide whether to set the step changing
!   factors to their nominal values.
! KSTEP  (*IVAMC) Number of steps taken since the start of integration.
! L      Used for temporary storage.  In *IVAIN, L is the initial value
!   of LDT, except L=1 if LDT=-1, and MAXINT >= 0.  (Used in *IVAA,BU
!   CR,DA,DB,IN,PR.)
! LAHAG  (*IVADB) Used to get proper offset into an diagnostic message.
! LAIAG  (*IVADB) Used to get proper offset into an diagnostic message.
! LDIS   (*IVAA) Count of steps since user flagged a discontinuity.
! LDT    (*IVASC) Used to keep track of state of difference table.
!   = -5  Used only on first step to indicate that an extra iteration
!         is desired to get a firm estimate on the error.
!   = -4  Set on initialization before there is any difference table.
!   = -3  Set just after predicting, interpolation is not allowed when
!         this value is set.
!   = -2  Set when difference table is to be updated to the end of the
!         current step, but no interpolation is to be done.  (For
!         dumping the solution.)
!   =  0  Calculations for current step are complete, it is not o.k. to
!         update the difference table.
!   =  1  Difference table has been updated to the end of the current
!         step, e.g. by doing an interpolation.
!   =  2  Set when doing a special interpolation during computation of
!         derivatives.  (For delay equations.)
! LEX    (*IVAMC) Indicates how to get values at next output point:
!   = -1  Extrapolate
!   =  0  Interpolate (The usual case.)
!   =  1  Integrate to the output point, integration is not continued.
! LGO    (*IVAIN) Used as an an assigned go to.  Result is to add in
!   extra correction term when LDT has been set to 2.
! LGROUP (formal) This is a part of KORD passed into *IVACR.  The first
!   location is the start of the information on the grouping of
!   equations for error control.
! LINC   (*IVAMC) Used to indicate state of step size selection.
!   = -10 After computed derivatives at base time, after computing other
!         extra derivatives for the noise test.
!   = -9  After computed second extra derivative for noise test.
!   = -8  After computed first extra derivative for noise test.
!   = -7  Dumping the solution and then doing a user initiated restart,
!         or getting ready to compute extra derivatives for the noise
!         test.
!   = -6  Dumping the solution before a restart.
!   = -5  Set on the first step, and also set when dumping the solution
!         after a discontinuity.
!   = -4  Repeat step with no change in the step size.
!   = -3  Set when the error tolerance is set improperly.
!   = -2  User has complete control of selecting the step size.
!   = -1  Step is being repeated.
!   =  0  Step size is not to be increased on this step.
!   = k>0 Step size can be increased by HINCC**k.
! LINCD  (*IVAMC) Value of smallest k for which HINCC**k >= 2.
!   (=-2 if user is specifying all step size changes.)
! LINCQ  (*IVAMC) Value of smallest k for which HINCC**k >= 4.
! LIOPT  (*IVAOP) Value of the last index in IOPT on the last call.
!   Used so *IVA can print IOPT in error messages.
! LL     (*IVACR) Temporary variable used when equations are grouped
!   for integration order control.
! LNOTM1 (*IVAIN) Logical variable = L /= -1.  If LNOTM1 is true,
!   storage in Y() is different in some way lost to antiquity.  Such
!   a case can only arise in the case of stiff equations.
! LOCF1  (*IVADB) Gives packed data needed for output of tables by the
!   message processor MESS.  See comments there under METABL for defs.
! LOCF2  (*IVADB) As for LOCF1 above.
! LOCM   (*IVAA) Parameter = 32*256, used to unpack integers stored
!   in MLOC for use in error message processing.
! LPRINT (formal, *IVADB) Defines how much printing is to be done in
!   *IVADB.  Let |LPRINT| = 10*N1 + N2     (N1,N2 digits)
!    N1=1   Do not print any variables external to the integrator.
!    N1=2   Print  tspecs, current y, past y, current f, all pertinent
!           contents of KORD, and TOL.
!    N1=3   Above + difference tables up to highest difference used.
!    N1=4   Same as N1=1 + all in storage allocated for differences.
!    N2=1   Do not print any variables internal to the integrator.
!    N2=2   Print all scalar variables in interpolation common block.
!    N2=3   Above + all scalar variables in main integ. common block.
!    N2=4   Same as N1=3 + all used in arrays XI,BETA,ALPHA, first
!           column of G, GS,RBQ,SIGMA
!    N2=5   Same as N1=4 + all used in arrays G,D,DS,V
! LSC    (*IVAMC) Indicates if starting or if noise may be present.
!   =k<0 -k steps have been taken for which noise appears to be limiting
!        the precision.
!   = 0  Usual case
!   = 1  Doing 1 derivative per step after initial part of start.
!   = 2  Used as flag that it is time to set LSC=0.
!   = 3  Third step, hold the order constant.
!   = 4  Second step, increase orders from 2 to 3.
!   = 5  First step, third time through the first step (if required).
!   = 6  First step, second time through.
!   = 7  First step, first time through.
!   = 8  Set on initialization.
! LTXT?? Names of this form are used in setting up data statements for
!   error messages.  These names are generated automatically by PMESS,
!   the program that makes up these messages.
! LX     (*IVAA) Used for temporary storage in computing TMARKA().
!        ( formal *IVADE)  An integer array containing extra
!   information, as follows.
!  LX(1) Points to a location in Y beyond those already in use.  Values
!        of Y requested are computed at TSPECS(1) - Y(LX(1)) and stored
!        starting at Y(LX(1)+1).  If this index is 0, no more extra Y
!        values are to be computed.
!  LX(2) Index of the first equation for which the Y's above are to be
!        computed.  Y(LX(1)+1) will correspond to this first equation
!        index.
!  LX(3) Index of the last equation for which the Y's above are to be
!        computed.  Thus the Y's stored starting at Y(LX(1)+1) will
!        require no more space than half the space ordinarily required
!        for the array Y(), and may require significantly less.
!  LX(4) Maximum number of times to integrate F to get Y.  This should
!        be > 0, and less than or equal to the order of the highest
!        order differential equation.  (= 0 is allowed, but probably
!        not what you want.  It would give a value only for F.)  Space
!        must be set aside for all integrals of F, even if not all are
!        requested.  For a first order system, all Y's are just the
!        first integrals of the corresponding F's.  For higher order
!        equations, the first Y associated with a given F is the d-th
!        integral of the corresponding F, where d is the order of the
!        equation, and the last Y corresponding to the F is the first
!        integral of that F.
!  LX(5) As for LX(4), but gives the index for the fewest number of
!        times to integrate F.  Ordinarily this should be > 0.  If 0 is
!        requested, an estimate for the value of F at the delay point is
!        computed.  This should not be 0 more than once, for equations
!        covering the same index, since later such requests would write
!        over the earlier results.
!  LX(5i+k) , k = 1, 2, ... 5.  Treated as for the cases above.  If i
!        different cases of delayed Y's are to be computed, then
!        LX(5i+1) must be 0.
! LX2    (*IVADE) Value of LX(5i+2), when working on the i-th delay.
! MACT   Used in the programs which call the error message program.
!   This array difines the actions to be taken by that program.  (In
!   (*IVA,A,DA,DE,G,IN,OP)
! MACT0  (*IVADB) Used to call the message program, see MACT.
! MACT?  As for MACT, in (*IVA,CR,DB)
! MACTFV (*IVADB) As for MACT0.
! MAXDIF (*IVASC) Maximum differentiations required for stiff equations.
! MAXINT (*IVASC) Maximum integrations required.  (= max. order of
!   differential equations if equations are not stiff.)
! MAXKQ  (*IVA, BU)e
! MAXKQD (*IVAMC) Largest integration order allowed for stiff equations.
! MAXKQI (*IVAMC) Largest integ. order allowed for nonstiff equations.
! ME???? Parameters defining constants used for interaction with the
!   error message program MESS.  See comments there for definitions.
!   (In *IVA,A,DA,DE,G,IN,OP)
! METHOD (*IVAMC) Defines kind of methods being used.
!   = -1  Only stiff equations are being integrated.
!   =  0  Only nonstiff equations are being integrated.
!   =  1  Both kinds of methods are required.
! MLOC   (*IVA,A,DE) Contains locations in MTEXT for error messages.  In
!   *IVAA this data is packed using MLOC??, see below.
! MLOC?? (*IVAA) Parameters constructed to aid in making up packed data
!   for processing error messages.  Low two digits give the value of
!   KORD1I to use for the error index and later processing, the next two
!   give the error severity level, and the rest point to text used for
!   the message.
! MODF2  (*IVADB) Used in constructing the same kind of packed data as
!   described for LOCF1 above.
! MULTJ  Local to DIVAOP for calls not using F.
! MTEXT  (*IVA,A,CR,IN,OP) Text for error messages.
! MTXT?? (*IVA,A,CR,DA,DB,DE,G,IN,OP) Equivalenced into MTEXT.
! N      Used for temporary storage.  (In *IVAHC,IN,PR)
! NDTF   (*IVASC) Location in F() where difference table starts.
! NE     (*IVAMC) Number of equations in the first group.  (=NTE if
!   option 18 is not used.)
! NEDDIG (*IVADB) Parameter = -MEDDIG.
! NEPTOL (*IVAMC) Used for temporary storage and to save the value of
!   ITOLEP for error messages.
! NEQ    (formal) Total number of equations being integrated.
! NG     (*IVAMC) Used in *ivag for the number of g's in the current
!   context.
! NGSTOP (*IVAG) Dimension 2 array equivalenced to IOP6, and IOP7.  To
!   get the number of interpolating and extrapolating G-Stops.
! NGTOT  (*IVAMC) NGTOT(1) gives the number of interpolating G-Stops,
!   and NGTOT(2) gives the number of extrapolating G-Stops.
! NKDKO  (*IVASC) If this is nonzero (option 17), it gives the location
!   in KORD() where a vector defining the order of each equation is
!   specified.
! NLX    (*IVADE) Temporary index used to keep track of interpolations
!   being done to get Y() values for a delay differential equation.
! NOISEQ (*IVAMC) max(2, order of equation for which (error estimate)/
!   (error requested) is a maximum).
! NOUTKO (*IVAMC) If nonzero, gives the index in KORD where information
!   on what equations are to be included in the diagnostic output is
!   given.   See option 10.
! NSTOP  (formal) In *IVAG.  Index of the G-stop, see IFLAG.
! NTE    (*IVASC) Total number of equations being integrated = NEQ.
! NTEXT  (formao *IVADB) Character variable containing heading text.
! NTOLF  (*IVAMC) First location in F() where tolerance specifying
!   accuracy desired is stored.
! NUMDT  (*IVASC) Maximum allowed number of differences available for
!   doing an integration.
! NXTCHK (*IVA) Equivalenced to INTCHK(1), which gives the next
!   available location in INTCHK for storing data on storage allocation.
! NY     (*IVAMC) Total order of the system.
! NYNY   (*IVASC) Location in Y() where the base value for Y() is saved.
! NYNYO  (*IVADE) Equivalenced to the saved value from common of NYNY.
! OUTPUT (formal) Name of subroutine to be called for the output of
!   data or for computing G-Stops.  Not used if option 14 is set.
! OVD10  (*IVAEV) (largest floating point number) / 10.
! OVTM75 (*IVAEV) (largest floating point number) ** (-.75)
! RBQ    (*IVAMC) Array containing data for the preliminary noise test.
! RD     (formal *IVACO) Array use to contain floating point data from
!   common.  Values are returned in locations 1 to 3 as follows.
!   1    EMAX =  Max. ratio of estimated error to requested error
!   2            Reserved for future use
!   3            Reserved for future use
! REF    (*IVACR) Array of length 3 used for translating error tolerance
!   type into the factor used for exponential averaging for that type.
! RND    (*IVACR) Usually the current estimated error.  Used in deciding
!   if noise is limiting precision.
! RNOISE (*IVACR) Value used in comparison with RBQ() for preliminary
!   noise test.
! ROBND  (*IVAMC) Used to influence the selection of integration order.
!   The larger ROBND, the harder it is to increase the order and the
!   easier it is to decrease it.
! RVC2   (*IVADB) Array used for output of variables DNOISE to SNOISE in
!   common block *IVAMC.  These are variables that don't require a great
!   deal of precision.
! S      (*IVACR) Estimate of (step size) * eigenvalue of Jacobian.
! SIGMA  (*IVAMC) The k-th entry of this array contains a factor that
!   gives the amount the k-th difference is expected to increase if the
!   step size in increased.  These numbers get bigger it there is a past
!   history of increasing the step size.
! SIGMAS (*IVAA) Saved value of SIGMA(k) from the last step, where k =
!   integration order for equation with index KEMAX.
! SNOISE (*IVAMC) Value used in comparison with RBQ() on equation with
!   largest value for (error estimate) / (error requested).
! T      (formal) in *IVAIN. T(1) contains the point to be interpolated
!   to, and T(2) is used in a check that |HI| <= |T(2)|.  When used by
!   other routines in this package, TSPECS is passed in for T.
! TB      (*IVADE) Base time for current interpolation.
! TC      (*IVADE) Original value of TN when getting past Y's for a
!   delay differential equation.
! TEMP   Used for temporary storage, in *IVAHC,PR
! TEMPA  (*IVACR) Array equivalenced to (TPS1,TPS2,TPS3,TPS4).
! TEMPAO (*IVACR) Array used to accumulate values in TEMPA.
! TG     (*IVAMC) TG(1) gives the last value of TSPECS(1) for which an
!   interpolatory G-Stop has been computed.  TG(2) is defined similarly
!   for extrapolatory G-Stops.
! TGSTOP (*IVAMC) TGSTOP(1) gives the value of TSPECS(1) where the last
!   0 for an interpolatory G-Stop was found.  TGSTOP(2) is defined
!   similarly for extrapolatory G-Stops.
! TMARK  (*IVAMC) Location of the next output point.
! TMARKA (*IVAA)  Array of length 2 equivalenced to TMARK (and TMARKX).
! TMARKX (*IVAMC) Location of the next output point to be found using
!   integration or extrapolation.  This variable must follow immediately
!   after TMARK in the common block.
! TN     (*IVASC) The value of TSPECS(1) at the conclusion of the last
!   step.
! TNEQ   (*IVADB) Array of dimension 1 equivalenced to TN so that an
!   array can be passed to *MESS.
! TOL    (formal) This is a part of F passed into *IVACR.  The first
!   location is the start of the information on the tolerances for error
!   control.
! TOLD   (*IVAG) Value of TSPECS(1) on one side of a zero.
! TOLG   (*IVAMC) Tolerance to pass to dzero when locating G-Stops.
! TOUT   (*IVAMC) Location of next output point defined by value of
!   TSPECS(3).  Such output is given with KORD(1) = 2.
! TP     (*IVA,A,DA,DE,HC) Used for temporary storage.
! TP1    (*IVAA,DA,HC,IN,PR) Used for temporary storage.
! TP2    (*IVAA,DA,HC,PR) Used for temporary storage.
! TP3    (*IVAA) Used for temporary storage.
! TPD    (*IVABU) Used for temporary storage.
! TPP    (*IVACR) Used for temporary storage.  Usually same as TPS3.
! TPS1   (*IVAA,CR) Used for temporary storage.  (In *IVACR is the
!   difference of order KQQ-2)
! TPS2   (*IVAA,CR) Used for temporary storage.  (In *IVACR is the
!   difference of order KQQ-1)
! TPS3   (*IVACR) Contains the difference of order KQQ.  This is the
!   last difference used in the corrector.
! TPS4   (*IVACR) Contains the difference of order KQQ+1.
! TPS5   (*IVACR) Temporary storage.
! TPS6   (*IVACR) Temporary storage.
! TPS7   (*IVACR) Temporary storage.
! TSAVE  (*IVAG) Value of TSPECS(1) before starting the search for a 0.
! TSPECS (formal *IVA,A,DB,DE,G)
!   TSPECS(1) is the current value of the independent variable.
!   TSPECS(2) is the current value of the step size.
!   TSPECS(3) is the increment to use between output points that give
!             output with KORD(1) = 2.
!   TSPECS(4) is the "final" output point.
! V      (*IVAMC) Array used in computing integration coefficients.
! XI     (*IVASC) XI(K) = TSPECS(1) - value of TSPECS(1) K steps
!   previous.
! W      (*IVAHC) Array used in computing integration coefficients.
! WDE    (formal, *IVADE)  Array used for working storage.  This storage
!   is used to save derivative values when iterating to get started.  To
!   be safe one should allow as much space as is allowed for differences
!   in F.  In most cases the start will not require this much space
!   however.  This array is also intended for the support of saving long
!   past histories.
! Y      (formal, *IVA,A,CR,DA,DB,DE,G,IN,PR) Array containing the
!   independent variable and all derivatives up to order one less than
!   the order of the differential equation.  Also use to save these
!   values at the beginning of the current step, the base values.
! YN     (formal, in *IVAPR)  Base values of y, these follow the
!   current values of the dependent variable, y, in Y().
!
!++S Default KDIM = 16
!++  Default KDIM = 20
!++  Default MAXORD = 2, MAXSTF = 1
!++  Default INTEGO, VAREQ, OUTPUT, DUMP, GSTOP, EXTRAP
!++  Default STIFF=.F., ARGM=.F., ERRSTO=.F.
!
!### History
!  * 2015-03-15 DIVA  Krogh  Removed extra call divabu after noise test
!  * 2015-03-15 DIVA  Krogh  Forced restart needs more reduction in h.
!  * 2010-02-20 DIVA  Krogh  Fixed calling DIVAOP with array other than F
!  * 2009-11-03 DIVA  Krogh  Added option 11, more variables initialized.
!  * 2009-10-30 DIVA  Krogh  Gave KSSTRT and ROBND initial values.
!  * 2009-10-30 DIVA  Krogh  Fixed reference to undefined location in F.
!  * 2009-10-21 DIVA  Krogh  Got rid of NaN in diag. print when LSC=3.
!  * 2009-10-15 DIVA  Krogh  A few changes on how noise is handled.
!  * 2002-11-12 DIVA  Krogh  Fixed problem integrating to final output pt
!  * 2002-08-29 DIVA  Krogh  Added test for invalid HMIN/HMAX.
!  * 2002-07-26 DIVA  Krogh  Added KOUTKO to fully support Option 10.
!  * 2002-05-14 DIVA  Krogh  Fix starting prob. for Option 18.
!  * 2002-05-13 DIVA  Krogh  Put exponent letter in  numbers missing them
!  * 2002-05-12 DIVA  Krogh  Added error message for bad option 5 usage.
!  * 2001-09-07 DIVA  Krogh  Changes to allow user tol on G-Stops.
!  * 2001-05-25 DIVA  Krogh  Minor change for making .f90 version.
!  * 2001-05-18 DIVA  Krogh  Less computing with no error test
!  * 2001-05-17 DIVA  Krogh  Fixed so with no error test can't start dump
!  * 2001-04-24 DIVA  Krogh  Inserted comments from ivacom.
!  * 2000-12-01 DIVA  Krogh  Removed (some of) unused C1, MAXSTF, METEXT.
!  * 1999-12-28 DIVA  Krogh  Saved S in DIVACR for output consistency.
!  * 1999-08-19 DIVA  Krogh  Removed superfluous test above label 3520.
!  * 1997-04-22 DIVA  Krogh  Got rid of assigned go to's. F=0 if diag.
!  * 1996-08-26 DIVA  Krogh  Initialize F to 0 if dumping solution.
!  * 1996-08-23 DIVA  Krogh  Print TN not TSPECS(1) in error messages.
!  * 1996-05-30 DIVA  Krogh  Changed DERIVS/OUTPUT to  DIVAF/DIVAO.
!  * 1996-04-27 DIVA  Krogh  Changes to use .C. and C%%.
!  * 1996-03-30 DIVA  Krogh  Added external statement.
!  * 1996-03-25 DIVA  Krogh  Introduced TEXT1 to comply with F77.
!  * 1996-02-27 DIVA  Krogh  Fixed so DUMP not affected by ignored eqs.
!  * 1995-12-18 DIVA  Krogh  Fixed so no solution dump on 0 length integ.
!  * 1995-11-09 DIVA  Krogh  Fixed so char. data at col. 72 is not ' '.
!  * 1995-06-19 DIVA  Krogh  Fixed prob. with discon. just after restart.
!  * 1995-05-09 DIVA  Krogh  Fixed G-Stop/discontinuity code interaction
!  * 1995-04-26 DIVA  Krogh  Use KQMAXS instead of KQMAXI when LDIS>1000.
!  * 1995-04-26 DIVA  Krogh  Keep current KQL on discontinutiy.
!  * 1994-12-16 DIVA  Krogh  Fixed option 12 with K12 < 0.
!  * 1994-11-11 DIVA  Krogh  Declared all vars.
!  * 1994-11-02 DIVA  Krogh  Changes to use M77CON
!  * 1994-09-08 DIVA  Krogh  Added CHGTYP code.
!  * 1994-07-11 DIVA  Krogh  Fix to get same state with/without var. eqs.
!  * 1994-03-07 DIVA  Krogh  Allow larger order in single precision.
!  * 1994-01-14 DIVA  Krogh  Minor change to allow changing TFINAL.
!  * 1993-04-27 DIVA  Krogh  Additions for Conversion to C.
!  * 1993-04-12 DIVA  Krogh  Converted to use slightly altered MESS.
!  * 1993-04-12 DIVA  Krogh  Fixed LSC so sol. saved when HMAX is small.
!  * 1992-10-13 DIVA  Krogh  Fixed G-Stop/discontinuity code interaction.
!  * 1992-09-21 DIVA  Krogh  Fixed bug in discontinuity code.
!  * 1992-09-09 DIVA  Krogh  Fixed bug - Var. Eqs. with discontinuities.
!  * 1992-08-07 DIVA  Krogh  Storage map printed only if option 10 /= 0
!  * 1992-07-16 DIVA  Krogh  Restored correct discontinuity code.
!  * 1992-06-16 DIVA  Krogh  Eliminate reuse of storage for option 12.
!  * 1992-04-08 DIVA  Krogh  Removed unused labels, 1020, 2120.
!  * 1992-03-30 DIVA  Krogh  Fixed bug in DIVAOP error message.
!  * 1992-03-12 DIVA  Krogh  Simplified DIVABU, more digits in B's.
!  * 1992-01-16 DIVA  Krogh  Fixed minor bug in error messages.
!  * 1991-12-03 DIVA  Krogh  Major change for improved error checks.
!  * 1991-06-17 DIVA  Krogh  Fixed bug in checking storage allocation.
!  * 1991-04-11 DIVA  Krogh  Fixed minor bug re. option 12 in DIVAOP.
!  * 1991-03-28 DIVA  Krogh  Removed check at label 650 for KORD2I<0.
!  * 1991-02-08 DIVA  Krogh  Changed some floats to generics
!  * 1990-11-08 DIVA  Krogh  Fixed bug on TSPECS on discon.
!  * 1990-09-14 DIVA  Krogh  Fixed bug when discon. and sol. save.
!  * 1990-09-13 DIVA  Krogh  Increased dimension of BETA by 1.
!  * 1990-09-13 DIVA  Krogh  Added one more poss. on rel. error test.
!  * 1990-09-11 DIVA  Krogh  Recent change messed up getting dump output.
!  * 1990-06-05 DIVA  Krogh  Fixed bug in noise test, comments in IVACOM.
!  * 1990-05-08 DIVA  Krogh  Fixed new bug when TMARK hit in DIVAG.
!  * 1990-04-17 DIVA  Krogh  Fixed minor problem in DIVAIN error msg.
!  * 1990-04-10 DIVA  Krogh  Fixed interaction between discon. & dump.
!  * 1990-03-23 DIVA  Krogh  Fixed bug on option "-2", see 1989-12-07.
!  * 1990-03-20 DIVA  Krogh  Fixed rarely occuring loop.
!  * 1990-01-29 DIVA  Krogh  Removed unneeded labels.
!  * 1989-12-14 DIVA  Krogh  Saved common block DIVAEV.
!  * 1989-12-07 DIVA  Krogh  Added option "2" to DIVAOP.
!  * 1989-11-09 DIVA  Krogh  Made GG a save var. in DIVAHC
!  * 1989-08-21 DIVA  Krogh  Fix out of bounds ref. to V in DIVABU
!  * 1989-07-26 DIVA  Krogh  Fix bug in initial dim. check
!  * 1989-07-21 DIVA  Krogh  Code for integrating discontinuities
!  * 1987-12-07 DIVA  Krogh  Initial code.

    subroutine DIVA(TSPECS, Y, F, KORD, NEQ, DIVAF, DIVAO, IDIMT,     &
                    IDIMY, IDIMF, IDIMK, IOPT)

      integer NEQ, IDIMT, IDIMY, IDIMF, IDIMK
      integer KORD(*), IOPT(*)
!--D Next line special: P=>D, X=>Q
      double precision TSPECS(*), Y(*)
      double precision F(*)
      external DIVAF, DIVAO
!
! *********************** Internal Variables ***************************
!
! Comments for variables used in this package can be found in the file
!   IVACOM.
!
! *********************** Type Declarations ****************************
!
      integer KDIM, MAXORD, MAXSTF
!++ Substitute for KDIM, MAXORD, MAXSTF below
      parameter (KDIM = 20, MAXORD = 2, MAXSTF = 1)
!--D Next line special: P=>D, X=>Q
      double precision TN
      double precision XI(KDIM)
!
!--D Next line special: P=>D, X=>Q
      double precision TG(2), TGSTOP(2), TMARK, TMARKX, TOUT, TOLG
      double precision ALPHA(KDIM), BETA(KDIM+1)
      double precision  D(MAXSTF+MAXORD,MAXORD), G(KDIM,MAXORD)
      double precision V(KDIM+MAXORD)
      double precision HC, HDEC, HINC, HINCC, HMAX, HMAXP9, HMIN
      double precision FDAT(11)
!
      double precision DS(MAXSTF+MAXORD, MAXORD), GS(KDIM)
      double precision SIGMA(KDIM), RBQ(KDIM), DNOISE
      double precision EAVE, EIMAX, EIMIN, EMAX, EREP, ROBND, SNOISE
!
!.    SPECIFICATION OF ENVIRONMENTAL CONSTANTS.
      double precision EEPS10, EEPS16, EROV10, EEPS2
      double precision EEPT75, EOVEP2, OVTM75, OVD10
      common / DIVAEV / EEPS2, EEPT75, EOVEP2, OVTM75, OVD10, EEPS10,   &
     &   EEPS16, EROV10
      save / DIVAEV /
      integer IOPST, KORDI, KQMAXD, KQMAXI, LDT, MAXDIF, MAXINT, NKDKO, &
     &   NTE, NYNY, NDTF, NUMDT
      common / DIVASC / TN, XI, IOPST, KORDI, KQMAXD, KQMAXI, LDT,      &
     &   MAXDIF, MAXINT, NKDKO, NTE, NYNY, NDTF, NUMDT
!
      integer ICF,ICS,IGFLG,IGTYPE(2),IGSTOP(2),ILGREP,INGS,IOP3,IOP4,  &
     &   IOP5,IOP6,IOP7,IOP8,IOP9,IOP10,IOP11,IOP12,IOP13,IOP14,IOP15,  &
     &   IOP16,IOP17,IOP18,IOP19,IOP20,IOP21,IOP22,IOP21S,ITOLEP,IY,    &
     &   KEMAX,KIS,KMARK,KORD1I,KORD2I,KPRED,KQDCON,KQICON,KQMAXS,      &
     &   KQMXDS,KQMXIL,KQMXIP,KQMXIS,KSC,KSOUT,KSSTRT,KSTEP,LEX,LINC,   &
     &   LINCD,LINCQ,LSC,MAXKQD,MAXKQI,METHOD,NE,NEPTOL,NG,NGTOT,       &
     &   NOISEQ,NOUTKO,NTOLF,NY,IDAT(6)
      common /DIVAMC/ TG,TGSTOP,TMARK,TMARKX,TOUT,TOLG,HC,HDEC,HINC,    &
     &   HINCC,HMAX,HMAXP9,HMIN,ALPHA,BETA,D,G,V,DS,GS,SIGMA,RBQ,DNOISE,&
     &   EAVE,EIMAX,EIMIN,EMAX,EREP,ROBND,SNOISE,FDAT,ICF,ICS,IGFLG,    &
     &   IGTYPE,IGSTOP,ILGREP,INGS,IOP3,IOP4,IOP5,IOP6,IOP7,IOP8,IOP9,  &
     &   IOP10,IOP11,IOP12,IOP13,IOP14,IOP15,IOP16,IOP17,IOP18,IOP19,   &
     &   IOP20,IOP21,IOP22,IOP21S,ITOLEP,IY,KEMAX,KIS,KMARK,KORD1I,     &
     &   KORD2I,KPRED,KQDCON,KQICON,KQMAXS,KQMXDS,KQMXIL,KQMXIP,KQMXIS, &
     &   KSC,KSOUT,KSSTRT,KSTEP,LEX,LINC,LINCD,LINCQ,LSC,MAXKQD,MAXKQI, &
     &   METHOD,NE,NEPTOL,NG,NGTOT,NOISEQ,NOUTKO,NTOLF,NY,IDAT
      save / DIVAMC / , / DIVASC /
!
      integer KGO, INTCHK(0:30), NXTCHK
      integer IHI, JL, J, ILOW, K, KQQ, JLIM, KFERR
!
      equivalence (INTCHK(1), NXTCHK)
      double precision CM1
      parameter (CM1 = (-1.D0))
      integer IOPIVA(2)
      save IOPIVA
!
!                      Declarations for error message processing.
!
      character TEXT1(1)*10
      integer  MENTXT,MEIDAT,MEMDA1,MECONT,MERET,MEEMES,METEXT,MEIVEC
      parameter (MENTXT =23)
      parameter (MEIDAT =24)
      parameter (MEMDA1 =27)
      parameter (MECONT =50)
      parameter (MERET  =51)
      parameter (MEEMES =52)
      parameter (METEXT =53)
      parameter (MEIVEC =57)
!
      integer MACT(16), MLOC(12), MACT1(4)
!
! ********* Error message text ***************
![Last 2 letters of Param. name]  [Text generating message.]
!AA DIVA$B
!AB The interval [1, 10**6], bounds the allowed values for NTE=$I.$E
!AC For option $I, the interval [$I, $I], bounds the allowed $C
!   values for the integration order which is set to $I.$E
!AD Option 16 must be used for error control.$E
!AE F($I) = $F, but it must be -1.0 when skipping the error check.$E
!AF For option $I, the interval [$I, $I] bounds the allowed $C
!   values for KORD($I)=$I, which is used to specify an $B
!AG output type for printing.$E
!AH output group for printing.$E
!AI equation group for variational equations.$E
!AJ order for a differential equation.$E
!AK equation group for diagnostic print.$E
!AL equation group for integration order control.$E
!AM equation group for error control.$E
!AN Option 5 argument must be <= 0 or > 4.$E
!   $
!AO KORD values for this option starting at KORD($M) are:$E
      integer LTXTAA,LTXTAB,LTXTAC,LTXTAD,LTXTAE,LTXTAF,LTXTAG,LTXTAH,  &
     & LTXTAI,LTXTAJ,LTXTAK,LTXTAL,LTXTAM,LTXTAN,LTXTAO
      parameter (LTXTAA=  1,LTXTAB=  7,LTXTAC= 71,LTXTAD=183,LTXTAE=226,&
     & LTXTAF=290,LTXTAG=400,LTXTAH=427,LTXTAI=455,LTXTAJ=498,          &
     & LTXTAK=535,LTXTAL=573,LTXTAM=620,LTXTAN=655,LTXTAO=  1)
      character MTXTAA(3) * (233)
      character MTXTAB(1) * (55)
      data MTXTAA/'DIVA$BThe interval [1, 10**6], bounds the allowed val&
     &ues for NTE=$I.$EFor option $I, the interval [$I, $I], bounds the$&
     & allowed values for the integration order which is set to $I.$EOpt&
     &ion 16 must be used for error control.$EF($I) = ','$F, but it must&
     & be -1.0 when skipping the error check.$EFor option $I, the interv&
     &al [$I, $I] bounds the allowed values for KORD($I)=$I, which is us&
     &ed to specify an $Boutput type for printing.$Eoutput group for pri&
     &nting.$Eequation gro','up for variational equations.$Eorder for a$&
     & differential equation.$Eequation group for diagnostic print.$Eequ&
     &ation group for integration order control.$Eequation group for err&
     &or control.$EOption 5 argument must be <= 0 or > 4.$E'/
      data MTXTAB/'KORD values for this option starting at KORD($M) are:&
     &$E'/

! End of automatically generated error message code.
!
!        for KGO =     1      2      3      4      5      6      7
      data MLOC / LTXTAI,LTXTAK,LTXTAL,LTXTAM,LTXTAJ,LTXTAG,LTXTAH,     &
     &            LTXTAB,LTXTAC,LTXTAD,LTXTAE, LTXTAN /
!           KGO        8      9     10     11      12
!
!                      1  2  3 4       5  6       7       8       9 10
      data MACT / MEEMES,38,24,0, MENTXT, 0, METEXT, MECONT, MEMDA1,0,  &
     &  METEXT, MEIDAT,0, MEIVEC,0, MECONT /
!           11      12 13     14 15     16
      data MACT1 / METEXT, MEIVEC, 0, MERET /
      data IOPIVA(1) / 1111 /
      data TEXT1 / 'IOPT()= $B' /
!
! ************** START OF EXECUTABLE CODE ******************
!
!     **** TEST IF CONTINUING AN INTEGRATION
      if (KORD(1) /= 0) go to 330
!     **** INITIALIZE VARIOUS SCALARS
      KSTEP = 0
      KQMXIS = 0
      KORD2I = -5
      KORD(2) = -1
      NTE = NEQ
      NE = NTE
      TOLG = 0.D0
!     **** SET UP OPTIONS
      if (IOPT(1) /= 0) call DIVAOP(IOPT, F)
      call DIVAOP(IOPIVA, F)
      if (IOPT(1) == 0) IOPIVA(2) = 1
!
      if ((NE <= 0) .or. (NE > 1000000)) then
         IDAT(1) = NE
         KGO = 8
         go to 650
      end if
!                         Set up diagnostic print on storage allocation.
      INTCHK(0) = 245
      if (IOP10 /= 0) INTCHK(0) = 247
!
!     **** CHECK TSPECS STORAGE ALLOCATION
      INTCHK(2) = IDIMT
      INTCHK(3) = 4
      NXTCHK = 4
      if (IOP5 /= 0) then
         INTCHK(4) = 5
         INTCHK(5) = 5
         if (IOP5 > 0) then
            INTCHK(6) = IOP5 - 4
            if (IOP5 < 5) then
               KGO = 12
               go to 600
            end if
         else
            IHI = -IOP5
            JL = 4
            do 15 IHI = IHI, IDIMK-3, 3
               J = abs(KORD(IHI))
               if (J == 0) go to 20
               if (abs(KORD(IHI + 2)) > 1) then
                  IDAT(2) = -1
                  IDAT(3) = 1
                  KGO = 6
                  go to 600
               end if
               if ((J <= JL) .or. (J > KORD(IHI+1))) then
                  KGO = 7
                  IDAT(2) = JL + 1
                  IDAT(3) = KORD(IHI+1)
                  go to 610
               end if
               JL = KORD(IHI+1)
   15       continue
            if (KORD(IHI) /= 0) IHI = IHI + 3
   20       INTCHK(6) = JL - 4
         end if
         NXTCHK = 7
      end if
   25 call OPTCHK(INTCHK, IOPT, 'DIVA / TSPECS$E')
      if (NXTCHK < 0) KORD2I = -4
!
!     **** CHECK KORD STORAGE ALLOCATION
      INTCHK(2) = IDIMK
      INTCHK(3) = NE + 3
      NXTCHK = 4
      if (IOP5 < 0) then
         INTCHK(4) = 5
         INTCHK(5) = -IOP5
         INTCHK(6) = IHI + IOP5
         NXTCHK = 7
      end if
!
!++  Code for VAREQ is active
      if (IOP18 /= 0) then
         NE = abs(KORD(IOP18))
         INTCHK(NXTCHK) = 18
         ILOW = IOP18
         KGO = 1
!.       **** CHECK OPTION FOR VALID INPUT
         go to 430
      end if
!++  End
   30 continue
      if (NKDKO /= 0) then
!                        **** STORAGE ALLOCATED FOR ODE ORDERS
         INTCHK(NXTCHK) = 17
         INTCHK(NXTCHK+1) = NKDKO
         INTCHK(NXTCHK+2) = NTE
         NXTCHK = NXTCHK + 3
      end if
!++  Code for STIFF is inactive
!      IF (IOPST /= 0) then
!         INTCHK(NXTCHK) = 17
!         INTCHK(NXTCHK+1) = IOPST
!         INTCHK(NXTCHK+2) = NTE
!         NXTCHK = NXTCHK + 3
!      end if
!++  End
!
! **** SET INITIAL INTEGRATION ORDERS, TEST ODE ORDERS ****
!
      MAXINT = 0
      MAXDIF = 0
      NY = 0
      do 80 K = 1, NTE
         if (NKDKO /= 0) KORDI = KORD(NKDKO + K - 1)
         NY = NY + abs(KORDI)
!++  Code for STIFF is inactive
!      IF (IOPST == 0) GO TO 60
!c.    **** CHECK FOR POTENTIAL STIFF EQUATION
!      JS = abs(KORD(IOPST+K-1)) - 1
!      IF ( JS ) 52,60,54
!c.    **** EQUATION IS NOT ACTIVE
!   52 KQQ = 0
!      GO TO 56
!c.    **** EQUATION USES IMPLICIT METHOD
!   54 KQQ = -1
!      IF (JS > abs(KORDI)) then
!        Set up an error message.
!      end if
!      MAXINT = max(MAXINT, abs(KORDI) - JS)
!   56 IF (KORDI >= 0) GO TO 70
!      KORDI = -1 - KORDI
!      JS = JS - 1
!      MAXDIF = max(MAXDIF, JS, 1)
!      GO TO 70
!++  End
!     **** EQUATION IS TO USE AN EXPLICIT METHOD
   60    KQQ = 1
         MAXINT = max(MAXINT, KORDI)
   70    if ((KORDI > MAXORD) .or. (KORDI <= 0)) then
!                    Set up error message.  KORDI is out of range.
            IDAT(1) = 17
            IDAT(2) = 1
            IDAT(3) = MAXORD
            if (NKDKO /= 0) then
               KGO = 5
               ILOW = NKDKO
               IHI = NKDKO + K - 1
               go to 640
            else
               KGO = 9
               IDAT(4) = KORDI
               go to 650
            end if
         end if
         KORD(K + 3) = KQQ
   80 continue
!     **** SET FLAGS WHICH DEPEND ON METHOD USED
!++  Code for STIFF is inactive
!      METHOD = 1
!      IF (MAXINT > 0) IF (MAXDIF) 85,90,85
!      METHOD = -1
!   85 CONTINUE
!      KPRED = 5
!      GO TO 100
!++  End
   90 METHOD = 0
      KPRED = 1
  100 continue
!
! ******* CHECK KORD FOR DIAGNOSTIC OUTPUT CONTROL *********
!
!++  Code for OUTPUT is active
      if (IOP10 > 0) then
         if (NOUTKO /= 0) then
            INTCHK(NXTCHK) = 10
            ILOW = NOUTKO
!.    **** Check option for valid input
            KGO = 2
            go to 430
         end if
      end if
!++  End
  110 continue
!
! ********** CHECK KORD FOR INTEGRATION ORDER CONTROL ******
!
!++  Code for INTEGO is active
      if (IOP19 /= 0) then
!.           **** Check option for valid input
         INTCHK(NXTCHK) = 19
         ILOW = IOP19
         JLIM = -30
         KGO = 3
         go to 430
      end if
!++  End
  120 continue
!
! ********** CHECK SET UP FOR ERROR TOLERANCES *************
!
      INTCHK(NXTCHK) = 16
      ILOW = IOP16
      JLIM = -5
      KGO = 4
      if (IOP16 /= 0) go to 430
!.                      **** IN CURRENT CODE, IOP16=0 IS AN ERROR
      KGO = 10
      go to 650
  150 continue
!     **** CHECK KORD STORAGE ALLOCATION
      call OPTCHK(INTCHK, IOPT, 'DIVA / KORD$E')
      if (NXTCHK < 0) KORD2I = -4
!
!     ******** DONE CHECKING KORD STORAGE ALLOCATION *******
!
!     **** CHECK  Y  STORAGE ALLOCATION
      INTCHK(2) = IDIMY
      INTCHK(3) = NY + NY
      NXTCHK = 4
      NYNY = NY + 1
      call OPTCHK(INTCHK, IOPT, 'DIVA / Y$E')
      if (NXTCHK < 0) KORD2I = -4
!
!     **** CHECK  F  STORAGE ALLOCATION
      INTCHK(2) = IDIMF
      INTCHK(3) = NTE
      NXTCHK = 4
      if (IOP16 /= 0) then
!                                Error tolerance info.
         INTCHK(4) = 16
         INTCHK(5) = NTOLF
         INTCHK(6) = IHI - IOP16 + 1
         NXTCHK = 7
      end if
      if (IOP12 > 0) then
         INTCHK(NXTCHK) = 12
         INTCHK(NXTCHK+1) = IOP12
         INTCHK(NXTCHK+2) = 4
         NXTCHK = NXTCHK + 3
      end if
      if (IOP21 > 0) then
         INTCHK(NXTCHK) = 21
         INTCHK(NXTCHK+1) = IOP21
         INTCHK(NXTCHK+2) = 1
         NXTCHK = NXTCHK + 3
      end if
!
!++  Code for ERRSTO is inactive
!      IF (IOP20 /= 0) then
!c.                                Space for saving error estimates
!         INTCHK(NXTCHK) = 20
!         INTCHK(NXTCHK) = IOP20
!         INTCHK(NXTCHK) = NTE
!         NXTCHK = NXTCHK + 3
!      end if
!++  Code for STIFF is inactive
!      if (IOP21 > 0) then
!c.                               Info. for stiff equations
!         INTCHK(NXTCHK) = 21
!         INTCHK(NXTCHK+1) = IOP21
!         INTCHK(NXTCHK+2) = IOP21S
!         NXTCHK = NXTCHK + 3
!      end if
!      MAXKQD = min(MAXKQI, 6)
!++  End
!                          Set aside space for the difference tables.
      INTCHK(NXTCHK) = 0
      INTCHK(NXTCHK+1) = -KDIM * NTE
      INTCHK(NXTCHK+2) = 0
      NXTCHK = NXTCHK + 3
      INTCHK(NXTCHK) = -5 * NTE
      call OPTCHK(INTCHK, IOPT, 'DIVA / F$E')
      if (NXTCHK < 0) then
         KORD2I = -4
      else if (KORD2I /= -4) then
         do 290 K = NXTCHK+1, INTCHK(NXTCHK)
            if (INTCHK(INTCHK(K)) == 0) then
               NDTF = INTCHK(INTCHK(K)+1)
               NUMDT = min(KDIM, (INTCHK(INTCHK(K)+2)-NDTF+1) / NTE)
               MAXKQI = NUMDT - 1
            else
!         Take a quick return if needed space was not specified by user.
               KORD2I = -4
            end if
  290    continue
      end if
      if (IOP9 + abs(IOP10) + IOP11 /= 0) then
! Insure user doesn't get in trouble with F not iniitalized.
         do 300 K = NDTF, NDTF + NTE*NUMDT - 1
            F(K) = 0.D0
  300    continue
      end if
  320 continue
      if ((KORD2I == -4) .or. (IOP10 /= 0)) then
         MACT1(3) = IOPIVA(2)
         call MESS(MACT1, TEXT1, IOPT)
         KORD1I = 24
         KORD(1) = 24
      end if
      TMARK = TSPECS(1)
      TMARKX = TSPECS(4) + TSPECS(2)
!
!     **** DONE WITH INITIALIZATION AND CHECKING INPUTS
      if (IOP13 + IOP14 + IOP15 /= 0) return
  330 call DIVAA(TSPECS, Y, F, KORD, DIVAF, DIVAO)
      return
!
! ************ LOOP TO CHECK OPTION SPECIFICATIONS *********
!
  430 JL = 0
      do 560 IHI = ILOW, IDIMK
         J = KORD(IHI)
         go to (460, 480, 490, 490), KGO
!     **** CHECK ON VARIATIONAL EQUATIONS
  460    continue
!++  Code for VAREQ is active
         if (J - NTE) 470, 565, 620
  470    if (J == 0) go to 560
         if (J <= JL) go to 620
!++  End
!     **** Check on diagnostic output option
  480    continue
!++  Code for OUTPUT is active
!.    **** CHECK IF DONE
         if (J >= NTE) go to 565
         if (J <= JL) go to 620
         go to 550
!++  End
  490    continue
!     **** Check integration order control (KGO=3) and
!     **** error tolerance equation grouping (KGO=4).
         if (J - NTE) 500, 565, 620
  500    if (J) 510, 530, 540
  510    if ((JL <= 0) .and. (IHI /= ILOW)) go to 620
         if (J < JLIM) then
!                         Output an error message.
            IDAT(2) = JLIM
            IDAT(3) = 0
            go to 630
         end if
  520    JL = -JL
         go to 560
  530    if (KGO == 3) go to 520
         KFERR = NTOLF + IHI - ILOW
         if (F(KFERR) == CM1) go to 510
!                         Set up error message, TOL must be -1.
            IDAT(1) = KFERR
            KGO = 11
            go to 650
  540    if (abs(JL) >= abs(J)) go to 620
  550    JL = J
  560    continue
  565 NXTCHK = NXTCHK + 3
      INTCHK(NXTCHK-2) = ILOW
      INTCHK(NXTCHK-1) = IHI - ILOW + 1
      go to (30, 110, 120, 150), KGO
!
!     **** AN ERROR HAS BEEN MADE
!                  Error in setting up TSPECS for extra output
  600 IHI = IHI + 2
  610 ILOW = -IOP5
      go to 630
!                  Error in KORD indices
  620 IDAT(2) = abs(JL) + 1
      IDAT(3) = NTE
!                  Set up for print of message about KORD
  630 IDAT(1) = INTCHK(NXTCHK)
  640 IDAT(4) = IHI
      IDAT(5) = KORD(IHI)
!
! ***************** Process Errors *************************************
!
  650 KORD2I = -4
      MACT(4) = LTXTAF
      if (KGO >= 8) MACT(4) = -1
      MACT(6) = MLOC(KGO)
!--D Next line special: P=>S, X=>D
      CALL DMESS(MACT, MTXTAA, IDAT, FDAT)
      if (KGO < 8) then
         MACT(10) = ILOW
         MACT(13) = ILOW
         MACT(15) = -min(IHI+2, IDIMK)
         CALL MESS(MACT(9), MTXTAB, KORD)
         if (KGO <= 4) go to 565
      end if
!              5   6   7    8    9   10   11  12
      go to (100, 25, 25, 320, 100, 150, 660, 25), KGO - 4
  660 KGO = 4
      go to 565
    end subroutine diva
!*************************************************************************

!*************************************************************************
!>
!  Main subroutine for variable order integration of ordinary
!  differential equations
!
!### History
!  * 1989-02-24 DIVAA  Krogh   Big error with BETA(2)=1+epsilon -- looped
!  * 1988-07-27 DIVAA  Krogh   Fixed to allow restart on a restart.
!  * 1988-03-07 DIVAA  Krogh   Initial code.

    subroutine DIVAA(TSPECS, Y, F, KORD, DIVAF, DIVAO)

      integer KORD(*)
!--D Next line special: P=>D, X=>Q
      double precision TSPECS(*), Y(*)
      double precision F(*)
      external DIVAF, DIVAO
!
      integer KDIM, MAXORD, MAXSTF
!++ Substitute for KDIM, MAXORD, MAXSTF below
      parameter (KDIM = 20, MAXORD = 2, MAXSTF = 1)
!--D Next line special: P=>D, X=>Q
      double precision TN
      double precision XI(KDIM)
!
!--D Next line special: P=>D, X=>Q
      double precision TG(2), TGSTOP(2), TMARK, TMARKX, TOUT, TOLG
      double precision ALPHA(KDIM), BETA(KDIM+1)
      double precision  D(MAXSTF+MAXORD,MAXORD), G(KDIM,MAXORD)
      double precision V(KDIM+MAXORD)
      double precision HC, HDEC, HINC, HINCC, HMAX, HMAXP9, HMIN
      double precision FDAT(11)
!
      double precision DS(MAXSTF+MAXORD, MAXORD), GS(KDIM)
      double precision SIGMA(KDIM), RBQ(KDIM), DNOISE
      double precision EAVE, EIMAX, EIMIN, EMAX, EREP, ROBND, SNOISE
!
!.    SPECIFICATION OF ENVIRONMENTAL CONSTANTS.
      double precision EEPS10, EEPS16, EROV10, EEPS2
      double precision EEPT75, EOVEP2, OVTM75, OVD10
      common / DIVAEV / EEPS2, EEPT75, EOVEP2, OVTM75, OVD10, EEPS10,   &
     &   EEPS16, EROV10
      save / DIVAEV /
      integer IOPST, KORDI, KQMAXD, KQMAXI, LDT, MAXDIF, MAXINT, NKDKO, &
     &   NTE, NYNY, NDTF, NUMDT
      common / DIVASC / TN, XI, IOPST, KORDI, KQMAXD, KQMAXI, LDT,      &
     &   MAXDIF, MAXINT, NKDKO, NTE, NYNY, NDTF, NUMDT
!
      integer ICF,ICS,IGFLG,IGTYPE(2),IGSTOP(2),ILGREP,INGS,IOP3,IOP4,  &
     &   IOP5,IOP6,IOP7,IOP8,IOP9,IOP10,IOP11,IOP12,IOP13,IOP14,IOP15,  &
     &   IOP16,IOP17,IOP18,IOP19,IOP20,IOP21,IOP22,IOP21S,ITOLEP,IY,    &
     &   KEMAX,KIS,KMARK,KORD1I,KORD2I,KPRED,KQDCON,KQICON,KQMAXS,      &
     &   KQMXDS,KQMXIL,KQMXIP,KQMXIS,KSC,KSOUT,KSSTRT,KSTEP,LEX,LINC,   &
     &   LINCD,LINCQ,LSC,MAXKQD,MAXKQI,METHOD,NE,NEPTOL,NG,NGTOT,       &
     &   NOISEQ,NOUTKO,NTOLF,NY,IDAT(6)
      common /DIVAMC/ TG,TGSTOP,TMARK,TMARKX,TOUT,TOLG,HC,HDEC,HINC,    &
     &   HINCC,HMAX,HMAXP9,HMIN,ALPHA,BETA,D,G,V,DS,GS,SIGMA,RBQ,DNOISE,&
     &   EAVE,EIMAX,EIMIN,EMAX,EREP,ROBND,SNOISE,FDAT,ICF,ICS,IGFLG,    &
     &   IGTYPE,IGSTOP,ILGREP,INGS,IOP3,IOP4,IOP5,IOP6,IOP7,IOP8,IOP9,  &
     &   IOP10,IOP11,IOP12,IOP13,IOP14,IOP15,IOP16,IOP17,IOP18,IOP19,   &
     &   IOP20,IOP21,IOP22,IOP21S,ITOLEP,IY,KEMAX,KIS,KMARK,KORD1I,     &
     &   KORD2I,KPRED,KQDCON,KQICON,KQMAXS,KQMXDS,KQMXIL,KQMXIP,KQMXIS, &
     &   KSC,KSOUT,KSSTRT,KSTEP,LEX,LINC,LINCD,LINCQ,LSC,MAXKQD,MAXKQI, &
     &   METHOD,NE,NEPTOL,NG,NGTOT,NOISEQ,NOUTKO,NTOLF,NY,IDAT
      save / DIVAMC / , / DIVASC /
!
      double precision CMP75, C0, C1M5, C1M3, C2P5M3, C8M3, CP0625, CP1
      double precision CP25, CP3, CP4, CP5, CP875, C1, C1P125, C1P3, C2
      double precision C6, C10, C16, C4096
      parameter (CMP75 = -.75D0)
      parameter (C0 = 0.D0)
      parameter (C1M5 = 1.D-5)
      parameter (C1M3 = 1.D-3)
      parameter (C2P5M3 = 2.5D-3)
      parameter (C8M3 = 8.D-3)
      parameter (CP0625 = .0625D0)
      parameter (CP1 = .1D0)
      parameter (CP25 = .25D0)
      parameter (CP3 = .3D0)
      parameter (CP4 = .4D0)
      parameter (CP5 = .5D0)
      parameter (CP875 = .875D0)
      parameter (C1 = 1.D0)
      parameter (C1P125 = 1.125D0)
      parameter (C1P3 = 1.3D0)
      parameter (C2 = 2.D0)
      parameter (C6 = 6.D0)
      parameter (C10 = 10.D0)
      parameter (C16 = 16.D0)
      parameter (C4096 = 4096.D0)
!
      integer LDIS, KEXIT, I, J, K, J1, J2, L, LX
      double precision TP, TP1, TP2, TP3, HH, DISADJ
      double precision SIGMAS, TPS1, TPS2, EIMINO, EXR
!--D Next line special: P=>D, X=>Q
      double precision  TMARKA(2), XP, XP1
      equivalence (G(1, 1), HH), (TMARKA(1), TMARK)
      equivalence (KEXIT, IOP17)
      save EIMINO, EXR, SIGMAS, TP, TP1, TP2, TP3, TPS1, TPS2, DISADJ,  &
     &   LDIS
!
!                      Declarations for error message processing.
!
      integer MACT(17), MLOC(8), MENTXT, MERET, MEEMES, METEXT
      parameter (MENTXT =23)
      parameter (MERET  =51)
      parameter (MEEMES =52)
      parameter (METEXT =53)
!
! ********* Error message text ***************
![Last 2 letters of Param. name]  [Text generating message.]
!AA DIVAA$B
!AB At: TN=$F, KSTEP=$I, with H=$F$E
!AC A previously reported error was fatal.$E
!AD Print points not properly ordered: TSPEC($I)=$F$E
!AE An error tolerance of 0 requires setting special flags.  $B
!AF Step size reduced too fast, doing a restart.  $B
!AG H is so small that TN + H = TN.  $B
!AH Error tolerance too small.  $B
!AI Step size at end of start < HMIN=$F, $B
!AJ Error estimates require a stepsize < HMIN=$F, $B
!AK (Estimated Error) / (Requested Error) for equation $I is $F.  $B
!AL Tolerance $I is F($I) = $F.$B
!AM Tolerance $I is F($I) * F($I) = $F * $F = $F.$B
!AN   Replacing F($I) with $F.$E
      integer LTXTAA,LTXTAB,LTXTAC,LTXTAD,LTXTAE,LTXTAF,LTXTAG,LTXTAH,  &
     & LTXTAI,LTXTAJ,LTXTAK,LTXTAL,LTXTAM,LTXTAN
      parameter (LTXTAA=  1,LTXTAB=  8,LTXTAC= 40,LTXTAD= 80,LTXTAE=129,&
     & LTXTAF=189,LTXTAG=237,LTXTAH=272,LTXTAI=302,LTXTAJ=342,          &
     & LTXTAK=390,LTXTAL=455,LTXTAM=484,LTXTAN=531)
      character MTXTAA(3) * (186)
!
      integer LOCM, MLOCAC, MLOCAD, MLOCAE, MLOCAF, MLOCAG, MLOCAH,     &
     &   MLOCAI, MLOCAJ

      parameter (LOCM = 32 * 256)
!                     KORD1I   Severity   Loc. message
      parameter (MLOCAC = 23 + 32 * (99 + 256 * LTXTAC))
      parameter (MLOCAD = 13 + 32 * (38 + 256 * LTXTAD))
      parameter (MLOCAE = 22 + 32 * (38 + 256 * LTXTAE))
      parameter (MLOCAF =  3 + 32 * (14 + 256 * LTXTAF))
      parameter (MLOCAG = 21 + 32 * (38 + 256 * LTXTAG))
      parameter (MLOCAH =  2 + 32 * (25 + 256 * LTXTAH))
      parameter (MLOCAI = 11 + 32 * (38 + 256 * LTXTAI))
      parameter (MLOCAJ = 12 + 32 * (38 + 256 * LTXTAJ))
!
      data MTXTAA/'DIVAA$BAt: TN=$F, KSTEP=$I, with H=$F$EA previously r&
     &eported error was fatal.$EPrint points not properly ordered: TSPEC&
     &($I)=$F$EAn error tolerance of 0 requires setting special flags. $&
     & ','$BStep size reduced too fast, doing a restart.  $BH is so smal&
     &l that TN + H = TN.  $BError tolerance too small.  $BStep size at$&
     & end of start < HMIN=$F, $BError estimates require a steps','ize <&
     & HMIN=$F, $B(Estimated Error) / (Requested Error) for equation $I$&
     & is $F.  $BTolerance $I is F($I) = $F.$BTolerance $I is F($I) * F(&
     &$I) = $F * $F = $F.$B  Replacing F($I) with $F.$E'/
      data MLOC / MLOCAC, MLOCAD, MLOCAE, MLOCAF, MLOCAG, MLOCAH,       &
     &   MLOCAI, MLOCAJ /
!
!                      1 2 3 4       5 6       7       8 9      10
      data MACT / MEEMES,0,0,0, MENTXT,0, METEXT, MENTXT,0, METEXT,     &
     &  MENTXT, 0, METEXT, MENTXT, LTXTAN, METEXT, MERET /
!           11 12      13      14      15      16     17
!
      data EXR, EIMINO / CP1, C8M3 /
      data LDIS / 0 /
!
! ************** START OF EXECUTABLE CODE ******************************
!
  660 if (KORD2I) 670, 1380, 1840
  670 if (KORD2I == -1) go to 2140
      if (KORD2I == -5) go to 720
      if (KORD2I+8) 1840, 710, 1840
!     **** SPECIAL OUTPUT CASE (EXTRAPOLATION OR GSTOP)
  680 go to (1220, 1190, 830, 690, 690, 2190), KEXIT
!     **** RESET TMARK BEFORE GOING WHERE DIRECTED BY KEXIT
  690 KEXIT = KEXIT - 2
      go to 1890
!     **** USER HAS REQUESTED A RESTART
  700 KORD2I = -5
      IGFLG = 0
      if ((KORD(2) >= 2) .and. (LSC < 3)) KORD2I = -8
      if ((KORD1I <= 3) .or. (KORD1I == 5)) go to 1890
      if (KORD2I == -5) go to 720
!                Set up for a discontinuity
  710 XP = TSPECS(1)
      if (KORD(2) == 3) then
!                Adjust for discontinuity in Y
         J = 0
         do 712 I = 1, NTE
            if (NKDKO /= 0) KORDI = KORD(NKDKO + I -1)
            K = 1
            J = J + KORDI
            XP1 = Y(J)
  711       Y(NYNY + J - K) = Y(NYNY + J - K) + XP1
            if (K < KORDI) then
               XP1 = Y(J-K) + (TN - XP) * XP1 / dble(K)
               K = K + 1
               go to 711
            end if
  712    continue
      else
         IGFLG = -3
      end if
      DISADJ = HH
      XP1 = (XP - TN) / XI(1)
      if (XP1 < CP25) then
         K = 1
         if (XP1 < CMP75) K = 2
         TSPECS(1) = TN - XI(K)
         TSPECS(2) = 2.D0*(TN - TSPECS(1))
         call DIVAIN(TSPECS(1), Y, F, KORD)
         do 713 J = 1, NY
            Y(NYNY + J - 1) = Y(J)
  713    continue
!          Move difference tables back one step
         TN = TSPECS(1)
  714    call DIVABU(F, KORD)
         if (K == 2) then
            KSC = max(KSC - 1, 1)
            do 715 K = max(1, KSC), IOP11-1
               BETA(K+1) = BETA(K) * (XI(K) / (XI(K+1) - XI(1)))
  715       continue
            K = 0
            go to 714
         end if
      end if
!      Take step to discontinuity.
      HH = (XP - TN)
      EREP = -abs(EREP)
      KIS = 1000
      LDIS = 1
      LSC = 0
      LINC = 0
      HINCC = C1P125
      KSSTRT = KSTEP + 2
      go to 1115
! ********
! INITIALIZE FOR STARTING AN INTEGRATION
! ********
  720 HH = TSPECS(2)
      LSC = 8
      LINC = -5
      LINCD = 64
      LDT = -4
      KQMAXI = 0
      EIMIN = CP3
      EAVE = C10
      XI(1) = C0
      IGFLG = 0
      KIS = 0
      KSSTRT = 0
      ROBND = 0.D0


! GO COMPUTE INITIAL DERIVATIVES
  730 KORD2I = -6
!++  Code for VAREQ is active
      ICS = 0
      ICF = NE
!++  End
  740 KORD1I = KPRED
      go to 1360
!   RETURN AFTER COMPUTING INITIAL (OR NOISE TEST) DERIVATIVES
  750 continue
!++  Code for VAREQ is active
      if (IOP18 == 0) go to 790
!.    **** SPECIAL LOGIC TO COMPUTE VARIATIONAL DERIVATIVES
  760 if (KORD1I == 3) go to 770
      KORD1I = 3
      KORD(3) = 0
  770 if (ICF == NTE) go to 790
      if (KORD2I == -6) then
         if (ICS == 1) go to 790
      end if
      ICS = ICF + 1
  780 KORD(3) = KORD(3) + 1
      ICF = IOP18 + KORD(3)
      ICF = abs(KORD(ICF))
      if (ICF) 1360, 780, 1360
!++  End
  790 ICS = 1
      ICF = NE
!++  Code for VAREQ is active
      if (KORD2I == 0) if (EREP) 2220, 2220, 1430
!++  End
      if (LINC + 5) 1490, 810, 1490
! END OF SPECIAL CODE FOR INITIALIZATION AT THE START
! ********
! UPDATE VARIABLES TO PREPARE FOR NEXT STEP
! ********
  800 LDT = 0
      EIMIN = C2P5M3+EIMIN*(C6*EAVE+EIMAX)/((EIMIN+C1)*(C6*EIMAX+EAVE))
  810 do 820 J = 1, NY
  820    Y(NYNY + J - 1) = Y(J)
      TN = TSPECS(1)
! ********
! TEST FOR VARIOUS TYPES OF OUTPUT
! ********
!     TEST FOR END OF STEP OUTPUT (OR IF DUMP OUTPUT TO BE TESTED FOR)
      if (IOP6 == 0) go to 840
!   SET UP FOR END OF STEP OUTPUT
      KORD1I = 6
      go to 1810
!     SET UP AFTER OTHER COMPUTATIONS HAVE BEEN MADE (TSPECS(1)/=TN)
  830 KORD1I = 6
      go to 860
!     TEST FOR AN OUTPUT POINT
  840 if (HH * (TMARK - TN)) 1760, 1760, 850
!     TEST FOR TOO MANY STEPS OUTPUT
  850 continue
      if (KSOUT > KSTEP) go to 890
      KORD1I = 4
  860 if (TSPECS(1) == TN) go to 1810
!     GO INTERPOLATE VALUES AT END OF LAST STEP
  870 TSPECS(1) = TN
      go to 1780
!     CONTINUE AFTER TOO MANY STEPS OUTPUT
  880 continue
      KSOUT = KSTEP + IOP4
  890 continue
!++  Code for DUMP is active
      if (IOP9 == 0) go to 920
!++  Code for DUMP & STIFF is inactive
!      KQMXDS=KQMAXD
!++  Code for DUMP is active
      KQMXIS = KQMAXI
      KIS = KIS + 1
      if (KIS /= 0) go to 920
!.   TIME TO DUMP THE SOLUTION
  900 KORD1I = 9
!.    SET TO UPDATE THE DIFFERENCE TABLE
      if (LDT == 1) go to 1810
! Note that doing this update can lead to very small differences in the
! results because of round off differences.
      LDT = -2
      go to 1780
!++  End
! RETURN AFTER DUMPING THE SOLUTION
  910 continue
!++  Code for DUMP is active
      KIS = 2
!.    TEST IF SOLUTION DUMP DUE TO RESTART, END, OR
!.    DROP IN INTEG. ORDER
      if (LINC < 0) if (LINC + 6) 1860, 1750, 1180
!++  End
! END OF TESTING FOR VARIOUS TYPES OF OUTPUT
!   TEST IF STEPSIZE MAY BE INCREASED OR IF TESTS SHOULD BE MADE
!   FOR DECREASING THE STEPSIZE (OR IF STARTING OR USER SELECTING H)
  920 KSTEP = KSTEP + 1
      if (LINC) 930, 980, 940
  930 continue
!     **** ONLY POSSIBLE VALUES AT THIS POINT ARE
!          LINC = -2 OR -5
      if (LINC + 5) 1120, 1110, 1120
! ********
! ERROR ESTIMATES INDICATE STEPSIZE CAN BE INCREASED
! ********
  940 HC = HINCC
      if (LINC > 1) HC = HC ** LINC
      HH = HC * HH
!     TEST IF NEW STEPSIZE IS TOO BIG
      if (abs(HH) > HMAX) if (HMAX) 970, 970, 960
  950 EAVE = EIMAX
      ROBND = CP3 + HINCC
      go to 1110
!     NEW STEPSIZE IS TOO BIG
  960 if (abs(XI(1)) >= HMAXP9) go to 970
      HH = sign(HMAX, HH)
      go to 950
!     RESTORE THE OLD STEPSIZE
  970 HH = XI(1)
      LINC = 0
      go to 1150
! END OF CODE FOR CASE WHEN ERROR ESTIMATES INDICATE STEPSIZE INCREASE
! ********
! TEST IF ESTIMATED ERRORS INDICATE STEPSIZE SHOULD BE DECREASED
! ********
  980 ROBND = C1P3
      if (EIMAX <= EAVE) go to 990
      EAVE = EAVE + CP4 * (EIMAX - EAVE)
      if ((EIMAX * EMAX) - C1M3) 1000, 1010, 1010
  990 EAVE = EIMAX
      if ((EIMAX * EMAX) >= EIMIN) go to 1010
 1000 ROBND = CP3 + (SIGMA(KQMAXS)/SIGMA(KQMAXS-1))
      go to 1180
!     TEST IF STEPSIZE SHOULD BE REDUCED
 1010 if (EMAX * EIMAX < EXR * EAVE) go to 1180
! ********
! ERROR ESTIMATES INDICATE STEPSIZE SHOULD BE REDUCED
! ********
      HC = HDEC
      if (EIMIN <= EIMINO) go to 1030
      EIMIN = EIMINO
      HC = CP875
 1030 HH = HC * XI(1)
      if (LSC - 1) 1040, 1080, 1090
 1040 if (abs(HH) >= HMIN) go to 1090
      if (abs(CP875*XI(1)) <= HMIN) if (LINC) 1050, 970, 970
      HH = sign(HMIN, HH)
      go to 1090
!     STEPSIZE IS TOO SMALL TO BE REDUCED
!     SET UP ERROR INDICATORS AND PREPARE FOR RETURN TO USER
 1050 KORD1I = 8
      go to 2240
!     PROCEED WITH CURRENT STEPSIZE DESPITE ERROR BEING TOO BIG
 1070 HH = XI(1)
      TSPECS(2) = HH
      EMAX = C0
      LINC = 0
      if (KORD1I - 2) 1420, 1430, 1430
!     SET LSC TO END STARTING PHASE
 1080 LSC = 2
!     CHECK IF REPEATING A STEP
 1090 if (LINC /= -1) go to 1110
!   WHEN REPEATING A STEP, BACK UP THE DIFFERENCES AND STEPSIZE INFO.
 1100 call DIVABU(F, KORD)
!   TEST IF NOISE TEST (LINC = -7) OR IF H IS NOT
!     BEING CHANGED (LINC = -4)
      if (LINC + 4) 1780, 1180, 1110
! ********
! STEPSIZE IS TO BE CHANGED
! ********
 1110 continue
! MODIFY STEPSIZE TO REDUCE ROUNDOFF ERROR IN ACCUMULATING INDEP. VAR.
      TP = C2 * abs(TN) + C4096 * abs(HH)
      TP = (TP + abs(HH)) - TP
      if (TP /= C0) HH = sign(TP, HH)
!     TEST IF NEW STEPSIZE SELECTED ACTUALLY GIVES A CHANGE
      if (HH == TSPECS(2)) go to 1140
 1115 TSPECS(2) = HH
      if (IOP8 == 0) go to 1140
!     SETUP TO TELL USER ABOUT STEPSIZE CHANGE (OR TO CHANGE STEPSIZE)
 1120 KORD1I = 8
      go to 860
!     RETURN AFTER TELLING USER ABOUT STEPSIZE CHANGE
 1130 HH = TSPECS(2)
 1140 if (HH /= XI(1)) KQICON = -1
 1150 HC = min(EOVEP2, abs(HH)) / EEPS2
! ********
! PREPARE FOR BEGINNING A NEW STEP
! ********
      if (LINC > 0) then
         LINC = min(LINC, LINCQ) + LINCQ
         go to 1190
      end if
 1180 LINC = LINCD
 1190 if (HC > abs(TN)) go to 1200
!     **** GIVE SINGULARITY DIAGNOSTIC
      KORD1I = 5
      go to 2240
 1200 TSPECS(1) = TN + HH
      if (LEX == 0) go to 1250
      if (HH * (TSPECS(1) - TMARKX) < C0) go to 1250
      TSPECS(1) = TMARKX
      HH = TMARKX - TN
      LINC = 64
      if (LEX > 0) go to 1240
      if ((LSC < 4) .and. (HH / XI(1) < CP3)) go to 1230
 1220 HH = CP875 * HH
      go to 1110
!     **** GIVE OUTPUT AT CURRENT TMARK (WITH EXTRAPOLATION)
 1230 KORD1I = -KMARK
      go to 1770
!     **** INTEGRATE TO TMARKX
 1240 KQICON = -1
!   TEST IF SUBROUTINE FOR COMPUTING INTEGRATION COEFF. SHOULD BE CALLED
 1250 continue
!++  Code for STIFF is inactive
!      IF ((KQMAXI <KQICON) .OR. (KQMAXD<KQDCON)) GO TO 1320
!++  Code for ~STIFF is active
      if (KQMAXI < KQICON) go to 1320
!++  End
!   GO COMPUTE COEFFICIENTS REQUIRED FOR THE INTEGRATION
!     TEST IF STARTING
      if (LSC < 7) go to 1310
 1260 KQMAXI = 2
!++  Code for STIFF is inactive
!      IF (METHOD) 1262,1270,1264
! 1262 KQMAXI=0
! 1264 KQMAXD=max(MAXDIF,2)
!      CALL DIVAHC
!c.  SET UP TO GO DO INITIALIZATION FOR CASE OF STIFF EQUATIONS
!      KORD1I=5
!      GO TO 1350
!++  End
!   INITIALIZE FOR EQUATIONS WHICH ARE NOT STIFF
 1270 KQMAXD = 0
      call DIVAHC
      J = NDTF
      do 1300 I = 1, NTE
!++  Code for STIFF is inactive
!         if (KORD(I + 3) <= 0) go to 1290
!++  End
         KORD(I + 3) = 1
!     INITIALIZE THE DIFFERENCE TABLE
         if (LDT == -4) F(J) = F(I)
         F(J + 1) = C0
         F(J + 2) = C0
 1290    continue
         J = J + NUMDT
 1300    continue
      if (LSC == 5) go to 1340
      LSC = 7
      LDT = 1
      go to 1330
!   INTEGRATION IS NOT BEING STARTED
 1310 K = KORD(KEMAX + 3)
      SIGMAS = SIGMA(K)
      call DIVAHC
!     **** ADJUST EAVE
      TPS1 = BETA(K)
      if (TPS1 > C1) TPS1 = CP5 * TPS1 + CP5
      EAVE = EAVE * TPS1 * (SIGMA(K) / SIGMAS)
!     TEST BELOW USED TO GET SAME RESULTS WITH/WITHOUT EXTRA EQUATIONS
      if (K > KQICON) LSC = max(LSC, -3)
! END OF SPECIAL LOGIC FOR CASE WHEN INTEG. COEFF. ROUTINE IS CALLED
 1320 continue
! ********
! PREDICT Y
! ********
 1330 continue
!++  Code for ~ARGM is active
      call DIVAPR(Y, Y(NYNY), F, KORD)
!++  Code for ARGM is inactive
!      CALL DIVAPE
!++  End
!     GO GET PREDICTED DERIVATIVES
 1340 KORD1I = KPRED
! ********
! CALL DIVAF  (OR RETURN)
! ********
 1350 KORD2I = 0
 1360 KORD(1) = KORD1I
      KORD(2) = 0
      if (IOP13 /= 0) return
      call DIVAF(TSPECS(1), Y, F, KORD(1))
!     TEST FOR SPECIAL USER RETURN
 1380 if (KORD(1) < 0) go to 2130
!     TEST FOR SPECIAL CASE
      if (KORD2I /= 0) go to 660
! ********
! TRANSFER CONTROL TO PROPER PLACE AFTER COMPUTING DERIVATIVES
! ********
      if (KORD1I - 2) 1400, 800, 1390
 1390 continue
!++  Code for VAREQ is active
      if (ICS - ICF) 1410, 1410, 760
!++  End
! ********
! PREPARE FOR CORRECTING, AND CORRECT Y
! ********
 1400 ITOLEP = 0
      ILGREP = 0
      IY = 1
      EIMAX = C1M5
      EMAX = C0
      KQMAXI = 2
      KQMAXS = 2
!++  Code for STIFF is inactive
!      IF (METHOD) 1404,1410,1406
! 1404 KQMAXI=0
! 1406 KQMAXD=2
!++  End
 1410 continue
!++  Code for ~ARGM is active
      call DIVACR(Y, F, KORD, F(NTOLF), KORD(IOP16))
!++  Code for ARGM is inactive
!      CALL DIVACE
!++  End
!     TEST IF ESTIMATED ERROR IS TOO BIG (OR IF DIAGNOSTIC CALLED FOR)
      if (EMAX > EREP) if (EREP) 2210, 2210, 1670
 1420 continue
!++  Code for VAREQ is active
      if (IOP18 /= 0) go to 760
!++  End
 1430 KORD1I = 2
!     TEST IF NOISE APPEARS TO LIMIT PRECISION
      if (EMAX < C0) go to 1470
!++  Code for ~STIFF is active
      if (LSC) 1450, 1360, 1610
!++  Code for STIFF is inactive
!      IF (LSC) 1450,1460,1610
!++  End
!     SET LSC=0 IF NOISE NO LONGER APPEARS TO LIMIT PRECISION
!     OR IF THE END OF THE STARTING PHASE HAS BEEN REACHED
 1450 LSC = 0
 1460 if (METHOD) 800, 1350, 1350
! ********
! NOISE APPEARS TO BE LIMITING PRECISION
! ********
 1470 continue
      if (LSC <= 0) LSC = max(LSC - 1, -KQMAXS)
      if (LSC == -1) go to 1460
      if (abs(EMAX) < EXR) go to 1590
      LINC = -7
      TPS2 = (C1 + BETA(NOISEQ - 1)) ** NOISEQ
      if (SNOISE < EEPS10 * TPS2) go to 1550
      TP = sign(EEPT75 * abs(TN) + OVTM75, HH)
      if (abs(TP) > abs(HH)) go to 1550
      TSPECS(1) = TN + TP
      KORD1I = 0
!     **** GO TO BACK UP THE DIFFERENCES AND GET F(TSPECS(1))
      go to 1100
!     **** SOLUTION HAS BEEN INTERPOLATED AND F COMPUTED
 1490 continue
      KORD1I = 0
      LINC = LINC - 1
      if (LINC + 9) 1510, 1520, 1500
 1500 TSPECS(1) = TN + (TP + TP)
      TP1 = F(KEMAX)
      TP2 = F(NDTF + NUMDT * KEMAX - NUMDT)
      go to 1780
!     **** COMPUTE 2-ND DIFFERENCE AT CLOSE SPACED T VALUES
 1510 TP2 = TP3
 1520 TP3 = F(KEMAX)
      TPS1 = abs((TP3 - TP1) - (TP1 - TP2))
      if ((C16 * TPS1 * TPS2) >= DNOISE) if (LINC + 9) 1550, 870, 1550
 1530 continue
      TPS2 = CP25 * SNOISE / RBQ(NOISEQ)
      do 1540 K = 2, NUMDT
         TPS1 = TPS1 + TPS1
         RBQ(K) = max(TPS1, TPS2 * RBQ(K))
 1540    continue
      LINC = 0

!FTK Next two lines added 2009-10-15
      if (abs(EMAX) < EREP) go to 1460
!FTK  LINC = -1  And then on 2015-03-14 commented out this line

      HH = CP875 * HH
      go to 1040
!     **** SET UP TO GIVE NOISE DIAGNOSTIC
 1550 KORD1I = 6
      go to 2240
!     **** AFTER GIVING NOISE DIAGNOSTIC
 1560 KORD1I = 2
      if (KORD(2) >= 0) then
        TPS1 = EEPS10
!FTK Next line added 2009-10-15
        if (TPS1 < .49D0 * RBQ(2)) go to 1530
      end if
!     **** SET NEW VALUE FOR OFFENDING TOL
      F(NTOLF + ITOLEP - 1) = FDAT(7)
      if (LINC + 7) 1180, 1570, 1180
 1570 LINC = 0
 1580 if (LSC) 1460, 1460, 1610
!     **** CHANGE HINCC AND ADJUST SIGMA( )
 1590 if (LSC /= -4) go to 1580
      if (HINCC == C1P125) go to 1580
      TPS1 = C1P125 / HINCC
      TPS2 = 1.0D0
      do 1600 K = 2, IOP11
         TPS2 = TPS2 * TPS1
         SIGMA(K) = SIGMA(K) * TPS2
 1600    continue
      EAVE = EAVE * TPS1 ** (1-KORD(KEMAX+3))
      LINCD = 6
      LINCQ = 12
      HINCC = C1P125
      go to 1460
!   END OF CODE FOR CASE WHEN NOISE APPEARS TO LIMIT PRECISION
! ********
! SPECIAL LOGIC FOR STARTING THE INTEGRATION
! ********
 1610 if (LSC == 1) go to 800
      LSC = LSC - 1
      if (LSC - 2) 1620, 1640, 1650
 1620 if (EIMAX <= (CP0625*EAVE*(SIGMA(KQMAXS)/SIGMAS)*(BETA(KQMAXS+  &
     &   1))**2)) go to 800
 1630 KSSTRT = KSTEP + 2
!   TEST IF STEPSIZE IS TOO SMALL BEFORE ENDING STARTING PHASE
      if (abs(HH) >= HMIN) go to 1450
!     GIVE DIAGNOSTIC FOR STEPSIZE TOO SMALL AT END OF START
      KORD1I = 7
      go to 2240
!   SET LSC TO DO ONE DERIVATIVE EVAL. PER STEP
 1640 LSC = 1
      go to 800
!     TEST IF FIRST TIME THROUGH THE FIRST STEP
 1650 if (LSC == 6) go to 1340
!     END STARTING PHASE IF CONVERGENCE OF CORRECTOR ITERATES TOO SLOW
      if (LDT == -5) go to 1660
      LSC = min(LSC, 4)
      go to 800
 1660 LDT = 0
      if (LSC - 4) 1260, 1630, 1260
! END OF SPECIAL LOGIC FOR STARTING THE INTEGRATION
! ********
! ESTIMATED ERROR IS TOO BIG
! ********
 1670 if (BETA(2) - C1) 1690, 1730, 1680
 1680 HC = C1 / BETA(2)
      if (BETA(2) >= C1P125) go to 1740
 1690 if (BETA(2) > CP1) go to 1730
!   REQUIRED STEPSIZE REDUCTION IS TOO RAPID -- GIVE A DIAGNOSTIC
      KORD1I = 4
      go to 2240
!
!     TEST KORD(2) AFTER ABOVE DIAGNOSTIC OR A DISCONTINUITY DIAGNOSTIC
 1700 continue
      if (KORD(2) == 0) go to 1730
!  TEST IF SOLUTION MUST BE DUMPED BEFORE A RESTART
 1710 LINC = -1
!++  Code for DUMP is active
      if (IOP9 == 0) go to 1750
      if (KIS == 2) go to 1750
      LINC = -6
!.    GO DUMP SOLUTION BEFORE REPEATING THE STEP
 1720 KQMAXI = KQMXIS
!++  Code for DUMP & STIFF is inactive
!      KQMAXD=KQMXDS
!++  Code for DUMP is active
      call DIVABU(F, KORD)
      go to 900
!++  End
!   SET UP TO REPEAT THE STEP
 1730 HC = CP5
 1740 LINC = -1
      if (LSC <= 3) go to 1030
!   RESTART THE INTEGRATION IF ERROR IS TOO BIG ON FIRST OR SECOND STEP
! LOOP TO SELECT A NEW INITIAL STEPSIZE
 1750 LSC = 7
 1755 HH = HH * CP5
      EMAX = EMAX * CP25
      if (EMAX >= CP3) go to 1755
      go to 1090
!   END OF SELECTING A NEW INITIAL STEPSIZE
! END OF LOGIC FOR CASE WHEN ESTIMATED ERROR IS TOO BIG
! ********
! INTEGRATION HAS REACHED AN OUTPUT POINT
! ********
 1760 if (KMARK == 0) go to 1920
      KORD1I = min(KMARK, 5)
      KORD(3) = KMARK
      if (TSPECS(1) == TMARK) go to 1790
 1770 TSPECS(1) = TMARK
 1780 call DIVAIN(TSPECS(1), Y, F, KORD)
 1790 continue
      if (KORD1I) 1800, 730, 1810
!   OUTPUT POINT IS OBTAINED BY EXTRAPOLATION
 1800 continue
!++  Code for EXTRAP is active
      KORD1I = -KORD1I
      KORD2I = -7
      KEXIT = 4
!.  TEST IF GSTOP-S ARE PRESENT
!++  Code for EXTRAP &  GSTOP is active
      if (NGTOT == 0) go to 1820
      IGFLG = 4
      KEXIT = 2
      KORD1I = 7
      if (IOP7) 740, 1820, 740
!++  End
! ********
! CALL DIVAO  (OR RETURN)
! ********
 1810 KORD2I = 1
 1820 KORD(1) = KORD1I
      KORD(2) = 1
      if (IOP14 /= 0) return
      call DIVAO(TSPECS(1), Y, F, KORD(1))
!     TEST FOR SPECIAL USER RETURN OR FOR A RESTART
!++  Code for ~DUMP is inactive
! 1840 IF (KORD(1)) 2130,700,1880
!++  Code for DUMP is active
 1840 if (KORD(1) > 0) go to 1880
 1850 if (IOP9 == 0) go to 1870
!.    **** GO DUMP THE SOLUTION
      LINC = -7
      ITOLEP = KORD(1)
      IDAT(1) = KORD(2)
      NEPTOL = KORD1I
      if (LSC /= 8) go to 900
 1860 LINC = min(0, LINCD)
      KORD1I = NEPTOL
      KORD(1) = ITOLEP
      KORD(2) = IDAT(1)
 1870 if (KORD(1)) 2130, 700, 2100
!++  End
 1880 if (KORD2I < 0) go to (2140, 1810, 1350, 2110, 720, 750, 680,  &
     &   710), -KORD2I
      if (KORD2I == 0) go to 1380
! ********
! TRANSFER CONTROL TO PROPER PLACE AFTER OUTPUT
! ********
 1890 if (KORD1I - 5) 1910, 1930, 1900
 1900 if (KORD1I - 8) 840, 1130, 910
 1910 if (KORD1I - 3) 1920, 1930, 880
!   GET NEW TOUT
 1920 TOUT = TSPECS(1) + TSPECS(3)
! GET NEW TMARK (NEXT INDEP. VAR. OUTPUT POINT)
 1930 XP = TMARK
      K = KMARK
      TMARK = TOUT
      KMARK = 2
      LEX = 0
      if (IOP5) 1940, 1980, 1970
 1940 I = -IOP5
 1950 I = I + 3
      J1 = KORD(I - 3)
      if (J1) 1950, 1980, 1960
 1960 J2 = KORD(I - 2)
      L = KORD(I - 1)
      go to 1990
 1970 J1 = 5
      J2 = IOP5
      L = 0
      if (J2 >= J1) go to 1990
 1980 J1 = 4
      J2 = 4
      L = IOP3
!
!     **** LOOP TO SET NEW TMARK (AND TMARKX)
 1990 do 2060 J = J1, J2
!        **** TEST IF EXTRAPOLATION NOT POSSIBLE
         if (L == 0) go to 2010
         LX = 2
         if (LEX) 2020, 2030, 2020
 2000    LEX = L
 2010    LX = 1
 2020    if (HH * (TSPECS(J) - TMARKA(LX))) 2030, 2060, 2060
 2030    if (J == 4) go to 2050
         if (HH * (TSPECS(J) - XP)) 2060, 2040, 2050
 2040    if ((K >= J) .or. (K == 3)) go to 2060
 2050    TMARKA(LX) = TSPECS(J)
         if (LX == 2) go to 2000
         KMARK = J
 2060    continue
      if (IOP5 < 0) go to 1950
      if (J1 /= 4) go to 1980
      if (KMARK == 4) KMARK = 3
!     **** TEST IF NEW TMARK IS ACCEPTABLE
      if (HH * (XP - TMARK)) 2070, 2080, 2090
 2070 if (KORD2I - 1) 670, 840, 670
 2080 if (K /= KMARK) go to 2070
!++  Code for DUMP is active
      if (KORD1I == 3) go to 1850
!++  Code for ~DUMP is inactive
!      IF (KORD1I == 3) GO TO 2100
!++  End
 2090 if (KORD1I == 13) go to 2190
! SETUP TO INDICATE ERROR IN SPECIFICATION OF OUTPUT POINTS
      KORD1I = 2
      IDAT(2) = KMARK
      if (KMARK <= 3) IDAT(2) = KMARK + 1
      FDAT(3) = TSPECS(IDAT(2))
      go to 2240
!     SET KORD1I=1 TO INDICATE THAT END OF INTEGRATION HAS BEEN REACHED
 2100 KORD1I = 1
! ********
! RETURN TO USER
! ********
 2110 KORD2I = -1
      KORD(1) = KORD1I
 2130 KORD(2) = -1
      return
! ********
! TRANSFER CONTROL TO PROPER PLACE AFTER RETURN TO USER
! ********
 2140 if (KORD1I - 2) 2150, 1560, 2160
 2150 KORD2I = 1
      go to 1930
 2160 if (KORD1I - 4) 1700, 2200, 2170
 2170 if (KORD1I - 13) 2180, 1930, 2190
 2180 if (abs(HH) >= HMIN) if (KORD1I - 11) 1030, 1450, 1030
      if (KORD(2) == 0) if (KORD1I - 11) 1070, 800, 1070
!   ERROR MESSAGES HAVE BEEN IGNORED -- COMPUTATION CAN NOT CONTINUE
 2190 KORD1I = 1
      go to 2240
!
!        AFTER A DISCONTINUITY RETURN
 2200 LINC = -4
      if (KORD(2)) 1710, 1730, 1100
! ********
! PROBLEM ENCOUNTERED WHEN CORRECTING
! ********
 2210 if (LDIS == 0) go to 2230
!           Extra checks when had a user specified discontinuity.
!++  Code for VAREQ is active
      if (IOP18 /= 0) go to 760
!++  End
 2220 KORD1I = 2
      LDIS = LDIS + 1
      TP = DISADJ / HH
      if (KIS >= 1000) then
         if (LDIS == 2) then
            if (KQMAXS <= 3) then
               LDIS = 0
               EREP = abs(EREP)
               TSPECS(2) = HH*min(min(TP, TP**2),                       &
     &            (CP25 * EXR/EMAX)**.333333333D0)
               go to 720
            end if
            LINC = -5
            if (IOP9 == 0) KIS = 1001
            go to 800
         end if
         if (IOP9 == 0) KIS = KIS + 1
         if (KQMAXS <= LDIS + 2) KIS = LDIS + 1
         LINC = min(LINC, LDIS-2)
      end if
      if (LDIS > 2*KQMAXS) then
         EREP = abs(EREP)
         LDIS = 0
         if (EMAX > EREP) go to 1670
         go to 1430
      end if
      if (TP >= HINCC**(LINC+2)) then
         if ((LDIS /= 3) .and. (TP > dble(KQMAXS))) LSC = 1
         EIMIN = CP5
         EAVE = EAVE * TP**8
      end if
      if (LSC == 2) go to 1630
      if (EMAX > EXR) go to 1730
      go to 1430
!
 2230 EREP = abs(EREP)
!++  Code for DUMP is active
      if (LINC < -3) go to 1720
!++  End
!     BAD TOL
      KORD1I = 3
! ********
! ERROR PROCESSING
! ********
 2240 FDAT(1) = TN
      FDAT(2) = HH
      IDAT(1) = KSTEP
      ITOLEP = max(NEPTOL, -NEPTOL - 1)
      J = 3
      if (KORD1I >= 7) then
         J = 4
         FDAT(3) = HMIN
      end if
      if (KORD1I <= 3) then
         if (KORD1I < 3) then
            K = 8
         else
            MACT(9) = LTXTAL
            FDAT(3) = C0
            IDAT(2) = ITOLEP
            IDAT(3) = ITOLEP + NTOLF - 1
            K = 11
         end if
      else
         MACT(9) = LTXTAK
         FDAT(J) = EMAX
         IDAT(2) = KEMAX
         IDAT(3) = ITOLEP
         IDAT(4) = ITOLEP + NTOLF - 1
         FDAT(J+1) = F(IDAT(4))
         K = 14
         if (KORD1I == 6) then
            K = 17
            IDAT(5) = IDAT(4)
            FDAT(7) = 32.D0 * abs(EMAX) * FDAT(J+1)
            FDAT(J+2) = FDAT(7)
         end if
         MACT(12) = LTXTAL
         if (NEPTOL < 0) then
            MACT(12) = LTXTAM
            IDAT(6) = IDAT(4)
            IDAT(5) = IDAT(4) + 1
            FDAT(J+2) = F(IDAT(5))
            FDAT(J+3) = FDAT(J+1) * FDAT(J+2)
         end if
      end if
! Set the location for the first part of the message that varies, set
! the error severity, and the index number, print the error and
! return or stop.
      L = MLOC(KORD1I)
      MACT(6) = L / LOCM
      MACT(2) = (L - MACT(6) * LOCM) / 32
      KORD1I = mod(L, 32)
      MACT(3) = KORD1I
      MACT(K) = MERET
!--D Next line special: P=>S, X=>D
      call DMESS(MACT, MTXTAA, IDAT, FDAT)
      MACT(K) = MENTXT
      go to 2110

    end subroutine divaa
!*************************************************************************

!*************************************************************************
!>
! THIS SUBROUTINE RESTORES THE DIFFERENCE TABLE TO ITS STATE
! AT THE BEGINNING OF THE CURRENT STEP.  IF THE INTEGRATION ORDER
! WAS INCREASED, IT IS REDUCED. THE COMMON ARRAY XI IS ALSO
! RESTORED TO ITS STATE AT THE BEGINNING OF THE STEP. IF THE
! STEPSIZE IS NOT BEING CHANGED, THE ARRAY V USED TO COMPUTE
! INTEGRATION COEFFICIENTS IS RESTORED.
!
!### History
!  * 1987-12-07 DIVABU Krogh   Initial code.

    subroutine DIVABU(F, KORD)

      integer KORD(*)
      double precision F(*)
!
      integer KDIM, MAXORD, MAXSTF
!++ Substitute for KDIM, MAXORD, MAXSTF below
      parameter (KDIM = 20, MAXORD = 2, MAXSTF = 1)
!--D Next line special: P=>D, X=>Q
      double precision TN
      double precision XI(KDIM)
!
!--D Next line special: P=>D, X=>Q
      double precision TG(2), TGSTOP(2), TMARK, TMARKX, TOUT, TOLG
      double precision ALPHA(KDIM), BETA(KDIM+1)
      double precision  D(MAXSTF+MAXORD,MAXORD), G(KDIM,MAXORD)
      double precision V(KDIM+MAXORD)
      double precision HC, HDEC, HINC, HINCC, HMAX, HMAXP9, HMIN
      double precision FDAT(11)
!
      double precision DS(MAXSTF+MAXORD, MAXORD), GS(KDIM)
      double precision SIGMA(KDIM), RBQ(KDIM), DNOISE
      double precision EAVE, EIMAX, EIMIN, EMAX, EREP, ROBND, SNOISE
!
      integer IOPST, KORDI, KQMAXD, KQMAXI, LDT, MAXDIF, MAXINT, NKDKO, &
     &   NTE, NYNY, NDTF, NUMDT
      common / DIVASC / TN, XI, IOPST, KORDI, KQMAXD, KQMAXI, LDT,      &
     &   MAXDIF, MAXINT, NKDKO, NTE, NYNY, NDTF, NUMDT
!
      integer ICF,ICS,IGFLG,IGTYPE(2),IGSTOP(2),ILGREP,INGS,IOP3,IOP4,  &
     &   IOP5,IOP6,IOP7,IOP8,IOP9,IOP10,IOP11,IOP12,IOP13,IOP14,IOP15,  &
     &   IOP16,IOP17,IOP18,IOP19,IOP20,IOP21,IOP22,IOP21S,ITOLEP,IY,    &
     &   KEMAX,KIS,KMARK,KORD1I,KORD2I,KPRED,KQDCON,KQICON,KQMAXS,      &
     &   KQMXDS,KQMXIL,KQMXIP,KQMXIS,KSC,KSOUT,KSSTRT,KSTEP,LEX,LINC,   &
     &   LINCD,LINCQ,LSC,MAXKQD,MAXKQI,METHOD,NE,NEPTOL,NG,NGTOT,       &
     &   NOISEQ,NOUTKO,NTOLF,NY,IDAT(6)
      common /DIVAMC/ TG,TGSTOP,TMARK,TMARKX,TOUT,TOLG,HC,HDEC,HINC,    &
     &   HINCC,HMAX,HMAXP9,HMIN,ALPHA,BETA,D,G,V,DS,GS,SIGMA,RBQ,DNOISE,&
     &   EAVE,EIMAX,EIMIN,EMAX,EREP,ROBND,SNOISE,FDAT,ICF,ICS,IGFLG,    &
     &   IGTYPE,IGSTOP,ILGREP,INGS,IOP3,IOP4,IOP5,IOP6,IOP7,IOP8,IOP9,  &
     &   IOP10,IOP11,IOP12,IOP13,IOP14,IOP15,IOP16,IOP17,IOP18,IOP19,   &
     &   IOP20,IOP21,IOP22,IOP21S,ITOLEP,IY,KEMAX,KIS,KMARK,KORD1I,     &
     &   KORD2I,KPRED,KQDCON,KQICON,KQMAXS,KQMXDS,KQMXIL,KQMXIP,KQMXIS, &
     &   KSC,KSOUT,KSSTRT,KSTEP,LEX,LINC,LINCD,LINCQ,LSC,MAXKQD,MAXKQI, &
     &   METHOD,NE,NEPTOL,NG,NGTOT,NOISEQ,NOUTKO,NTOLF,NY,IDAT
      save / DIVAMC / , / DIVASC /
!
      integer I, L, KQQ, J, K
      double precision TPD, C0, C2
      parameter (C0 = 0.D0)
      parameter (C2 = 2.D0)
! ********* START OF EXECUTABLE CODE **********
!
! ********
! BEGIN LOOP TO BACK UP DIFFERENCE TABLES
! ********
      L = NDTF - 1
      do 2410 I = 1, NTE
         KQQ = KORD(I + 3)
!++  Code for STIFF is inactive
!         IF (KQQ) 2302,2400,2310
!c.           EQUATION IS STIFF
! 2302    IF (LINC>=0) GO TO 2310
!         IF (F(L+1+I)) 2306,2308,2304
!c.     ORDER WAS INCREASED, AND THUS MUST BE DECREASED (KQQ<0)
! 2304    KQQ=KQQ+1
!         KORD(I+3) = KQQ
!         GO TO 2308
!c.     ORDER WAS DECREASED
! 2306    KQQ=KQQ-1
! 2308    KQQ=max(2,-KQQ)
!         GO TO 2350
!++  End
!     EQUATION IS NOT STIFF
 2310    if (KQQ > 2) then
            if (F(L + KQQ) == C0) then
!                 ORDER WAS INCREASED, AND THUS MUST BE DECREASED
               KQQ = KQQ - 1
               KORD(I + 3) = KQQ
            end if
         end if
         J = min(KQQ, KSC)
         KQMAXI = max(KQMAXI, KQQ)
         if (KQQ /= 1) F(L + KQQ + 1) = 0.D0
!           BACK UP FOR BACKWARD DIFFERENCES
         do 2360 K = 1, J
            F(L + K) = F(L + K) - F(L + K + 1)
 2360    continue
         if (KQQ > KSC) then
!           BACK UP FOR MODIFIED DIVIDED DIFFERENCES
            do 2390 K = J+1, KQQ
               F(L + K) = (F(L+K) - F(L+K+1)) / BETA(K)
 2390       continue
         end if
 2400    F(L + KQQ + 1) = F(L + KQQ + 1) / BETA(KQQ + 1)
         L = L + NUMDT
 2410 continue
! END OF LOOP TO BACK UP DIFFERENCE TABLES
! ********
! BACK UP XI TO BEGINNING OF THE STEP
! ********
      I = KSC + 1
      if (I - IOP11 - 1) 2420, 2440, 2450
 2420 TPD = XI(1)
!                Check below needed when starting?
      if (TPD == XI(2)) go to 2450
      do 2430 K = I, IOP11
 2430    XI(K - 1) = XI(K) - TPD
 2440 XI(IOP11) = C2 * XI(IOP11 - 1)
      if (IOP11 /= 2) XI(IOP11) = XI(IOP11) - XI(IOP11 - 2)
 2450 KQICON = -1
      ICF = NE
      ICS = 1
      LDT = 1
      return
    end subroutine divabu
!*************************************************************************

!*************************************************************************
!>
! THIS SUBROUTINE RETURNS THE FOLLOWING DATA FROM COMMON
!
! ID(1) = KEMAX  =  INDEX OF EQUATION WITH LARGEST ERROR ESTIMATE
! ID(2) = KSTEP  =  CURRENT STEP NUMBER
! ID(3) = NUMDT  =  NUMBER OF DIFFERENCES USED FOR EACH EQUATION
! ID(4) =           RESERVED FOR FUTURE USE
! ID(5) =           RESERVED FOR FUTURE USE
! RD(1) = EMAX   =  MAX. RATIO OF ESTIMATED ERROR TO REQUESTED ERROR
! RD(2) =           RESERVED FOR FUTURE USE
! RD(3) =           RESERVED FOR FUTURE USE
!
!### History
!  * 1987-12-07 DIVACO Krogh   Initial code.

    subroutine DIVACO(ID, RD)

      integer ID(5)
      double precision RD(3)
!
      integer KDIM, MAXORD, MAXSTF
!++ Substitute for KDIM, MAXORD, MAXSTF below
      parameter (KDIM = 20, MAXORD = 2, MAXSTF = 1)
!--D Next line special: P=>D, X=>Q
      double precision TN
      double precision XI(KDIM)
!
!--D Next line special: P=>D, X=>Q
      double precision TG(2), TGSTOP(2), TMARK, TMARKX, TOUT, TOLG
      double precision ALPHA(KDIM), BETA(KDIM+1)
      double precision  D(MAXSTF+MAXORD,MAXORD), G(KDIM,MAXORD)
      double precision V(KDIM+MAXORD)
      double precision HC, HDEC, HINC, HINCC, HMAX, HMAXP9, HMIN
      double precision FDAT(11)
!
      double precision DS(MAXSTF+MAXORD, MAXORD), GS(KDIM)
      double precision SIGMA(KDIM), RBQ(KDIM), DNOISE
      double precision EAVE, EIMAX, EIMIN, EMAX, EREP, ROBND, SNOISE
!
      integer IOPST, KORDI, KQMAXD, KQMAXI, LDT, MAXDIF, MAXINT, NKDKO, &
     &   NTE, NYNY, NDTF, NUMDT
      common / DIVASC / TN, XI, IOPST, KORDI, KQMAXD, KQMAXI, LDT,      &
     &   MAXDIF, MAXINT, NKDKO, NTE, NYNY, NDTF, NUMDT
!
      integer ICF,ICS,IGFLG,IGTYPE(2),IGSTOP(2),ILGREP,INGS,IOP3,IOP4,  &
     &   IOP5,IOP6,IOP7,IOP8,IOP9,IOP10,IOP11,IOP12,IOP13,IOP14,IOP15,  &
     &   IOP16,IOP17,IOP18,IOP19,IOP20,IOP21,IOP22,IOP21S,ITOLEP,IY,    &
     &   KEMAX,KIS,KMARK,KORD1I,KORD2I,KPRED,KQDCON,KQICON,KQMAXS,      &
     &   KQMXDS,KQMXIL,KQMXIP,KQMXIS,KSC,KSOUT,KSSTRT,KSTEP,LEX,LINC,   &
     &   LINCD,LINCQ,LSC,MAXKQD,MAXKQI,METHOD,NE,NEPTOL,NG,NGTOT,       &
     &   NOISEQ,NOUTKO,NTOLF,NY,IDAT(6)
      common /DIVAMC/ TG,TGSTOP,TMARK,TMARKX,TOUT,TOLG,HC,HDEC,HINC,    &
     &   HINCC,HMAX,HMAXP9,HMIN,ALPHA,BETA,D,G,V,DS,GS,SIGMA,RBQ,DNOISE,&
     &   EAVE,EIMAX,EIMIN,EMAX,EREP,ROBND,SNOISE,FDAT,ICF,ICS,IGFLG,    &
     &   IGTYPE,IGSTOP,ILGREP,INGS,IOP3,IOP4,IOP5,IOP6,IOP7,IOP8,IOP9,  &
     &   IOP10,IOP11,IOP12,IOP13,IOP14,IOP15,IOP16,IOP17,IOP18,IOP19,   &
     &   IOP20,IOP21,IOP22,IOP21S,ITOLEP,IY,KEMAX,KIS,KMARK,KORD1I,     &
     &   KORD2I,KPRED,KQDCON,KQICON,KQMAXS,KQMXDS,KQMXIL,KQMXIP,KQMXIS, &
     &   KSC,KSOUT,KSSTRT,KSTEP,LEX,LINC,LINCD,LINCQ,LSC,MAXKQD,MAXKQI, &
     &   METHOD,NE,NEPTOL,NG,NGTOT,NOISEQ,NOUTKO,NTOLF,NY,IDAT
      save / DIVAMC / , / DIVASC /
!
      ID(1) = KEMAX
      ID(2) = KSTEP
      ID(3) = NUMDT
      RD(1) = EMAX
      return
    end subroutine DIVACO
!*************************************************************************

!*************************************************************************
!>
! THIS SUBROUTINE
!   1. CORRECTS Y FOR EQUATIONS WHICH ARE NOT STIFF
!   2. ESTIMATES ERRORS
!   3. SELECTS INTEGRATION ORDERS
!   4. TESTS IF NOISE LIMITS THE PRECISION
!
!     Y = VECTOR OF PREDICTED VALUES ON ENTRY, AND OF CORRECTED
!         VALUES WHEN THE RETURN IS MADE.
! LGROUP= VECTOR INDICATING HOW ERROR TOLERANCES ARE TO BE GROUPED
!         (AND POSSIBLY HOW INTEGRATION ORDERS ARE TO BE GROUPED).
!   TOL = VECTOR CONTAINING ERROR TOLERANCES (AND POSSIBLY RELATIVE
!         ERROR FACTORS).
!     F = VECTOR GIVING PREDICTED DERIVATIVE VALUES AND DIFF. TABLES.
!    KD = VECTOR GIVING ORDERS OF THE DIFFERENTIAL EQUATIONS
!         (IF EQUATIONS HAVE DIFFERENT ORDERS).
!    KQ = VECTOR OF INTEGRATION ORDERS.
!
!### History
!  * 1988-08-25 DIVACR Krogh   Fix bug in relative error test.
!  * 1988-01-15 DIVACR Krogh   Initial code.

    subroutine DIVACR(Y, F, KORD, TOL, LGROUP)

      integer LGROUP(*), KORD(*)
!--D Next line special: P=>D, X=>Q
      double precision Y(*)
      double precision TOL(*), F(*)
!
      integer KDIM, MAXORD, MAXSTF
!++ Substitute for KDIM, MAXORD, MAXSTF below
      parameter (KDIM = 20, MAXORD = 2, MAXSTF = 1)
!--D Next line special: P=>D, X=>Q
      double precision TN
      double precision XI(KDIM)
!
!--D Next line special: P=>D, X=>Q
      double precision TG(2), TGSTOP(2), TMARK, TMARKX, TOUT, TOLG
      double precision ALPHA(KDIM), BETA(KDIM+1)
      double precision  D(MAXSTF+MAXORD,MAXORD), G(KDIM,MAXORD)
      double precision V(KDIM+MAXORD)
      double precision HC, HDEC, HINC, HINCC, HMAX, HMAXP9, HMIN
      double precision FDAT(11)
!
      double precision DS(MAXSTF+MAXORD, MAXORD), GS(KDIM)
      double precision SIGMA(KDIM), RBQ(KDIM), DNOISE
      double precision EAVE, EIMAX, EIMIN, EMAX, EREP, ROBND, SNOISE
!
      integer IOPST, KORDI, KQMAXD, KQMAXI, LDT, MAXDIF, MAXINT, NKDKO, &
     &   NTE, NYNY, NDTF, NUMDT
      common / DIVASC / TN, XI, IOPST, KORDI, KQMAXD, KQMAXI, LDT,      &
     &   MAXDIF, MAXINT, NKDKO, NTE, NYNY, NDTF, NUMDT
!
      integer ICF,ICS,IGFLG,IGTYPE(2),IGSTOP(2),ILGREP,INGS,IOP3,IOP4,  &
     &   IOP5,IOP6,IOP7,IOP8,IOP9,IOP10,IOP11,IOP12,IOP13,IOP14,IOP15,  &
     &   IOP16,IOP17,IOP18,IOP19,IOP20,IOP21,IOP22,IOP21S,ITOLEP,IY,    &
     &   KEMAX,KIS,KMARK,KORD1I,KORD2I,KPRED,KQDCON,KQICON,KQMAXS,      &
     &   KQMXDS,KQMXIL,KQMXIP,KQMXIS,KSC,KSOUT,KSSTRT,KSTEP,LEX,LINC,   &
     &   LINCD,LINCQ,LSC,MAXKQD,MAXKQI,METHOD,NE,NEPTOL,NG,NGTOT,       &
     &   NOISEQ,NOUTKO,NTOLF,NY,IDAT(6)
      common /DIVAMC/ TG,TGSTOP,TMARK,TMARKX,TOUT,TOLG,HC,HDEC,HINC,    &
     &   HINCC,HMAX,HMAXP9,HMIN,ALPHA,BETA,D,G,V,DS,GS,SIGMA,RBQ,DNOISE,&
     &   EAVE,EIMAX,EIMIN,EMAX,EREP,ROBND,SNOISE,FDAT,ICF,ICS,IGFLG,    &
     &   IGTYPE,IGSTOP,ILGREP,INGS,IOP3,IOP4,IOP5,IOP6,IOP7,IOP8,IOP9,  &
     &   IOP10,IOP11,IOP12,IOP13,IOP14,IOP15,IOP16,IOP17,IOP18,IOP19,   &
     &   IOP20,IOP21,IOP22,IOP21S,ITOLEP,IY,KEMAX,KIS,KMARK,KORD1I,     &
     &   KORD2I,KPRED,KQDCON,KQICON,KQMAXS,KQMXDS,KQMXIL,KQMXIP,KQMXIS, &
     &   KSC,KSOUT,KSSTRT,KSTEP,LEX,LINC,LINCD,LINCQ,LSC,MAXKQD,MAXKQI, &
     &   METHOD,NE,NEPTOL,NG,NGTOT,NOISEQ,NOUTKO,NTOLF,NY,IDAT
      save / DIVAMC / , / DIVASC /
!.    SPECIFICATION OF ENVIRONMENTAL CONSTANTS.
      double precision EEPS10, EEPS16, EROV10, EEPS2, EEPT75, EOVEP2
      double precision OVTM75, OVD10
      common / DIVAEV / EEPS2, EEPT75, EOVEP2, OVTM75, OVD10, EEPS10,   &
     &   EEPS16, EROV10
      save / DIVAEV /
!
      integer L, I, KQL, KQN, KQD, JLGREP, J, K, ILGROR, ITOLOR, JLGROR,&
     &   IORD, KOUTKO, KQLORD, LL, LKQMAX
      double precision CM8, CM2, CMP5, C0, CQ3125, CP1, CP125, CP25, CP5
      double precision CP75, CP8, CP9375, C1, C1P4, C2, C4, C10, C20
      double precision C1000, C40
      parameter (CM8 = -8.D0)
      parameter (CM2 = -2.D0)
      parameter (CMP5 = -.5D0)
      parameter (C0 = 0.D0)
      parameter (CQ3125 = .03125D0)
      parameter (CP1 = .1D0)
      parameter (CP125 = .125D0)
      parameter (CP25 = .25D0)
      parameter (CP5 = .5D0)
      parameter (CP75 = .75D0)
      parameter (CP8 = .8D0)
      parameter (CP9375 = .9375D0)
      parameter (C1 = 1.D0)
      parameter (C1P4 = 1.4D0)
      parameter (C2 = 2.D0)
      parameter (C4 = 4.D0)
      parameter (C10 = 10.D0)
      parameter (C20 = 20.D0)
      parameter (C40 = 40.D0)
      parameter (C1000 = 1000.D0)
      double precision TPP, HH, E, EI, EPS, ERCOEF, RND, RNOISE, S
      double precision TP2, TPS1, TPS2, TPS3, TPS4, TPS5, TPS6, TPS7
      double precision REF(4)
      double precision EIBND(KDIM-1)
!++  Code for INTEGO is active
      double precision TEMPA(4), TEMPAO(4)
!++  End
      save KOUTKO, LKQMAX
      equivalence (TPS1,TEMPA(1)), (TPS2,TEMPA(2)), (TPS3,TEMPA(3)),    &
     &   (TPS4, TEMPA(4))
      equivalence (G(1, 1), HH)
      integer MACT1(2), MACT2(12)
!             Parameters for Interface to MESS and DMESS
      integer MERET, METEXT, METABL
      parameter (MERET  =51)
      parameter (METEXT =53)
      parameter (METABL =55)
! ********* Error message text ***************
![Last 2 letters of Param. name]  [Text generating message.]
!AA KSTEP=$(I6) T=$(E15.8) H=$(E12.5) LSC=$(I3) $C
!   EIMIN=$(E8.2) EAVE=$G KSC=$(I2) SIGMA($J)=$G $C
!   RQ=$(E11.5)$G$E
!   $
!AB I$HKQ$HLI$HE$HEI$HEPS$HF$H$H$H$H
!   HIGH ORDER PREDICTED DIFFERENCES$HRNOISE$HSTIFF$HBETA$E
      integer LTXTAA,LTXTAB
      parameter (LTXTAA=  1,LTXTAB=  1)
      character MTXTAA(1) * (104)
      character MTXTAB(1) * (88)
      data MTXTAA/'KSTEP=$(I6) T=$(E15.8) H=$(E12.5) LSC=$(I3) EIMIN=$(E&
     &8.2) EAVE=$G KSC=$(I2) SIGMA($J)=$G RQ=$(E11.5)$G$E'/
      data MTXTAB/'I$HKQ$HLI$HE$HEI$HEPS$HF$H$H$H$HHIGH ORDER PREDICTED$&
     & DIFFERENCES$HRNOISE$HSTIFF$HBETA$E'/
!
      data MACT1 / METEXT, MERET /
! (rr=repeat, t=3/5 for I/E format)  wwddtrr  wwddtrr  wwddtrr
      data MACT2 / METABL, 1, 0, 14, 0400201, 0300202, 0801503,         &
     &    1507501, 1103504, 0901501, 1002501, 1205501 /
!         wwddtrr  wwddtrr  wwddtrr  wwddtrr  wwddtrr
!          End of stuff for interface to message processor
!
      data REF(1), REF(2), REF(3), REF(4) / C1, CP9375, CP75, CP5 /
!++ Save data by elements if ~.C.
!++ Of next 20 lines, only the first KDIM-1 are active
      data EIBND(1) / .1D0 /
      data EIBND(2) / .1D0 /
      data EIBND(3) / .14D0 /
      data EIBND(4) / .19D0 /
      data EIBND(5) / .26D0 /
      data EIBND(6) / .36D0 /
      data EIBND(7) / .50D0 /
      data EIBND(8) / .69D0 /
      data EIBND(9) / .94D0 /
      data EIBND(10) / C1 /
      data EIBND(11) / C1 /
      data EIBND(12) / C1 /
      data EIBND(13) / C1 /
      data EIBND(14) / C1 /
      data EIBND(15) / C1 /
      data EIBND(16) / C1 /
      data EIBND(17) / C1 /
      data EIBND(18) / C1 /
      data EIBND(19) / C1 /
!     data EIBND(20) / C1 /
!
!++  Code for ARGM is inactive
!      RETURN
!      ENTRY DIVACE
!++  End
! ********
! START OF CODE
! ********
      L = NDTF - 1
      if (ICS /= 1) L = L + (ICS - 1) * NUMDT
      do 3340 I = ICS, ICF
         if (NKDKO /= 0) KORDI = KORD(NKDKO + I - 1)
         IY = IY + abs(KORDI)
         KQL = KORD(I + 3)
         KQN = abs(KQL)
         KQD = max(2, KQN)
! ********
! OBTAIN ERROR TOLERANCE SPECIFIED BY THE USER
! ********
         if (I <= ILGREP) if (KQL) 2600, 3310, 2610
         ITOLEP = abs(ITOLEP) + 1
         EPS = TOL(ITOLEP)
         ILGREP = LGROUP(ITOLEP)
!   TEST IF SIMPLE ABSOLUTE ERROR TEST IS BEING USED
         if (ILGREP > 0) go to 2580
         JLGREP = ILGREP
!     GET OLD RELATIVE ERROR FACTOR
         TPS6 = TOL(ITOLEP + 1)
         ILGREP = LGROUP(ITOLEP + 1)
         ITOLEP = -ITOLEP - 1
!
         if (JLGREP + 1) 2540, 2570, 2510
!   NO CHECK ON THE ERROR ESTIMATE IS TO BE MADE
 2510    if (EPS + C1) 2520, 2590, 2520
!   ERROR TOLERANCE IS SPECIFIED IMPROPERLY
 2520    KEMAX = I
         NEPTOL = ITOLEP
         LINC = -3
         EREP = -abs(EREP)
         return
!   COMPUTE NEW RELATIVE ERROR FACTOR
 2540    continue
         TPS1 = C0
         do 2550 J = I, ILGREP
            TPS1 = TPS1 + abs(F(J))
 2550       continue
         TPS1 = abs(HH) * TPS1 / dble(ILGREP - I + 1)
         if (LSC <= 2) go to 2560
!     ON FIRST 3 STEPS INCREASE TPS6 WHEN COMPUTING REL. ERROR FACTOR
         TPS6 = max(C4 * TPS1, TPS6)
!     ON 1-ST TIME THROUGH THE FIRST STEP, REL. ERR. FAC. IS NOT STORED
         if (LSC == 7) go to 2570
 2560    continue
         TPS6 = max(TPS1, TPS6)
!   STORE NEW RELATIVE ERROR FACTOR
         TOL(-ITOLEP) = TPS6 * REF(-JLGREP - 1)
!   COMPUTE ABSOLUTE ERROR TOLERANCE
 2570    EPS = EPS * TPS6
 2580    if (EPS <= C0) go to 2520
 2590    if (KQL) 2600, 3330, 2610
! END OF OBTAINING ERROR TOLERANCE
! ********
! OBTAIN INFORMATION USED FOR ERROR ESTIMATION, ORDER SELECTION, ETC.
! ********
! EQUATION IS STIFF
 2600    continue
!++  Code for STIFF is inactive
!      JS=abs(KORD(NJSKO+I-1))-1
!      JSI=JS
!      TPP=C0
!      TPS4=F(L+KQD+2)
!      TPS3=F(L+KQD+1)
!      TPS2=F(L+KQD)
!      TPS1=F(L+KQD-1)
!      IF (KQD==2) TPS1=Y(IY-1)
!      E=ABS(TPS3)+ABS(TPS4)
!      EI=E+ABS(TPS2)
!      RND=EI
!      IF (KORDI>=0) GO TO 2604
!c.    EQUATION IS IMPLICIT
!      JSI=JSI-1
!      IF (JSI/=0) GO TO 2604
!      IF (KORDI==-1) GO TO 2602
!      ERCOEF=GS(KQN+1)
!      GO TO 2606
! 2602 ERCOEF=.5D0*DS(KQD,1)
!      JSI=1
!      GO TO 2606
!c.    END OF SPECIAL CODE FOR IMPLICIT EQUATIONS
! 2604 ERCOEF = DS(KQD,JSI)
! 2606 ERCOEF = ABS(ERCOEF) / EPS
!      IF (LSC<=2)  GO TO 2710
!      IF (LSC-5) 2650,2710,2710
!c.  END OF CODE FOR STIFF EQUATIONS
!++  End
!
! EQUATION IS NOT STIFF
 2610    TPP = F(I) - F(L + 1)
         TPS3 = TPP
         TPS4 = TPP - F(L + KQD + 1)
         TPS2 = TPP + F(L + KQD)
         TPS1 = TPP + F(L + KQD - 1)
         E = abs(TPS3) + abs(TPS4)
         RND = E
         EI = E + abs(TPS2)
         ERCOEF = abs(GS(KQN + 1)) / EPS
         if (KQL >= 4) go to 2710
!   TEST IF STARTING OR IF INTEGRATION ORDER IS ONE
         if (LSC <= 2) if (KQL - 2) 2660, 2710, 2710
! ********
! LOGIC ASSOCIATED WITH STARTING THE INTEGRATION
! ********
         TPS4 = C0
         if (LSC - 4) 2650, 2640, 2620
! FIRST STEP
 2620    E = E * CQ3125
         TPS3 = C0
         F(L + 4) = C0
         S = C0
!   TEST IF FIRST TIME THROUGH THE FIRST STEP
         if (LSC == 7) go to 2690
!   COMPUTE S=ESTIMATE OF H * EIGENVALUE OF JACOBIAN = 2*(F(A)-F(B))/
!   (F(B)-F(C)) WHERE F(A)=CURRENT F(I), AND F(B) AND F(C) PRECEDING
!   VALUES OR ESTIMATES OF F(I)
         TPP = F(I) - F(L + 5)
         TPS4 = TPP
         E = C2 * abs(TPS4)
         if (S /= C0) S = (TPS4 + TPS4) / S
         if (S + CP125) 2630, 2700, 2700
!     SET LDT=-5  TO INDICATE POSSIBLE PROBLEMS DUE TO INSTABILITY
 2630    LDT = -5
         go to 2690
!   ADJUST CORRECTION MADE ON SECOND STEP
 2640    TPP = CP8 * TPP
!   ADJUST ESTIMATED ERRORS ON SECOND AND THIRD STEPS
 2650    E = abs(TPS3)
         RND = C4 * E
         go to 2710
! END OF SPECIAL LOGIC FOR STARTING
! ********
! INTEGRATION ORDER =1 IS TREATED AS A SPECIAL CASE
! ********
 2660    TPP = TPP + F(L + 2)
         if (BETA(2) >= C1P4) EI = EI * C1000
!   ESTIMATE NEW VALUE FOR S
         S = F(L + 4)
         if (S == C0) go to 2680
         S = max(CM8, C2 * BETA(2) * (TPS1 - TPS2 - F(L + 5)) / S)
         if (S >= CMP5) go to 2670
!   MODIFY TPP (TO GET BETTER STABILITY CHARACTERISTICS)
         TPP = TPP * max(CP25, (CM2 - C2 * S) / (S * S))
 2670    TPS4 = TPS4 * abs(S)
 2680    E = CP25 * (E + abs(TPS4))
         EI = EI + abs(TPS4 * S)
!     STORE INFORMATION REQUIRED TO ESTIMATE S ON NEXT STEP
 2690    F(L + 4) = TPP
 2700    F(L + 5) = F(I)
! END OF SPECIAL CODE FOR INTEGRATION ORDER =1
! ********
! CODE FOR NOISE TEST AND GETTING ERROR ESTIMATE
! ********
 2710    E = E * ERCOEF
         RNOISE = C0
         if (EPS < C0) go to 2810
         TPS5 = abs(F(L + 2)) + abs(F(I))
         if (TPS5 == C0) go to 2760
 2720    RNOISE = RND / TPS5
         if (RNOISE > RBQ(KQD)) if (RNOISE - C1) 2760, 2750, 2750
!   NOISE IS APPARENTLY SLOWING CONVERGENCE OF THE DIFFERENCES
!     REDUCE EI
         EI = RND
         TPS5 = abs(EEPS2 * Y(IY - 1)) / EPS
         if (TPS5 < abs(E)) if (LSC) 2730, 2730, 2760
         E = TPS5
         RNOISE = C0
 2730    E = -abs(E)
         if (EIMIN > CP1) EI = (C10 * EIMIN) * EI
!     COMPUTE REDUCTION TO BE MADE IN EI
         if (RNOISE > (C20 * RBQ(KQD))) go to 2760
         K = -6 - LSC
 2740    if (K <= 0) go to 2760
!     REDUCE EI WHEN NOISE APPARENTLY LIMITS PRECISION
         K = K - 1
         EI = CP5 * EI
         if (EI > EIMIN) go to 2740
         go to 2760
 2750    TPS4 = 1.1D0 * RND
         TPS3 = RND
 2760    continue
!   TEST FOR STIFFNESS GOES HERE WHEN IMPLEMENTED
! *       INGREDIENTS OF TEST MAY INCLUDE --
! *       RNOISE, WHETHER (ABS(TPS4)>ABS(TPS3)),
! *       WHETHER EMAX IS INCREASING, RESULT OF TEST ON
! *       PREVIOUS STEPS, ETC.
!
! ********
! COMPUTE ERROR ESTIMATES AND INFORMATION FOR SELECTING THE STEPSIZE
! ********
         if (E >= abs(EMAX)) go to 2770
         if (-E <= abs(EMAX)) go to 2780
         SNOISE = RNOISE
         DNOISE = RND
         NOISEQ = KQD
!   STORE PARAMETERS ASSOCIATED WITH LARGEST VALUE OF E
 2770    EMAX = E
         KEMAX = I
         NEPTOL = ITOLEP
!   DETERMINE HOW MUCH STEPSIZE CAN BE INCREASED
 2780    EI = EI * ERCOEF * SIGMA(KQD)
         EIMAX = max(EIMAX, EI)
         if (LINC <= 0) go to 2810
         K = 0
 2790    if (EI >= min(EIMIN, EIBND(KQN))) go to 2800
         K = K + 1
         if (K == LINC) go to 2810
         EI = EI * SIGMA(KQD)
         go to 2790
 2800    LINC = K
! END OF COMPUTING ERROR ESTIMATES
 2810    continue
!++  Code for ERRSTO is inactive
!      IF (IOP20 == 0) GO TO 780
!c.********
!c.STORE ERROR ESTIMATE (OPTIONAL)
!c.********
!      F(IOP20+I-1)=TPS3*GS(KQN+1)
!c.END OF STORING ERROR ESTIMATE
!++  Code for INTEGO | ERRSTO is active
         if (IOP19 == 0) go to 3090
!.********
!.EQUATIONS ARE GROUPED TO USE SAME INTEGRATION METHOD (OPTIONAL)
!.********
!++  Code for INTEGO is active
         if (I > 1) if (I - ILGROR) 2900, 2900, 2830
         ITOLOR = IOP19
 2830    JLGROR = KORD(ITOLOR)
         ITOLOR = ITOLOR + 1
         if (JLGROR > 0) go to 2870
         ILGROR = KORD(ITOLOR)
         ITOLOR = ITOLOR + 1
         if (JLGROR + 1) 2840, 2850, 2890
 2840    if (JLGROR < -2) if (KQD + JLGROR) 2850, 2880, 2880
!.INITIALIZE FOR ACCUMULATING VARIABLES USED IN ORDER SELECTION
 2850    IORD = I
         KQLORD = KQL
         do 2860 K = 1, 4
 2860       TEMPAO(K) = abs(TEMPA(K))
         go to 2930
!.ORDERS IN CURRENT GROUP CAN BE DIFFERENT
 2870    ILGROR = JLGROR
         go to 3090
!.ORDER IS NOT GOING TO BE CHANGED
 2880    JLGROR = 0
 2890    if (KQL) 3240, 3270, 3270
!.TAKE ACTION FOR EQUATION WHICH IS NOT THE FIRST IN THE GROUP
 2900    if (JLGROR) 2910, 2890, 3090
!.ACCUMULATE VARIABLES USED IN ORDER SELECTION
 2910    do 2920 K = 1, 4
 2920       TEMPAO(K) = TEMPAO(K) + abs(TEMPA(K))
!.    TEST IF THIS IS LAST EQUATION IN THE GROUP
 2930    if (I /= ILGROR) if (KQL) 3310, 3290, 3290
!.SET UP TO GO SELECT INTEGRATION ORDER
         KQL = 0
         do 2940 K = 1, 4
 2940       TEMPA(K) = TEMPAO(K)
         go to 3090
!.INTEGRATION ORDER HAS BEEN SELECTED
!++  Code for INTEGO | STIFF is active
 2950    continue
!++  Code for INTEGO is active
         KQL = KQLORD
         if (KQN - abs(KQL)) 2960, 2980, 3020
!.  TEST IF ORDER CAN BE DECREASED
 2960    if (JLGROR >= -2) if (KQL) 3010, 3040, 3040
!.    INTEGRATION ORDER WAS SELECTED OUTSIDE PERMITTED RANGE
 2970    KQN = abs(KQL)
!.    INTEGRATION ORDER IS NOT GOING TO BE CHANGED
 2980    if ((KQL /= 1) .or. (LSC > 0)) if (KQL) 3030, 3040, 3040
!.    SET  4-TH ENTRY IN DIFFERENCE TABLES SO THAT STANDARD ADAMS
!.    METHOD IS USED WHEN KQL=1
 2990    do 3000 K = IORD, I
 3000       F(NDTF + K*NUMDT - NUMDT + 3) = C0
         go to 3270
!.  ORDER FOR STIFF EQUATION WAS REDUCED
 3010    continue
!++  Code for INTEGO & STIFF is inactive
!      IF (KQN<JSI) GO TO 990
!      TPP=-C1
!      GO TO 1090
!c.  TEST IF ORDER CAN BE INCREASED
!++  Code for INTEGO is active
 3020    if (JLGROR == -2) go to 2970
!++  Code for INTEGO & STIFF is inactive
!      IF (KQL>=0) GO TO 1140
!      IF ((JSI/=0).AND.(KQN>(MAXKQD+JSI))) GO TO 990
!      TPP=C1
!c.  STORE RESULTS FOR STIFF EQUATIONS
!++  Code for INTEGO is active
 3030    continue
!++  Code for INTEGO & STIFF is inactive
!      DO 3035 K=IORD,I
!      KORD(K+3) = -KQN
! 3035 F(NDTF+K*NUMDT-NUMDT)=TPP
!      GO TO 3245
!c.  STORE RESULTS FOR EQUATIONS WHICH ARE NOT STIFF
!++  Code for INTEGO is active
 3040    LL = NDTF + NUMDT * IORD - NUMDT
         do 3080 J = IORD, I
            KORD(J + 3) = KQN
            if (KQN - KQL) 3050, 3070, 3060
 3050       F(LL + KQD - 1) = F(LL + KQD - 1) + (F(J) - F(LL))
            go to 3080
 3060       F(LL + KQN) = F(LL + KQD)
 3070       F(LL + KQD) = C0
 3080       LL = LL + NUMDT
         if (KQN - 1) 3270, 2990, 3270
!++  End
!.********
!.SELECT INTEGRATION ORDER
!.********
 3090    if (LSC <= 0) go to 3120
!. SPECIAL ORDER SELECTION WHEN STARTING
         if (LSC - 3) 3110, 3210, 3100
 3100    if (LSC == 5) if (S + .125D0) 3160, 3130, 3130
         if (LSC - 6) 3130, 3130, 3210
 3110    if (C40 * min(abs(TPS4), abs(TPS3)) > abs(TPS2)) then
            if (EPS /= -C1) LSC = 2
         end if
         if (abs(TPS4) < abs(TPS3)) if (C4 * abs(TPS4) - abs(TPS2))  &
     &      3130, 3130, 3210
!.  CHECK IF ORDER CAN BE INCREASED OR SHOULD BE DECREASED
 3120    TPS5 = ROBND * abs(TPS4)
         TPS6 = ROBND * (TPS5 + abs(TPS3))
         TPS7 = abs(TPS1) + abs(TPS2)
         if (TPS5 >= abs(TPS3)) go to 3140
         if (TPS6 >= TPS7) go to 3210
 3130    if (KQN >= MAXKQI) go to 3210
!.    INCREASE THE INTEGRATION ORDER
         KQN = KQN + 1
!++  Code for INTEGO | STIFF is active
         if (KQL) 3230, 2950, 3250
!++  Code for ~(INTEGO | STIFF) is inactive
!      GO TO 3250
!++  End
!.  CHECK IF ORDER SHOULD BE DECREASED
 3140    if (TPS6 < TPS7) go to 3210
         if (TPS5 < abs(TPS3 - TPS4)) go to 3210
         if ((TPS3==TPS4) .and. (LSC<=0)) go to 3210
         if (KQN - 2) 3210, 3160, 3180
 3160    KQN = 1
!++  Code for INTEGO | STIFF is active
         if (KQL) 3220, 2950, 3170
!++  End
!.    WHEN ORDER IS REDUCED TO 1 WITH ADAMS METHOD SET F(L+4)=0
 3170    F(L + 4) = C0
         go to 3260
!.    DECREASE THE INTEGRATION ORDER
 3180    KQN = KQN - 1
!++  Code for INTEGO | STIFF is active
         if (KQL) 3220, 2950, 3200
!++  End
 3200    F(L+KQD) = F(L+KQD) + TPP
         go to 3260
!   NO CHANGE IN INTEGRATION ORDER IS BEING MADE
 3210    continue
!++  Code for INTEGO is active
         if (KQL) 3240, 2950, 3270
!++  Code for ~INTEGO is inactive
!         TPS1 = EEPS10
!      GO TO 1530
!++  End
! END OF SELECTING INTEGRATION ORDER
! ********
! COMPUTE MAXIMUM INTEGRATION ORDERS AND SET NEW ONES (IF ANY)
! ********
! EQUATION IS STIFF
!     ORDER WAS DECREASED
!++  Code for INTEGO | STIFF is active
 3220    continue
!++  Code for STIFF is inactive
!      IF (KQN<JSI) GO TO 3236
!      F(L+1)=-C1
!      GO TO 3233
!c.    ORDER WAS INCREASED
!++  Code for INTEGO |  STIFF  is active
 3230    continue
!++  Code for STIFF is inactive
!      IF ((JSI/=0).AND.(KQN>(MAXKQD+JSI))) GO TO 3236
!      F(L+1)=C1
! 3233 KORD(I+3) = -KQN
!      GO TO 3245
!      ORDER WAS SET TO AN UNACCEPTABLE VALUE
! 3236 KQN=abs(KQL)
!      ORDER IS NOT BEING CHANGED
!++  Code for STIFF |  INTEGO is active
 3240    continue
!++  Code for STIFF is inactive
!      F(L+1)=C0
! 3245 IF (JSI/=0) KQMAXD=max(KQN,KQMAXD)
!      IF (JS<abs(KORDI)) KQMAXI=max(KQN,KQMAXI)
!      GO TO 3290
!++  End
! EQUATION IS NOT STIFF
!     ORDER INCREASED
 3250    F(L + KQN + 1) = -F(L + KQD + 1)
         if (LSC > 0) F(L + KQN + 1) = F(L + 1) - F(I)
!     ORDER CHANGED
 3260    KORD(I + 3) = KQN
 3270    KQMAXI = max(KQN, KQMAXI)
         if (EPS > C0) KQMAXS = max(KQN, KQMAXS)
         F(L + KQD + 1) = C0
 3290    continue
         if (KQN > KIS) go to 3310
!.********
!.DETERMINE IF TIME TO STORE SOLUTION (OPTIONAL)
!.********
         if (KIS >= 1000) then
            TP2 = max(1.5D0, dble(KQN) * C2 ** (1001 - KIS)) * abs(TPS4)
 3295       if (TP2 > abs(F(L+KQN))) then
               if (KQN <= KQL) then
                  KQN = KQN - 1
                  if (KQN > 1) go to 3295
                  KQN = 1
               end if
            end if
            KORD(I+3) = KQN
            if (I == 1) LKQMAX = 0
            LKQMAX = max(KQN, LKQMAX)
            KQMAXI = LKQMAX
            if (KIS == 1000) then
               if (I == KEMAX) EMAX = dble(8 + KQN**2) * abs(EMAX)
               go to 3325
            end if
!++  Code for DUMP is active
         else if ((E /= C0) .and. (EPS > C0)) then
            if (IOP9 > 0) if((abs(E)*dble(KIS-KQN+2)**(KQN+1))-1.D-2)&
     &         3310, 3310, 3300
 3300       KIS = -1
!++  End
         end if
 3310    continue
! ********
! CORRECT
! ********
         do 3320 K = 1, KORDI
!++  Code for ~{p,x} is active
            Y(IY - K) = Y(IY - K) + G(KQL + 1, K) * TPP
!++  Code for {p,x} is inactive
!c--D Next line special: P=>D, X=>Q
!            Y(IY - K) = Y(IY - K) + dble(G(KQL + 1, K)) * dble(TPP)
!++  END
 3320    continue
! END OF CORRECTING
 3325 continue
!++  Code for OUTPUT is active
      if (IOP10 > 0) then
         if (I == 1) then
            IDAT(1) = KSTEP
            IDAT(2) = LSC
            IDAT(3) = KSC
            IDAT(4) = IOP11
            FDAT(1) = TN
            FDAT(2) = HH
            FDAT(3) = EIMIN
            FDAT(4) = EAVE
            FDAT(5) = SIGMA(IOP11)
            FDAT(6) = ROBND
            MACT2(3) = NTE
!--D Next line special: P=>S, X=>D
            call DMESS(MACT1, MTXTAA, IDAT, FDAT)
            KOUTKO = NOUTKO
         end if
         if (KOUTKO /= 0) then
            if (KORD(KOUTKO) > 0) then
               if (I < KORD(KOUTKO)) go to 3328
               KOUTKO = KOUTKO + 1
            else
               if (I >= abs(KORD(KOUTKO))) KOUTKO = KOUTKO + 1
            end if
         end if
         IDAT(1) = I
         IDAT(2) = KQL
         IDAT(3) = LINC
         FDAT(1) = E
         FDAT(2) = EI
         FDAT(3) = EPS
         FDAT(4) = F(I)
         FDAT(5) = TPS1
         FDAT(6) = TPS2
         FDAT(7) = TPS3
         FDAT(8) = TPS4
         FDAT(9) = RNOISE
         FDAT(10) = 0.D0
         if (KQL == 1) FDAT(10) = S
         FDAT(11) = BETA(KQD)
!--D Next line special: P=>S, X=>D
         call DMESS(MACT2, MTXTAB, IDAT, FDAT)
 3328    if (I == NTE) IOP10 = IOP10 - 1
      end if
!++  End
 3330    L = L + NUMDT
 3340    continue
      return
    end subroutine divacr
!*************************************************************************

!*************************************************************************
!>
! SUBROUTINE TO COMPUTE COEFFICIENTS REQUIRED FOR INTEGRATING
! ORDINARY DIFFERENTIAL EQUATIONS
!
!### History
!  * 1988-05-20 DIVAHC Krogh   Initial code.

    subroutine DIVAHC

      integer KDIM, MAXORD, MAXSTF
!++ Substitute for KDIM, MAXORD, MAXSTF below
      parameter (KDIM = 20, MAXORD = 2, MAXSTF = 1)
!--D Next line special: P=>D, X=>Q
      double precision TN
      double precision XI(KDIM)
!
!--D Next line special: P=>D, X=>Q
      double precision TG(2), TGSTOP(2), TMARK, TMARKX, TOUT, TOLG
      double precision ALPHA(KDIM), BETA(KDIM+1)
      double precision  D(MAXSTF+MAXORD,MAXORD), G(KDIM,MAXORD)
      double precision V(KDIM+MAXORD)
      double precision HC, HDEC, HINC, HINCC, HMAX, HMAXP9, HMIN
      double precision FDAT(11)
!
      double precision DS(MAXSTF+MAXORD, MAXORD), GS(KDIM)
      double precision SIGMA(KDIM), RBQ(KDIM), DNOISE
      double precision EAVE, EIMAX, EIMIN, EMAX, EREP, ROBND, SNOISE
!
      integer IOPST, KORDI, KQMAXD, KQMAXI, LDT, MAXDIF, MAXINT, NKDKO, &
     &   NTE, NYNY, NDTF, NUMDT
      common / DIVASC / TN, XI, IOPST, KORDI, KQMAXD, KQMAXI, LDT,      &
     &   MAXDIF, MAXINT, NKDKO, NTE, NYNY, NDTF, NUMDT
!
      integer ICF,ICS,IGFLG,IGTYPE(2),IGSTOP(2),ILGREP,INGS,IOP3,IOP4,  &
     &   IOP5,IOP6,IOP7,IOP8,IOP9,IOP10,IOP11,IOP12,IOP13,IOP14,IOP15,  &
     &   IOP16,IOP17,IOP18,IOP19,IOP20,IOP21,IOP22,IOP21S,ITOLEP,IY,    &
     &   KEMAX,KIS,KMARK,KORD1I,KORD2I,KPRED,KQDCON,KQICON,KQMAXS,      &
     &   KQMXDS,KQMXIL,KQMXIP,KQMXIS,KSC,KSOUT,KSSTRT,KSTEP,LEX,LINC,   &
     &   LINCD,LINCQ,LSC,MAXKQD,MAXKQI,METHOD,NE,NEPTOL,NG,NGTOT,       &
     &   NOISEQ,NOUTKO,NTOLF,NY,IDAT(6)
      common /DIVAMC/ TG,TGSTOP,TMARK,TMARKX,TOUT,TOLG,HC,HDEC,HINC,    &
     &   HINCC,HMAX,HMAXP9,HMIN,ALPHA,BETA,D,G,V,DS,GS,SIGMA,RBQ,DNOISE,&
     &   EAVE,EIMAX,EIMIN,EMAX,EREP,ROBND,SNOISE,FDAT,ICF,ICS,IGFLG,    &
     &   IGTYPE,IGSTOP,ILGREP,INGS,IOP3,IOP4,IOP5,IOP6,IOP7,IOP8,IOP9,  &
     &   IOP10,IOP11,IOP12,IOP13,IOP14,IOP15,IOP16,IOP17,IOP18,IOP19,   &
     &   IOP20,IOP21,IOP22,IOP21S,ITOLEP,IY,KEMAX,KIS,KMARK,KORD1I,     &
     &   KORD2I,KPRED,KQDCON,KQICON,KQMAXS,KQMXDS,KQMXIL,KQMXIP,KQMXIS, &
     &   KSC,KSOUT,KSSTRT,KSTEP,LEX,LINC,LINCD,LINCQ,LSC,MAXKQD,MAXKQI, &
     &   METHOD,NE,NEPTOL,NG,NGTOT,NOISEQ,NOUTKO,NTOLF,NY,IDAT
      save / DIVAMC / , / DIVASC /
!.    SPECIFICATION OF ENVIRONMENTAL CONSTANTS.
      double precision EEPS10, EEPS16, EROV10, EEPS2, EEPT75, EOVEP2
      double precision OVTM75, OVD10
      common / DIVAEV / EEPS2, EEPT75, EOVEP2, OVTM75, OVD10, EEPS10,   &
     &   EEPS16, EROV10
      save / DIVAEV /
!                 K - 1 + 1 / K  is equivalent to max(1, K-1)
      double precision GG(MAXORD - 1 + 1/MAXORD), B(KDIM+MAXORD),       &
     &   W(KDIM+MAXORD)
      integer  K, N, J
      double precision C0, CP1, CRBQI, CP5, CP5625, C1, C1P125
      parameter (C0 = 0.D0)
      parameter (CP1 = .1D0)
      parameter (CRBQI = .421875D0)
      parameter (CP5 = .5D0)
      parameter (CP5625 = .5625D0)
      parameter (C1 = 1.D0)
      parameter (C1P125 = 1.125D0)
!++  Code for STIFF is inactive
!      INTEGER          GODIF
!++  End
      double precision TP1, TP2, HH, TEMP, TP
      equivalence (G(1, 1), HH)
!
      save GG, W
!
!  B(K)= 1/(K*(K+1))
!++ Save data by elements if ~.C.
!++ Of next 23 lines, only the first KDIM+MAXORD are active
      data B(1)  / 5.000000000000000000000000000000000000000D-1 /
      data B(2)  / 1.666666666666666666666666666666666666667D-1 /
      data B(3)  / 8.333333333333333333333333333333333333333D-2 /
      data B(4)  / 5.000000000000000000000000000000000000000D-2 /
      data B(5)  / 3.333333333333333333333333333333333333333D-2 /
      data B(6)  / 2.380952380952380952380952380952380952381D-2 /
      data B(7)  / 1.785714285714285714285714285714285714286D-2 /
      data B(8)  / 1.388888888888888888888888888888888888889D-2 /
      data B(9)  / 1.111111111111111111111111111111111111111D-2 /
      data B(10) / 9.090909090909090909090909090909090909091D-3 /
      data B(11) / 7.575757575757575757575757575757575757576D-3 /
      data B(12) / 6.410256410256410256410256410256410256410D-3 /
      data B(13) / 5.494505494505494505494505494505494505495D-3 /
      data B(14) / 4.761904761904761904761904761904761904762D-3 /
      data B(15) / 4.166666666666666666666666666666666666667D-3 /
      data B(16) / 3.676470588235294117647058823529411764706D-3 /
      data B(17) / 3.267973856209150326797385620915032679739D-3 /
      data B(18) / 2.923976608187134502923976608187134502924D-3 /
      data B(19) / 2.631578947368421052631578947368421052632D-3 /
      data B(20) / 2.380952380952380952380952380952380952381D-3 /
      data B(21) / 2.164502164502164502164502164502164502165D-3 /
      data B(22) / 1.976284584980237154150197628458498023715D-3 /
!     data B(23) / 1.811594202898550724637681159420289855072D-3 /
!
! ********
! START OF CODE
! ********
!     SET STEP NUMBER OF METHOD
!++  Code for STIFF is inactive
!      IOP11 = MIN(max(KQMAXI,KQMAXD) + 1), KDIM)
!++  Code for ~STIFF is active
      IOP11 = MIN(KQMAXI + 1, KDIM)
!++  End
!     TEST IF STEPSIZE WAS CHANGED
      if (KQICON >= 0) go to 3510
! ********
! STEPSIZE JUST CHANGED
! ********
!     SET CONSTANTS DEPENDING ON NEW STEPSIZE
      KQMXIL = KQMAXI
      TP1 = HH
      GG(1) = TP1 * TP1
      G(1, 2) = GG(1) * CP5
      if (MAXINT <= 2) go to 3450
      do 3440 K = 3, MAXINT
         GG(K - 1) = G(1, K - 1) * TP1
         G(1, K) = GG(K - 1) / dble(K)
 3440    continue
!     SET CONSTANTS INDICATING STEP CHANGE
 3450 KQICON = 0
!++  Code for STIFF is inactive
!      KQDCON=0
!++  End
      KQMXIP = 1
      KSC = 1
      if (LSC < 7) go to 3490
!     SPECIAL SET-UP OF CONSTANTS ON THE VERY FIRST STEP
      HINCC = C1P125
      LINCD = 6
      LINCQ = 12
      if (HINC > C0) go to 3460
      LINCD = -2
      LINC = -2
      ROBND = C1
 3460 SIGMA(1) = 1.0D0
      BETA(1) = C1
      do 3470 N = 1, IOP11
!++  Code for STIFF is inactive
!      D(1,N)=C0
!++  End
         XI(N) = TP1
         ALPHA(N) = C1
         BETA(N + 1) = C1
         SIGMA(N + 1) = dble(N + 1) * SIGMA(N) * HINCC
 3470    continue
      TEMP = EEPS16
      RBQ(1) = C1
      RBQ(2) = CP1
      TP = CRBQI
!     **** IN THE LOOP BELOW RBQ(K) IS COMPUTED TO BE
!          APPROXIMATELY (3/4 ** ((K-1) ** 2 - 1) / 10
!          .5625 = (3/4) ** 2    TP = (3/4) ** (2*K -3)
      do 3480 K = 3, KDIM
         TEMP = TEMP + TEMP
         RBQ(K) = max(TEMP, RBQ(K - 1) * TP)
 3480    TP = TP * CP5625
      go to 3560
!     SET-UP AFTER THE FIRST STEP
 3490 TP2 = XI(1)
      XI(1) = TP1
      BETA(2) = TP1 / TP2
      K = 2
      if (HINCC == HINC) go to 3540
      if ((LSC /= 0) .or. ((KSTEP-KSSTRT-KQMAXS) < 10)) go to 3540
      HINCC = C1
      LINCD = 0
 3500 LINCD = LINCD + 1
      HINCC = HINCC * HINC
      if (HINCC < 2.D0) go to 3500
      LINC = (LINC * (LINCD + LINCD)) / LINCQ
      LINCQ = LINCD + LINCD
      HINCC = HINC
      go to 3540
! END OF LOGIC FOR CASE WHEN STEPSIZE JUST CHANGED
!     TEST IF MAXIMUM INTEGRATION ORDER DID NOT INCREASE
 3510 if (KQMAXI > KQMXIL) then
! ********
! INTEGRATION ORDER WAS INCREASED -- GET NEW V'S
! ********
         KQMXIL = KQMAXI
         KQMXIP = KQMXIL + MAXINT
         K = KQMXIP
         V(K) = B(K)
         if (KQICON == 1) go to 3530
!     if (KQICON == K) KQICON = KQICON - 1 --- Removed 1999-08-19
         do 3520 N = 2, KQICON
            K = K - 1
 3520       V(K) = V(K) - ALPHA(N) * V(K + 1)
! END OF GETTING NEW V'S
      else
         IOP11 = max(IOP11, KQMXIL+1)
      end if
 3530 if (IOP11 <= KSC) go to 3560
! ********
! COMPUTE PARAMETERS WHICH ARE STILL CHANGING AS A RESULT OF
! A CHANGE IN THE STEPSIZE
! ********
      TP2 = XI(KSC)
!     UPDATE CONSTANT STEP COUNTER
      KSC = KSC + 1
      K = KSC
      BETA(K) = C1
 3540 continue
      TEMP = HINCC
!
!   LOOP TO COMPUTE NEW VALUES OF PARAMETERS
      do 3550 N = K, IOP11
         TP1 = TP2 + HH
         TP2 = XI(N)
         XI(N) = TP1
         ALPHA(N) = HH / TP1
         BETA(N + 1) = BETA(N) * (TP1 / TP2)
         TEMP = max(TEMP, dble(N) * (ALPHA(N) * HINCC))
         SIGMA(N) = SIGMA(N - 1) * TEMP
 3550    continue
      if (IOP11 /= KDIM) XI(IOP11 + 1) = TP2 + HH
! END OF CODE FOR COMPUTING PARAMETERS WHICH ARE STILL CHANGING
!
 3560 if (KQICON >= KQMXIP) go to 3690
! ********
! COMPUTE INTEGRATION COEFFICIENTS WHICH ARE STILL CHANGING
! ********
      KQMXIL = max(KQMAXI, KQMXIL)
      KQMXIP = KQMXIL + MAXINT
      J = KQMXIP - KQICON
      N = KQICON + 1
      KQICON = N
      if (N /= 1) go to 3580
! INITIALIZE V AND W
      do 3570 K = 1, J
         V(K) = B(K)
 3570    W(K) = V(K)
      go to 3600
! UPDATE V AND INITIALIZE W
 3580 if (N == KDIM) go to 3690
      do 3590 K = 1, J
         V(K) = V(K) - ALPHA(N) * V(K + 1)
 3590    W(K) = V(K)
! SET TRANSFER FOR LOOP BELOW DEPENDING ON VALUE OF MAXINT
 3600 continue
      go to 3660
!
 3640 J = J - 1
! INNER LOOP FOR COMPUTING INTEGRATION COEFFICIENTS
      do 3650 K = 1, J
 3650    W(K) = W(K) - ALPHA(N) * W(K + 1)
!     STORE INTEGRATION COEFFICIENTS
 3660 G(N + 1, 1) = HH * W(1)
      GS(N + 1) = G(N + 1, 1) - G(N, 1)
!++  Code for MAXORD >= 2 is active
      if (MAXINT >= 2) then
         G(N + 1, 2) = GG(1) * W(2)
!++  Code for MAXORD >= 3 is inactive
!        if (MAXINT > 2) then
!           DO 3665 K=3,MAXINT
!3665          G(N+1,K)=GG(K-1)*W(K)
!        end if
!++  Code for MAXORD >= 2 is active
      end if
!++  End
      N = N + 1
      if (N <= KQMXIL) go to 3640
! END OF COMPUTING INTEGRATION COEFFICIENTS
!
 3690 continue
!++  Code for STIFF is inactive
!      IF (KQDCON>KQMAXD) GO TO 4662
!c.********
!c.COMPUTE DIFFERENTIATION COEFFICIENTS WHICH ARE STILL CHANGING
!c.********
!c.SET TRANSFER FOR LOOP BELOW, DEPENDING ON VALUE OF MAXDIF
!++  Code for STIFF & MAXORD >= 2 is inactive
!      IF (MAXDIF-2) 3693,3692,3691
! 3691 ASSIGN 3696 TO GODIF
!      GO TO 3694
! 3692 ASSIGN 3698 TO GODIF
!      GO TO 3694
! 3693 ASSIGN 3699 TO GODIF
!++  Code for STIFF is inactive
! 3694 KQDCON=KQDCON+1
!c.LOOP FOR COMPUTING DIFFERENTIATION COEFFICIENTS
!      DO 3699 N=KQDCON,KQMAXD
!      DS(N+1,2)=C1/XI(N)
!      D(N+1,1)=DS(N+1,2)+D(N,1)
!      DS(N+1,1)=DS(N+1,2)/D(N+1,1)
!++  Code for STIFF & MAXORD >= 2 is inactive
!      GO TO GODIF, (3696,3698,3699)
! 3696 CONTINUE
!++  Code for STIFF & MAXORD >= 3 is inactive
!      DO 3697 K=3,MAXDIF
!      DS(N+1,K)=D(N,K-2) * (K-1)/XI(N)
! 3697 D(N+1,K-1)=DS(N+1,K) + D(N,K-1)
!++  Code for STIFF is inactive
! 3698 CONTINUE
!++  Code for STIFF & MAXORD >= 2 is inactive
!      D(N+1,MAXDIF)=D(N,MAXDIF) + D(N,MAXDIF-1) * (MAXDIF)/XI(N)
!++  Code for STIFF is inactive
! 3699 CONTINUE
!++  End
!
! END OF COMPUTING DIFFERENTIATION COEFFICIENTS
      return
    end subroutine divahc
!*************************************************************************

!*************************************************************************
!>
!  SUBROUTINE TO DO INTERPOLATION FOR VARIABLE ORDER INTEG. ROUTINE
!
!### History
!  * 1988-01-14 DIVAIN Krogh   Initial code.

    subroutine DIVAIN(T, Y, F, KORD)

      integer KORD(*)
!--D Next line special: P=>D, X=>Q
      double precision T(*), Y(*)
      double precision F(*)
      integer KDIM, MAXORD
!++ Substitute for KDIM, MAXORD below
      parameter (KDIM = 20, MAXORD = 2)
!--D Next line special: P=>D, X=>Q
      double precision TN
      double precision XI(KDIM)
!
      integer IOPST, KORDI, KQMAXD, KQMAXI, LDT, MAXDIF, MAXINT, NKDKO, &
     &   NTE, NYNY, NDTF, NUMDT
      common / DIVASC / TN, XI, IOPST, KORDI, KQMAXD, KQMAXI, LDT,      &
     &   MAXDIF, MAXINT, NKDKO, NTE, NYNY, NDTF, NUMDT
      save / DIVASC /
      integer I, ICI, IDT, INTERP, INTEG, INTEGZ, IY, IYI, IYN, IYNI, J,&
     &    K, KQMXI, KQMXS, KQQ, L, N
      double precision C0, C1, C2
      parameter (C0 = 0.D0)
      parameter (C1 = 1.D0)
      parameter (C2 = 2.D0)
      double precision C(KDIM+MAXORD-1), ETA(KDIM)
      double precision GAMMA(KDIM)
      double precision TP1, HI
      double precision CSUM(KDIM+MAXORD-1)
!--D Next line special: P=>D, X=>Q
      double precision XP1
      logical LNOTM1
!
!              Stuff for processing error messages
      integer IDAT(1)
      double precision FDAT(6)
      integer MENTXT, MERET, MEEMES, METEXT
      parameter (MENTXT =23)
      parameter (MERET  =51)
      parameter (MEEMES =52)
      parameter (METEXT =53)
      integer MACT(8)
! ********* Error message text ***************
![Last 2 letters of Param. name]  [Text generating message.]
!AA DIVAIN$B
!AB Interpolating at T(1)=$F with $B
!AC TN=$F, T(2)=$F and H=$F.  T(1) must be in [$F, $F].$E
!AD internal variable LDT = $I.  Interpolation not allowed now.$E
      integer LTXTAA,LTXTAB,LTXTAC,LTXTAD
      parameter (LTXTAA=  1,LTXTAB=  9,LTXTAC= 41,LTXTAD= 94)
      character MTXTAA(1) * (154)
      data MTXTAA/'DIVAIN$BInterpolating at T(1)=$F with $BTN=$F, T(2)=$&
     &F and H=$F.  T(1) must be in [$F, $F].$Einternal variable LDT = $I&
     &.  Interpolation not allowed now.$E'/
!
!                      1 2 3 4       5 6       7      8
      data MACT / MEEMES,0,0,0, MENTXT,0, METEXT, MERET /
!
!++  Code for ARGM is inactive
!      ENTRY DIVAIE
!++  End
! ********
! START OF CODE -- CHECK ON STATE OF DIFFERENCE TABLE
! ********
      L = LDT
      if (L) 3710, 3730, 3780
 3710 if (L + 2) 4170, 3730, 3720
 3720 if (MAXINT >= 0) L = 1
      go to 3840
! ********
! UPDATE DIFFERENCE TABLE TO START OF NEXT STEP
! ********
 3730 K = NDTF
      do 3770 I = 1, NTE
         KQQ = KORD(I + 3)
         if (KQQ <= 0) go to 3760
! EQUATION IS NOT STIFF
         TP1 = F(I) - F(K)
! LOOP TO DO UPDATING
         N = K + max(abs(KQQ), 2)
         do 3750 J = K, N
 3750       F(J) = F(J) + TP1
 3760    continue
         K = K + NUMDT
 3770    continue
      LDT = 1
      if (L /= 0) return
! END OF UPDATING DIFFERENCE TABLE
! ********
! INITIALIZE FOR COMPUTATION OF COEFFICIENTS
! ********
 3780 INTERP = 0
      HI = T(1) - TN
      GAMMA(1) = HI / XI(1)
      if (GAMMA(1)) 3790, 3800, 3810
 3790 if (GAMMA(1) >= -C1) go to 3820
      INTERP = 1
      if (abs(HI) - abs(T(2))) 3820, 3820, 4180
 3800 INTERP = 2 - KQMAXI
      go to 3820
 3810 if (GAMMA(1) > C2) if (LDT - 2) 4180, 3820, 3820
 3820 KQMXI = KQMAXI + INTERP - 1
!++  Code for STIFF is inactive
!      KQMXS=max(KQMXI,KQMAXD)
!++  Code for ~STIFF is active
      KQMXS = KQMXI
!++  End
      do 3830 N = 2, KQMXS
 3830    GAMMA(N) = (HI + XI(N-1)) / XI(N)
 3840 LNOTM1 = L /= -1
      INTEG = MAXINT
      if (INTEG <= 0) if (INTEG + MAXDIF) 4160, 3950, 3950
! ********
! COMPUTE INTEGRATION COEFFICIENTS
! ********
!     INITIAL SET-UP
!         COMPUTE INITIAL C VALUES
      do 3850 N = 1, INTEG
         C(N) = HI / dble(N)
 3850 continue
      I = INTEG + 1
      INTEG = INTEG + KQMXI
      do 3860 N = I, INTEG
         C(N) = C(N - 1) * (dble(N - MAXINT) / dble(N))
 3860    continue
!         COMPUTE ETA'S
      do 3870 N = 1, KQMXI
 3870    ETA(N) = HI / XI(N)
!         COMPUTE C(K)'S TO CORRESPOND TO G(K-MAXINT+1,MAXINT),
!         K=MAXINT, MAXINT+1,..., MAXINT+KQMXI-1
      I = INTEG
 3880 J = INTEG
      INTEG = J - 1
      if (INTEG <= MAXINT) go to 3900
      do 3890 N = J, I
 3890    C(N) = ETA(N - INTEG) * C(N) + C(N - 1)
      go to 3880
 3900 do 3910 N = J, I
 3910    C(N) = ETA(N - INTEG) * C(N)
!         END OF COMPUTING  G(---,MAXINT)
      INTEGZ = 0
      go to 3940
!         COMPUTE C(K)-S TO CORRESPOND TO G(K-INTEG+1,INTEG),
!         K=INTEG+1,INTEG+2,..., INTEG+KQMXI
 3920 do 3930 N = 1, KQMXI
 3930    C(INTEG+N) = GAMMA(N)*C(INTEG+N-1) - ETA(N)*C(INTEG+N)
 3940 ICI = INTEG - 1
      go to 4020
! END OF COMPUTING INTEGRATION COEFFICIENTS
! ********
! COMPUTE COEFFICIENTS FOR INTERPOLATION
! ********
 3950 C(1) = C1
      ICI = 0
      do 3960 N = 1, KQMXS
 3960    C(N + 1) = GAMMA(N) * C(N)
      if (INTEG + 1) 3970, 3990, 4010
! END OF COMPUTING INTERPOLATION COEFFICIENTS
!
!     SET-UP TO COMPUTE DIFFERENTIATION COEFFICIENTS REQUIRED
!     IN ORDER TO GET COEFFICIENTS ACTUALLY USED
 3970 INTEG = 0
      ICI = 1
 3980 INTEG = INTEG - 1
      if (INTEG == MAXINT) ICI = 0
! ********
! COMPUTE DIFFERENTIATION COEFFICIENTS
! ********
 3990 INTERP = max(INTERP, 0)
      TP1 = dble(-INTEG)
      C(1) = TP1 * C(1) / XI(-INTEG)
      J = KQMAXD + INTEG
      do 4000 N = 1, J
 4000    C(N + 1) = (TP1*C(N)) / XI(N - INTEG) + GAMMA(N - INTEG) * C(N)
!     C(N) NOW CORRESPONDS TO THE DIFFERENTIAL COEFFICIENT
!          D(N-INTEG,-INTEG)
 4010 INTEGZ = INTEG
      if (ICI /= 0) go to 3980
! END OF COMPUTING DIFFERENTIATION COEFFICIENTS
! ********
! BEGINNING OF LOOP TO DO
!         INTEGRATION       (INTEG>0)
!         INTERPOLATION     (INTEG==0)
!         DIFFERENTIATION   (INTEG<0)
! TO THE POINT INDICATED BY T.
! ********
!     SET UP INITIAL INDICES
 4020 if (NYNY < 0) then
         IY = -NYNY
         IYNI = NYNY + ICI + 1
         if (LDT == 2) then
            CSUM(ICI+1) = C(ICI+1)
            do 4025 J = ICI+2, INTEG+KQMXI
            CSUM(J) = CSUM(J-1) + C(J)
 4025       continue
         end if
      else
         IY = 1
         IYNI = NYNY + ICI - 1
      end if
      IDT = NDTF - INTEGZ
      do 4140 I = 1, NTE
         if (NKDKO /= 0) KORDI = KORD(NKDKO + I - 1)
         IY = IY + abs(KORDI)
         KQQ = KORD(I + 3)
!         GET INDEX OF HIGHEST ORDER DIFFERENCE TO BE USED
         K = max(abs(KQQ) + INTERP, 2)
         IYI = -INTEG
         if (KQQ) 4030, 4130, 4040
! EQUATION IS STIFF
 4030    continue
!++  Code for STIFF is inactive
!      JS=abs(KORD(NJSKO+I-1))-1
!      IYI=IYI-JS
!      IF(LNOTM1) IF (IYI) 4034,4032,4130
!      IF (KORDI<0) IYI=IYI+1
!      IYI=IYI+MAXINT-abs(KORDI)
!      IF (IYI) 4034,4130,4130
!c.      IF EQUATION IS IMPLICIT DO NOT COMPUTE AN F
! 4032 IF (KORDI<0) GO TO 4130
!c.      TEST IF INTEG TOO BIG FOR THIS EQUATION
! 4034 IF (abs(KORDI)<-IYI) GO TO 4130
!      IYI=IYI+IY
!      IYN=IYI+IYNI
!c. COMPUTE INNER PRODUCT FOR STIFF EQUATIONS
!      IF (INTEGZ==0) GO TO ???
!c.    DIFFERENTIATING
!      TP1 = C0
!      DO 4036 J = K+INTEGZ, 1, -1
!         TP1 = TP1 + C(J) * F(IDT+J-1)
! 4036 CONTINUE
!c.    TEST WHETHER TO STORE RESULT IN Y OR F
!      IF (IYI-IY) 4080, 4090, 4080
!c.    INTEGRATING OR INTERPOLATING
!      TP1 = C0
!      DO 4037 J = ICI + K, ICI + 2, -1
!         TP1 = TP1 + C(J) * F(IDT+J-ICI-1)
! 4037 CONTINUE
!      IF (INTEG==0) GO TO 4120
!      TP1=TP1 + C(ICI+1)*Y(IYN+1)
!++  End
         go to 4100
! END OF SPECIAL CODE FOR STIFF EQUATIONS
!
! EQUATION IS NOT STIFF
 4040    if (LNOTM1) if (IYI) 4050, 4060, 4130
         IYI = IYI + MAXINT - KORDI
         if (IYI >= 0) go to 4130
!       TEST IF INTEG TOO BIG FOR THIS EQUATION
 4050    if (KORDI < -IYI) go to 4130
 4060    IYI = IYI + IY
         IYN = IYI + IYNI
!  COMPUTE INNER PRODUCT FOR EQUATION WHICH IS NOT STIFF
         XP1 = C0
         if (LDT == 2) then
            if (KQQ /= KQMAXI) XP1 = CSUM(K+INTEGZ+ICI) *             &
     &         F(IDT+INTEGZ+NUMDT-1)
         end if
         do 4070 J = K + INTEGZ + ICI, ICI + 1, -1
            XP1 = XP1 + C(J) * F(IDT - ICI - 1 + J)
 4070       continue
         if (INTEG) 4080, 4090, 4100
! STORE FINAL RESULT IN Y WHEN DIFFERENTIATING
 4080    continue
         Y(IYI) = XP1
         go to 4130
! STORE INTERPOLATED VALUE IN F (OR STIFF DIFFERENTIATION)
 4090    F(I) = XP1
         go to 4130
! PICK UP EXTRA STUFF TO ADD TO INNER PRODUCT WHEN INTEGRATING
 4100    K = ICI
         if (K == 0) go to 4120
 4110    continue
         XP1 = C(K) * (XP1 + Y(IYN))
         IYN = IYN - 1
         K = K - 1
         if (K /= 0) go to 4110
! STORE FINAL RESULT IN Y WHEN INTEGRATING (OR STIFF INTERPOLATION)
 4120    Y(IYI) = XP1 + Y(IYN)
 4130    continue
         IDT = IDT + NUMDT
 4140    continue
!
      INTEG = INTEG - 1
      if (INTEG >= -MAXDIF) if (INTEG) 3990, 3950, 3920
 4160 return
! ********
! ERROR PROCESSING
! ********
 4170 MACT(2) = 68
      MACT(3) = 11
      MACT(6) = LTXTAD
      IDAT(1) = LDT
      go to 4190
 4180 MACT(2) = 28
      MACT(3) = 1
      MACT(6) = LTXTAC
      FDAT(2) = TN
      FDAT(3) = T(2)
      FDAT(4) = XI(1)
      FDAT(5) = TN - T(2)
      FDAT(6) = TN + C2 * XI(1)
      if (XI(1) < 0) then
         FDAT(5) = FDAT(6)
         FDAT(6) = TN - T(2)
      end if
 4190 FDAT(1) = T(1)
!--D Next line special: P=>S, X=>D
      call DMESS(MACT, MTXTAA, IDAT, FDAT)
      if (MACT(2) < 50) go to 3820
      return
    end subroutine divain
!*************************************************************************

!*************************************************************************
!>
!  SUBROUTINE TO SET UP OPTIONS FOR DIFFERENTIAL EQUATION PACKAGE -IVA
!
!### History
!  * 1987-12-07 DIVAOP Krogh   Initial code.

    subroutine DIVAOP(IOPT, FOPT)
      double precision FOPT(*)
      integer IOPT(*)
!
      integer KDIM, MAXORD, MAXSTF
!++ Substitute for KDIM, MAXORD, MAXSTF below
      parameter (KDIM = 20, MAXORD = 2, MAXSTF = 1)
!--D Next line special: P=>D, X=>Q
      double precision TN
      double precision XI(KDIM)
!
!--D Next line special: P=>D, X=>Q
      double precision TG(2), TGSTOP(2), TMARK, TMARKX, TOUT, TOLG
      double precision ALPHA(KDIM), BETA(KDIM+1)
      double precision  D(MAXSTF+MAXORD,MAXORD), G(KDIM,MAXORD)
      double precision V(KDIM+MAXORD)
      double precision HC, HDEC, HINC, HINCC, HMAX, HMAXP9, HMIN
      double precision FDAT(11)
!
      double precision DS(MAXSTF+MAXORD, MAXORD), GS(KDIM)
      double precision SIGMA(KDIM), RBQ(KDIM), DNOISE
      double precision EAVE, EIMAX, EIMIN, EMAX, EREP, ROBND, SNOISE
!
      integer IOPST, KORDI, KQMAXD, KQMAXI, LDT, MAXDIF, MAXINT, NKDKO, &
     &   NTE, NYNY, NDTF, NUMDT
      common / DIVASC / TN, XI, IOPST, KORDI, KQMAXD, KQMAXI, LDT,      &
     &   MAXDIF, MAXINT, NKDKO, NTE, NYNY, NDTF, NUMDT
!
      integer ICF,ICS,IGFLG,IGTYPE(2),IGSTOP(2),ILGREP,INGS,IOP3,IOP4,  &
     &   IOP5,IOP6,IOP7,IOP8,IOP9,IOP10,IOP11,IOP12,IOP13,IOP14,IOP15,  &
     &   IOP16,IOP17,IOP18,IOP19,IOP20,IOP21,IOP22,IOP21S,ITOLEP,IY,    &
     &   KEMAX,KIS,KMARK,KORD1I,KORD2I,KPRED,KQDCON,KQICON,KQMAXS,      &
     &   KQMXDS,KQMXIL,KQMXIP,KQMXIS,KSC,KSOUT,KSSTRT,KSTEP,LEX,LINC,   &
     &   LINCD,LINCQ,LSC,MAXKQD,MAXKQI,METHOD,NE,NEPTOL,NG,NGTOT,       &
     &   NOISEQ,NOUTKO,NTOLF,NY,IDAT(6)
      common /DIVAMC/ TG,TGSTOP,TMARK,TMARKX,TOUT,TOLG,HC,HDEC,HINC,    &
     &   HINCC,HMAX,HMAXP9,HMIN,ALPHA,BETA,D,G,V,DS,GS,SIGMA,RBQ,DNOISE,&
     &   EAVE,EIMAX,EIMIN,EMAX,EREP,ROBND,SNOISE,FDAT,ICF,ICS,IGFLG,    &
     &   IGTYPE,IGSTOP,ILGREP,INGS,IOP3,IOP4,IOP5,IOP6,IOP7,IOP8,IOP9,  &
     &   IOP10,IOP11,IOP12,IOP13,IOP14,IOP15,IOP16,IOP17,IOP18,IOP19,   &
     &   IOP20,IOP21,IOP22,IOP21S,ITOLEP,IY,KEMAX,KIS,KMARK,KORD1I,     &
     &   KORD2I,KPRED,KQDCON,KQICON,KQMAXS,KQMXDS,KQMXIL,KQMXIP,KQMXIS, &
     &   KSC,KSOUT,KSSTRT,KSTEP,LEX,LINC,LINCD,LINCQ,LSC,MAXKQD,MAXKQI, &
     &   METHOD,NE,NEPTOL,NG,NGTOT,NOISEQ,NOUTKO,NTOLF,NY,IDAT
      save / DIVAMC / , / DIVASC /
!.    SPECIFICATION OF ENVIRONMENTAL CONSTANTS.
      double precision EEPS10, EEPS16, EROV10, EEPS2, EEPT75, EOVEP2
      double precision OVTM75, OVD10
      common / DIVAEV / EEPS2, EEPT75, EOVEP2, OVTM75, OVD10, EEPS10,   &
     &   EEPS16, EROV10
      save / DIVAEV /
!
      integer IOPTS(23), INCOP(22), IOPTC(23), I, IA, J, K, LIOPT, MULTJ
      double precision CMP75, C0, CP25, CP3, CP5, CP625, CP75, CP875,   &
     &   CP9, C1, C1P125, C2, C4, C10, C16
      parameter (CMP75 = (-.75D0))
      parameter (C0 = 0.D0)
      parameter (CP25 = .25D0)
      parameter (CP3 = .3D0)
      parameter (CP5 = .5D0)
      parameter (CP625 = .625D0)
      parameter (CP75 = .75D0)
      parameter (CP875 = .875D0)
      parameter (CP9 = .9D0)
      parameter (C1 = 1.D0)
      parameter (C1P125 = 1.125D0)
      parameter (C2 = 2.D0)
      parameter (C4 = 4.D0)
      parameter (C10 = 10.D0)
      parameter (C16 = 16.D0)
      !external D1MACH
      !double precision D1MACH
      equivalence (IOPTC(3), IOP3)
      save IOPTS, LIOPT
!
!                      Declarations for error message processing.
!
      integer MECONT, MERET, MEEMES, MEIVEC
      parameter (MECONT =50)
      parameter (MERET  =51)
      parameter (MEEMES =52)
      parameter (MEIVEC =57)
      integer MACT(7), MACT1(5)
!
! ********* Error message text ***************
![Last 2 letters of Param. name]  [Text generating message.]
!AA DIVAOP$B
!AB Error in IOPT() specifications: IOPT =$E
!AC HMIN = $F is > HMAX = $F.$E
      integer LTXTAA,LTXTAB,LTXTAC
      parameter (LTXTAA= 1,LTXTAB= 9,LTXTAC=49)
      character MTXTAA(1) * (76)
      data MTXTAA/'DIVAOP$BError in IOPT() specifications: IOPT =$EHMIN$&
     & = $F is > HMAX = $F.$E'/
! **** End of text generated by pmess
!                      1   2  3   4       5  6      7
      data MACT / MEEMES, 88, 24, 0, MEIVEC, 0, MERET /
      data MACT1 / MEEMES, 28, 24, LTXTAC, MERET /
!
!                      IOP4       IOP17
      data IOPTS / 3*0, 500000, 12*0, 1, 6*0 /
!
!                  1  2    3  8  9 10 11 12  13  16   17 21 22
      data INCOP / 1, 3, 5*2, 1, 2, 3, 1, 2, 3*1, 3, 4*2, 2, 2 /
!
! ********* START OF EXECUTABLE CODE ***********************
!
      MULTJ = 1
      K = 1
 4200 I = IOPT(K)
      IA = abs(I)
! 1 and 6 lines below need 21 changed if more options are added.
      if (IA <= 21) if (I) 4220, 4520, 4280
      if (IA /= 1111) go to 4490
      if (I < 0) then
        MULTJ = -1
        K = K + 1
        go to 4200
      end if
      IOPT(2) = LIOPT
!
!     ****  INITIALIZE FOR STARTING A NEW INTEGRATION
      do 4210 J = 3, 23
 4210    IOPTC(J) = IOPTS(J)
      KSOUT = IOPTS(4)
      KMARK = 1 - IOPTS(1)
      KORDI = IOPTS(17)
      NKDKO = max(-KORDI, 0)
      IOPST = IOPTS(22)
      go to 4260
!
!     **** SET A NOMINAL VALUE
 4220 IOPTS(IA) = 0
      if (IA == 12) go to 4420
      if (IA - 2) 4400, 4240, 4230
 4230 if (IA == 4) IOPTS(4) = 500000
      if (IA == 21) TOLG = 0.D0
      go to 4390
!
!     **** SET ALL OPTIONS TO THEIR NOMINAL VALUES
 4240 IA = 1
      IOPTS(1) = 0
      do 4250 J = 3, 22
         IOPTS(J) = 0
         IOPTC(J) = 0
 4250 continue
      IOPTS(4) = 500000
      IOPTC(4) = IOPTS(4)
      IOPTS(17) = 1
      TOLG = 0.D0
 4260 NGTOT = IOPTS(7) + max(IOPTS(6), 0)
      if (IOPTS(12) == 0) go to 4420
 4270 return
!
!     **** SET SPECIFIED OPTION
 4280 J = IOPT(K + 1)
      if (INCOP(IA) - 2) 4290, 4330, 4300
!     **** OPTION INVOLVES NO EXTRA PARAMETERS
 4290 IOPTS(IA) = 1
      if (IA - 2) 4400, 4400, 4390
!     **** TAKE CARE OF SECOND EXTRA PARAMETER
 4300 if (IA /= 10) go to 4310
      NOUTKO = IOPT(K + 2)
      if (NOUTKO) 4500, 4350, 4350
 4310 if (IA /= 16) go to 4320
      NTOLF = IOPT(K + 2)
      if (NTOLF) 4500, 4500, 4350
 4320 if (J == 3) then
        if (KMARK /= 3) then
           if (XI(1)*(FOPT(IOPT(K+2)) - TMARK) >= C0) go to 4400
        end if
      end if
      TMARK = FOPT(IOPT(K+2))
      KMARK = J
      go to 4400
!     **** TAKE CARE OF FIRST EXTRA PARAMETER
 4330 continue
      if (IA == 12) go to 4410
      if (IA == 4) KSOUT = J
      if (IA == 21) TOLG = FOPT(J)
 4350 IOPTS(IA) = J * MULTJ
      if (abs(IA - 7) > 1) go to 4360
!     **** SET SPECIAL PARAMETERS FOR GSTOP-S
      IGFLG = 0
      NGTOT = IOPTS(7) + max(IOPTS(6), 0)
!     **** TEST FOR ERROR
      if (J > 500) go to 4500
 4360 if (J > 0) go to 4390
      if ((IA == 5) .or. (IA == 17)) go to 4390
      if (J + 1) 4500, 4370, 4380
 4370 if (IA == 7) go to 4500
 4380 if ((IA == 4) .or. (IA == 11) .or. (IA >= 16)) go to 4500
!     **** STORE SAVED VALUE IN COMMON
 4390 IOPTC(IA) = IOPTS(IA)
!
!     **** INCREMENT K TO GET NEXT OPTION
 4400 K = K + INCOP(IA)
      go to 4200
!
! ******* SET UP INFORMATION FOR CHANGING STEPSIZE *********
!
!     **** TEST IF VALUES ARE ALREADY SET
 4410 if (IOPTS(12) /= 0) go to 4430
!     **** SET NOMINAL VALUES FOR VARIABLES ONLY SET ONCE
 4420 EREP = CP3
!     **** SET NOMINAL VALUES FOR STEPSIZE CONTROL AND ENV. CONSTANTS
      EEPS2 = D1MACH(4)
      EEPS16 = C16 * EEPS2
      EEPS10 = CP625 * EEPS16
      EEPT75 = EEPS2 ** CP75
      EEPS2 = EEPS2 + EEPS2
      OVD10 = D1MACH(2)
      EROV10 = C10 / OVD10
      EOVEP2 = OVD10 * EEPS2
      OVTM75 = OVD10 ** CMP75
      OVD10 = OVD10 / C10
      HINC = C2
      HDEC = CP5
      HMIN = EROV10
      HMAX = OVD10
      if (I /= 12) go to 4470
 4430 IOPTS(12) = J * MULTJ
      if (J) 4450, 4470, 4460
!     **** SET UP TO GIVE USER COMPLETE STEPSIZE CONTROL
 4450 EREP = C1 / EROV10
      HINC = -C2
      IOP8 = 1
!## Recent code 12/16/94
      LINCD = -2
      LINC = -2
      ROBND = C1
!## Recent code 9/6/2001
      TOLG = 0.D0
!## End of recent code
      go to 4480
!     **** SET USER VALUES FOR STEPSIZE CONTROL
 4460 if (FOPT(J) /= C0) HINC = max(C1P125, min(FOPT(J),C4))
      if (FOPT(J + 1) /= C0) HDEC = min(CP875, max(FOPT(J + 1), CP25))
      if (FOPT(J + 2) /= C0) HMIN = FOPT(J + 2)
      if (FOPT(J + 3) /= C0) HMAX = FOPT(J + 3)
      if ((HMIN > HMAX) .or. (HMAX <= 0.D0)) then
         call DMESS(MACT1, MTXTAA, IOPT, FOPT(J+2))
         KORD2I = -4
      end if
 4470 KQICON = -1
 4480 HMAXP9 = HMAX * CP9
      if (I - 1111) 4400, 4270, 4400
!
! ***************** ERROR  IN  IOPT ************************
!
 4490 IA = 1
 4500 MACT(6) = K + INCOP(IA) - 1
      call MESS(MACT, MTXTAA, IOPT)
      KORD1I = 24
      KORD2I = -4
! Usual return with no error is here.
 4520 LIOPT = K
      return
    end subroutine divaop
!*************************************************************************

!*************************************************************************
!>
! THIS SUBROUTINE
!  1. UPDATES THE DIFFERENCE TABLE FROM THE PREVIOUS STEP (IF NOT
!     DONE ALREADY).
!  2. PREDICTS WHAT THE VALUES OF THE DEPENDENT VARIABLES, Y, AND
!     THE DIFFERENCE TABLE, DT, WILL BE AT THE END OF THE CURRENT STEP.
!
!   Y = VECTOR OF PREDICTED VALUES COMPUTED BY THIS SUBROUTINE.
!   YN= VECTOR OF VALUES OF Y COMPUTED ON THE LAST STEP.
!   F = VECTOR OF DERIVATIVE VALUES.
!   DT= ARRAY CONTAINING DIFFERENCE TABLES.
!   KD= VECTOR GIVING ORDERS OF THE DIFFERENTIAL EQUATIONS (IF
!       EQUATIONS HAVE DIFFERENT ORDERS).
!   KQ= VECTOR OF INTEGRATION ORDERS.
!
!### History
!  * 1988-01-13 DIVAPR Krogh   Initial code.

    subroutine DIVAPR(Y, YN, F, KORD)

      integer KORD(*)
!--D Next line special: P=>D, X=>Q
      double precision Y(*), YN(*)
      double precision F(*)
!
      integer KDIM, MAXORD, MAXSTF
!++ Substitute for KDIM, MAXORD, MAXSTF below
      parameter (KDIM = 20, MAXORD = 2, MAXSTF = 1)
!--D Next line special: P=>D, X=>Q
      double precision TN
      double precision XI(KDIM)
!
!--D Next line special: P=>D, X=>Q
      double precision TG(2), TGSTOP(2), TMARK, TMARKX, TOUT, TOLG
      double precision ALPHA(KDIM), BETA(KDIM+1)
      double precision  D(MAXSTF+MAXORD,MAXORD), G(KDIM,MAXORD)
      double precision V(KDIM+MAXORD)
      double precision HC, HDEC, HINC, HINCC, HMAX, HMAXP9, HMIN
      double precision FDAT(11)
!
      double precision DS(MAXSTF+MAXORD, MAXORD), GS(KDIM)
      double precision SIGMA(KDIM), RBQ(KDIM), DNOISE
      double precision EAVE, EIMAX, EIMIN, EMAX, EREP, ROBND, SNOISE
!
!.    SPECIFICATION OF ENVIRONMENTAL CONSTANTS.
      double precision EEPS10, EEPS16, EROV10, EEPS2
      double precision EEPT75, EOVEP2, OVTM75, OVD10
      common / DIVAEV / EEPS2, EEPT75, EOVEP2, OVTM75, OVD10, EEPS10,   &
     &   EEPS16, EROV10
      save / DIVAEV /
      integer IOPST, KORDI, KQMAXD, KQMAXI, LDT, MAXDIF, MAXINT, NKDKO, &
     &   NTE, NYNY, NDTF, NUMDT
      common / DIVASC / TN, XI, IOPST, KORDI, KQMAXD, KQMAXI, LDT,      &
     &   MAXDIF, MAXINT, NKDKO, NTE, NYNY, NDTF, NUMDT
!
      integer ICF,ICS,IGFLG,IGTYPE(2),IGSTOP(2),ILGREP,INGS,IOP3,IOP4,  &
     &   IOP5,IOP6,IOP7,IOP8,IOP9,IOP10,IOP11,IOP12,IOP13,IOP14,IOP15,  &
     &   IOP16,IOP17,IOP18,IOP19,IOP20,IOP21,IOP22,IOP21S,ITOLEP,IY,    &
     &   KEMAX,KIS,KMARK,KORD1I,KORD2I,KPRED,KQDCON,KQICON,KQMAXS,      &
     &   KQMXDS,KQMXIL,KQMXIP,KQMXIS,KSC,KSOUT,KSSTRT,KSTEP,LEX,LINC,   &
     &   LINCD,LINCQ,LSC,MAXKQD,MAXKQI,METHOD,NE,NEPTOL,NG,NGTOT,       &
     &   NOISEQ,NOUTKO,NTOLF,NY,IDAT(6)
      common /DIVAMC/ TG,TGSTOP,TMARK,TMARKX,TOUT,TOLG,HC,HDEC,HINC,    &
     &   HINCC,HMAX,HMAXP9,HMIN,ALPHA,BETA,D,G,V,DS,GS,SIGMA,RBQ,DNOISE,&
     &   EAVE,EIMAX,EIMIN,EMAX,EREP,ROBND,SNOISE,FDAT,ICF,ICS,IGFLG,    &
     &   IGTYPE,IGSTOP,ILGREP,INGS,IOP3,IOP4,IOP5,IOP6,IOP7,IOP8,IOP9,  &
     &   IOP10,IOP11,IOP12,IOP13,IOP14,IOP15,IOP16,IOP17,IOP18,IOP19,   &
     &   IOP20,IOP21,IOP22,IOP21S,ITOLEP,IY,KEMAX,KIS,KMARK,KORD1I,     &
     &   KORD2I,KPRED,KQDCON,KQICON,KQMAXS,KQMXDS,KQMXIL,KQMXIP,KQMXIS, &
     &   KSC,KSOUT,KSSTRT,KSTEP,LEX,LINC,LINCD,LINCQ,LSC,MAXKQD,MAXKQI, &
     &   METHOD,NE,NEPTOL,NG,NGTOT,NOISEQ,NOUTKO,NTOLF,NY,IDAT
      save / DIVAMC / , / DIVASC /
      double precision C0
      parameter (C0 = 0.D0)
!
      integer I, INTEG, INTEGS, J, K, KQQ, L, N
      double precision TEMP(KDIM)
      double precision TP1
!--D Next line special: P=>D, X=>Q
      double precision XP
      data INTEGS / -1 /
!
!++  Code for ARGM is inactive
!      RETURN
!      ENTRY DIVAPE
!++  End
! ********
! START OF CODE
! ********
      IY = 0
      L = NDTF - 1
      do 4680 I = 1, NTE
         INTEG = KORDI
         if (NKDKO /= 0) INTEG = KORD(NKDKO + I - 1)
         KQQ = KORD(I + 3)
         K = max(abs(KQQ), 2)
         if (KQQ) 4530, 4520, 4540
 4520    IY = IY + abs(INTEG)
         go to 4670
! ********
! EQUATION IS STIFF, OR IMPLICIT
! ********
 4530    continue
!++  Code for STIFF is inactive
!      KQQ=-KQQ
!      N=KQQ-1
!      JS=abs(KORD(NJSKO+I-1))-1
!      IMPLIC=INTEG
!      INTEG=abs(IMPLIC)-JS
!c.    SET INTEGS FOR STIFF EQUATIONS
!      INTEGS=0
!      IF (K-KSC) 160,160,140
!c.END OF SET-UP FOR STIFF EQUATIONS
!++  End
! ********
! EQUATION IS NOT STIFF
! ********
 4540    N = KQQ
         if (LDT /= 0) if (K - KSC) 4570, 4570, 4550
!     DIFFERENCE TABLE HAS NOT BEEN UPDATED
         TP1 = F(I) - F(L + 1)
         if (K - KSC) 4610, 4610, 4590
! END OF SET-UP FOR EQUATIONS WHICH ARE NOT STIFF
! ********
! GET PREDICTED DIFFERENCES FROM UPDATED DIFFERENCE TABLE
! ********
 4550    F(L + K + 1) = F(L + K + 1) * BETA(K + 1)
         TEMP(K) = F(L + K) * BETA(K)
         F(L + K) = TEMP(K)
! LOOP FOR MODIFIED DIVIDED DIFFERENCES
 4560    K = K - 1
         if (K <= KSC) go to 4580
         TEMP(K) = F(L + K) * BETA(K)
         F(L + K) = TEMP(K) + F(L + K + 1)
         go to 4560
! CODE FOR BACKWARD DIFFERENCES
 4570    F(L + K + 1) = F(L + K + 1)
         TEMP(K) = F(L + K)
         K = K - 1
!
 4580    TEMP(K) = F(L + K)
         F(L + K) = TEMP(K) + F(L + K + 1)
         K = K - 1
         if (K /= 0) go to 4580
         go to 4630
! ********
! UPDATE DIFFERENCE TABLE AND GET PREDICTED DIFFERENCES
! ********
! CODE FOR MODIFIED DIVIDED DIFFERENCES
 4590    F(L + K + 1) = (F(L+K+1) + TP1) * BETA(K + 1)
         TEMP(K) = (F(L + K) + TP1) * BETA(K)
         F(L + K) = TEMP(K)
 4600    K = K - 1
         if (K <= KSC) go to 4620
         TEMP(K) = (F(L + K) + TP1) * BETA(K)
         F(L + K) = TEMP(K) + F(L + K + 1)
         go to 4600
! CODE FOR BACKWARD DIFFERENCES
 4610    F(L + K + 1) = (F(L+K+1) + TP1)
         TEMP(K) = F(L + K) + TP1
         F(L + K) = TEMP(K)
         K = K - 1
!
 4620    TEMP(K) = F(L + K) + TP1
         F(L + K) = TEMP(K) + F(L + K + 1)
         K = K - 1
         if (K /= 0) go to 4620
! ********
! COMPUTE Y-S OBTAINED USING INTEGRATION
! ********
!     TEST IF NEXT Y TO BE OBTAINED BY INTERPOLATION
 4630    continue
!++  Code for STIFF is inactive
!      IF (INTEG==0) GO TO 4662
!++  End
         IY = IY + 1
!     FORM INNER PRODUCT
         XP = C0
         do 4650 J = INTEGS + N + 1, INTEGS + 2, -1
!++  Code for ~{p,x} is active
            XP = XP + G(J, INTEG) * TEMP(J)
!++  Code for {p,x} is inactive
!c--D Next line special: P=>D, X=>Q
!            XP = XP + dble(G(J, INTEG)) * dble(TEMP(J))
!++  END
 4650    continue
         K = INTEG + INTEGS
         do 4660 J = K, 1, -1
!++  Code for ~{p,x} is active
            XP = XP + G(1, J) * YN(IY + J)
!++  Code for {p,x} is inactive
!c--D Next line special: P=>D, X=>Q
!            XP = XP + dble(G(1, J)) * dble(YN(IY + J))
!++  END
 4660    continue
         Y(IY) = YN(IY) + XP
         INTEG = INTEG - 1
         if (K) 4670, 4670, 4630
! END OF COMPUTING Y-S OBTAINED BY INTEGRATION
! ********
! COMPUTE Y-S OBTAINED USING INTERPOLATION AND DIFFERENTIATION
! ********
!++  Code for STIFF is inactive
!c.    RESTORE INTEGS FOR EQUATIONS WHICH ARE NOT STIFF
! 4662 INTEGS=-1
!      IY=IY+1
!c.    COMPUTE Y USING INTERPOLATION
!      Y(IY)=YN(IY) + F(L+2)
!      IF (KQQ==1) Y(IY)=YN(IY)
! 4663 INTEG=INTEG+1
!      IF (INTEG==JS) IF (IMPLIC) 4680,4680,4664
!c.    COMPUTE INTEG-TH DERIVATIVE
!      XP = C0
! 4664 DO 4666 J = KQQ+1, INTEG+1, -1
!         XP = XP + D(J, INTEG) * TEMP(J)
! 4666 CONTINUE
!      IF (INTEG==JS) GO TO 4667
!      IY=IY+1
!      Y(IY)=XP
!      GO TO 4663
!c.STORE PREDICTED VALUE FOR F
! 4667 CONTINUE
!      F(L+NUMDT)=XP
!++  End
 4670    L = L + NUMDT
 4680    continue
      LDT = -3
      return
    end subroutine divapr
!*************************************************************************

!*************************************************************************
!>
! SUBROUTINE TO GIVE DEBUG PRINT FOR VARIABLE ORDER INTEGRATOR
!
!  LET ABS(LPRINT)=  10*N1 + N2     (N1,N2 DIGITS)
!    N1=1   DO NOT PRINT ANY VARIABLES EXTERNAL TO THE INTEGRATOR
!    N1=2   PRINT  TSPECS, CURRENT Y, PAST Y, CURRENT F,
!           ALL PERTINENT CONTENTS OF KORD, AND TOL.
!    N1=3   ABOVE + DIFFERENCE TABLES UP TO HIGHEST DIFFERENCE USED
!    N1=4   SAME AS N1=1 + ALL IN STORAGE ALLOCATED FOR DIFFERENCES
!
!    N2=1   DO NOT PRINT ANY VARIABLES INTERNAL TO THE INTEGRATOR
!    N2=2   PRINT ALL SCALAR VARIABLES IN INTERPOLATION COMMON BLOCK
!    N2=3   ABOVE + ALL SCALAR VARIABLES IN MAIN INTEG. COMMON BLOCK
!    N2=4   SAME AS N1=3 + ALL USED IN ARRAYS XI,BETA,ALPHA, FIRST
!           COLUMN OF G, GS,RBQ,SIGMA
!    N2=5   SAME AS N1=4 + ALL USED IN ARRAYS G,D,DS,V
!
!### History
!  * 2009-11-04 DIVADB Krogh Included TOLG, initilized the unitialized.
!  * 2000-12-01 DIVADB Krogh Removed unused parameter METXTF.
!  * 1996-07-02 DIVADB Krogh Transpose flag for matrix output in C.
!  * 1996-03-25 DIVADB Krogh Introduced TEXT1-TEXT4 to comply with F77.
!  * 1996-01-19 DIVADB Krogh Changed NTEXT to TEXT to agree with doc.
!  * 1995-04-26 DIVADB Krogh Fixed print of V & G's for high order eqs.
!  * 1994-11-11 DIVADB Krogh Declared all vars.
!  * 1994-10-20 DIVADB Krogh Changes to use M77CON
!  * 1994-09-12 DIVADB Krogh Added CHGTYP code.
!  * 1994-03-07 DIVADB Krogh Allow larger order in single precision.
!  * 1993-05-03 DIVADB Krogh Additions for Conversion to C.
!  * 1993-04-14 DIVADB Krogh Changes for new MESS usage.
!  * 1992-04-08 DIVADB Krogh Unused labels 10 and 60 removed.
!  * 1992-03-10 DIVADB Krogh Fixed value for KDIM in single p. version.
!  * 1992-02-17 DIVADB Krogh Made tabs depend on # digits output.
!  * 1991-11-04 DIVADB Krogh Switched to use MESS, DMESS
!  * 1990-03-08 DIVADB Krogh Unused stiff vars. set to 0.
!  * 1989-07-21 DIVADB Krogh Code for integrating discontinuities
!  * 1988-06-07 DIVADB Krogh Dim. of IVC2 and DVC2 upped by 1 (old bug)
!  * 1987-12-07 DIVADB Krogh Initial code.

    subroutine DIVADB(LPRINT, TSPECS, Y, F, KORD, TEXT)

      integer LPRINT, KORD(*)
      character TEXT*(*)
      character TEXT1(1)*11, TEXT2(1)*4, TEXT3(1)*5, TEXT4(1)*4
      integer IVC1(12), IVC2(65), J, K, L, N1, N2
      double precision  DVC2(7), RVC2(8), EVC(8)
!--D Next line special: P=>D, X=>Q
      double precision TSPECS(*), Y(*), TNEQ(1), DVC1(7)
      double precision F(*)
!
!++S Default KDIM = 16
!++  Default KDIM = 20
!++  Default MAXORD = 2, MAXSTF = 1
!++  Default STIFF=.F.
      integer KDIM, MAXORD, MAXSTF
!++ Substitute for KDIM, MAXORD, MAXSTF below
      parameter (KDIM = 20, MAXORD = 2, MAXSTF = 1)
!--D Next line special: P=>D, X=>Q
      double precision TN
      double precision XI(KDIM)
!
!--D Next line special: P=>D, X=>Q
      double precision TG(2), TGSTOP(2), TMARK, TMARKX, TOUT, TOLG
      double precision ALPHA(KDIM), BETA(KDIM+1)
      double precision  D(MAXSTF+MAXORD,MAXORD), G(KDIM,MAXORD)
      double precision V(KDIM+MAXORD)
      double precision HC, HDEC, HINC, HINCC, HMAX, HMAXP9, HMIN
      double precision FDAT(11)
!
      double precision DS(MAXSTF+MAXORD, MAXORD), GS(KDIM)
      double precision SIGMA(KDIM), RBQ(KDIM), DNOISE
      double precision EAVE, EIMAX, EIMIN, EMAX, EREP, ROBND, SNOISE
!
!.    SPECIFICATION OF ENVIRONMENTAL CONSTANTS.
      double precision EEPS10, EEPS16, EROV10, EEPS2
      double precision EEPT75, EOVEP2, OVTM75, OVD10
      common / DIVAEV / EEPS2, EEPT75, EOVEP2, OVTM75, OVD10, EEPS10,   &
     &   EEPS16, EROV10
      save / DIVAEV /
      integer IOPST, KORDI, KQMAXD, KQMAXI, LDT, MAXDIF, MAXINT, NKDKO, &
     &   NTE, NYNY, NDTF, NUMDT
      common / DIVASC / TN, XI, IOPST, KORDI, KQMAXD, KQMAXI, LDT,      &
     &   MAXDIF, MAXINT, NKDKO, NTE, NYNY, NDTF, NUMDT
!
      integer ICF,ICS,IGFLG,IGTYPE(2),IGSTOP(2),ILGREP,INGS,IOP3,IOP4,  &
     &   IOP5,IOP6,IOP7,IOP8,IOP9,IOP10,IOP11,IOP12,IOP13,IOP14,IOP15,  &
     &   IOP16,IOP17,IOP18,IOP19,IOP20,IOP21,IOP22,IOP21S,ITOLEP,IY,    &
     &   KEMAX,KIS,KMARK,KORD1I,KORD2I,KPRED,KQDCON,KQICON,KQMAXS,      &
     &   KQMXDS,KQMXIL,KQMXIP,KQMXIS,KSC,KSOUT,KSSTRT,KSTEP,LEX,LINC,   &
     &   LINCD,LINCQ,LSC,MAXKQD,MAXKQI,METHOD,NE,NEPTOL,NG,NGTOT,       &
     &   NOISEQ,NOUTKO,NTOLF,NY,IDAT(6)
      common /DIVAMC/ TG,TGSTOP,TMARK,TMARKX,TOUT,TOLG,HC,HDEC,HINC,    &
     &   HINCC,HMAX,HMAXP9,HMIN,ALPHA,BETA,D,G,V,DS,GS,SIGMA,RBQ,DNOISE,&
     &   EAVE,EIMAX,EIMIN,EMAX,EREP,ROBND,SNOISE,FDAT,ICF,ICS,IGFLG,    &
     &   IGTYPE,IGSTOP,ILGREP,INGS,IOP3,IOP4,IOP5,IOP6,IOP7,IOP8,IOP9,  &
     &   IOP10,IOP11,IOP12,IOP13,IOP14,IOP15,IOP16,IOP17,IOP18,IOP19,   &
     &   IOP20,IOP21,IOP22,IOP21S,ITOLEP,IY,KEMAX,KIS,KMARK,KORD1I,     &
     &   KORD2I,KPRED,KQDCON,KQICON,KQMAXS,KQMXDS,KQMXIL,KQMXIP,KQMXIS, &
     &   KSC,KSOUT,KSSTRT,KSTEP,LEX,LINC,LINCD,LINCQ,LSC,MAXKQD,MAXKQI, &
     &   METHOD,NE,NEPTOL,NG,NGTOT,NOISEQ,NOUTKO,NTOLF,NY,IDAT
      save / DIVAMC / , / DIVASC /
      equivalence (IVC1(1), IOPST), (IVC2(1), ICF)
      equivalence (TNEQ, TN)
      equivalence (RVC2(1), DNOISE), (DVC1(1), TG), (DVC2(1), HC),      &
     &   (EVC(1), EEPS2)
!
!                      Declarations for error message processing.
      integer MEDDIG, NEDDIG, METDIG, METABS, MERET, METEXT,            &
     &   METABL, MEIVEC, MEFVEC, MEFMAT
      parameter (MEDDIG =12)
      parameter (NEDDIG =-MEDDIG)
      parameter (METDIG =22)
      parameter (METABS =32)
      parameter (MERET  =51)
      parameter (METEXT =53)
      parameter (METABL =55)
      parameter (MEIVEC =57)
      parameter (MEFVEC =61)
      parameter (MEFMAT =62)
      integer MACT0(3), MACT1(2), MACT2(7), MACT3(7), MACT4(8),         &
     &   MACT5(11), MACT6(3), MACT7(14), MACTFV(4)
      integer KPI, KPE
!                      wddtrr        wwddtrr
      parameter (KPI = 400301, KPE = 1305501)
!
! ********* Error message text ***************
![Last 2 letters of Param. name]  [Text generating message.]
!AA $NKORD:    $B
!AB Int. Ord.: $B
!   $
!AC D.E. Ord.: $B
!   $
!AD Meth.Type: $B
!   $
!AE Tolerance Groups: $B
!AF Tolerances: $B
!   $
!AG $NDifferences$B
!AH Eq. $#
!AI Ord. $#
!   $
!AJ $NTN=$F$E
!   $
!AK $NIOPST=$I$TKORDI=$I$TKQMAXD=$I$TKQMAXI=$I$TLDT=$I$T$C
!   MAXDIF=$I$TMAXINT=$I$TNKDKO=$I$TNTE=$I$TNYNY=$I$TNDTF=$I$C
!   $TNUMDT=$I$E
!   $
!AL $NICF=$I$TICS=$I$TIGFLG=$I$TIGTYPE(1)=$I$TIGTYPE(2)=$I$T$C
!   IGSTOP(1)=$I$TIGSTOP(2)=$I$TILGREP=$I$TINGS=$I$TIOP3=$I$T$C
!   IOP4=$I$TIOP5=$I$TIOP6=$I$TIOP7=$I$TIOP8=$I$TIOP9=$I$T$C
!   IOP10=$I$TIOP11=$I$TIOP12=$I$TIOP13=$I$TIOP14=$I$TIOP15=$I$T$C
!   IOP16=$I$TIOP17=$I$TIOP18=$I$TIOP19=$I$TIOP20=$I$TIOP21=$I$T$C
!   IOP22=$I$TIOP21S=$I$TITOLEP=$I$TIY=$I$TKEMAX=$I$TKIS=$I$T$C
!   KMARK=$I$TKORD1I=$I$TKORD2I=$I$TKPRED=$I$TKQDCON=$I$T$C
!   KQICON=$I$TKQMAXS=$I$TKQMXDS=$I$TKQMXIL=$I$TKQMXIP=$I$T$C
!   KQMXIS=$I$TKSC=$I$TKSOUT=$I$TKSSTRT=$I$TKSTEP=$I$TLEX=$I$T$C
!   LINC=$I$TLINCD=$I$TLINCQ=$I$TLSC=$I$TMAXKQD=$I$TMAXKQI=$I$T$C
!   METHOD=$I$TNE=$I$TNEPTOL=$I$TNG=$I$TNGTOT=$I$TNOISEQ=$I$T$C
!   NOUTKO=$I$TNTOLF=$I$TNY=$I$E
!AM $NDNOISE=$F$TEAVE=$F$TEIMAX=$F$TEIMIN=$F$TEMAX=$F$T$C
!   EREP=$F$TROBND=$F$E
!   $
!AN $NTG(1)=$F$TTG(2)=$F$TTGSTOP(1)=$F$TTGSTOP(2)=$F$C
!   $TTMARK=$F$TTMARKX=$F$TTOUT=$F$E
!   $
!AO HC=$F$THDEC=$F$THINC=$F$THINCC=$F$THMAX=$F$T$C
!   HMAXP9=$F$THMIN=$F$T$N$E
!   $
!AP K$HXI(K)$HBETA(K)$HALPHA(K)$HG(K,1)$HRBQ(K)$HSIGMA(K)$H
!   GS(K)$HV(K)$HG(K,2..MAXINT)$E
!   $
!AQ $NEEPS2=$F$TEEPT75=$F$TEOVEP2=$F$TOVTM75=$F$TOVD10=$F$T$C
!   EEPS10=$F$TEEPS16=$F$TEROV10=$F$E
      integer LTXTAA,LTXTAB,LTXTAC,LTXTAD,LTXTAE,LTXTAF,LTXTAG,LTXTAH,  &
     & LTXTAI,LTXTAJ,LTXTAK,LTXTAL,LTXTAM,LTXTAN,LTXTAO,LTXTAP,LTXTAQ
      parameter (LTXTAA=  1,LTXTAB= 14,LTXTAC=  1,LTXTAD=  1,LTXTAE=  1,&
     & LTXTAF= 21,LTXTAG=  1,LTXTAH= 16,LTXTAI= 22,LTXTAJ=  1,          &
     & LTXTAK=  1,LTXTAL=  1,LTXTAM=655,LTXTAN=  1,LTXTAO=  1,          &
     & LTXTAP=  1,LTXTAQ=  1)
      character MTXTAA(1) * (26)
      character MTXTAB(1) * (13)
      character MTXTAC(1) * (13)
      character MTXTAD(1) * (34)
      character MTXTAE(1) * (28)
      character MTXTAF(1) * (9)
      character MTXTAG(1) * (120)
      character MTXTAH(3) * (242)
      character MTXTAI(1) * (80)
      character MTXTAJ(1) * (68)
      character MTXTAK(1) * (84)
      character MTXTAL(1) * (88)
      data MTXTAA/'$NKORD:    $BInt. Ord.: $B'/
      data MTXTAB/'D.E. Ord.: $B'/
      data MTXTAC/'Meth.Type: $B'/
      data MTXTAD/'Tolerance Groups: $BTolerances: $B'/
      data MTXTAE/'$NDifferences$BEq. $#Ord. $#'/
      data MTXTAF/'$NTN=$F$E'/
      data MTXTAG/'$NIOPST=$I$TKORDI=$I$TKQMAXD=$I$TKQMAXI=$I$TLDT=$I$TM&
     &AXDIF=$I$TMAXINT=$I$TNKDKO=$I$TNTE=$I$TNYNY=$I$TNDTF=$I$TNUMDT=$I$&
     &E'/
      data MTXTAH/'$NICF=$I$TICS=$I$TIGFLG=$I$TIGTYPE(1)=$I$TIGTYPE(2)=$&
     &I$TIGSTOP(1)=$I$TIGSTOP(2)=$I$TILGREP=$I$TINGS=$I$TIOP3=$I$TIOP4=$&
     &I$TIOP5=$I$TIOP6=$I$TIOP7=$I$TIOP8=$I$TIOP9=$I$TIOP10=$I$TIOP11=$I&
     &$TIOP12=$I$TIOP13=$I$TIOP14=$I$TIOP15=$I$TIOP16=$I$TIOP17','=$I$TI&
     &OP18=$I$TIOP19=$I$TIOP20=$I$TIOP21=$I$TIOP22=$I$TIOP21S=$I$TITOLEP&
     &=$I$TIY=$I$TKEMAX=$I$TKIS=$I$TKMARK=$I$TKORD1I=$I$TKORD2I=$I$TKPRE&
     &D=$I$TKQDCON=$I$TKQICON=$I$TKQMAXS=$I$TKQMXDS=$I$TKQMXIL=$I$TKQMXI&
     &P=$I$TKQMXIS=$I$TKSC=$I$TKSOUT=$I$TKSS','TRT=$I$TKSTEP=$I$TLEX=$I$&
     &TLINC=$I$TLINCD=$I$TLINCQ=$I$TLSC=$I$TMAXKQD=$I$TMAXKQI=$I$TMETHOD&
     &=$I$TNE=$I$TNEPTOL=$I$TNG=$I$TNGTOT=$I$TNOISEQ=$I$TNOUTKO=$I$TNTOL&
     &F=$I$TNY=$I$E$NDNOISE=$F$TEAVE=$F$TEIMAX=$F$TEIMIN=$F$TEMAX=$F$TER&
     &EP=$F$TROBND=$F$E  '/
      data MTXTAI/'$NTG(1)=$F$TTG(2)=$F$TTGSTOP(1)=$F$TTGSTOP(2)=$F$TTMA&
     &RK=$F$TTMARKX=$F$TTOUT=$F$E'/
      data MTXTAJ/'HC=$F$THDEC=$F$THINC=$F$THINCC=$F$THMAX=$F$THMAXP9=$F&
     &$THMIN=$F$T$N$E'/
      data MTXTAK/'K$HXI(K)$HBETA(K)$HALPHA(K)$HG(K,1)$HRBQ(K)$HSIGMA(K)&
     &$HGS(K)$HV(K)$HG(K,2..MAXINT)$E'/
      data MTXTAL/'$NEEPS2=$F$TEEPT75=$F$TEOVEP2=$F$TOVTM75=$F$TOVD10=$F&
     &$TEEPS10=$F$TEEPS16=$F$TEROV10=$F$E'/
!
      data MACT0 / METABS, 10, MERET /
      data MACT1 / METEXT, MERET /
!                       1       2  3       4       5  6      7
      data MACT2 / METEXT, MEIVEC, 3, METEXT, MEIVEC, 0, MERET /
!                       1       2  3       4       5  6      7
      data MACT3 / METEXT, MEIVEC, 0, METEXT, MEFVEC, 0, MERET /
!                       1       2  3  4  5       6       7
      data MACT4 / METEXT, MEFMAT, 0, 0, 0, LTXTAI, LTXTAH, MERET /
!                       1   2       3       4   5       6  7       8
      data MACT5 / METABS, 12, METEXT, METABS, 18, METDIG, 5, METEXT,   &
     &             METABS, 0, MERET /
!                       9 10     11
      data MACT6 / NEDDIG, 0, MERET /
!                         2 3 4   5 6 7 8 9  10  11  12 13 14
      data MACT7 / METABL,0,0,0,KPI,0,0,0,0,KPE,KPE,KPE, 0, 0 /
!                        1       2  3      4
      data MACTFV / METEXT, MEFVEC, 3, MERET /
!
      data TEXT1 / '$NTSPECS:$B' /
      data TEXT2 / 'Y:$B' /
      data TEXT3 / 'YN:$B' /
      data TEXT4 / 'F:$B' /
!
! ********
! START OF CODE -- PRINT TEXT AND SET INDEX FOR F
! ********
!    Getting variables that are not yet assigned some values.
!++  Code for ~STIFF is active
      KQDCON = 0
      KQMXDS = 0
      MAXKQD = 0
!++  End
      GS(1) = 1.D0
      if (IOP6 == 0) then
        IGTYPE(1) = 0
        IGSTOP(1) = 0
        TG(1) = 0.D0
        TGSTOP(1) = 0.D0
      end if
      if (IOP7 == 0) then
        IGTYPE(2) = 0
        IGSTOP(2) = 0
        TG(2) = 0.D0
        TGSTOP(2) = 0.D0
      end if
      if (IOP6 + IOP7 == 0) then
        INGS = 0
        NG = 0
      end if
      if (IOP10 == 0) NOUTKO = 0
      J = 0
      call MESSFT(MACT0, TEXT)
!
      N1 = LPRINT / 10
      N2 = LPRINT - 10 * N1
      if (N1 <= 1) go to 80
! ********
! PRINT ALL EXTERNAL VARIABLES EXCEPT FOR THE DIFFERENCES
! ********
      MACTFV(3) = max(IOP5, 4)
!--D Next line special: P=>D, X=>Q
      call DMESS(MACTFV, TEXT1, KORD, TSPECS)
      MACTFV(3) = NY
!--D Next line special: P=>D, X=>Q
      call DMESS(MACTFV, TEXT2, KORD, Y)
!--D Next line special: P=>D, X=>Q
      call DMESS(MACTFV, TEXT3, KORD, Y(NYNY))
      MACTFV(3) = NTE
!--D Next line special: P=>S, X=>D
      call DMESS(MACTFV, TEXT4, KORD, F)
      MACT2(6) = NTE
      call MESS(MACT2, MTXTAA, KORD)
      if (NKDKO > 0) call MESS(MACT2(4), MTXTAB, KORD(NKDKO))
      if (IOPST > 0) call MESS(MACT2(4), MTXTAC, KORD(IOPST))
! WRITE TOL
      K = IOP16
   70 if (KORD(K) < 0) K = K + 1
      K = K + 1
      if (KORD(K - 1) < NTE) go to 70
      MACT3(3) = K - IOP16
      MACT3(6) = MACT3(3)
!--D Next line special: P=>S, X=>D
      call DMESS(MACT3, MTXTAD, KORD(IOP16), F(NTOLF))
      if (N1 == 2) go to 80
! ********
! WRITE THE DIFFERENCE TABLES
! ********
      K = NUMDT
      if (N1 == 3) K = KQMAXS
      MACT4(3) = NUMDT
      MACT4(4) = -K
      MACT4(5) = NTE
!--D Next line special: P=>S, X=>D
      call DMESS(MACT4, MTXTAE, KORD, F(NDTF))
!
   80 if (N2 <= 1) return
! ********
! WRITE SCALARS IN COMMON
! ********
!--D Next line special: P=>D, X=>Q
      call DMESS(MACT1, MTXTAF, KORD, TNEQ)
!
! ===== COMMON 1  -- INTEGER
!
      call MESS(MACT1, MTXTAG, IVC1)
      if (N2 == 2) return
      call MESS(MACT6, MTXTAA, IDAT)
      MACT5(10) = MACT6(2) + 14
!
! ===== COMMON 2  -- INTEGER AND FLOATING POINT
!
!--D Next line special: P=>S, X=>D
      call DMESS(MACT5, MTXTAH, IVC2, RVC2)
!--D Next line special: P=>D, X=>Q
      call DMESS(MACT1, MTXTAI, IVC2, DVC1)
!--D Next line special: P=>S, X=>D
      call DMESS(MACT1, MTXTAJ, IVC2, DVC2)
      if (N2 == 3) return
!         wddtrr              wddtrr
      J = 101000 * MACT6(2) + 800501
      MACT7(2) = 1
      MACT7(3) = KQMAXS
      MACT7(4) = 8
      do 90 K = 6, 9
         MACT7(K) = J
   90 continue
      if (N2 > 0) then
         MACT7(4) = 8 + MAXINT
         MACT7(13) = J
         L = min(MAXINT, 4)
         MACT7(14) = J + L - 2
      end if
      do 100 K = 1, MACT7(3)
         FDAT(1) = XI(K)
         FDAT(2) = BETA(K)
         FDAT(3) = ALPHA(K)
         FDAT(4) = G(K, 1)
         FDAT(5) = RBQ(K)
         FDAT(6) = SIGMA(K)
         FDAT(7) = GS(K)
         if (N2 >= 4) then
            FDAT(8) = V(K)
            do 95 J = 2, L
               FDAT(7+J) = G(K, J)
   95       continue
         end if
!--D Next line special: P=>S, X=>D
         call DMESS(MACT7, MTXTAK, IDAT, FDAT)
  100 continue
!++  Code for STIFF is inactive
!     if (MAXDIF <= 0) return
!        Need to define MACT8 and set values
!c--D Next line special: P=>S, X=>D
!     call DMESS(MACT8, 'D$B', IDAT, D)
!c--D Next line special: P=>S, X=>D
!     call DMESS(MACT8, 'DS$B', IDAT, DS)
!++  End
!
!--D Next line special: P=>S, X=>D
      call DMESS(MACT1, MTXTAL, IDAT, EVC)
      return
    end subroutine divadb
!*************************************************************************

!*************************************************************************
!>
!  SUBROUTINE TO LOCATE OUTPUT POINTS AT ZEROS OF ARBITRARY
!  FUNCTIONS  **** GSTOPS **** FOR DIFFERENTIAL EQUATION
!  INTEGRATOR -ODE (OR -IVA).
!
!### History
!  * 2001-09-07 DIVAG  Krogh  Changes to allow user tol on G-Stops.
!  * 1995-06-20 DIVAG  Krogh  Fixed problem introduced with last change.
!  * 1995-05-09 DIVAG  Krogh  Fixed G-Stop/discontinuity code interaction
!  * 1994-11-11 DIVAG  Krogh  Declared all vars.
!  * 1994-10-20 DIVAG  Krogh  Changes to use M77CON
!  * 1994-09-12 DIVAG  Krogh  Added CHGTYP code.
!  * 1994-08-17 DIVAG  Krogh  Modified internal comment.
!  * 1994-03-07 DIVAG  Krogh  Allow larger order in single precision.
!  * 1993-04-27 DIVAG  Krogh  Additions for Conversion to C.
!  * 1992-10-12 DIVAG  Krogh  Fixed G-Stop/discontinuity code interaction
!  * 1992-09-17 DIVAG  Krogh  Slight change in check for sign change.
!  * 1992-04-08 DIVAG  Krogh  Unused labels 140,150,230, and 250 removed.
!  * 1992-03-10 DIVAG  Krogh  Fixed value for KDIM in single p. version.
!  * 1990-01-29 DIVAG  Krogh  Added arg to call to DERMN.
!  * 1988-03-04 DIVAG  Krogh  Initial code.

    subroutine DIVAG(TSPECS, Y, F, KORD, IFLAG, NSTOP, GNEW, GT)

      integer KORD(*), IFLAG, NSTOP
!--D Next line special: P=>D, X=>Q
      double precision TSPECS(*), Y(*), TOLD, TSAVE
      double precision F(*), GNEW(*), GT(*), GOLD
      save GOLD, TOLD, TSAVE
!
!++SP Default KDIM = 16
!++  Default KDIM = 20
!++  Default MAXORD = 2, MAXSTF = 1
      integer KDIM, MAXORD, MAXSTF
!++ Substitute for KDIM, MAXORD, MAXSTF below
      parameter (KDIM = 20, MAXORD = 2, MAXSTF = 1)
!--D Next line special: P=>D, X=>Q
      double precision TN
      double precision XI(KDIM)
!
!--D Next line special: P=>D, X=>Q
      double precision TG(2), TGSTOP(2), TMARK, TMARKX, TOUT, TOLG
      double precision ALPHA(KDIM), BETA(KDIM+1)
      double precision  D(MAXSTF+MAXORD,MAXORD), G(KDIM,MAXORD)
      double precision V(KDIM+MAXORD)
      double precision HC, HDEC, HINC, HINCC, HMAX, HMAXP9, HMIN
      double precision FDAT(11)
!
      double precision DS(MAXSTF+MAXORD, MAXORD), GS(KDIM)
      double precision SIGMA(KDIM), RBQ(KDIM), DNOISE
      double precision EAVE, EIMAX, EIMIN, EMAX, EREP, ROBND, SNOISE
!
!.    SPECIFICATION OF ENVIRONMENTAL CONSTANTS.
      double precision EEPS10, EEPS16, EROV10, EEPS2
      double precision EEPT75, EOVEP2, OVTM75, OVD10
      common / DIVAEV / EEPS2, EEPT75, EOVEP2, OVTM75, OVD10, EEPS10,   &
     &   EEPS16, EROV10
      save / DIVAEV /
      integer IOPST, KORDI, KQMAXD, KQMAXI, LDT, MAXDIF, MAXINT, NKDKO, &
     &   NTE, NYNY, NDTF, NUMDT
      common / DIVASC / TN, XI, IOPST, KORDI, KQMAXD, KQMAXI, LDT,      &
     &   MAXDIF, MAXINT, NKDKO, NTE, NYNY, NDTF, NUMDT
!
      integer ICF,ICS,IGFLG,IGTYPE(2),IGSTOP(2),ILGREP,INGS,IOP3,IOP4,  &
     &   IOP5,IOP6,IOP7,IOP8,IOP9,IOP10,IOP11,IOP12,IOP13,IOP14,IOP15,  &
     &   IOP16,IOP17,IOP18,IOP19,IOP20,IOP21,IOP22,IOP21S,ITOLEP,IY,    &
     &   KEMAX,KIS,KMARK,KORD1I,KORD2I,KPRED,KQDCON,KQICON,KQMAXS,      &
     &   KQMXDS,KQMXIL,KQMXIP,KQMXIS,KSC,KSOUT,KSSTRT,KSTEP,LEX,LINC,   &
     &   LINCD,LINCQ,LSC,MAXKQD,MAXKQI,METHOD,NE,NEPTOL,NG,NGTOT,       &
     &   NOISEQ,NOUTKO,NTOLF,NY,IDAT(6)
      common /DIVAMC/ TG,TGSTOP,TMARK,TMARKX,TOUT,TOLG,HC,HDEC,HINC,    &
     &   HINCC,HMAX,HMAXP9,HMIN,ALPHA,BETA,D,G,V,DS,GS,SIGMA,RBQ,DNOISE,&
     &   EAVE,EIMAX,EIMIN,EMAX,EREP,ROBND,SNOISE,FDAT,ICF,ICS,IGFLG,    &
     &   IGTYPE,IGSTOP,ILGREP,INGS,IOP3,IOP4,IOP5,IOP6,IOP7,IOP8,IOP9,  &
     &   IOP10,IOP11,IOP12,IOP13,IOP14,IOP15,IOP16,IOP17,IOP18,IOP19,   &
     &   IOP20,IOP21,IOP22,IOP21S,ITOLEP,IY,KEMAX,KIS,KMARK,KORD1I,     &
     &   KORD2I,KPRED,KQDCON,KQICON,KQMAXS,KQMXDS,KQMXIL,KQMXIP,KQMXIS, &
     &   KSC,KSOUT,KSSTRT,KSTEP,LEX,LINC,LINCD,LINCQ,LSC,MAXKQD,MAXKQI, &
     &   METHOD,NE,NEPTOL,NG,NGTOT,NOISEQ,NOUTKO,NTOLF,NY,IDAT
      save / DIVAMC / , / DIVASC /
!
      integer I, IG, IGFLGS, IZFLAG, KEXIT, NGSTOP(2)
      double precision HH
      equivalence (G(1,1), HH), (NGSTOP(1), IOP6), (KEXIT, IOP17),      &
     &   (IZFLAG, IY), (IGFLGS, ITOLEP)
!
!                      Declarations for error message processing.
!
      integer MEMDA1, MERET, MEEMES
      parameter (MEMDA1 =27)
      parameter (MERET  =51)
      parameter (MEEMES =52)
      integer MACT(7)
! ********* Error message text ***************
![Last 2 letters of Param. name]  [Text generating message.]
!AA DIVAG$B
!AB Call with bad values of KORD.  KORD(1)=$I, KORD(2)=$I, when $C
!   TSPECS(1)=$F and KSTEP=$M.$E
      integer LTXTAA,LTXTAB
      parameter (LTXTAA= 1,LTXTAB= 8)
      character MTXTAA(1) * (95)
      data MTXTAA/'DIVAG$BCall with bad values of KORD.  KORD(1)=$I, KOR&
     &D(2)=$I, when TSPECS(1)=$F and KSTEP=$M.$E'/
      data MACT / MEMDA1, 0, MEEMES, 68, 24, 0, MERET /
!
! ****************** START OF EXECUTABLE CODE **************
!
      IFLAG = 1
      IG = KORD(2)
      if ((IG /= 0) .and. (IG /= 1)) go to 500
      if (IGFLG - 3) 10, 80, 70
   10 if (IGFLG - 1) 20, 210, 60
!
! ******************** INITIAL POINT ***********************
!
   20 if (IGFLG == -3) then
         IGFLG = 5
         return
      end if
      IGTYPE(IG + 1) = 0
      if ((NGSTOP(IG + 1) <= 0) .or. (IG + IGFLG == -1)) go to 30
      IGFLG = IG - 2
      go to 40
   30 IGFLG = 5
   40 NG = NGSTOP(2 - IG)
      do 50 I = 1, NG
   50    GT(I) = GNEW(I)
      go to 480
!
!     **** USER HAS BEEN TOLD THAT A GSTOP WAS FOUND
!          TEST IF CALLED FROM OUTPUT WHEN CALL SHOULD
!          BE FROM DERIVS
   60 if ((IGTYPE(1) == 0) .and. (IG /= 0)) go to 420
!     **** PROTECT AGAINST NOISEY G NEAR THE ZERO
      if (GNEW(INGS) == 0.D0) GNEW(INGS) = GT(INGS)
!
! ****** TEST FOR CHANGE IN THE SIGN OF A G ****************
!
   70 NG = NGSTOP(2 - IG)
      INGS = 0
   80 INGS = INGS + 1
      if (INGS > NG) if (IGFLG - 4) 400, 380, 480
      if (GNEW(INGS)) 90, 100, 110
   90 if (GT(INGS)) 120, 120, 130
  100 if (GT(INGS)) 130, 80, 130
  110 if (GT(INGS)) 130, 120, 120
  120 GT(INGS) = GNEW(INGS)
      go to 80
!
! ********* A SIGN CHANGE HAS BEEN FOUND *******************
!
  130 NSTOP = INGS
      if (IG == 0) NSTOP = -INGS
      if (IGFLG /= 5) go to 200
!     **** USUAL CASE -- TEST IF OUTPUT POINT PRECEDES THE
!          SIGN CHANGE, OR IF PREDICTING, CORRECTING, OR
!          NOT THROUGH THE FIRST STEP.
!     **** TEST IF AN INTERPOLATED G WAS WHAT CHANGED SIGN
      if (IG /= 0) go to 180
!     **** BACK UP DIFFERENCES AND OTHER STUFF TO BEGINNING
!          OF THE STEP
      call DIVABU(F, KORD)
!     **** TEST IF CORRECTING
      if (KORD1I == 2) go to 170
!     **** TEST IF THROUGH THE FIRST STEP
      if (LSC < 4) go to 180
!     **** IF FIRST DERIVATIVE EVALUATION OF THE FIRST
!          STEP, FIND THE GSTOP, AND USE IT TO GET A NEW
!          INITIAL STEPSIZE
      if (LSC == 7) go to 200
!     **** SET NEW STEPSIZE AFTER SIGN CHANGE WHEN STARTING
  160 HH = TSPECS(1) - TN
!     **** SET KEXIT TO TRY NEW STEPSIZE
  170 KEXIT = 1
      TSPECS(1) = TN
      go to 460
!     **** SET KEXIT FOR USUAL CASE
  180 KEXIT = IG + 2
!     **** TEST IF SIGN CHANGE IN G PRECEDES NEXT OUTPUT PT.
      if (HH * (TSPECS(1) - TMARK)) 200, 200, 190
!     **** SET UP TO EVALUATE G AT OUTPUT POINT
  190 IGFLG = 4
      TSPECS(1) = TMARK
      NSTOP = 0
      go to 240
!
! ***************** FIND THE ZERO OF G *********************
!
!     **** INITIALIZE ZERO FINDER
  200 TOLD = TG(2 - IG)
      GOLD = GT(INGS)
      TSAVE = TSPECS(1)
      IGFLGS = IGFLG
      IGFLG = 1
      IZFLAG = 0
      go to 220
!     **** TEST IF ZERO ALREADY FOUND
  210 if (IZFLAG - 1) 350, 220, 310
  220 continue
      call DZERO(TSPECS(1), GNEW(INGS), TOLD, GOLD, IZFLAG, TOLG)
!     **** TEST FOR CONVERGENCE
      if (IZFLAG /= 1) go to 260
!     **** INTERPOLATE NEW Y, AND GO COMPUTE G AGAIN
  240 call DIVAIN(TSPECS(1), Y, F, KORD)
      IFLAG = 4
      KORD2I = IG - 3
      return
!     **** CONVERGENCE -- CHOOSE TOLD TO GIVE A CHANGE
!          IN SIGN
  260 if (GNEW(INGS) == 0.D0) go to 290
      if (TSPECS(1) - TOLD) 270, 300, 280
  270 if (HH) 290, 300, 300
  280 if (HH) 300, 300, 290
  290 TOLD = TSPECS(1)
!     **** CHECK IF SIGN CHANGE DUE TO NOISE
  300 TSPECS(1) = TOLD + XI(1)
      go to 240
  310 TSPECS(1) = TOLD
      if (GNEW(INGS)) 320, 340, 330
  320 if (GT(INGS)) 340, 360, 360
  330 if (GT(INGS)) 360, 360, 340
!     **** ZERO WAS EVIDENTLY DUE TO NOISE
  340 TSPECS(1) = TSAVE
      IZFLAG = 0
      go to 370
  350 IGFLG = IGFLGS
!     SET KORD2I TO INITIAL VALUE TO AVOID LOOP
      KORD2I = IG
      go to 80
!     **** SAVE INFORMATION ABOUT THE STOP
  360 IGFLG = 3
      TGSTOP(2 - IG) = TSPECS(1)
      IGTYPE(2 - IG) = IZFLAG + 3
      IGSTOP(2 - IG) = NSTOP
  370 NSTOP = 0
      go to 240
!
! ************** AFTER SEARCH FOR A SIGN CHANGE ************
!
!     **** NO SIGN CHANGE AT A T OUTPUT POINT
!     TEST IF CALLED FROM OUTPUT
  380 if (IG /= 0) go to 390
!     SET UP FOR CALL TO OUTPUT
      KORD1I = 7
      KORD2I = -2
      IFLAG = 3
      go to 480
!     **** ADJUST KEXIT AND SET UP TO GIVE OUTPUT
  390 KEXIT = KEXIT + 2
      KORD1I = min(KMARK, 5)
      KORD(3) = KMARK
      KORD(1) = KORD1I
      IGFLG = 5
      IFLAG = 2
      go to 470
!     **** TEST IF USER HAS BEEN TOLD OF GSTOP
  400 if (IGFLG == 2) go to 450
!     **** A GSTOP HAS BEEN FOUND
!     TEST IF STARTING
      if (LSC == 7) go to 160
      IFLAG = IGTYPE(2 - IG)
      NSTOP = IGSTOP(2 - IG)
      INGS = abs(NSTOP)
      if (INGS == 0) go to 410
      GT(INGS) = -GT(INGS)
      if (IG == 0) go to 430
      IGFLG = 2
!     If interpolated GSTOP was found set to check again in case of
!     multiple stops at exactly the same point.
      if (IGTYPE(1) /= 0) go to 440
!     **** TELL USER OF AN EXTRAPOLATED GSTOP
  410 IFLAG = IGTYPE(2)
      NSTOP = IGSTOP(2)
      INGS = abs(NSTOP)
!     **** SET SO DERIVS IS CALLED WITH KORD(1) = KPRED
  420 KORD1I = KPRED
      KORD2I = -3
      return
!     **** AN EXTRAPOLATED GSTOP WAS FOUND, SET UP TO CHECK
!          INTERPOLATED STOPS (IF ANY)
  430 NG = NGSTOP(1)
      INGS = 0
      IFLAG = 3
      NSTOP = IGSTOP(2)
!     **** SET SO OUTPUT IS CALLED WITH KORD(1) = 7
  440 KORD1I = 7
      KORD2I = -2
      go to 490
!     **** CHECK THAT AN EXTRAPOLATED G-STOP IS NOT MISSED
  450 if ((IG == 0) .or. (IGTYPE(2) == 0)) go to 460
!     SET TO CHECK FOR INTERPOLATED G-S.
      TG(1) = TSPECS(1)
      IGTYPE(1) = 0
      TSPECS(1) = TGSTOP(2)
      INGS = 0
      IGFLG = 3
      go to 240
!     **** SET SO INTEGRATOR GOES TO PLACE DIRECTED BY KEXIT
  460 NSTOP = 0
      IGFLG = 5
      IFLAG = 3
  470 KORD2I = -7
!     **** STORE INFO. ON LAST G COMPUTED
  480 IGTYPE(2 - IG) = 0
  490 TG(2 - IG) = TSPECS(1)
      return
!
! ********************** ERROR PROCESSING ******************
!
  500 MACT(2) = KSTEP
!--D Next line special: P=>S, X=>D
      call DMESS(MACT, MTXTAA, KORD, TSPECS)
      IFLAG = 8
      KEXIT = 6
      go to 470
    end subroutine divag
!*************************************************************************

!*************************************************************************
!>
! Processes Messages -- Actions are controlled by MACT().  See
! comment is subroutine MESS.  This program is for the extra
! argument of type real.
!
! BUF    In common CMESSC, see MESS.
! DOLS   In common for intitialization, not used here.  See MESS.
! EUNIT  In common for intitialization, not used here.  See MESS.
! FDAT   Formal argument -- gives floating point data to print.
! FBIG   Largest magnitude of floating point number to output.
! FSMA   Smalles magnitude floating point number to be printed.
! FMTF   In common CMESSC, format for printing floating point number.
! FMTG   In common CMESSC, user format to use in place of FMTF.
! FMTSP  Format for printing sparse vectors.
! FOUT   Floating point number to be output.
! FSMA   Smallest postitive floating point number.
! ICOL   In common CMESSI, see MESS.
! ID     Number of decimal digits for floating point format statement.
! IDAT   Integer data -- passed to MESS.
! IVAR   In common CMESSI, see MESS.
! IWF    In common CMESSI, see MESS.
! IWG    In common CMESSI, see MESS.
! J      Temporary index.
! K      Temporary index.
! KSMA   Number of leading 0's in "F" format.  If < 0, -KSMA gives the
!    number of extra digits to the left of the decimal point.
!    KSMA depends on abs(smallest number printed).
! KBIG   Number of extra digits before the decimal point required by the
!    largest number to be printed if "F" format is used.
! KDF    In common CMESSI, see MESS.
! KDFDEF In common CMESSI, see MESS.
! KDIAG  In common CMESSI, not used here, see MESS.
! KEXE   Extra space required for E format.
! KF     In common CMESSI, see MESS.
! KLINE  In common CMESSI, see MESS.
! KSCRN  In common CMESSI, see MESS.
! KRES1  In common CMESSI, see MESS.
! KSPEC  In common CMESSI, see MESS.
! LASTER In common CMESSI, not used here, see MESS.
! LASTI  In common CMESSI, see MESS.
! LBUF   In common CMESSI, see MESS.
! LDFDEF Value of KDFDEF for this routine.  (Saved)
! LENBUF In common CMESSI, see MESS.
! LENLIN In common CMESSI, not used here, see MESS.
! LENTRY In common CMESSI, see MESS.
! LHEAD  In common CMESSI, not used here, see MESS.
! LINERR In common CMESSI, not used here, see MESS.
! LINMSG In common CMESSI, not used here, see MESS.
! LOCBEG In common CMESSI, see MESS.
! LPRINT In common CMESSI, not used here, see MESS.
! LSTOP  In common CMESSI, not used here, see MESS.
! LSTRT  In common CMESSI, see MESS.
! MACT   Formal argument, see MESS.
! MDAT   In common CMESSI, not used here, see MESS.
! MEMDA5 In common CMESSI, see MESS.
! MESS   Program called for most of the message processing.
! MPT    In common CMESSI, see MESS.
! MUNIT  In common CMESSI, not used here, see MESS.
! NCOL   In common CMESSI, see MESS.
! NDIM   In common CMESSI, see MESS.
! NEG    1 if any number is < 0, else it is 0.
! NFDAT  In common CMESSI, see MESS.
! NIDAT  In common CMESSI, not used here, see MESS.
! NMDAT  In common CMESSI, not used here, see MESS.
! NTEXT  In common CMESSI, not used here, see MESS.
! OUNIT  In common CMESSI, not used here, see MESS.
! D1MACH External func. giving floating pt. info. about the environment.
! SUNIT  In common CMESSI, not used here, see MESS.
! TEXT   Formal argument, passed to MESS, see there.
! XARGOK In common CMESSI, see MESS.
!
!%% void dmessxc(long int);
!
!### History
!  * 2009-09-27 DMESS Krogh  Same as below, in another place.
!  * 2009-07-23 DMESS Krogh  Changed ,1x to :1x in write to FMTF.
!  * 2008-06-13 DMESS Krogh  Changed -0's to 0.
!  * 2007-09-08 DMESS Krogh  Fixed definitions of MEVLAS.
!  * 2006-10-08 DMESS Krogh  Another try, see 2005-05-26
!  * 2006-10-08 DMESS Krogh  Fixed minor problem in matrix/vector output.
!  * 2006-10-01 DMESS Krogh  Print NaN's and infity (at least for g77).
!  * 2006-07-01 DMESS Krogh  messxc => dmessxc (and not static) (for C)
!  * 2006-04-07 DMESS Krogh  Major rewrite of code for F.Pt. format.
!  * 2006-04-04 DMESS Krogh  Fixes in C code for vectors & matrices.
!  * 2006-04-02 DMESS Krogh  Added code for output of sparse vector.
!  * 2005-07-10 DMESS Krogh  Small adjustment for last correction.
!  * 2005-05-26 DMESS Krogh  Fixed "*****" output in boundary case.
!  * 2002-05-16 DMESS Krogh  Added way for user to get error count.
!  * 2002-03-27 DMESS Krogh  Fixed crash when number is -INF.
!  * 2001-06-08 DMESS Krogh  Eliminated Hollerith in formats.
!  * 2001-05-25 DMESS Krogh  Added a couple of commas in formats.
!  * 1997-06-17 DMESS Krogh  In C code made messxc, static.
!  * 1996-07-12 DMESS Krogh  Changes to use .C. and C%%.
!  * 1996-03-30 DMESS Krogh  Added external statement.
!  * 1994-10-20 DMESS Krogh  Changes to use M77CON
!  * 1994-09-21 DMESS Krogh  Added CHGTYP code.
!  * 1994-09-08 DMESS Krogh  Added new matrix/vector capabilities.
!  * 1994-08-17 DMESS Krogh  Removed duplicate save statement.
!  * 1994-04-19 DMESS Krogh  Removed blank line from DMESS.
!  * 1993-05-14 DMESS Krogh  Changed TEXT to array of character strings.
!  * 1993-04-14 DMESS Krogh  Fixes for conversion to C. (C%% comments.)
!  * 1992-07-12 DMESS Krogh  Fixed so negative KDFDEF works.
!  * 1992-05-27 DMESS Krogh  Initialized LDFDEF in a data statement.
!  * 1992-05-14 DMESS Krogh  Put common blocks in save statement.
!  * 1992-04-28 DMESS Krogh  Corrected minor error in floating pt. format
!  * 1992-02-28 DMESS Krogh  Initial Code.

    subroutine DMESS (MACT, TEXT, IDAT, FDAT)

      !external         D1MACH
      integer          MACT(*), IDAT(*)
      double precision FDAT(*)
      character        TEXT(*)*(*)
      character        FMTSP*29
      integer          ICOL, ID, J, K, KSMA, KBIG, KEXE, LDFDEF, NEG
      double precision FBIG, FOUT, FSMA !, D1MACH
      save LDFDEF, FMTSP
      save /CMESSI/, /CMESSC/
!++ CODE for .C. is inactive
!      integer  kciwid, kccwid, kcrwid, lbeg, lend, lfprec, lgprec
!      common /MESSCC/ kciwid,kccwid,kcrwid,lbeg,lend,lfprec,lgprec
!++ END
!
! ************************** Data from common block ********************
!
! For comments on these variables, see the listing for MESS.
!
      integer   LENBUF, MEVBAS, MEVLAS
      parameter (LENBUF=250)
      parameter (MEVBAS=10)
      parameter (MEVLAS=33)
      logical          XARG, GOTFMT, XARGOK
      integer          ERRCNT, EUNIT, ICHAR0, IRC, IVAR(MEVBAS:MEVLAS), &
     &   IMAG, INC, ITEXT, IWF, IWG, KDF, KDFDEF, KDI, KDIAG, KDJ,      &
     &   KLINE, KSCRN, KSHIFT, KSPEC, KT, MAXERR, LASTI, LBUF, LENLIN,  &
     &   LENOUT, LENTRY, LENTXT, LHEAD, LINERR, LINMSG, LOCBEG, LPRINT, &
     &   LSTOP, LSTRT, LTEXT, MAXWID(2), MDAT(5), MPT, MUNIT, NCOL,     &
     &   NDIM, NFDAT, NIDAT, NMDAT, NROW, NTEXT, OUNIT, SUNIT, TABSPA
!
      character BUF*(LENBUF), DOLS*72, FMTC*7, FMTF*20, FMTG*15,        &
     &  FMTI*7, FMTIM(2)*7, FMTJ*7, FMTR*7, FMTT*15
      common /CMESSI/ SUNIT, LHEAD, KDFDEF, LINMSG, LINERR, MUNIT,      &
     &   EUNIT, KSCRN, KDIAG, MAXERR, LSTOP, LPRINT, KDF, NTEXT, NIDAT, &
     &   NFDAT, NMDAT, MDAT, TABSPA, ERRCNT, ICHAR0, IMAG, INC, IRC,    &
     &   ITEXT, IWF, IWG, KDI, KDJ, KLINE, KSHIFT, KSPEC, KT, LASTI,    &
     &   LBUF, LENLIN, LENOUT, LENTRY, LENTXT, LOCBEG, LSTRT, LTEXT,    &
     &   MAXWID, MPT, NROW, NCOL, NDIM, OUNIT, GOTFMT, XARG, XARGOK
      common /CMESSC / BUF, DOLS, FMTF, FMTG, FMTI, FMTJ, FMTT, FMTIM
      equivalence (IVAR(MEVBAS), SUNIT)
      equivalence (FMTIM(1), FMTR), (FMTIM(2), FMTC)
!
      data LDFDEF / 0 /
!

! ************************* Start of Executable Code *******************
!
      XARGOK = .true.
      if (LDFDEF == 0) then
         LDFDEF = 1 - int(log10(d1mach(3)))
      end if
      KDFDEF = LDFDEF
      KDF = KDFDEF
   10 call MESS (MACT, TEXT, IDAT)
!             4    5    6    7    8    9   10   11
      go to (20, 100, 100, 200, 300, 400, 100, 500), LENTRY-3
      XARGOK = .false.
      LDFDEF = KDFDEF
      return
!                                      Print from FDAT
   20 J = LBUF + 1
      FOUT = FDAT(NFDAT)
      NFDAT = NFDAT + 1
      if (KSPEC >= 8) then
         LBUF = LBUF + IWG
!%% messcc.lend = cmessi.lbuf;
!%% cmessc.buf[messcc.lend] = ' ';
!%% if ((j > 1) && (cmessc.buf[j-2] >= '0') &&
!%%    (cmessc.buf[j-2] <= '9')) cmessc.buf[j++ - 1] = ' ';
!%% sprintf(&cmessc.buf[j-1], cmessc.fmtg, cmessi.iwg,
!%%    messcc.lgprec, fout);
!%% if (cmessc.buf[messcc.lend] != 0) {messcc.lbeg=j; dmessxc(kexe);}
         write (BUF(J:LBUF), FMTG) FOUT
         go to 10
      end if
      if (FOUT <= 0.D0) then
        if (FOUT == 0.D0) then
          FDAT(NFDAT-1) = 0.D0
          FOUT = 0.D0
        else
          NEG = 1
        end if
      else if (FOUT > 0.D0) then
        NEG = 0
      else
!               Must be a Nan
        NEG = 0
        FBIG = 1.0
        FSMA = 1.0
        IWF = 2
        go to 40
      end if

      FBIG = abs(FOUT)
      FSMA = FBIG
      IWF = 2
!                                      Get the format.
   40 continue
      if (KDF == 0) KDF = KDFDEF
      KEXE = 0
      if (FBIG /= 0.D0) then
        if (FSMA == 0.D0) FSMA = 1.D0
        FBIG = FBIG * (1.D0 + .5D0 * .1D0**abs(KDF))
        IWF = IWF + NEG
        if (KDF < 0) then
          KSMA = 0
        else
          KSMA = -log10(FSMA)
          if (FSMA < 1.D0) KSMA = KSMA + 1
        end if
        KBIG = log10(FBIG)
        if (FBIG < 1.D0) then
          KBIG = KBIG - 1
          if (FBIG > 1.D0 - .1D0**abs(KDF-1)) KBIG = 0
        end if
!         This is to get ininities out (at least with g77)
        if ((KBIG < -1000) .or. (KBIG > 1000)) KBIG = 8
        if ((KSMA < -1000) .or. (KSMA > 1000)) KSMA = 8
        if (max(KBIG, 0) + max(KSMA,0) >= 4) then
!               Want to use an "E" format
          KEXE = 3 + max(0, int(log10(dble(max(KBIG,abs(KSMA))+1.D-5))))
          if (KDF < 0) then
            ID = -KDF
          else
            ID = KDF - 1
          end if
          IWF = IWF + ID + KEXE
!++ CODE for ~.C. is active
          if (LENTRY == 10) IWF = IWF - 1
          write (FMTF, '(''(1P,99(E'',I2,''.'',I2,''E'',I1,'':1X))'')') &
     &      IWF, ID, KEXE - 2
!++ CODE for .C. is inactive
!c WATCOM C and others (?) man need an extra 1 here??
!%%    strcpy(cmessc.fmtf, "%*.*E ");
!      lfprec = id
!++ END
          go to 60
        end if
      else
        KSMA = 1
        KBIG = -1
      end if
!               Want to use an "F" format
      if (KDF < 0) then
        ID = -KDF
      else
        ID = KDF + KSMA - 1
      end if
!++ CODE for ~.C. is active
      IWF = IWF + ID + max(KBIG, -1)
      write (FMTF, '(''(0P,99(F'',I2,''.'',I2,'':1X))'')') IWF,ID
!++ CODE for .C. is inactive
!      IWF = IWF + ID + max(KBIG, 0)
!%%    strcpy(cmessc.fmtf, "%*.*f ");
!      lfprec = id
!++ END
   60 if (LENTRY /= 4) then
        IWF = IWF + 1
        if (LENTRY /= 10) go to 10
!               Take care of setup for sparse vector
        IMAG = 0
        do 70 J = LOCBEG, LASTI
          IMAG = max(abs(IMAG), IDAT(J))
   70   continue
        call MESSFI
!  Format forms:     12345678901234567890   123456789012345678  1234567
!                    (1P,99(Edd.ddEd:1X))   (0P,99(Fxx.xx:1X))  (99Idd)
!++ CODE for ~.C. is active
        if (FMTF(8:8) == 'F') then
          FMTSP=                                                        &
     &      '(99(' // FMTI(4:6) // ''') '',0P,' // FMTF(8:18)
        else
          FMTSP=                                                        &
     &     '(99(' // FMTI(4:6) // ''')'' ,1P,' // FMTF(8:20)
        end if
!++ CODE for .C. is inactive
!c Using cmessc.fmtf in place of fmtsp
!%%      if (cmessc.fmtf[4] == 'f') {
!%%      strcpy(cmessc.fmtf, "%*ld) %*.*f ");
!%%      kexe = 0;
!%%      }
!%%      else strcpy(cmessc.fmtf, "%*ld) %*.*E ");
!%%      cmessi.iwf++;
!++ END
        IWF = IWF + KDI + 1
        go to 10
      end if
!
      LBUF = LBUF + IWF
!%% messcc.lend = cmessi.lbuf;
!%% cmessc.buf[messcc.lend] = ' ';
!%% if ((j > 1) && (cmessc.buf[j-2] >= '0') &&
!%%    (cmessc.buf[j-2] <= '9')) cmessc.buf[j++ - 1] = ' ';
!%% sprintf(&cmessc.buf[j-1], cmessc.fmtf, cmessi.iwf,
!%%    messcc.lfprec, fout);
!%% if (cmessc.buf[messcc.lend] != 0) {messcc.lbeg=j; dmessxc(kexe);}
      write (BUF(J:LBUF),FMTF) FOUT
      go to 10
!                                     Get format for a vector or matrix
  100 ICOL = 1
      if (FDAT(LOCBEG) < 0.D0) then
        NEG = 1
      else if (FDAT(LOCBEG) >= 0.D0) then
        NEG = 0
      else
!               Must be a Nan
        NEG = 0
        FBIG = 1.0
        FSMA = 1.0
        go to 110
      end if

      FBIG = abs(FDAT(LOCBEG))
      FSMA = FBIG
  110 do 120 J = LOCBEG, LASTI, INC
        if (FDAT(J) <= 0.D0) then
          if (FDAT(J) == 0.D0) then
            FDAT(J) = 0.D0
          else
            NEG = 1
          end if
        end if
        FBIG = max(FBIG, abs(FDAT(J)))
        if (FSMA == 0.D0) then
          FSMA = abs(FDAT(J))
        else if (FDAT(J) /= 0.D0) then
          FSMA = min(FSMA, abs(FDAT(J)))
        end if
  120 continue
      if (NCOL /= 0) then
         ICOL = ICOL + 1
         LOCBEG = LOCBEG + NDIM
         LASTI = LASTI + NDIM
         if (ICOL <= NCOL) go to 110
      end if
      IWF = 2
      go to 40
!                                    Floating point vector output
  200 continue
!%% messcc.lend = cmessi.lstrt-1;
!%% neg = 0;
!%% for (j=cmessi.mpt; j<cmessi.mpt+cmessi.kline; j++){
!%%   messcc.lbeg = messcc.lend;
!%%   messcc.lend = messcc.lend + cmessi.iwf;
!%%   if (kexe) {
!%%     if (kexe == 5)
!%%        neg = ((fabs(fdat[cmessi.inc*j-1]) < 1.e100)) ? -1: 0;
!%%     else if (kexe ==3) neg = 1;
!%%   }
!%%   sprintf(&cmessc.buf[messcc.lbeg], cmessc.fmtf,
!%%    cmessi.iwf+neg, messcc.lfprec, fdat[cmessi.inc*j-1]);
!%%   if ((kexe == 3) || ((kexe == 5) && neg)) dmessxc(kexe); }
      write(BUF(LSTRT:LBUF),FMTF)(FDAT(K),K=MPT,MPT+INC*(KLINE-1),INC)


!      print '(/A/)', BUF(1:LBUF)

      MPT = MPT + KLINE * INC
      go to 10
!                                    Floating point matrix output
  300 continue
!%% messcc.lend = cmessi.lstrt-1;
!%% neg = 0;
!%% for (j=cmessi.mpt; j<=cmessi.lasti; j+=cmessi.ndim){
!%%    messcc.lbeg = messcc.lend;
!%%    messcc.lend = messcc.lend + cmessi.iwf;
!%%   if (kexe) {
!%%     if (kexe == 5) neg = ((fabs(fdat[j-1]) < 1.e100)) ? -1: 0;
!%%     else if (kexe ==3) neg = 1;
!%%   }
!C%%    if ((messcc.lbeg > 1) && (cmessc.buf[messcc.lbeg-1] >= '0') &&
!C%%       (cmessc.buf[messcc.lbeg-1] <= '9'))
!C%%       cmessc.buf[messcc.lbeg++] = ' ';
!%%    sprintf(&cmessc.buf[messcc.lbeg],
!%%       cmessc.fmtf, cmessi.iwf+neg, messcc.lfprec, fdat[j-1]);
!%%    if ((kexe == 3) || ((kexe == 5) && neg)) dmessxc(kexe); }
      write (BUF(LSTRT:LBUF), FMTF) (FDAT(K), K = MPT, LASTI, NDIM)
      go to 10
!                                    Table output
  400 continue
!%% messcc.lend = cmessi.lstrt-1;
!%% neg=0;
!%% for (j=cmessi.mpt; j<cmessi.mpt+cmessi.kline; j++){
!%%    messcc.lbeg = messcc.lend;
!%%    messcc.lend = messcc.lend + cmessi.iwf;
!C%%    if ((messcc.lbeg > 1) && (cmessc.buf[messcc.lbeg-1] >= '0') &&
!C%%       (cmessc.buf[messcc.lbeg-1] <= '9'))
!C%%       cmessc.buf[messcc.lbeg++] = ' ';
!%%    sprintf(&cmessc.buf[messcc.lbeg], cmessc.fmtf,
!%%       cmessi.iwf+neg, messcc.lfprec, fdat[j-1]); }
      write (BUF(LSTRT:LBUF), FMTF) (FDAT(K), K = MPT, MPT+KLINE-1)
      go to 10


!                                   Sparse vector output
  500 continue
!%%  messcc.lend = -1;
!%%  neg = 0;
!%%  for (j=cmessi.mpt; j<cmessi.mpt+cmessi.kline; j++) {
!%%    messcc.lbeg = messcc.lend + 1;
!%%     messcc.lend = messcc.lend + cmessi.iwf;
!%%   if (kexe) {
!%%     if (kexe == 5) neg = ((fabs(fdat[j-1]) < 1.e100)) ? -1: 0;
!%%     else if (kexe == 3) neg = 1;
!%%   }
!%%   sprintf(&cmessc.buf[messcc.lbeg], cmessc.fmtf, cmessi.kdi,
!%%     idat[j-1],cmessi.iwf-cmessi.kdi-2+neg,messcc.lfprec,fdat[j-1]);
!%%     if ((kexe == 3) || ((kexe == 5) && neg)) dmessxc(kexe); }
      write (BUF(1:LBUF), FMTSP) (IDAT(K),FDAT(K),K=MPT,MPT+KLINE-1)
      MPT = MPT + KLINE
      go to 10

    end subroutine dmess
!*************************************************************************

!%%  void dmessxc(long int kexe)
!%%{
!%%  /* Adjust for lack of control on digits in exponent */
!%%  char c;
!%% if (cmessc.fmtf[4] == 'f') return;
!%% if (kexe == 4) return;
!%% if (kexe == 3) { // Should only be one digit in the exponent
!%%   cmessc.buf[messcc.lend-1] = cmessc.buf[messcc.lend];
!%%   cmessc.buf[messcc.lend] = ' ';
!%% }
!%% else { // Should be at least 3 digits in the exponent.
!%%   c =cmessc.buf[messcc.lend-4];
!%%   if ((c < '0') || (c > '9')) {
!%%     cmessc.buf[messcc.lend-1] = cmessc.buf[messcc.lend-2];
!%%     cmessc.buf[messcc.lend-2] = cmessc.buf[messcc.lend-3];
!%%     cmessc.buf[messcc.lend-3] = '0';
!%%     cmessc.buf[messcc.lend] = ' ';
!%%   }
!%% }
!%% return;
!%%} /* end of function */

!*************************************************************************
!>
! SUBROUTINE TO FIND A BOUNDED ZERO
!
! analysis and coding by Fred T.Krogh at the Jet Propulsion
! Laboratory, Pasadena, Calif.  April 25, 1972.
! Modified for portability, April 1984 by Krogh.
! Algorithmic changes, vars. added to save stmt., Sept. 1987 by Krogh
!
! Parameters in the calling sequence are defined as follows:
!
!  X1  = independent variable
!  F1  = dependent variable --  initially   F1=F(X1).
!        When MODE=1 (or 5) the user is to compute F(X1) given X1
!  X2  = second value of independent variable
!  F2  = F(X2) on the initial entry.  When MODE = 2-4, F2=F(X2) and
!        F1*F2 <= 0.
!  MODE  is a parameter used for communication between this
!        subroutine and the user. (The user should set MODE
!        only to initialize it to 0 before the first call)
!      =1  compute F(X1) and call $ZERO
!      =2  F(X1) is approximately 0, the iteration is finished
!          and the error criterion is satisfied.
!      =3  same as MODE=2, except the error criterion can
!          not be satisfied.
!      =4  apparently the function has a discontinuity
!          between X1 and X2 -- No zero can be found
!      =5  F1*F2 was greater than zero on the first call, and an attempt
!          to bound the zero on both sides have failed.
!      =6  fatal error -- $ZERO was called after mode was set >=2.
!          If $ZERO is called again, the program will be stopped.
!          (Unless MODE is set to 0)
!      <0  If MODE is set <0 and $ZERO is called, no action is taken
!          except that print is turned on for -MODE calls to $ZERO.
!          This print gives all values of X and F used in the iteration.
!  TOL    is the error tolerance
!     TOL>0  Iterate until values of X1 and X2 are known
!              for which abs(X1-X2) <= tol and F1*F2 <= 0.
!     TOL<0  Iterate until a value of X1 is found for which
!              abs(F1) <= abs(TOL).
!     TOL  = 0  Iterate until the zero is determined as
!              precisely as possible.  MODE = 3 is impossible
!              in this case.
!
! Usage is as follows (of course, variations are possible.)
!         Somehow one has available X1, F1, X2, and F2 such
!         that F1 = F(X1), F2 = F(X2) and F1*F2 <= 0.
!         In addition, one should assign a value to TOL.
!     MODE = 0
!***  In the statement below, $ is replaced by an 'S' for single
!***  precision and a 'D' for double.
! XXX call $ZERO(X1,F1,X2,F2,MODE,TOL)
!     go to  (N1,N2,N3,N4,N5,N6), MODE
!  N1 COMPUTE  F1=F(X1)
!     go to XXX
!
!  N4 continue
!  N5 continue
!  N6 stop
!  N3 If you want to -- print results to note that error
!                       is too big.
!  N2 zero is found, do whatever you want to with it.
!
! End of comments explaining usage.
!
! ************************* Usage of internal variables ****************
!
! C0     Parameter = 0.
! C1     Parameter = 1.
! C1P25  Parameter = 1.25
! C2     Parameter = 2.
! C4     Parameter = 4.
! CP01   Parameter = 0.01
! CP125  Parameter = 1.25
! CP25   Parameter = 0.25
! CP5    Parameter = 0.5
! CP75   Parameter = 0.75
! CP99   Parameter = 0.99
! D1MACH Gets constants associated with floating point arithmetic.
! DFDXXX = (XXMXL/FFMFL) * (est. deriv. of f w.r.t. x at x = XX).  All
!   derivatives are base on a second degree polynonial that interpolates
!   the last three points generated.
! DFDXXX = (XXMXO/FFMFL) * (est. deriv. of f w.r.t. x at x = X0).
! DIV    If artificial subdivision of the interval is used, determines
!   the amount of the sudivision.  (-XXMXOO * DIV / (1. + DIV))
! DMESS  Prints error messages.
! DXDFFF = (FFMFL/XXMXL) * (est. deriv. of x w.r.t. f at f = FF).
! DXDFFO = (FFMFO/XXMXL) * (est. deriv. of x w.r.t. f at f = FO).
! F1     (formal arg.) The last value of F computed, on return the value
!   of F(X1).
! F2     (formal arg.) The other initial value provided for F.  Set to
!   the value of F(X2) on returns.
! FDAT   Temporary storage for floating point values for messages.
! FF     Value of F1 after F is computed.
! FFDFO  FF / FO
! FFMFB  FFMFL + FLMFB = FF - value of FF 2 iterations back.
! FFMFL  FF - FL
! FL     Value of FF from the previous iteration.
! FLMFB  Value of FFMFL from the previous iteration
! FO     F(XO)
! I      Comments for LCHNG define how I is set.
! IDAT   Temporary storage for integer values for messages.
! J      This is 1 if FF <= 0., and is 2 if FF > 0.
! KNKP   KNKP(J) (see J above) is decreased by 3 when there are signs of
!   decent convergence.  It is counted up when convergence is slow.
! KS     =-1  initially,
!        = 0  whenever F changes sign, otherwise
!        = number of times the sign of F has remained the same
! KTYP   = 1 if interpolation was used to get the last iterate, = 0 if
!   an artificial subdivision was used.
! LCHG  the J-th continuation in the data statement for LCHG below gives
! new states for the case when the state number is J-1.  State 0 is the
! initial state.  The I-th entry on a row gives the state for case on I
! as follows:  (TP is the ratio (new f) / (last f of same sign)
!    I = 1   TP < 0.01
!    I = 2   .01 <= TP < 1
!    I = 3   TP = 1
!    I = 4   1 < TP <= 4
!    I = 5   TP > 4.
! States are as follows:
!    0   initial state, or big increase, or small increase in state 0
!    1   after big decrease, perhaps followed by small decreases
!    2   after one or more small decreases from state 0
!    3   one or more small increases from state 2
!    4   one or more small decreases from state 3
!    5   decision made that noise is a problem on this side of zero.
! LINIT  - the maximum number of iterations that can be taken with out
!   getting a sign change in F.
! LMODE  The value of MODE the last time in this routine.
! LNLP
! LTXTxx Names of this form are used in setting up data statements for
!   error messages.  These names are generated automatically by PMESS,
!   the program that makes up these messages.
! MACT   This array difines the actions to be taken by the error message
!   program.  See comments in MESS for details.  MODE is set to MACT(3)
!   on exit.
! MACT1  As for MACT except used for the diagnostic print.
! MExxxx Parameters defining constants used for interaction with the
!   error message program MESS.  See comments there for definitions.
! MLOC   Contains locations in MTXTAA for error messages.
! MODE   (formal) See comments above.
! MTXTAA Text for error messages.
! MTXTAB Text for diagnostic message.
! MTXTAC Text for diagnostic message.
! NP     If > 0, gives number of iterations till diagnostic print stops.
! QFM
! QXM
! RND    Largest relative difference between succesive floating point
!   numbers.
! SMALL  .5 / (RND * largest floating point number)
! TOL    (Formal) See description above.
! TOLX   Actually tolerance required for accuracy in X.  Usually =
!   max(TOL, XRND).  It can be cut by a factor of 2 for use in setting
!   bounds on an acceptable interval.
! TP     Ordinarily the ratio (FF / prev. FF of the same sign.
! TP1    Used for temporary storage.
! X1     (Formal) Value of x where F is to be computed, and value
!   returned for the zero after convergence.
! X2     (Formal) Initially other value of x where F is given.  After
!   convergence gives the other closest x which gives an F of opposite
!   sign from that given by x1.
! XL     Value of XX from the previous iteration.
! XLMXB  Value of XXMXL from the previous iteration.
! XO     Value of x on opposite side of the zero from the current x.
! XRND   Best accuracy that one could hope for based on the finite
!   precision of floating point numbers.
! XX     Current x, the last value of X1 where F was computed.
! XXMXL  XX - XL
! XXMXO  XX - XO = length of interval in which 0 lies.
! XXMXOL Value of XXMXO from a previous iteration.
!
!### History
!  * 2010-04-14 DZERO  Krogh  No discontinuity message if |F1-F2| small.
!  * 2010-04-12 DZERO  Krogh  Fixed KNKP to get discontinuity diagnostic.
!  * 2010-02-20 DZERO  Krogh  $G => $F for print out of iterations
!  * 2008-03-01 DZERO  Krogh  Minor change in diagnostic print.
!  * 2000-12-01 DZERO  Krogh  Removed unused variable C1P01.
!  * 1998-11-01 DZERO  Krogh  Set so errors stop less easily.
!  * 1998-11-01 DZERO  Krogh  For "mangle", INDIC replaced with MACT(3).
!  * 1996-03-30 DZERO  Krogh  Added external statement.
!  * 1995-11-09 DZERO  Krogh  Fixed so char. data at col. 72 is not ' '.
!  * 1994-11-11 DZERO  Krogh  Declared all vars.
!  * 1994-10-20 DZERO  Krogh  Changes to use M77CON
!  * 1994-09-08 DZERO  Krogh  Added CHGTYP code.
!  * 1993-04-27 DZERO  Krogh  Additions for Conversion to C.
!  * 1993-04-13 DZERO  Krogh  Minor change for new MESS.
!  * 1992-04-08 DZERO  Krogh  Unused label 400 removed.
!  * 1992-01-09 DZERO  Krogh  Moved calc. of XXMXO up (for error msg.)
!  * 1991-11-26 DZERO  Krogh  Converted to new error processor.
!  * 1988-08-14 DZERO  Krogh  Labels runumbered.
!  * 1988-03-07 DZERO  Krogh  Initial code.

    subroutine DZERO(X1, F1, X2, F2, MODE, TOL)

      integer MODE
      double precision X1, X2, F1, F2, TOL

      !external D1MACH
      integer LINIT, KS, KTYP, J, I
      parameter (LINIT = -40)
      integer KNKP(2), LCHG(30), LMODE, LNLP(2), NP
      double precision XX, XO, XL, FF, FO, FL, FFDFO
      double precision DIV, QFM, QXM, TP, TP1, XXMXO, XXMXOL
      double precision RND, XRND, SMALL, TOLX
      double precision XXMXL, XLMXB, FFMFL, FFMFB, FLMFB
      double precision DXDFFF, DXDFFO, DFDXXX, DFDXXO
      double precision C0, C1, C2, C4, CP125, CP25, CP5, CP75, C1P25
      double precision C8, CP01, CP99, CP001, C1P031
      !double precision D1MACH
!
      parameter (C0 = 0.D0, C1 = 1.D0, C2 = 2.D0, C4 = 4.D0)
      parameter (C8 = 8.D0)
      parameter (CP125 = 0.125D0, CP25 = 0.25D0, CP75 = 0.75D0)
      parameter (CP5 = 0.5D0)
      parameter (C1P25 = 1.25D0)
      parameter (CP01 = 0.01D0)
      parameter (CP001 = 0.001D0)
      parameter (CP99 = 0.99D0)
      parameter (C1P031 = 1.03125D0)
!
!                      Declarations for error message processing.
!
      integer MERET, MEEMES, METEXT
      double precision FDAT(4)
      integer MACT(5), MACT1(2), MLOC(4), IDAT(2)
      save DIV, FL, FLMFB, FO, KNKP, KS, KTYP, LCHG, LMODE,             &
     &   LNLP, MACT, NP, RND, SMALL, XL, XLMXB, XO, XX, XXMXOL
      parameter (MERET  =51)
      parameter (MEEMES =52)
      parameter (METEXT =53)
! ********* Error message text ***************
![Last 2 letters of Param. name]  [Text generating message.]
!AA DZERO$B
!AB Best bound for zero is [$F, $F], but tolerance is $F.$E
!AC Apparent discontinuity in function near X = $F.$E
!AD Can not find a sign change: X1=$F, X2=$F, F1=$F, F2=$F$E
!AE Called with MODE = $I.$E
!   $
!AF In DZERO -- X1=$F F1=$F KTYP=$I DIV=$G KS=$I$E
!   $
!AG             X2=$F F2=$F$E
      integer LTXTAA,LTXTAB,LTXTAC,LTXTAD,LTXTAE,LTXTAF,LTXTAG
      parameter (LTXTAA=  1,LTXTAB=  8,LTXTAC= 63,LTXTAD=112,LTXTAE=169,&
     & LTXTAF=  1,LTXTAG=  1)
      character MTXTAA(1) * (193)
      character MTXTAB(1) * (46)
      character MTXTAC(1) * (25)
      data MTXTAA/'DZERO$BBest bound for zero is [$F, $F], but tolerance&
     & is $F.$EApparent discontinuity in function near X = $F.$ECan not$&
     & find a sign change: X1=$F, X2=$F, F1=$F, F2=$F$ECalled with MODE$&
     & = $I.$E'/
      data MTXTAB/'In DZERO -- X1=$F F1=$F KTYP=$I DIV=$F KS=$I$E'/
      data MTXTAC/'            X2=$F F2=$F$E'/
! **** End of automatically generated text
!                      1  2  3  4      5
      data MACT / MEEMES, 0, 0, 0, MERET /
      data MACT1 / METEXT, MERET /
      data MLOC / LTXTAB, LTXTAC, LTXTAD, LTXTAE /
!
      data RND / C0 /
      data KS, KTYP, LMODE, DIV / 0, 0, 2, C0 /
      data LCHG /                                                       &
     &   1, 2, 0, 0, 0,                                                 &
     &   1, 1, 4, 5, 0,                                                 &
     &   1, 2, 3, 3, 0,                                                 &
     &   1, 4, 4, 3, 0,                                                 &
     &   1, 4, 5, 5, 0,                                                 &
     &   1, 5, 5, 5, 0 /
      data NP / 0 /

!
! INITIALIZE
!
      if (MODE < 0) then
         NP = -1 - MODE
         return
      end if
      if (NP > 0) then
         NP = NP - 1
         FDAT(1) = X1
         FDAT(2) = F1
         FDAT(3) = DIV
         IDAT(1) = KTYP
         IDAT(2) = KS
         call DMESS(MACT1, MTXTAB, IDAT, FDAT)
         if (MODE /= 0) if (LMODE - 1) 70, 80, 450
         FDAT(1) = X2
         FDAT(2) = F2
         call DMESS(MACT1, MTXTAC, IDAT, FDAT)
      else if (MODE /= 0) then
         if (LMODE - 1) 70, 80, 450
      end if
!
      if (RND == C0) then
         RND = D1MACH(4)
         SMALL = CP5 / (RND * D1MACH(2))
      end if
      XL = X2
      FL = F2
   30 TP = C1
      MODE = 1
      MACT(3) = 2
      XXMXOL = C0
      KNKP(1) = 0
      KNKP(2) = 0
      LNLP(1) = 0
      LNLP(2) = 0
      KS = -1
      XX = X1
      FF = F1
      if (FL) 40, 75, 50
   40 if (FF) 60, 230, 100
   50 if (FF) 100, 230, 60
   60 LMODE = 0
!             Take care of points on same side of zero.
   70 FF = F1
      XX = X1
      TP = FF / FL
      if (TP < C0) go to 30
      LMODE = LMODE - 1
      if (LMODE < LINIT) then
         MACT(3) = 5
         FDAT(1) = XX
         FDAT(2) = XL
         FDAT(3) = FF
         FDAT(4) = FL
         go to 250
      end if
      if (TP > C1) then
         FF = FL
         XX = XL
         FL = F1
         XL = X1
      end if
      if (abs(FF) >= C8 * abs(FL-FF)) then
         TP = C8
      else
         TP = max(-CP25*dble(LMODE), FF / (FL - FF))
      end if
      FL = FF
      XO = XL
      XL = XX
      if (XX == XO) XO = C1P031 * XX + sign(CP001, XX)
      XX = XX + TP * (XX - XO)
      X1 = XX
      MODE = 1
      return
!
   75 X1 = XL
      F1 = FL
      go to 250
! END OF INITIALIZATION
!
!
! ENTRY AFTER COMPUTING F FOR THE LAST ITERATE
   80 FF = F1
      TP = FF / FL
      if (TP) 90, 230, 110
   90 TP = FF / FO
      KS = 0
  100 FO = FL
      XO = XL
      go to 120
  110 KS = KS + 1
  120 J = 1
      if (FF > C0) J = 2
      if (TP - C1) 150, 140, 130
  130 I = 4
      if (TP > C4) I = 5
      go to 160
  140 I = 3
      go to 160
  150 I = 2
      if (TP < CP01) I = 1
      if (TP < CP99) go to 170
  160 KNKP(J) = KNKP(J) + 1
      go to 180
  170 KNKP(J) = 0
  180 XXMXO = XX - XO
      LNLP(J) = LCHG(5*LNLP(J) + I)
      if (LNLP(J) >= 4) then
         if (LNLP(3 - J) >= 4) go to 210
      end if
! XXMXO GIVES THE LENGTH OF THE INTERVAL INSIDE WHICH
! THE ZERO IS KNOWN TO LIE.
      if (C2 * abs(XXMXO) < abs(XXMXOL)) then
         KNKP(J) = max(0, KNKP(1) - 3)
      end if
      XXMXOL = XXMXO
      XRND = RND * (abs(XX) + abs(XO) + SMALL)
!
! TEST FOR CONVERGENCE
      if (TOL) 190, 200, 200
  190 continue
      if (abs(FF) <= abs(TOL)) go to 220
  200 continue
      TOLX = max(TOL, XRND)
      if (abs(XXMXO) > TOLX) go to 310
!
! CONVERGENCE -- PREPARE FOR FINAL EXIT
  210 if ((abs(XXMXO) > TOL) .and. (TOL /= C0)) then
         MACT(3) = 3
         FDAT(3) = TOL
         if (XXMXO > 0) then
            FDAT(2) = XX
            FDAT(1) = XO
         else
            FDAT(1) = XX
            FDAT(2) = XO
         end if
      end if
! SET FINAL VALUES FOR X1,F1,X2,AND F2
  220 continue
      if (abs(FF) <= abs(FO)) go to 240
      F1 = FO
      X1 = XO
  230 FO = FF
      XO = XX
  240 X2 = XO
      F2 = FO
! TEST FOR DISCONTINUITY
      if ((KNKP(1) > 5) .or. (KNKP(2) > 5)) then
        if (abs(F1 - F2) > RND * max(X1, 1.D0)) then
          MACT(3) = 4
          FDAT(1) = XX
        end if
      end if
  250 MODE = MACT(3)
      if (MACT(3) - 2) 420, 420, 430
! END OF CODE FOR FINAL EXIT
!
! F NOT DECREASING (OR THE FIRST ITERATE)
! PREPARE TO DIVIDE THE INTERVAL
  260 TP = C1
      if (KS) 370, 280, 270
  270 if (KTYP == 0) go to 290
  280 DIV = C2
  290 continue
      DIV = max(DIV, FFDFO)
! KTYP=0 IF AND ONLY IF THE INTERVAL WAS DIVIDED (USING DIV)
! ON THE LAST ITERATION
      if (KTYP == 0) DIV = DIV * (C1P25 / (C1P25 - TP))
! DIVIDE THE INTERVAL AS SPECIFIED BY DIV
  300 TP1 = -XXMXO * (DIV/(DIV+C1))
      KTYP = 0
      go to 410
!
  310 continue
      XXMXL = XX - XL
      FFMFL = FF - FL
      FFDFO = abs(FF / FO)
      TOLX = CP5 * TOLX
      if (TP >= C1) go to 260
! DIVIDE THE INTERVAL IF F HAS HAD THE SAME SIGN FOR
! FOUR OR MORE TIMES IN SUCCESSION
      if (KS - 4) 320, 340, 290
  320 continue
      if (FLMFB == C0) go to 340
! BEGINNING OF CODE TO DETERMINE IF INVERSE QUADRATIC
! INTERPOLATION IS TO BE USED.
      FFMFB = FFMFL + FLMFB
      if (FFMFB == C0) go to 330
      QFM = C1 - (FFMFL / FLMFB) * (XLMXB / XXMXL)
      QXM = C1 - (XXMXL / XLMXB) * (FLMFB / FFMFL)
      DXDFFF = C1 + (FFMFL / FFMFB) * QFM
      DXDFFO = DXDFFF + C2 * ((FO - FF) / FFMFB) * QFM
      TP1 = XXMXL + XLMXB
      DFDXXX = C1 + (XXMXL / TP1) * QXM
      DFDXXO = DFDXXX + C2 * ((XO - XX) / TP1) * QXM
      TP1 = DXDFFF * DFDXXX
      if ((TP1 <= CP25) .or. (TP1 >= C4)) go to 330
      TP1 = DXDFFO * DFDXXO
      if ((TP1 > CP25) .and. (TP1 < C4)) go to 380
!
! DERIVATIVES DO NOT MATCH WELL ENOUGH
  330 continue
      if (KS == 0) if (FFDFO - C1) 350, 370, 360
  340 continue
      if ((KTYP == 0) .and. (TP >= CP75)) go to 290
      continue
      TP = C1 - TP
      if (TP <= FFDFO) go to 280
      FFDFO = FFDFO / TP
      DIV = CP125
      go to 290
  350 continue
      DIV = CP5 * max(max(CP25, FFDFO), TP / (C1P25 - min(TP, C1)))
      go to 300
  360 continue
      DIV = min(C4, CP5 * FFDFO)
      go to 300
! INTERPOLATE WITH SECANT METHOD
  370 TP1 = -XXMXL
      go to 390
!
! DERIVATIVES MATCH UP PRETTY WELL.
  380 continue
! INTERPOLATE USING THE INVERSE QUADRATIC
      TP1 = XXMXL * (QFM * (FL / FFMFB) - C1)
  390 TP1 = (FF/FFMFL) * TP1
      KTYP = 1
!
! EXIT TO GET F(X)
  410 continue
      FL = FF
      FLMFB = FFMFL
      XLMXB = XXMXL
      XL = XX
! COMPUTE X1, INSURING THAT IT IS NOT TOO CLOSE TO THE
! ENDS OF THE INTERVAL
      XX = min(max(XL + TP1, min(XL, XO) + TOLX), max(XL, XO) - TOLX)
      X1 = XX
  420 LMODE = MODE
      return
!
  430 MACT(2) = 11*MACT(3)  - 19
  440 MACT(4) = MLOC(MACT(3)-2)
      call DMESS(MACT, MTXTAA, IDAT, FDAT)
      go to 420
!
! A CALL TO THE SUBROUTINE HAS BEEN MADE WITH MODE/=1
  450 IDAT(1) = MODE
      MACT(3) = 6
      MODE = 6
      if (LMODE /= 6) go to 430
      MACT(2) = 99
      go to 440
    end subroutine dzero
!*************************************************************************

!++ CODE for .C. is inactive
!%% static FILE *c_handle[2], *scratch_file;
!%% static char *c_fname[2]={"MESSF-xx", "MESSF-xx"};
!%% char *ctmp;
!++ END

!*************************************************************************
!>
! Processes Messages -- Actions are controlled by MACT().
!
! This routine is intended for use primarily by other library routines.
! Users of library routines may want to use values of MACT from MERET-
! MESUNI, and may have an interest in using it to print messages
! from their own software.
! This routine has companion routines that are called with the same
! three arguments, plus one additional argument.  This argument is
! referred to here as FDAT since actions specified here can result
! in returns to print data from FDAT.  The name FDAT is used because
! this other routine will ordinarily print floating point data, but
! it could also print other kinds of data, e.g. logical.  At present
! only SMESS and DMESS are defined which are for single and double
! precision floating point data.
! MACT is a vector specifying sequentially the actions desired.
! Some of these actions require more than one location, in which
! case the next action follows the last datum required by the
! previous action.  Internal variables together with default
! values in parentheses which are used to keep track of locations
! are as follows:
!  NTEXT  (1)   The next text output starts at TEXT(NTEXT).
!  NIDAT  (1)   The next output from IDAT starts at IDAT(NIDAT).
!  NFDAT  (1)   The next output from FDAT starts at FDAT(NFDAT).
!  NMDAT  (1)   The next output from MDAT starts at MDAT(NMDAT), where
!               MDAT is defined by actions MEMDA1-MEMDA5 below, and
!               NMDAT is set to one at end of every text output.
! An action which uses data pointed to by one of the above will cause
! the pointer to be incremented to one past the last location used.  An
! exception is NMDAT which when it reaches 5 is not incremented and the
! value pointed to is incremented instead.
! Actions are encoded by values starting in MACT(1) as follows.
! (Arguments required are given in parentheses at the start of
! description.  These arguments follow the action index.  The next
! action follows the last argument for the preceding action.  Action
! indices have been selected so that it is easy to add new functionality
! without affecting codes using an earlier version.  Where bounds are
! indicated for an argument, if the argument is outside the bounds it is
! treated as if it had the value for the bound violated.)
! MESUNI=10  (0 <= K10 <= 99) Set the unit to use for a scratch
!            file.  The default unit for a scratch file is 30.  If a
!            scratch file is needed, (only needed here if a table
!            exceeds the line length), and unit 30 can not be opened as
!            a new scratch file, then units 29, 28, ..., will be tried
!            until an acceptable unit is found.  Library routines may
!            use this file but must be sure that the use does not
!            conflict with the printing of tables here, or the use by
!            any other library routines.  If K10 is 0, a scratch unit is
!            assumed not to be available, and tables with long lines
!            will be printed with each line on multiple lines.
! MEHEAD=11  (0 <= K11 <= 1) Defines the print that surrounds an
!            error message.  K11=0 gives nothing, and 1 gives the first
!            4 characters in TEXT repeated 18 times.  If this is not
!            used, one gets 72 $'s.  (To get a blank line use 1 with
!            TEXT = '    '.)
! MEDDIG=12  (-50 <= K12 <= 50) Set default digits to print for
!            floating point.  If K12 > 0 then K12 significant digits
!            will be printed, if K12 < 0, then -K12 digits will be
!            printed after the decimal point, and if K12 = 0, the
!            default will be used, which is the full machine precision.
!            Setting or getting this value will only work properly if
!            the action is taken by calling SMESS or DMESS as
!            appropriate.
! MEMLIN=13  (39 <= K13 <= 500) Set message line length to K13.
!            (Default is 128.)
! MEELIN=14  (39 <= K14 <= 500) Set error message line length to
!            K14. (Default is 79)
! MEMUNI=15  (-99 <= K15 <= 99) Messages go to unit K15.  If K15 = 0
!            (default), 'print' is used.  If K15 < 0, messages go to
!            both 'print' and to unit abs(K15).  If a write can not be
!            done to unit abs(K15), this unit will be opened with file
!            name MESS_Fxx.tmp, where xx is the value of abs(K15).
! MEEUNI=16  (-99 <= K16 <= 99) As for MEMUNI, except for Error
!            Messages.
! MESCRN=17  (0 <= K17 <= 100000000) Set number of lines to print to
!            standard output before pausing for "go" from user.  Default
!            is 0, which never stops.
! MEDIAG=18  (0 <= K18 <= 1000000000) Set the diagnostic level
!            desired.  This routine makes no use of K18.  It merely
!            serves as a place to set it and to answer inquiries on its
!            value.  It is intended to be set by users of library
!            software.  Library packages that make use of this number
!            are expected to use it as described below.  If K18 = 0 (the
!            default), no diagnostic print is being requested.  Else m =
!            mod(K18, 256) determines whether a package will do
!            diagnostic printing.  Associated with a library package is
!            a number L which must be a power of 2 < 129, and which
!            should be mentioned in the documentation for the package.
!            If the bit logical or(m,L) = L then diagnostic output for
!            the routine with the associated value of L is activated.
!            The value of L should have been selected by the following
!            somewhat vague rules.  Let base 2 log(L) = 2*i + j, where j
!            is 0 or 1.  Select i = level of the library package, where
!            the level is 0 if no other library routine that is likely
!            to be used with the package could reasonably be expected to
!            want any embedded diagnostics, and otherwise is
!            min(4, I+1), where I is the maximum level for any library
!            routine which is likely to be used with the package.
!            Select j = 0 if the user is relatively unlikely to want
!            diagnostics, and j = 1, if this is a routine for which
!            considering its level the user is relatively likely to want
!            diagnostic output.  The next 8 bits, mod(K18/256, 256), may
!            be used by the library routine to select the actual output
!            that is to be given.  These bits may be ignored,  but if
!            they are used, the lowest order bits should correspond to
!            less voluminous output that is more likely to be requested.
!            Finally, K18 / (2**16) may be used to give a count on how
!            many times to print the diagnostics that are given.  This
!            count may be interpreted by library routines in slightly
!            different ways, but when used it should serve to turn off
!            all output after a certain limit is reached.  By
!            convention, if this is 0 there is no upper bound on the
!            count.
! MEMAXE=19  (0 <= K19 <= 1000000000) Set the maximum error value.
!            When retrieving this value, it is the maximum value seen
!            for 10000*s + 1000*p + i, where s, p, and i are the stop
!            and print levels, and the index on the last error message
!            processed, respectively.  See MEEMES below.
! MESTOP=20  (0 <= K20 <= 8) Set the stop level for error messages.
!            If an error message has a stop index > min(K20, 8), the
!            program is stopped after processing the message.  The
!            default value is K20=3.
! MEPRNT=21  (0 <= K21 <= 8) Set the print level for error messages.
!            If an error message has a print index > K21, or the message
!            is going to stop when finished, information in an error
!            message is processed, else all the actions including
!            printing are skipped.  (MESTOP controls stopping.)  The
!            default value is MEPRNT = 3.
! An action index of -i, for i < METDIG, will return in the location
! ordinarily used for Ki the current default value for the internal
! variable set by Ki.  In the case of MESUNI, if the scratch unit has
! not been opened, it will be opened before returning the unit number.
!
! METDIG=22  (-50 <= K22 <= 50) As for MEDDIG, except the value here
!            is temporary, lasting until the return, or next use of this
!            action.  If 0, the internal value for K12 is used instead.
! MENTXT=23  (1 <= K23 <= 10000000) Set value of NTEXT to K23.
! MEIDAT=24  (1 <= K24 <= 1000000000) Set value of NIDAT to K24.
! MEFDAT=25  (1 <= K25 <= 1000000000) Set value of NFDAT to K25.
! MEMDAT=26  (1 <= K26 <= 5) Set value of NMDAT to K26.
! MEMDA1=27  (K27) set MDAT(1) to K27.  See description of NMDAT above.
! MEMDA2=28  (K28) set MDAT(2) to K28.
! MEMDA3=29  (K29) set MDAT(3) to K29.
! MEMDA4=30  (K30) set MDAT(4) to K30.
! MEMDA5=31  (K31) set MDAT(5) to K31.
! METABS=32  (1 <= K32 <= 100) set spacing for tabs to K32.
! MECONT=50  Exit, but no print of current print buffer.  The error or
!            diagnostic message is to be continued immediately.
! MERET=51   All done with diagnostic or error message, complete
!            processing and return, or for some error messages stop.
! MEEMES=52  (K52, L52, M52) Start an error message with severity level
!            K52,index for the error of L52, and message text starting
!            at TEXT(M52).  If M52 is 0, message text starts at
!            TEXT(NTEXT), and if M52 < 0, no message text is
!            printed as part of this action.  Library routines should
!            set K52 = 10*s + p, where s is the stop level desired, and
!            p the print level, and should have 10 > p >= s >= 0.
!            We offer the following guidelines as a yardstick for
!            setting the value of s.
!   = 9  User has ignored warning that program was going to be stopped.
!   = 8  Program has no way to continue.
!   = 7  User has given no indication of knowing that functionality of
!        results is reduced.  (E.g. not enough space for some result.)
!   = 6  Program could continue but with reduced functionality.
!   = 5  Results far worse than user expected to want.
!   = 4  User has given no indication of knowing that results do not
!        meet requested or expected accuracy.
!   = 3  Warning is given that program will be stopped without some
!        kind of response from the calling program.
!   = 2  Program is not delivering requested or expected accuracy.
!   = 1  Some kind of problem that user could correct with more care in
!        coding or in problem formulation.
!   = 0  Message is for information of uncritical nature.
!            Print levels might be counted down so that warnings given
!            several times are no longer given, or be counted up so
!            that a warning is only given after a certain threshold is
!            reached.  Levels should be selected with the understanding
!            that the default is to print only levels above 3.
! METEXT=53  Print TEXT, starting at TEXT(NTEXT).  Print ends
!            with the last character preceding the first '$'.  Special
!            actions are determined by the character following the '$'.
!            Except as noted, the '$' and the single character which
!            follows are not printed.  In the text below, "to continue",
!            means to continue print of TEXT with the next character
!            until the next "$".  Except for the one case noted, NTEXT
!            is set to point to the second character after the "$".
!            Note, letters must be in upper case.  Possibilities are:
!      B  Break text, but don't start a new line.
!      E  End of text and line.
!      R  Break text, don't change the value of NTEXT.  Thus next
!         text Repeats the current.
!      N  Start a New line, and continue.
!      I  Print IDAT(NIDAT), set NIDAT=NIDAT+1, and continue.
!      J  As for I above, except use the last integer format
!         defined by a "$(", see below.
!      F  Print FDAT(NFDAT), set NFDAT=NFDAT+1, and continue.
!      G  As for F above, except use the last floating format
!         defined by a "$(", see below.
!      M  Print MDAT(NMDAT), set NMDAT=NMDAT+1, and continue.
!      H  Marks terminator for column and row Headings, see table,
!         vector, and matrix output below.  This causes enough blanks to
!         be generated to keep column headings centered over their
!         columns.  After the blanks are generated, text is continued
!         until the next '$'.  This is not to be used except inside
!         column or row headings.  The last row or column should be
!         terminated with a '$E' or if appropriate, a '$#' for a row or
!         column label.
!      (  Starts the definition of a format for integer or floating
!         point output.  The format may not contain a "P" field, and
!         must require no more than 12 characters for floating point
!         (e.g. "(nnEww.ddEe)", where each of the lower case letters
!         represents a single digit), and no more than 7 characters for
!         integer output.  Following the ")" that ends the format, if
!         the next character is not a "$" then "$J" or "$G" type output
!         is done, see above.  In either case processing of TEXT then
!         continues.
!      T  Tab.
!      #  Used in matrix row or column labels this prints the current
!         row or column index, respectively, ends the text for the
!         current row or column, and resets the text pointer to where
!         it started.
!      $  a single '$' is printed, continue till the next '$'.
!      -  Start a negative number for skipping.
!     0-9 Digits for skipping.
!      C  Only used by PMESS which deletes it and the preceding '$'.
!         Used at the end of a line to indicate continued text.
!   other Don't use this -- the '$' is ignored, but new features may
!         change the action.  (E.g. $P might be added to get a prompt.)
! ME????=54  Not used.
! METABL=55  (K55, L55, M55, N55)  Note this action automatically
!            returns when done, further locations in MACT are not
!            examined.  This action prints a heading and/or data that
!            follows a heading.  If K55 is 1, then the heading text
!            starting in TEXT(NTEXT) is printed.  This text
!            should contain embedded "$H"'s to terminate columns of the
!            heading.  If there is no heading on a column, use " $H".
!            Note the leading blank.  If the heading is to continue
!            over k columns, begin the text with "$H" repeated k-1
!            times with no other embedded characters.  The very last
!            column must be terminated with "$E" rather than "$H".
!            After tabular data is printed, K55 is incremented by 1,
!            and compared with L55.  If K55 > L55, K55 is reset to 1,
!            and if the data that was to be printed had lines that were
!            too long, data saved in the scratch file is printed using
!            the headings for the columns that would not fit on the
!            first pass.  Note that only one line of tabular data can
!            be printed on one call to this subroutine.
!            M55 gives the number of columns of data associated with the
!            heading.
!            N55 is a vector containing M55 entries.  The k-th integer
!            in N55 defines the printing action for the k-th column
!            of the table.  Let such an integer have a value defined by
!            rr + 100 * (t + 10 * (dd + 100 * ww)), i.e. wwddtrr, where
!            0 <= rr,dd,ww < 100, and 0 <= t < 10.
!      rr    The number of items to print.
!      t     The type of output.
!            1  Print text starting at TEXT(NTEXT), rr = 01.
!            2  Print the value of K55, rr = 01.
!            3  Print integers starting at IDAT(NIDAT).
!            4  Print starting at FDAT(NFDAT), using an F format.
!            5  Print starting at FDAT(NFDAT), using an E format.
!            6  Print starting at FDAT(NFDAT), using an G format.
!      dd    Number of digits after the decimal point.
!      ww    The total number of column positions used by the column,
!            including the space used to separate this column from the
!            preceding one.  This must be big enough so that the column
!            headings will fit without overlap.
! MEIVEC=57  (K57) Print IDAT as a vector with K57 entries.  The vector
!            output starts on the current line even if the current line
!            contains text.  This is useful for labeling the vector.
!            The vector starts at IDAT(NIDAT).
!            If K57 < 0,  indices printed in the labels for the vector
!            start at at NIDAT, and entries from NIDAT to -K57 are
!            printed.
! MEIMAT=58  (K58, L58, M58, I58, J58) Print IDAT as a matrix with K58
!            declared rows, L58 actual rows, and M58 columns.  If K58<0,
!            instead of using 1 for the increment between rows, and K58
!            for the increment between columns, -K58 is used for the
!            increment between rows, and 1 is used for the increment
!            between columns.  If L58<0, the number of actual rows is
!            mod(-L58, 100000), and the starting row index is -L58 /
!            100000.  Similarly for M58<0. TEXT(I58) starts the text for
!            printing row labels.  If I58 < 0, no row labels are
!            printed.  If I58 = 0, it is as if it pointed to text
!            containing "Row $E".  Any "$" in a row or column label must
!            be followed by "H" or "E" which terminates the text for the
!            label.  In the case of $H, text for the next label follows
!            immediately, in the case of $E the current row index is
!            printed in place of the $E and the next label uses the same
!            text.  J58 is treated similarly to I58, except for column
!            labels, and with "Row $E" replaced with "Col $E".  The
!            matrix starts at IDAT(NIDAT), and NIDAT points one past the
!            end of the matrix when finished.
! MEJVEC=59  (K59) As for MEIVEC, except use format set with $(.
! MEJMAT=60  (K60, L60, M60, I60, J60) As for MEIMAT, except use the
!            format set with $(.
! MEFVEC=61  (K61) As for MEIVEC, except print FDAT as a vector with
!            K61 entries.  The vector starts at FDAT(NFDAT).
! MEFMAT=62  (K62, L62, M62, I62, J62) As for action MEIMAT, but
!            instead print FDAT, and use NFDAT in place of NIDAT.
! MEGVEC=63  (K63) As for MEFVEC, except use format set with $(.
! MEGMAT=64  (K64, L64, M64, I64, J64) As for MEIMAT, except use the
!            format set with $(.
! MEIVCI=65  (K65, L65) As for MEIVEC, except the vector entries have a
!            spacing of K65, and there are L65 entries in the vector.
! MEJVCI=66  (K66) As for MEIVCI, except use format set with $(.
! MEFVCI=67  (K67, L67) As for MEFVEC, except the vector entries have a
!            spacing of K67, and there are L67 entries in the vector.
! MEGVCI=68  (K68) As for MEFVCI, except use format set with $(.
! MEFSPV=69  (K69) Output IDAT, FDAT as a sparse vector.
!
!
! ************************** Internal Variables ************************
!
! BUF    Character string holding characters to be output.
! C      Used for temp. storage of a character.
! DOLS   A heading/trailing string of characters, default = $'s.
! ERMSG  Default part of error message.
! ERRCNT Used to keep a count of error messages.
! EUNIT  Unit number for output of error messages.
! FDAT   Formal array, containing floating point data to output.  Only
!   appears external to this subroutine.
! FIRST  Set = .true. initially, then .false. after MESS is called.
! FMTC   Format for integer output of matrix column headings.
! FMTF   Format for floating point or other output.
! FMTG   Format set by user for floating point or other output.
! FMTI   Character string holding format for integer output.
! FMTIM  Equivalenced to FMTR, FMTC.
! FMTJ   Format set by user for integer output.
! FMTR   Value of FMTI for format of row indices.
! FMTT   Format to be stored in FMTJ or FMTG.
! GETW   Set true if still need to get width for a format.
! GOTFMT Set .true. if format has been set by user.  When printing
!   tables, set true when heading has been output.
! I      Index of current action from MACT.
! ICHAR0 Value of ICHAR('0')
! ICOL   Current column index in matrix output.
! IDAT   Formal array, containing integer data to output.
! IMAG   Magnitude of integer to output, with negative sign if integer
!   is < 0.
! INC    Increment between successive elements in a vector or in the
!    column of a matrix.
! INCM   Array giving amount of space used by the options.
! INERR  0 if not processing an error message, 1 if printing an error
!   message, -1 if in an error message that is not being printed, and >1
!   if printing an error message that stops.  Set to -2 when the error
!   message is supposed to stop.
! IOUT   Integer to be output.
! IRC    = 1 for rows, = 2 for columns when determining labels for
!   matrix output.
! IROW   Row index for matrix output.  Also used in table output to
!    count lines for printing line index on read from scratch unit.
! IROW1  Starting row index for matrix output.
! ITEXT  Index of the element of TEXT use for the next text output.
! ITXTSV Saved value of NTEXT when doing matrix output.
! IVAR   Integer array, that is equivalenced to a number of integer
!   variables that can be set by the user.
! IWF    Width to be used in a floating pt. (or other) format.
! IWG    Value like IWF for user set format.
! J      Used as a temporary index.
! JJ      Used as a temporary index.
! K      Used as a temporary index.
! K      Used as a temporary index.
! K1     Used as a temporary index.
! K2     Used as a temporary index.
! KDF    Current number of digits to print for floating point.  See
!   description of MACT(*) = METDIG above.
! KDFDEF Current default for KDF, see description of MACT(*) = MEDDIG.
! KDI    Number of digits used to print last integer.
! KDIAG  Not directly referenced.  Holds place in IVAR for reference
!   from outside.  See comments above.
! KDILAB Length for printing index in vector output.
! KDJ    As for KDI, except for format set by user.
! KK     Temporary index.
! KLINE  Count of number of things to print on a line.  (In table print
!   is the number to print for one spec field.)
! KNT    In vector output gives the current index for output.
! KOLWID Length to use for printing a column heading.  Equivalenced to
!   MAXWID(2).
! KP     Index from error action input for the print action.
! KRES1  Holds place in common block for future use.
! KS     Index from error action input for the stop action.
! KSCRN  Number of lines to "print" before pausing.
! KSHIFT Amount to shift column heading before printing.
! KSPEC  Defines action after looking for character in TEXT.  (Also
!   used as a temporary index.)
!   1.  $B   Break the text here continue on same line.
!   2.  $E   Break text, print what is in BUF.
!   3.  $R   Break text, continue on same line, NTEXT set to repeat the
!            current text.
!   4.  $N   Print BUF, continue with following text.
!   5.  $I   Print IDAT(NIDAT), continue TEXT.
!   6.  $F   Print FDAT(NFDAT), continue TEXT.
!   7.  $M   Print MDAT(NMDAT), continue TEXT.
!   8.  $J   As for $I, but with user format.
!   9.  $G   As for $F, but with user format.
!  10.  $(   Set a user format.
!  11.  $T   Tab.
!  12.       Set when done with an action.
!  13.       Set when done with boiler plate text for an error message.
!   0. Other Ignore the "$", continue with TEXT.
! KT     Used for logic in output of headings.
!        = 1 Output table headings.
!        = 2 Get row/column widths for matrix output.
!        = 3 Output column headings for matrix output.
! LASKNT In vector output value index for last element to print.
! LASTI  Last index for matrix output, or for finding values that
!   determine format.
! LBUF   Position of characters in BUF, usually the last to print.
! LBUF1  Start of text to shift when shifting text in BUF to the left.
! LBUF2  End of text to shift when shifting text in BUF to the left.
! LENBUF Parameter giving the number of character in BUF.
! LENLIN Gives number of character in output lines.
! LENOUT Length of output for table or vector/matrix output.
! LENTXT Length of character array elements in TEXT.
! LENTRY Tells what to do on entry (and sometimes other places.)
!   = 1  The value on first entry.
!   = 2  A previous entry is to be continued.
!   = 3  A non printing error message is to be continued
!   = 4  Just done output from inside a METEXT action.
!   = 5  Got "maximum" value for entries in a vector.
!   = 6  Got "maximum" value for entries in a matrix.
!   = 7  Vector (either print or get format for label indices.)
!   = 8  Matrix (either print or get format for label indices.)
!   = 9  Output of data in a table.
!   =10  Get "maximum" valur for entries in a sparse vector.
!   =11  Output a sparse vector.
! LHEAD  If = 0 no print of DOLS, else DOLS printed in error messages.
! LINERR Gives LENLIN for error messages.
! LINMSG Gives LENLIN for diagnostic messages.
! LINSTR Space at start of line for label in vector and matrix output.
!   Equivalenced to MAXWID(1).
! LNERR  Parameter giving the default value for LENERR, only in XMESS.
! LNMSG  Parameter giving the default value for LINMSG, only in XMESS.
! LOCBEG Index of first item in vector and matrix output.
! LPRINT For error messages with a print level <= LPRINT nothing is
!   printed (unless the message would result in a stop).
! LSTOP  As for LPRINT, except with the stop level, and stopping.
! LSTRT  Starting character position in BUF for storing next characters.
! LTEXT  Length of heading text in TEXT.
! M      Index for the current action.
! MACT   Formal integer array specifying the actions, see above.
! MAXERR Value save in IVAR for user to get indication of last (really
!   the maximum error seen so far = 1000 * (10*stop + print) + index.
! MAXWID Equivalenced to  LINSTR and KOLWID.
! MBNDHI Upper bounds for inputs to IVAR.
! MBNDLO Lower bounds for inputs to IVAR.
! MDAT   Array where user can store integers in IVAR for later output.
!   Also used to store column indices when tables are saved on scratch
!   unit.
!
! The following parameter names starting with ME define actions
!   which the user can request.  They have been documented in the
!   description above, except for the ones defined just below.
! MEGBAS is 1 less than the smallest action that involves something
!   other than just storing or retrieving a value.
! MEMAXI is the largest action which a user can request.
! MEVBAS is the smallest action index, used to set the starting index in
!   IVAR.
! MEVLAS is the largest index for a variable in IVAR.
! MECONT,  MEDDI,  MEELI,  MEEME, MEEUNI, MEFDAT, MEFMAT, MEFSPV,
! MEFVCI, MEFVEC, MEGBAS, MEGMAT, MEGVCI, MEGVEC, MEHEAD, MEIDAT,
! MEIMAT, MEIVCI, MEIVEC, MEJMAT, MEJVCI, MEJVEC, MEMAXE, MEMAXI,
! MEMDA1, MEMDA2, MEMDA3, MEMDA4, MEMDA5, MEMDAT, MEMLIN, MEMUNI,
! MENTXT, MEPRNT, MESCRN, MERES1, MERES2, MERES3,  MERET, MESTOP,
! MESUNI,  METAB, METDIG, METEXT
! MPT  Current pointer to data for matrix or vector output.
! MTEXT  Equivalenced to MTEXTR and MTEXTC.
! MTEXTC TEXT(MTEXTC) starts text for printing column labels.
! MTEXTR TEXT(MTEXTR) starts text for printing row labels.
! MUNIT  Output unit used for messages that aren't in an error message.
! NCOL   Number of columns for matrix output, 0 for vector output,
!   count of column left for table output.
! NDIM   Distance between columns for matrix output.
! NFDAT  Index of next item in FDAT to print.
! NIDAT  Index of next item in IDAT to print.
! NLINE  Maximum number of data items to print on a line for matrix and
!   vector output.  If the scratch file is used for table output,
!   NLINE gives the original end of the buffer.
! NMDAT  Pointer to next thing to print from MDAT.
! NROCO  Equivalenced to (NROW, NCOL).  Used in matrix output.
! NROW   Number of rows for matrix output.    When printing tables,
!   MDAT(NROW) gives place where line was split.  (=0 if not split)
! NSKIP  The amount to skip ahead on the next floating or integer
!   output.
! NTEXT  Index inside an element of TEXT for the next text output.
! NTEXTR Value of NTEXT to use if get a $R.
! NTXTSV Saved value of NTEXT when doing matrix output.
! OUNIT  Index of the current output unit.
! SC     Parameter for special character used to introduce actions.
!   Default value is '$'.  If this is changed the "$"'s in comments
!   should be changed to the new value of the character.  (Note that
!   SC = '\' is not portable.)
! SCRNAM Name of file constructed for error output or message output.
! SUNIT  Index for the scratch unit, -1 if not yet assigned.
! TEXT   Formal argument giving the character string from which all text
!   is taken.
! UMESS  Name of subroutine called that does nothing, but which may be
!   modified by the user to cause different actions to be taken.
!   The usual version of MESS has the call to UMESS commented out.
! XARG   If .true., output data is not integer, and a return is made to
!   print data from FDAT.
! XARGOK Set .true. if call is from program that will print data from
!   FDAT.
!
!++ CODE for .C. is inactive
!      integer  kciwid, kccwid, kcrwid, lbeg, lend, lfprec, lgprec
!      common /MESSCC/ kciwid,kccwid,kcrwid,lbeg,lend,lfprec,lgprec
!%%    long int kc;
!++ END
!
!### History
!  * 2010-02-22 MESS  Krogh  Moved NSKIP=0 to start of code.
!  * 2009-10-30 MESS  Krogh  Defined DSCRN.
!  * 2009-02-28 MESS  Krogh  Added FMTT = ' ' for NAG compiler.
!  * 2009-02-28 MESS  Krogh  Fixed "f" format for C code.
!  * 2007-09-08 MESS  Krogh  Fixed definitions of MEVLAS.
!  * 2006-07-27 MESS  Krogh  Fixed boundary case in printing long text.
!  * 2006-03-20 MESS  Krogh  Added code for output of sparse vector.
!  * 2005-04-07 MESS  Krogh  Declared LFLGDB integer in MESSMH.
!  * 2004-12-15 MESS  Krogh  Added " - 1" at end of line on label 410.
!  * 2002-05-17 MESS  Krogh  Added way for user to get error count.
!  * 2001-12-28 MESS  Krogh  Added NSKIP for more flexible output values.
!  * 2000-12-30 MESS  Krogh  Fixed some types/casts in C code.
!  * 1997-12-12 MESS  Krogh  Prefixed 0P edit descriptor to F format.
!  * 1996-07-11 MESS  Krogh  Transpose matrix output for C.
!  * 1996-06-27 MESS  Krogh  fprintf(stdout, => printf( & memset now used
!  * 1996-06-18 MESS  Krogh  "Saved" NTEXTR.
!  * 1996-05-15 MESS  Krogh  Changes to use .C. and C%%.
!  * 1996-03-30 MESS  Krogh  Added external statement.
!  * 1996-01-24 MESS  Krogh  Fixed minor bug introduced with "$ " stuff.
!  * 1996-01-23 MESS  Krogh  Minor changes for C conversion.
!  * 1995-11-10 MESS  Krogh  Add code to change "$ " to " " in headings.
!  * 1995-08-11 MESS  Krogh  Made code default to not using UMESS.
!  * 1995-01-20 MESS  Krogh  Fixed unusual case in matrix output.
!  * 1994-12-15 MESS  Krogh  Removed block data for Cray T3D.
!  * 1994-11-11 MESS  Krogh  Declared all vars.
!  * 1994-09-14 MESS  Krogh  Fixed to get 1 more "$" in C output.
!  * 1994-09-08 MESS  Krogh  Added new matrix/vector capabilities.
!  * 1994-08-22 MESS  Krogh  Fix for conversion to C for new converter.
!  * 1994-07-05 MESS  Krogh  Fixed bug, KDI and FMTI could be inconsist.
!  * 1994-05-20 MESS  Krogh  Changes to MESSFT so line 1 can go to file.
!  * 1994-05-20 MESS  Krogh  Changes to setting output unit.
!  * 1994-05-09 MESS  Krogh  Integer vectors had overflow & space probs.
!  * 1993-05-19 MESS  Krogh  Changed TEXT to array of character strings.
!  * 1993-04-14 MESS  Krogh  Fixes for conversion to C. (C%% comments.)
!  * 1993-03-10 MESS  Krogh  Broke into smaller pieces.
!  * 1992-12-02 MESS  Krogh  Added save statement to block data subpr.
!  * 1992-07-13 MESS  Krogh  Add checks in heading set up.
!  * 1992-07-12 MESS  Krogh  Fixed so $$ prints a single $ in TEXT.
!  * 1992-07-12 MESS  Krogh  Set out of bound inputs to limit values.
!  * 1992-07-12 MESS  Krogh  Fixed so output works to alternate files.
!  * 1992-07-12 MESS  Krogh  Added integer declarations for parameters.
!  * 1992-06-24 MESS  Krogh  More blanks allowed on break of long lines.
!  * 1992-06-10 MESS  Krogh  Minor fix to vector output.
!  * 1992-05-27 MESS  Krogh  Fixed bug on line width setting.
!  * 1992-05-14 MESS  Krogh  Put common blocks in save statement.
!  * 1992-05-11 MESS  Krogh  Added label to assigned go to & a comment.
!  * 1992-04-08 MESS  Krogh  Unused labels 60, 220 and 320 removed.
!  * 1992-03-20 MESS  Krogh  Changed status on open to SCRATCH.
!  * 1992-03-10 MESS  Krogh  1 line below label 690 changed max to min.
!  * 1992-02-05 MESS  Krogh  Fixed bugs in printing matrix labels.
!  * 1992-01-29 MESS  Krogh  Added UMESS and multiple print option.
!  * 1991-12-09 MESS  Krogh  Fine tuning of vector output.
!  * 1991-10-10 MESS  Krogh  Insure no stop if stop level = 9.
!  * 1991-06-26 MESS  Krogh  Initial Code.

    subroutine MESS(MACT, TEXT, IDAT)

      integer LNMSG, LNERR
      parameter (LNMSG=128)
      parameter (LNERR=79)
!
! ************** Parameters Defining Actions (See Above) ***************
!
      integer   MESUNI, MEHEAD, MEDDIG, MEMLIN, MEELIN, MEMUNI, MEEUNI, &
     &  MESCRN, MEDIAG, MEMAXE, MESTOP, MEPRNT, METDIG, MENTXT, MEIDAT, &
     &  MEFDAT, MEMDAT, MEMDA1, MEMDA2, MEMDA3, MEMDA4, MEMDA5, METABS, &
     &  MEERRS, MECONT, MERET , MEEMES, METEXT, METABL, MERES3, MEIVCI, &
     &  MEIVEC, MEIMAT, MEJVCI, MEJVEC, MEJMAT, MEFVCI, MEFVEC, MEFMAT, &
     &  MEGVCI, MEGVEC, MEGMAT, MEMAXI, MEGBAS, MEVBAS, MEVLAS, MEFSPV
! Parameters for changing the environment.
      parameter (MESUNI=10,MEHEAD=11,MEDDIG=12,MEMLIN=13,MEELIN=14,     &
     & MEMUNI=15,MEEUNI=16,MESCRN=17,MEDIAG=18,MEMAXE=19,MESTOP=20,     &
     & MEPRNT=21,METDIG=22,MENTXT=23,MEIDAT=24,MEFDAT=25,MEMDAT=26,     &
     & MEMDA1=27,MEMDA2=28,MEMDA3=29,MEMDA4=30,MEMDA5=31,METABS=32,     &
     & MEERRS=33)
! Parameters for actions.
      parameter (MECONT=50, MERET=51,MEEMES=52,METEXT=53,MEFSPV=54,     &
     & METABL=55,MERES3=56,MEIVEC=57,MEIMAT=58,MEJVEC=59,MEJMAT=60,     &
     & MEFVEC=61,MEFMAT=62,MEGVEC=63,MEGMAT=64,MEIVCI=65,MEJVCI=66,     &
     & MEFVCI=67,MEGVCI=68)
! Parameter derived from those above.
      parameter (MEMAXI=68,MEGBAS=49,MEVBAS=10,MEVLAS=33)
!
! ************************** Variable Declarations *********************
!
      !external MESSGS
      integer    MACT(*), IDAT(*)
      character  TEXT(*)*(*)
!
      integer    I, ICOL, INCM(MECONT:MEIMAT), INERR, IOUT, IROW, IROW1,&
     &    ITEXTR, ITXTSV, J, JJ, K, K1, K2, KDILAB, KK, KNT, KOLWID, KP,&
     &    KS, LASKNT, LBUF1, LBUF2, LENBUF, LINSTR, M,                  &
     &    MBNDHI(MEVBAS:MEVLAS), MBNDLO(MEVBAS:MEVLAS), MTEXT(2),       &
     &    MTEXTC, MTEXTR, NLINE, NROCO(2), NSKIP, NTEXTR, NTXTSV
      !integer MESSGS
      logical   GETW, FIRST
      character ERMSG*63, ERMSG1*27
      character SC, C
      parameter (SC='$')
      save  FIRST, I, ICOL, INERR, IROW, IROW1, ITXTSV, KDILAB, KNT,    &
     &   LASKNT, M, MTEXT, NLINE, NSKIP, NTEXTR, NTXTSV
      save /CMESSI/, /CMESSC/
      equivalence (MTEXT(1), MTEXTR), (MTEXT(2), MTEXTC)
!
! ************************** Data from common block ********************
!
      parameter (LENBUF=250)
      logical          XARG, GOTFMT, XARGOK
      integer          ERRCNT, EUNIT, ICHAR0, IRC, IVAR(MEVBAS:MEVLAS), &
     &   IMAG, INC, ITEXT, IWF, IWG, KDF, KDFDEF, KDI, KDIAG, KDJ,      &
     &   KLINE, KSCRN, KSHIFT, KSPEC, KT, MAXERR, LASTI, LBUF, LENLIN,  &
     &   LENOUT, LENTRY, LENTXT, LHEAD, LINERR, LINMSG, LOCBEG, LPRINT, &
     &   LSTOP, LSTRT, LTEXT, MAXWID(2), MDAT(5), MPT, MUNIT, NCOL,     &
     &   NDIM, NFDAT, NIDAT, NMDAT, NROW, NTEXT, OUNIT, SUNIT, TABSPA
!
      character BUF*(LENBUF), DOLS*72, FMTC*7, FMTF*20, FMTG*15,        &
     &  FMTI*7, FMTIM(2)*7, FMTJ*7, FMTR*7, FMTT*15
      common /CMESSI/ SUNIT, LHEAD, KDFDEF, LINMSG, LINERR, MUNIT,      &
     &   EUNIT, KSCRN, KDIAG, MAXERR, LSTOP, LPRINT, KDF, NTEXT, NIDAT, &
     &   NFDAT, NMDAT, MDAT, TABSPA, ERRCNT, ICHAR0, IMAG, INC, IRC,    &
     &   ITEXT, IWF, IWG, KDI, KDJ, KLINE, KSHIFT, KSPEC, KT, LASTI,    &
     &   LBUF, LENLIN, LENOUT, LENTRY, LENTXT, LOCBEG, LSTRT, LTEXT,    &
     &   MAXWID, MPT, NROW, NCOL, NDIM, OUNIT, GOTFMT, XARG, XARGOK
      common /CMESSC / BUF, DOLS, FMTF, FMTG, FMTI, FMTJ, FMTT, FMTIM
      equivalence (IVAR(MEVBAS), SUNIT)
      equivalence (FMTIM(1), FMTR), (FMTIM(2), FMTC)
!
      equivalence (NROCO, NROW)
      equivalence (MAXWID(1), LINSTR), (MAXWID(2), KOLWID)
! ************************** End of stuff from common block ************
!
      data INERR, FIRST / 0, .true. /
      data ERMSG /                                                      &
     &' reports error: Stop level = x, Print level = y, Error index = '/
      data ERMSG1 / ': Print level = y, Index = ' /
!                 50  51, 52  53  54 55 56 57 58
      data INCM /  1,  1,  4,  1,  2, 0, 0, 2, 6 /
      data MBNDLO /  0, 0, -50,  39,  39, -99, -99,         0,          &
     &       0,          0, 0, 0, -50,        1,          1,            &
     &       1, 1, -1000000000, -1000000000, -1000000000, -1000000000,  &
     &       -1000000000,   1,  0 /
      data MBNDHI / 99, 1,  50, 500, 500,  99,  99, 100000000,          &
     &  1000000000, 1000000000, 8, 8,  50, 10000000, 1000000000,        &
     &  1000000000, 5, 1000000000, 1000000000, 1000000000, 1000000000,  &
     &  1000000000, 100, 1000000000 /
!
! ************************* Start of Executable Code *******************
!
!
      NSKIP = 0
      if (FIRST) then
         FIRST = .false.
! Initialize common block
         SUNIT = -1
         LHEAD = 1
         LINMSG = LNMSG
         LINERR = LNERR
         MUNIT = 0
         EUNIT = 0
         KSCRN = 0
         MAXERR = 0
         TABSPA = 6
         LSTOP = 3
         LPRINT = 3
         ERRCNT = 0
         ICHAR0 = ICHAR('0')
         KDI = 1
         KDJ = 6
         LENLIN = LNMSG
         LENTRY = 1
         OUNIT = 0
!++ CODE for ~.C. is active
         DOLS(1:40) = '$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$'
         DOLS(41:72) ='$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$'
         FMTI = '(99I01)'
         FMTJ = '(99I06)'
         FMTG = '(1P,99Exx.xx)  '
!++ CODE for .C. is inactive
!%%    memset(cmessc.dols,'$',72);
!      FMTI = '%*d'
!      FMTJ = '%*d\0'
!      FMTG = '%*.*E\0'
!++ END
      else
!               1  2  3   4    5    6    7    8   9   10   11
         go to (5,10,20,850,1160,1620,1130,1530,960,1210,1220), LENTRY
      end if
!                             First entry for a message
    5 LBUF = 0
!                             Usual continuation entry
   10 I = 1
      NTEXT = 1
      ITEXT = 1
      LENTXT = len(TEXT(1))
      NIDAT = 1
      NFDAT = 1
      NMDAT = 1
      go to 120
!                     Continuation entry when have a non printing error
! Skip all actions -- Inside non-printing error message.
   20 I = 1
   30 K = MACT(I)
      if (K <= MERET) then
         if (K == MERET) go to 120
         if (K == MECONT) return
         if (K <= -MEGBAS) go to 180
         I = I + 2
      else
         if (K > MEIMAT) then
            if (K > MEMAXI) go to 180
            K = MEIVEC + mod(K - MEIVEC, 2)
         end if
         I = I + INCM(K)
      end if
      go to 30
!
! Print BUF
   40 call MESSPR
!                             Usual place to end an action request.
  100 I = I + INCM(M)
!                             Pick up the next action request
  120 M = MACT(I)
      if (M > MEGBAS) go to 140
      I = I + 2
      if (abs(M) > MEVLAS) go to 180
      if (M > 0) then
         IVAR(M) = MACT(I-1)
         if (IVAR(M) < MBNDLO(M)) then
            IVAR(M) = MBNDLO(M)
         else if (IVAR(M) > MBNDHI(M)) then
            IVAR(M) = MBNDHI(M)
         end if
!            MEHEAD, MEDDIG, MEMLIN, MEELIN, MEMUNI, MEEUNI
         go to (122,    124,    126,    126,    128,    128), M - MESUNI
         if (M /= MENTXT) go to 120
         ITEXT = (NTEXT-1) / LENTXT
         NTEXT = NTEXT - LENTXT*ITEXT
         ITEXT = ITEXT + 1
         go to 120
  122    if (LHEAD /= 0) then
         end if
         go to 120
  124    KDF = KDFDEF
         go to 120
  126    LENLIN = LINMSG
         go to 120
  128    if (IVAR(M) /= 0) then
!%%          k = labs(cmessi.ounit);
!%%          c_fname[m-15][6] = k / 10 + '0';
!%%          c_fname[m-15][7] = k % 10 + '0';
!%%          if (strcmp(&c_fname[16-m][6], &c_fname[m-15][6]))
!%%             c_handle[m-15] = fopen(c_fname[m-15],"w");
!%%          else
!%%             c_handle[m-15] = c_handle[16-m];
            K = abs(IVAR(M))
         end if
         OUNIT = MUNIT
         go to 120
      end if
      if (M == -MESUNI) then
!%%      if (cmessi.sunit == -1L) {
!%%          scratch_file = tmpfile();
!%%          cmessi.sunit = 1L;}
         if (SUNIT <= 0) SUNIT = MESSGS()
      end if
!
      MACT(I-1) = IVAR(-M)
      go to 120
!  ME ..    CONT  RET EMES ETXT  FSPV TABL
  140 go to (170, 200, 310, 400, 1200, 910, 180), M-MEGBAS
      if (M <= MEGVCI) go to 1000
      go to 180
!
! Action MECONT -- Continue message on next entry
  170 LENTRY = 2
      return
!
! Some kind of error in message specification.
  180 continue
!++ CODE for ~.C. is active
      BUF(1:57) =                                                       &
     &   'Actions in MESS terminated due to error in usage of MESS.'
!++ CODE for .C. is inactive
!%%   memcpy(cmessc.buf,
!%%   "Actions in MESS terminated due to error in usage of MESS.",57);
!++ END
      LBUF = 57
!
! Action MERET -- Finish a message.
  200 LENTRY = 1
      J = INERR
      INERR = 0
      if (J >= 2) INERR = -2
      if (J > 0) go to 330
!                       Finish print before exit.
      call MESSPR
      return
!
! Action MEEMES -- Start an error message
  310 LENTRY = 3
      ERRCNT = ERRCNT + 1
!++  Code for UMESS is inactive
!      call UMESS(TEXT, MACT(I+1), IVAR)
!++  End
      IMAG = max( 0, min(999, MACT(I+2)))
      K = MACT(I+1)
      MAXERR = max(MAXERR, 1000*K + IMAG)
      KS = K / 10
      KP = K - 10 * KS
      if (KS <= min(LSTOP, 8)) then
         if (KP <= LPRINT) then
            INERR = -1
            go to 20
         end if
         INERR = 1
      else
         INERR = 2
      end if
      OUNIT = EUNIT
      LENLIN = LINERR
!                        Output a blank line.
      BUF(1:1) = ' '
      LBUF = 1
  330 call MESSPR
!                        Put out line of $'s
      if (LHEAD /= 0) then
         LBUF = min(len(DOLS), LENLIN)
!++ CODE for ~.C. is active
         BUF(1:LBUF) = DOLS(1:LBUF)
         if (INERR<0) BUF(5:37)=' Fatal error -- Program stopped. '
!++ CODE for .C. is inactive
!%%      memcpy(cmessc.buf, cmessc.dols, cmessi.lbuf);
!%%      if (inerr < 0L)
!%%      memcpy(&cmessc.buf[4]," Fatal error -- Program stopped. ",34);
!++ END
         call MESSPR
      end if
      if (INERR <= 0) then
!                                 Just finished an error message
         if (INERR /= 0) stop
         OUNIT = MUNIT
         LENLIN = LINMSG
         return
      end if
!                     Just starting an error message get program name
      NTEXTR = 0
      go to 410
!                     Got the program name in BUF.
  370 LBUF = min(LBUF, 40)
      if (KS == 0) then
         ERMSG1(17:17) = char(KP + ICHAR0)
!%%       memcpy(&cmessc.buf[cmessi.lbuf], ermsg1, strlen(ermsg1));
         BUF(LBUF+1:LBUF+len(ERMSG1)) = ERMSG1
         LBUF = LBUF + len(ERMSG1)
      else
         ERMSG(30:30) = char(KS + ICHAR0)
         ERMSG(47:47) = char(KP + ICHAR0)
!%%       memcpy(&cmessc.buf[cmessi.lbuf], ermsg, strlen(ermsg));
         BUF(LBUF+1:LBUF+len(ERMSG)) = ERMSG
         LBUF = LBUF + len(ERMSG)
      end if
      LSTRT = LBUF + 1
      call MESSFI
      LBUF = LBUF + KDI
!%%   sprintf(&cmessc.buf[cmessi.lstrt-1L], "%*ld",
!%%           (int)messcc.kciwid, cmessi.imag);
      write (BUF(LSTRT:LBUF), FMTI) IMAG
!          Finish up the start error message action.
      if (MACT(I+3) < 0) go to 40
      if (MACT(I+3) /= 0) then
         ITEXT = (MACT(I+3)-1) / LENTXT
         NTEXT = MACT(I+3) - LENTXT*ITEXT
         ITEXT = ITEXT + 1
      end if
      KSPEC = 13
      go to 480
!                  Take care of any left over print from error header
  390 if (LBUF /= 0) call MESSPR
!
! Action METEXT -- Print string from TEXT
  400 LENTRY = 4
      NTEXTR = NTEXT
      ITEXTR = ITEXT
!                  Continue with print from TEXT
! K     take at most K-1 chars., but if 0 take max number
! K1    is last loc. used from TEXT if LENTXT is BIG.
! NEXT  is first character location in TEXT(ITEXT)
! K2    is last character location in TEXT(ITEXT)
! LSTRT is first character position in BUF
! LBUF  is last used character position in BUF

  410 LSTRT = LBUF + 1
      K2 = min(LENTXT, NTEXT + (LENBUF - LSTRT))
!%%       if ((ctmp=memchr(TEXT(cmessi.itext-1L,cmessi.ntext-1), SC,
!%%          k2 - cmessi.ntext + 1)) == NULL)
!%%             k = 0;
!%%       else
!%%             k = ctmp - TEXT(cmessi.itext-1L,cmessi.ntext-1) + 1;
      K = index(TEXT(ITEXT)(NTEXT:K2), SC)
      if (K == 0) then
! Want to take all that we can.
         LBUF = LSTRT + K2 - NTEXT
!%%       memcpy(&cmessc.buf[cmessi.lstrt-1L], TEXT(cmessi.itext-1L,
!%%         cmessi.ntext-1), k2 - cmessi.ntext + 1L);
         BUF(LSTRT:LBUF) = TEXT(ITEXT)(NTEXT:K2)
         if (K2 == LENTXT) then
           ITEXT = ITEXT + 1
           NTEXT = 1
           if (LBUF <= LENLIN) go to 410
         else
           NTEXT = K2 + 1
         end if
         KSPEC = 12
         if (ITEXT - ITEXTR < 4000) go to 480
         KSPEC = 2
         go to 430
      end if
      LBUF = LBUF + K - 1
!%%   if (k >= 2) memcpy(&cmessc.buf[cmessi.lstrt-1],
!%%     TEXT(cmessi.itext-1L, cmessi.ntext-1), k - 1L);
      if (K >= 2) BUF(LSTRT:LBUF) = TEXT(ITEXT)(NTEXT:NTEXT+K-2)
!        Jump to location below if get $ after computing an NSKIP.
  415 continue
      NTEXT = NTEXT + K + 1
      if (NTEXT > LENTXT) then
         ITEXT = ITEXT + 1
         if (NTEXT == LENTXT + 1) then
            C = TEXT(ITEXT-1)(LENTXT:LENTXT)
            NTEXT = 1
         else
            C = TEXT(ITEXT)(1:1)
            NTEXT = 2
         end if
      else
         C = TEXT(ITEXT)(NTEXT-1:NTEXT-1)
      end if
      if (C == ' ') then
!                Special code to take care of " " following "$".
         NTEXT = NTEXT - 1
         if (NTEXT == 0) then
            ITEXT = ITEXT - 1
            NTEXT = LENTXT
         end if
         go to 410
      end if
      if (NTEXTR == 0) then
         if (LENTRY == 3) go to 370
         go to 1510
      end if
      KSPEC = index('BERNIMFJG(T', C)
  430 if (LBUF > LENLIN) go to 480
!              1   2   3   4   5   6   7   8   9  10  11  12, 13
!              B   E   R   N   I   M   F   J   G   (   T done end err
      go to (455,480,450,460,700,680,900,700,900,600,690,410,390), KSPEC
!               No match  -- Check for setting NSKIP
      if (((C >= '0') .and. (C <= '9')) .or. (C == '-')) then
         NSKIP = 0
         K1 = 1
         if (C /= '-') go to 436
         K1 = -1
  433    C = TEXT(ITEXT)(NTEXT:NTEXT)
         NTEXT = NTEXT + 1
         if (NTEXT >= LENTXT) then
            ITEXT = ITEXT + 1
            NTEXT = 1
         end if
  436    if ((C >= '0') .and. (C <= '9')) then
            NSKIP = 10 * NSKIP + K1 * (ICHAR(C) - ICHAR0)
            go to 433
         end if
         if (C == '$') then
            K = 0
            go to 415
         end if
      end if
!
! Continue with the text.
  440 LBUF = LBUF + 1
      BUF(LBUF:LBUF) = C
      go to 410
!                        Reset NTEXT for $R
  450 NTEXT = NTEXTR
      ITEXT = ITEXTR
!                             Done with METEXT action.
  455 NMDAT = 1
      go to 100
!           At this point want to output all in BUF
  460 do 470 LBUF = LBUF, 1, -1
         if (BUF(LBUF:LBUF) /= ' ') go to 480
  470 continue
  480 LBUF2 = LBUF
      if (LBUF2 == 0) then
         LBUF = 1
         BUF(1:1) = ' '
      else if (LBUF > LENLIN) then
         do 485 K = LENLIN+1, LENLIN/3, -1
            if (BUF(K:K) == ' ') then
               LBUF = K - 1
               go to 490
            end if
  485    continue
         LBUF = LENLIN
      end if
  490 LBUF1 = LBUF
      call MESSPR
      if (LBUF1 >= LBUF2) then
!                       The entire buffer has been printed.
         if (KSPEC <= 2) go to 455
         if (KSPEC /= 4) go to 430
         go to 410
      end if
!                       Remove trailing blanks
      do 510 LBUF1 = LBUF1+1, LBUF2
         if (BUF(LBUF1:LBUF1) /= ' ') go to 520
  510 continue
!                       Shift the contents of the buffer.
  520 LBUF = LBUF2-LBUF1+1
      LSTRT = 1
  530 if (LBUF >= LBUF1) then
!                              Take care of overlap.
         K = 2*LBUF1 - LSTRT
!%%memcpy(&cmessc.buf[cmessi.lstrt-1],&cmessc.buf[lbuf1-1],k-lbuf1);
         BUF(LSTRT:LBUF1-1) = BUF(LBUF1:K-1)
         LSTRT = LBUF1
         LBUF1 = K
         go to 530
      end if
!%% if (cmessi.lbuf>=cmessi.lstrt) memcpy(&cmessc.buf[cmessi.lstrt-1],
!%%       &cmessc.buf[lbuf1-1L], lbuf2-lbuf1+1);
      if (LBUF >= LSTRT) BUF(LSTRT:LBUF) = BUF(LBUF1:LBUF2)
      go to 430
!
! Get information on user format
  600 KSPEC = 8
!              I,   i,   F,   f,   E,   e,   G,   g
      go to (604, 604, 601, 601, 602, 602, 602, 602),                   &
     &   index('IiFfEeGg',TEXT(ITEXT)(NTEXT:NTEXT))
      go to 180
  601 continue
!++ CODE for ~.C. is active
      FMTG='(0P,99F'
!++ CODE for .C. is inactive
!%%   strcpy(cmessc.fmtg, "%*.*f\0");
!%%   messcc.lgprec = 0;
!++ END
      go to 603
  602 continue
!++ CODE for ~.C. is active
      FMTG='(1P,99'//TEXT(ITEXT)(NTEXT:NTEXT)
!++ CODE for .C. is inactive
!%%   strcpy(cmessc.fmtg, "%*.*E\0");
!      FMTG(5:5) = TEXT(ITEXT)(NTEXT:NTEXT)
!%%   messcc.lgprec = 0;
!++ END
  603 KSPEC = 9
  604 IMAG = 0
      GETW = .true.
      K = NTEXT
      FMTT = ' '
  606 continue
         NTEXT = NTEXT + 1
         if (NTEXT > LENTXT) then
            ITEXT = ITEXT + 1
            NTEXT = 1
         end if
!++ CODE for ~.C. is active
         FMTT(NTEXT-K:NTEXT-K) = TEXT(ITEXT)(NTEXT:NTEXT)
!++ END
         JJ = ichar(TEXT(ITEXT)(NTEXT:NTEXT)) - ICHAR0
         if (GETW) then
            if ((JJ >= 0) .and. (JJ <= 9)) then
               IMAG = 10*IMAG + JJ
            else
               if (TEXT(ITEXT)(NTEXT:NTEXT) == ')')  go to 610
               if (TEXT(ITEXT)(NTEXT:NTEXT) /= '.')  go to 180
               GETW = .false.
            end if
         else
            if (TEXT(ITEXT)(NTEXT:NTEXT) == ')') go to 610
            if ((JJ < 0) .or. (JJ > 9)) go to 180
!++ CODE for .C. is inactive
!%%         messcc.lgprec = 10*messcc.lgprec + jj;
!++ END
         end if
      go to 606
!
  610 NTEXT = NTEXT + 1
      if (NTEXT > LENTXT) then
         ITEXT = ITEXT + 1
         NTEXT = 1
      end if
!++ CODE for ~.C. is active
      if (KSPEC == 8) then
         KDJ = IMAG
         FMTJ(5:7) = FMTT
      else
         IWG = IMAG
         FMTG(8:15) = FMTT
      end if
!++ CODE for .C. is inactive
!%%   if (cmessi.kspec == 8)
!%%       cmessi.kdj = cmessi.imag;
!%%   else
!%%       cmessi.iwg = cmessi.imag;
!++ END
      if (TEXT(ITEXT)(NTEXT:NTEXT) == SC) go to 410
      if (KSPEC == 8) go to 700
      if (XARGOK) return
      go to 440
!
!                         Print from MDAT
  680 IOUT = MDAT(NMDAT)
      if (NMDAT >= 6) then
         MDAT(NMDAT) = MDAT(NMDAT) + 1
      else
         NMDAT = NMDAT + 1
      end if
      go to 720
!
!                         Process a tab
  690 LSTRT = LBUF + 1
      LBUF = min(LBUF + TABSPA - mod(LBUF, TABSPA), LENLIN+1)
!%%  for (kc=cmessi.lstrt-1; kc<cmessi.lbuf; kc++) cmessc.buf[kc]=' ';
      BUF(LSTRT:LBUF) = ' '
      go to 850
!                         Print from IDAT
  700 NIDAT = NIDAT + NSKIP
      NSKIP = 0
      IOUT = IDAT(NIDAT)
      NIDAT = NIDAT + 1
  720 LSTRT = LBUF + 1
      IMAG = IOUT
      if (KSPEC >= 8) then
         LBUF = LBUF + KDJ
!%%   sprintf(&cmessc.buf[cmessi.lstrt-1],"%*ld",(int)cmessi.kdj, iout);
      write (BUF(LSTRT:LBUF), FMTJ) IOUT
         go to 850
      end if
!
!                Get format for integer output.
      call MESSFI
      LBUF = LBUF + KDI
!%% sprintf(&cmessc.buf[cmessi.lstrt-1],"%*ld",(int)messcc.kciwid,iout);
      write (BUF(LSTRT:LBUF), FMTI) IOUT
!                         Entry here to check line after numeric output.
  850 if (LBUF <= LENLIN) go to 410
      KSPEC = 12
      go to 480
!
!                          Take care of output for extra argument.
  900 if (XARGOK) return
      go to 180
!
! Action METABL -- Start a table
  910 GOTFMT = MACT(I+1) /= 1
      if (.not. GOTFMT) then
         IROW = 0
         KOLWID = 0
      end if
      LENTRY = 9
      if (LBUF /= 0) call MESSPR
  920 continue
!%%   memset(cmessc.buf,' ',LENBUF);
      BUF = ' '
      NROW = 1
      NCOL = MACT(I+3)
      ICOL = I + 3
  940 ICOL = ICOL + 1
      JJ = MACT(ICOL)
      KLINE = mod(JJ, 100)
      LENOUT = JJ / 100000
      NCOL = NCOL - max(KLINE, 1)
      if (GOTFMT) then
!                                 Print the data
         LSTRT = LBUF + 1
         LBUF = min(LBUF + KLINE * LENOUT, LENBUF)
         JJ = JJ / 100
         KK = mod(JJ, 10)
!              Text,   I   I',   F    E    G
         go to (948, 941, 941, 943, 945, 944), KK
         go to 180
!                             Integer output
  941    continue
!++ CODE for ~.C. is active
         KDI = LENOUT
         FMTI(5:5) = char(LENOUT / 10 + ichar0)
         FMTI(6:6) = char(mod(LENOUT, 10) + ichar0)
!++ END
         if (KK == 3) then
!%%         sprintf(&cmessc.buf[cmessi.lstrt-1], "%*ld",
!%%            (int)cmessi.lenout, mact[i]);
            write (BUF(LSTRT:LBUF), FMTI) MACT(I+1)
            go to 960
         end if
!                            Regular integer output
         NIDAT = NIDAT + NSKIP
         NSKIP = 0
!++ CODE for ~.C. is active
         write (BUF(LSTRT:LBUF), FMTI) (IDAT(K), K = NIDAT,             &
     &      NIDAT+KLINE-1)
         NIDAT = NIDAT + KLINE
!++ CODE for .C. is inactive
!%%  kk = cmessi.nidat;
!%%  for (cmessi.nidat=kk; cmessi.nidat<kk+cmessi.kline; cmessi.nidat++)
!%%     sprintf(&cmessc.buf[cmessi.lstrt+cmessi.lenout*(cmessi.nidat
!%%       - kk) - 1], "%*ld", (int)cmessi.lenout, idat[cmessi.nidat-1]);
!++ END
         go to 960
!                           Various floating point output
  943    continue
!++ CODE for ~.C. is active
         FMTF = '(0P,99F  .  )'
!++ END
         go to 946
  944    continue
!++ CODE for ~.C. is active
         FMTF = '(1P,99G  .  )'
!++ END
         go to 946
  945    continue
!++ CODE for ~.C. is active
         FMTF = '(1P,99E  .  )'
!++ END
  946    JJ = mod(JJ/10, 100)
!++ CODE for ~.C. is active
         FMTF(8:8) = char(ICHAR0 + LENOUT / 10)
         FMTF(9:9) = char(ICHAR0 + mod(LENOUT, 10))
         FMTF(11:11) = char(ICHAR0 + JJ / 10)
         FMTF(12:12) = char(ICHAR0 + mod(JJ, 10))
!++ CODE for .C. is inactive
!%%      strcpy(cmessc.fmtf, "%*.*E\0");
!        IWF = LENOUT
!        lfprec = JJ
!++ END
         if (.not. XARGOK) go to 180
         MPT = NFDAT
         NFDAT = NFDAT + KLINE
         return
!                           Text output
  948    K1 = NTEXT + LBUF - LSTRT
!%%    memcpy(&cmessc.buf[cmessi.lstrt-1], TEXT(cmessi.itext-1,
!%%       cmessi.ntext -1), k1 - cmessi.ntext);
         BUF(LSTRT:LBUF) = TEXT(ITEXT)(NTEXT:K1-1)
         NTEXT = K1
      else
!                                 Print the heading
         KT = 1
         call MESSMH(TEXT)
         if (KT < 0) go to 180
      end if
  960 if ((LBUF <= MDAT(NROW)) .and. (NCOL > 0)) go to 940
      if (NROW == 1) then
         JJ = LBUF
         LBUF = MDAT(1)
         call MESSPR
         LBUF = JJ
      else
         if (IROW == 0) then
            if (NROW == 2) then
!++ CODE for ~.C. is active
               if (SUNIT <= 0) SUNIT = MESSGS()
               rewind(SUNIT)
!++ CODE for .C. is inactive
!%%        if (cmessi.sunit == -1) {
!%%           scratch_file = tmpfile();
!%%           cmessi.sunit = 1;}
!%%        rewind(scratch_file);
!++ END
            end if
         end if
!%%       fwrite(&cmessc.buf[4], cmessi.mdat[cmessi.nrow-1]-4, 1,
!%%          scratch_file);
         write(SUNIT) BUF(5:MDAT(NROW))
      end if
      if (LBUF > MDAT(NROW)) then
!%%  memcpy(&cmessc.buf[4], &cmessc.buf[cmessi.mdat[cmessi.nrow-1]],
!%%     cmessi.lbuf - cmessi.mdat[cmessi.nrow-1]);
         BUF(5:LBUF - MDAT(NROW) + 4) = BUF(MDAT(NROW)+1:LBUF)
         LBUF = LBUF - MDAT(NROW) + 4
         NROW = NROW + 1
         if (.not. GOTFMT) then
            if (NROW > 5) go to 180
            MDAT(NROW) = LBUF
         end if
         if (NCOL == 0) go to 960
         go to 940
      end if
      LBUF = 0
      if (.not. GOTFMT) then
         GOTFMT = .true.
         IROW = IROW - 1
         go to 920
      end if
      MACT(I+1) = MACT(I+1) + 1
      if (MACT(I+1) <= MACT(I+2)) go to 999
      MACT(I+1) = 1
      if (NROW == 1) go to 999
!%%    fputc(EOF, scratch_file);
      endfile SUNIT
      KK = 1
  994 KK = KK + 1
      if (KK > NROW) go to 999
!%%   rewind(scratch_file);
      rewind(SUNIT)
      IROW = -1
      K = KK
  995 LBUF = 5
      IROW = IROW + 1
      if (IROW /= 0) then
!%%      sprintf(cmessc.buf, "%4ld",  irow%10000);
         write(BUF(1:4), '(I4)') mod(IROW, 10000)
      else
!%%    memset(cmessc.buf,' ',4);
         BUF(1:4) = ' '
      end if
      do 996 J = 2, K
         if (J == K) LBUF = MDAT(KK)
!%%       if (fread(&cmessc.buf[4], cmessi.lbuf-4, 1,
!%%         scratch_file) == 0) goto L_994;
         read(SUNIT, END = 994) BUF(5:LBUF)
  996 continue
      K = NROW
      call MESSPR
      go to 995
  999 LENTRY = 1
      return
!
!                          Get started with vector or matrix output
 1000 INC = 1
      LOCBEG = NIDAT
      if (M > MEGMAT) then
! Have a user set increment between entries of a vector.
        M = MEIVEC + 2 * (M - MEIVCI)
        I = I + 1
        INC = MACT(I)
      end if
      XARG = M > MEJMAT
      if (XARG) then
         M = M - 4
         LOCBEG = NFDAT
         if (.not. XARGOK) go to 40
      end if
      GOTFMT = M > MEIMAT
      if (GOTFMT) M = M - 2
      LOCBEG = LOCBEG + NSKIP
      NSKIP = 0
      MPT = LOCBEG
      if (M == MEIMAT) go to 1300
!                           Take care of setup for vector output
      KNT = 0
      LASKNT = MACT(I+1)
      if (LASKNT <= 0) then
         LASKNT = -LASKNT
         KNT = LOCBEG - 1
         if (LASKNT <= KNT) go to 40
      end if
      IMAG = LASKNT
      LASTI = LOCBEG + INC * (LASKNT - 1 - KNT)
      NCOL = 0
!                          Get format for label output.
      call MESSFI
!%%   messcc.kcrwid = messcc.kciwid;
      FMTR = FMTI
      KDILAB = KDI+1
      LINSTR = 2*KDILAB+2
      if (XARG) then
         if (.not. GOTFMT) go to 1150
         IWF = IWG
         FMTF = FMTG
!++ CODE for .C. is inactive
!%%      cmessi.iwf = cmessi.iwg;
!%%      messcc.lfprec = messcc.lgprec;
!++ END
         go to 1160
      end if
      call MESSFD(IDAT)
!                                          After integer format
      LENOUT = KDI
      NIDAT = LASTI + 1
!                                          Common code continues here
 1080 NLINE = (LENLIN - LINSTR + 1) / LENOUT
      if (LBUF == 0) go to 1090
      K = max(LINSTR, LBUF+1)
      if (((K-LINSTR)/LENOUT + (LENLIN-K+1)/LENOUT) < NLINE) K = K + &
     &   LENOUT - mod(K-LINSTR, LENOUT)
      KLINE = (LENLIN - K + 1) / LENOUT
      if (KLINE < min(LASKNT-KNT, NLINE/2)) go to 1085
      LINSTR = K - LENOUT * ((K - LINSTR) / LENOUT)
      if (KLINE >= LASKNT-KNT)  then
         KLINE = LASKNT - KNT
         K = LBUF + 1
      end if
      KNT = KNT + KLINE
!%%    for (kc=cmessi.lbuf; kc < k; kc++) cmessc.buf[kc] = ' ';
      BUF(LBUF+1:K) = ' '
      LBUF = K
      go to 1110
 1085 call MESSPR
 1090 continue
!++ CODE for ~.C. is active
      BUF = ' '
      write (BUF(1:KDILAB), FMTR) KNT+1
!++ CODE for .C. is inactive
!%%   memset(cmessc.buf,' ',LENBUF);
!%%   sprintf(cmessc.buf, "%*ld", (int)messcc.kcrwid, knt+1);
!++ END
      BUF(KDILAB:KDILAB) = '-'
      KLINE = min(NLINE, LASKNT - KNT)
      KNT = KNT + KLINE
!%%    sprintf(&cmessc.buf[kdilab], "%*ld", (int)messcc.kcrwid, knt);
      write (BUF(KDILAB+1:2*KDILAB), FMTR) KNT
!%%    cmessc.buf[kdilab*2L-1] = ':';
!%%    for (kc=kdilab*2L; kc < *linstr-1; kc++) cmessc.buf[kc] = ' ';
      BUF(2*KDILAB:LINSTR-1) = ':'
      LBUF = LINSTR
 1110 LSTRT = LBUF
      LBUF = LBUF + LENOUT * KLINE - 1
      if (XARG) return
!                                    Integer output
!++ CODE for ~.C. is active
      write (BUF(LSTRT:LBUF), FMTI) (IDAT(K), K = MPT,                  &
     &    MPT+INC*(KLINE-1), INC)
!++ CODE for .C. is inactive
!%%   for (k=cmessi.mpt; k<=cmessi.mpt+cmessi.kline-1; k++)
!%%  sprintf(&cmessc.buf[cmessi.lstrt+messcc.kciwid*(k-cmessi.mpt)-1],
!%%      "%*ld", (int)messcc.kciwid, idat[cmessi.inc*k-1]);
!++ END
      MPT = MPT + KLINE * INC
!
!                                     Entry here after vector output.
 1130 if (MPT <= LASTI) go to 1085
      go to 40
!                                          Get other format
 1150 LENTRY = 5
      return
!                                          After other format
 1160 LENOUT = IWF
      LENTRY = 7
      NFDAT = LASTI + 1
      go to 1080


!                         Sparse vector output.
 1200 XARG = .true.
      if (.not. XARGOK) go to 40
      GOTFMT = .false.
      MPT = 1
      LOCBEG = 1
      INC = 1
      LASKNT = MACT(I+1)
      LASTI = LASKNT
      LENTRY = 10
      return

!                Entry after getting format for sparse data output.
 1210 LENOUT = IWF
      LENTRY = 11
      NLINE = LENLIN / IWF

 1220 call MESSPR
      KLINE = min(LASKNT - MPT + 1, NLINE)
      if (KLINE <= 0) go to 40
      LBUF = LENOUT * KLINE
      return

!
!                           Take care of setup for matrix output
 1300 continue
      NDIM = MACT(I+1)
      if (NDIM <= 0) then
         if (NDIM == 0) go to 40
         INC = -NDIM
         NDIM = 1
      end if
      ICOL = 1
      IROW1 = 1
      NROW = MACT(I+2)
      if (NROW <= 0) then
         if (NROW == 0) go to 40
         IROW1 = -NROW / 100000
         NROW = -NROW - 99999 * IROW1 - 1
      end if
      NCOL = MACT(I+3)
      if (NCOL <= 0) then
         if (NCOL == 0) go to 40
         ICOL = -NCOL / 100000
         NCOL = -NCOL - 99999 * IROW1 - 1
      end if
      NTXTSV = NTEXT
      ITXTSV = ITEXT
      IRC = 1
!                        Compute widths for row and column labels
 1320 MAXWID(IRC) = 0
      MTEXT(IRC) = MACT(I+IRC+3)
      IMAG = NROCO(IRC)
      KLINE = IMAG
 1330 NTEXT = MTEXT(IRC)
      if (NTEXT >= 0) then
         if (NTEXT == 0) then
            LTEXT = 5
         else
!                        Go get row/column widths
            KT = 2
            call MESSMH(TEXT)
            if (KT < 0) then
               MTEXT(IRC) = 0
               go to 1330
            end if
         end if
         call MESSFI
         MAXWID(IRC) = max(MAXWID(IRC), LTEXT + KDI+1)
!%%      if (cmessi.irc == 1)
!%%         messcc.kcrwid = cmessi.kdi;
!%%      else
!%%         messcc.kccwid = cmessi.kdi;
         FMTIM(IRC) = FMTI
      end if
      IRC = IRC + 1
      if (IRC == 2) go to 1320
!                 Widths for Row and column titles have been computed.
      KSHIFT = 1
      LASTI = LOCBEG + INC * (NROW - IROW1)
      if (XARG) then
         if (.not. GOTFMT) go to 1610
!++ CODE for ~.C. is active
         IWF = IWG
         FMTF = FMTG
!++ CODE for .C. is inactive
!%%      cmessi.iwf = cmessi.iwg;
!%%      messcc.lfprec = messcc.lgprec;
!++ END
         go to 1620
      end if
      call MESSFD(IDAT)
!
      If (KDI >= KOLWID) then
         LENOUT = KDI
      else
         KSHIFT = (KOLWID - KDI + 2) /2
         LENOUT = KOLWID
!++ CODE for ~.C. is active
         KDI = KOLWID
         FMTI(5:5) = char(ICHAR0 + KOLWID / 10)
         FMTI(6:6) = char(ICHAR0 + mod(KOLWID, 10))
!++ CODE for .C. is inactive
!%%  messcc.kciwid = *kolwid;
!++ END
      end if
      NIDAT = NIDAT + NDIM*NCOL
!                              Continue with commmon code
 1390 NLINE = (LENLIN - LINSTR) / LENOUT
      if (LBUF <= LINSTR) go to 1420
 1400 call MESSPR
 1420 IROW = IROW1
      KLINE = min(NLINE, NCOL-ICOL+1)
!                       Output column labels (if any)
      if (MTEXTC < 0) go to 1480
      NTEXT = MTEXTC
      IMAG = ICOL
      KT = 3
      call MESSMH(TEXT)
      if (KT < 0) go to 180
!                       Return from output of column labels.
      MTEXTC = NTEXT
 1480 ICOL = ICOL + KLINE
 1490 call MESSPR
!
!                      Output row labels (if any)
      if (MTEXTR < 0) go to 1520
      if (MTEXTR == 0) then
!%%       memcpy(&cmessc.buf[cmessi.lbuf],"Row ", 4);
         BUF(LBUF+1:LBUF+4) = 'Row '
         LBUF = LBUF + 4
         go to 1515
      end if
      NTEXT = MTEXTR
      ITEXT = (NTEXT-1) / LENTXT
      NTEXT = NTEXT - ITEXT * LENTXT
      ITEXT = ITEXT + 1

!                     Go get text for row label
      NTEXTR = 0
      go to 410
!                     Return from getting text for row label
 1510 if (C /= '#') then
         MTEXTR = NTEXT + LENTXT * (ITEXT-1)
!%%    for (kc=cmessi.lbuf; kc < *linstr; kc++) cmessc.buf[kc] = ' ';
         BUF(LBUF+1:LINSTR) = ' '
         go to 1520
      end if
 1515 continue
!%%   sprintf(&cmessc.buf[cmessi.lbuf],"%*ld",(int)messcc.kcrwid,irow);
!%%    for (kc=cmessi.lbuf+messcc.kcrwid;
!%%       kc < *linstr; kc++) cmessc.buf[kc] = ' ';
      write (BUF(LBUF+1:LINSTR), FMTR) IROW
 1520 LSTRT = LINSTR + 1
      LBUF = LINSTR + LENOUT*KLINE
      LASTI = MPT + NDIM * KLINE - 1
      if (XARG) return
!                                    Integer output
!%% for (k=cmessi.mpt; k<=cmessi.lasti; k+=cmessi.ndim)
!%%  sprintf(&cmessc.buf[cmessi.lstrt + messcc.kciwid*(k-cmessi.mpt)/
!%%     cmessi.ndim - 1], "%*ld", (int)messcc.kciwid, idat[k-1]);
         write (BUF(LSTRT:LBUF), FMTI) (IDAT(K), K=MPT,LASTI,NDIM)
!
!                                     Entry here after matrix output.
 1530 MPT = MPT + INC
      IROW = IROW + 1
!
      if (IROW <= NROW) go to 1490
      if (ICOL > NCOL) then
         NTEXT = NTXTSV
         ITEXT = ITXTSV
         go to 40
      end if
      MPT = NDIM*(ICOL-1) + 1
      MTEXTR = MACT(I+4)
      call MESSPR
      LBUF = 1
      BUF(1:1) = ' '
      go to 1400
!                                Need to get format for matrix print.
 1610 LENTRY = 6
      return
!                                Entry after got format for matrix print
 1620 If (IWF >= KOLWID) then
         LENOUT = IWF
      else
         KSHIFT = (KOLWID - IWF + 2) /2
         LENOUT = KOLWID
!%%      cmessi.iwf = *kolwid;
!%%      strcpy(cmessc.fmtf, "%*.*E\0");
         write (FMTF(7:8), '(I2)') KOLWID
      end if
      NFDAT = NFDAT + NDIM*NCOL
      LENTRY = 8
      go to 1390
    end subroutine mess
!*************************************************************************

!*************************************************************************
!>
! Get the format for data to be printed in vectors and arrays.
!
! ************** Variable only used here *******************************
!
! K      Temporary index.
! J      Temporary index.
! IDAT   Input array to MESS
! IMAX   Used when computing largest integer in array.
! IMIN   Used when computing smallest integer in array.

    subroutine MESSFD(IDAT)

      integer J, K, IDAT(*), IMAX, IMIN
!
! For comments on other variables, see the listing for MESS.
      integer   LENBUF, MEVBAS, MEVLAS
      parameter (LENBUF=250)
      parameter (MEVBAS=10)
      parameter (MEVLAS=33)
      logical          XARG, GOTFMT, XARGOK
      integer          ERRCNT, EUNIT, ICHAR0, IRC, IVAR(MEVBAS:MEVLAS), &
     &   IMAG, INC, ITEXT, IWF, IWG, KDF, KDFDEF, KDI, KDIAG, KDJ,      &
     &   KLINE, KSCRN, KSHIFT, KSPEC, KT, MAXERR, LASTI, LBUF, LENLIN,  &
     &   LENOUT, LENTRY, LENTXT, LHEAD, LINERR, LINMSG, LOCBEG, LPRINT, &
     &   LSTOP, LSTRT, LTEXT, MAXWID(2), MDAT(5), MPT, MUNIT, NCOL,     &
     &   NDIM, NFDAT, NIDAT, NMDAT, NROW, NTEXT, OUNIT, SUNIT, TABSPA
!
      character BUF*(LENBUF), DOLS*72, FMTC*7, FMTF*20, FMTG*15,        &
     &  FMTI*7, FMTIM(2)*7, FMTJ*7, FMTR*7, FMTT*15
      common /CMESSI/ SUNIT, LHEAD, KDFDEF, LINMSG, LINERR, MUNIT,      &
     &   EUNIT, KSCRN, KDIAG, MAXERR, LSTOP, LPRINT, KDF, NTEXT, NIDAT, &
     &   NFDAT, NMDAT, MDAT, TABSPA, ERRCNT, ICHAR0, IMAG, INC, IRC,    &
     &   ITEXT, IWF, IWG, KDI, KDJ, KLINE, KSHIFT, KSPEC, KT, LASTI,    &
     &   LBUF, LENLIN, LENOUT, LENTRY, LENTXT, LOCBEG, LSTRT, LTEXT,    &
     &   MAXWID, MPT, NROW, NCOL, NDIM, OUNIT, GOTFMT, XARG, XARGOK
      common /CMESSC / BUF, DOLS, FMTF, FMTG, FMTI, FMTJ, FMTT, FMTIM
      equivalence (IVAR(MEVBAS), SUNIT)
      equivalence (FMTIM(1), FMTR), (FMTIM(2), FMTC)
!
      save /CMESSI/, /CMESSC/
!
      if (GOTFMT) then
         KDI = KDJ
!%%      messcc.kciwid = cmessi.kdj;
         FMTI = FMTJ
         return
      end if
      K = 1
      IMAX = 1
      IMIN = 0
   10 do 20 J = LOCBEG, LASTI, INC
         IMAX = max(IMAX, IDAT(J))
         IMIN = MIN(IMIN, IDAT(J))
   20 continue
      if (NCOL /= 0) then
         K = K + 1
         LOCBEG = LOCBEG + NDIM
         LASTI = LASTI + NDIM
         if (K <= NCOL) go to 10
      end if
      IMAG = IMAX
      if ((IMAG/10) + IMIN < 0) IMAG = IMIN
      KDI = -KDI
      call MESSFI
      return
    end subroutine messfd
!*************************************************************************

!*************************************************************************
!>
! Get the format for the integer IMAG.
!
! ************** Variable only used here *******************************
!
! I, K, KD are used in determining number of characters needed to
!          represent IMAG.

    subroutine MESSFI()

      integer I, K, KD
!
! For comments on other variables, see the listing for MESS.
      integer   LENBUF, MEVBAS, MEVLAS
      parameter (LENBUF=250)
      parameter (MEVBAS=10)
      parameter (MEVLAS=33)
      logical          XARG, GOTFMT, XARGOK
      integer          ERRCNT, EUNIT, ICHAR0, IRC, IVAR(MEVBAS:MEVLAS), &
     &   IMAG, INC, ITEXT, IWF, IWG, KDF, KDFDEF, KDI, KDIAG, KDJ,      &
     &   KLINE, KSCRN, KSHIFT, KSPEC, KT, MAXERR, LASTI, LBUF, LENLIN,  &
     &   LENOUT, LENTRY, LENTXT, LHEAD, LINERR, LINMSG, LOCBEG, LPRINT, &
     &   LSTOP, LSTRT, LTEXT, MAXWID(2), MDAT(5), MPT, MUNIT, NCOL,     &
     &   NDIM, NFDAT, NIDAT, NMDAT, NROW, NTEXT, OUNIT, SUNIT, TABSPA
!
      character BUF*(LENBUF), DOLS*72, FMTC*7, FMTF*20, FMTG*15,        &
     &  FMTI*7, FMTIM(2)*7, FMTJ*7, FMTR*7, FMTT*15
      common /CMESSI/ SUNIT, LHEAD, KDFDEF, LINMSG, LINERR, MUNIT,      &
     &   EUNIT, KSCRN, KDIAG, MAXERR, LSTOP, LPRINT, KDF, NTEXT, NIDAT, &
     &   NFDAT, NMDAT, MDAT, TABSPA, ERRCNT, ICHAR0, IMAG, INC, IRC,    &
     &   ITEXT, IWF, IWG, KDI, KDJ, KLINE, KSHIFT, KSPEC, KT, LASTI,    &
     &   LBUF, LENLIN, LENOUT, LENTRY, LENTXT, LOCBEG, LSTRT, LTEXT,    &
     &   MAXWID, MPT, NROW, NCOL, NDIM, OUNIT, GOTFMT, XARG, XARGOK
      common /CMESSC / BUF, DOLS, FMTF, FMTG, FMTI, FMTJ, FMTT, FMTIM
      equivalence (IVAR(MEVBAS), SUNIT)
      equivalence (FMTIM(1), FMTR), (FMTIM(2), FMTC)
!
      save /CMESSI/, /CMESSC/
!
      KD = 1
      if (KDI < 0) then
!              KDI < 0 to flag need for extra space -- avoids overflows
         KDI = -KDI
         KD = 2
      end if
      K = 1
      if (IMAG < 0) then
         IMAG = -IMAG
         KD = KD + 1
      end if
      I = IMAG / 10
      if (I /= 0) then
   10    K = 10 * K
         KD = KD + 1
         if (I >= K) go to 10
      end if
      if (KD /= KDI) then
         KDI = KD
!++ CODE for ~.C. is active
         FMTI(5:5) = char(ICHAR0 + KDI / 10)
         FMTI(6:6) = char(ICHAR0 + mod(KDI, 10))
!++ CODE for .C. is inactive
!%%      messcc.kciwid = cmessi.kdi;
!++ END
      end if
      return
    end subroutine messfi
!*************************************************************************

!*************************************************************************
!>
! Get a scratch unit assigned.

    integer function MESSGS()
    integer J

      MESSGS = 31
   10 MESSGS = MESSGS - 1
      if (MESSGS == 0) stop 'Could not assign scratch unit in MESS.'
      open (MESSGS, STATUS='SCRATCH', ACCESS='SEQUENTIAL', &
            FORM='UNFORMATTED', IOSTAT=J)
      if (J /= 0) go to 10

    END function MESSGS
!*************************************************************************

!*************************************************************************
!>
! Processing of multiple headings:
!
! J     Used as a temporary index.
! K     Used as a temporary index.
! KB    Number of initial blanks
! KK    Used as a temporary index.
! KT    Used for logic in output of headings.  Set <0 on exit if there
!       is an input error.
!     KT = 1 Output table headings.  (Set to -1 on fatal error.)
!     KT = 2 Get row/column widths for matrix output.  (Set to -2 if
!            error results in no headings.)
!     KT = 3 Output column headings.  (Set to -1 on fatal error.)
! L     Used as a temporary index.
! LFLGDB 2 in usual case, 3 if see a "$ ".
! LSTRDB Value of LSTRT when see a "$ ".
! LTXTDB Value of LTEXT when see a "$ ".
! TEXT  Original input character vector.

    subroutine MESSMH(TEXT)

      integer J, K, KB, KK, L, LFLGDB, LSTRDB, LTXTDB
      character*(*)  TEXT(*)
      character SC, C
      parameter (SC='$')
! For comments on other variables, see the listing for MESS.
      integer   KOLWID, LINSTR, LENBUF, MEVBAS, MEVLAS
      parameter (LENBUF=250)
      parameter (MEVBAS=10)
      parameter (MEVLAS=33)
      logical          XARG, GOTFMT, XARGOK
      integer          ERRCNT, EUNIT, ICHAR0, IRC, IVAR(MEVBAS:MEVLAS), &
     &   IMAG, INC, ITEXT, IWF, IWG, KDF, KDFDEF, KDI, KDIAG, KDJ,      &
     &   KLINE, KSCRN, KSHIFT, KSPEC, KT, MAXERR, LASTI, LBUF, LENLIN,  &
     &   LENOUT, LENTRY, LENTXT, LHEAD, LINERR, LINMSG, LOCBEG, LPRINT, &
     &   LSTOP, LSTRT, LTEXT, MAXWID(2), MDAT(5), MPT, MUNIT, NCOL,     &
     &   NDIM, NFDAT, NIDAT, NMDAT, NROW, NTEXT, OUNIT, SUNIT, TABSPA
!
      character BUF*(LENBUF), DOLS*72, FMTC*7, FMTF*20, FMTG*15,        &
     &  FMTI*7, FMTIM(2)*7, FMTJ*7, FMTR*7, FMTT*15
      common /CMESSI/ SUNIT, LHEAD, KDFDEF, LINMSG, LINERR, MUNIT,      &
     &   EUNIT, KSCRN, KDIAG, MAXERR, LSTOP, LPRINT, KDF, NTEXT, NIDAT, &
     &   NFDAT, NMDAT, MDAT, TABSPA, ERRCNT, ICHAR0, IMAG, INC, IRC,    &
     &   ITEXT, IWF, IWG, KDI, KDJ, KLINE, KSHIFT, KSPEC, KT, LASTI,    &
     &   LBUF, LENLIN, LENOUT, LENTRY, LENTXT, LOCBEG, LSTRT, LTEXT,    &
     &   MAXWID, MPT, NROW, NCOL, NDIM, OUNIT, GOTFMT, XARG, XARGOK
      common /CMESSC / BUF, DOLS, FMTF, FMTG, FMTI, FMTJ, FMTT, FMTIM
      equivalence (IVAR(MEVBAS), SUNIT)
      equivalence (FMTIM(1), FMTR), (FMTIM(2), FMTC)
!
      equivalence (MAXWID(1), LINSTR), (MAXWID(2), KOLWID)

      save /CMESSI/, /CMESSC/
!++ CODE for .C. is inactive
!%%      long int kc;
!++ END
      save LFLGDB
      data LFLGDB / 2 /
!
      if (NTEXT /= 0) then
         ITEXT = (NTEXT-1) / LENTXT
         NTEXT = NTEXT - ITEXT * LENTXT
         ITEXT = ITEXT + 1
      end if
      do 300 J = 1, max(1,KLINE)
         if (NTEXT == 0) then
            K = KOLWID
            go to 210
         end if
         LFLGDB = 2
         LTEXT = 0
  110    continue
!%%       ctmp=memchr(TEXT(cmessi.itext-1L,cmessi.ntext-1), SC,
!%%          cmessi.lentxt - cmessi.ntext + 1);
!%%       if (ctmp == NULL)
!%%             l = 0;
!%%       else
!%%             l = ctmp - TEXT(cmessi.itext-1L,cmessi.ntext-1) + 1;
         L = index(TEXT(ITEXT)(NTEXT:LENTXT), SC)
         if (L == 0) then
            LTEXT = LTEXT + LENTXT - NTEXT + 1
            if (LTEXT < 80) then
               ITEXT = ITEXT + 1
               NTEXT = 1
               go to 110
            end if
            LTEXT = 0
            if (KT == 3) go to 310
            go to 160
         end if
         NTEXT = NTEXT + L + 1
         LTEXT = L + LTEXT - 1
         if (NTEXT > LENTXT) then
            ITEXT = ITEXT + 1
            if (NTEXT == LENTXT + 1) then
               C = TEXT(ITEXT-1)(LENTXT:LENTXT)
               NTEXT = 1
            else
               C = TEXT(ITEXT)(1:1)
               NTEXT = 2
            end if
         else
            C = TEXT(ITEXT)(NTEXT-1:NTEXT-1)
         end if
         if (C == 'H') go to (180, 190, 200), KT
         if (C == 'E') go to (180, 310, 200), KT
         if (C == '#') go to (140, 310, 200), KT
         if (C == ' ') then
!  Special code to set for removing the "$" preceding a blank.
            LSTRDB = LSTRT
            LFLGDB = 3
            LTXTDB = LTEXT
            LTEXT = LTEXT + 1
            go to 110
         end if
         if (KT /= 1) go to 160
  140    LTEXT = LTEXT + 2
         go to 110
  160    KT = -KT
         go to 310
!
  180    KOLWID = KOLWID + LENOUT
         if (LTEXT == 0) go to 300
         KB = KOLWID-LTEXT
         if (KB < 0) stop                                            &
     &   'Stopped in MESS -- Column width too small in a heading.'
         if (XARG)  KB = 1 + KB/2
         LSTRT = LBUF + KB + 1
         LBUF = LBUF + KOLWID
         if (LBUF <= LENLIN) MDAT(NROW) = LBUF
         KOLWID = 0
         go to 220
!
!                                  Set up column widths
  190    MAXWID(IRC) = max(MAXWID(IRC), LTEXT)
         go to 300
!
!                                  Output matrix column
  200    K = KOLWID
         if (C /= '#') K = LTEXT
  210    KB = LENOUT - KOLWID
         if (J == 1) then
!                        Special setup for the first column.
            if (XARG) KB = (KB + 1) / 2
            KB = KB + KSHIFT + LINSTR - LBUF
         end if
         KB = KB + KOLWID - K
         LSTRT = LBUF + KB + 1
         LBUF = LSTRT + K - 1
!                                  Set initial blanks
  220    continue
!%%      if (kb > 0) for (kc=cmessi.lstrt-kb-1; kc<cmessi.lstrt-1; kc++)
!%%         cmessc.buf[kc] = ' ';
         if (KB > 0) BUF(LSTRT-KB:LSTRT-1) = ' '
!                                  Move characters
         if (NTEXT == 0) then
!%%       memcpy(&cmessc.buf[cmessi.lstrt-1],"Col ", 4);
            BUF(LSTRT:LSTRT+3) = 'Col '
            C = '#'
            LSTRT = LSTRT+4
         else
            K = NTEXT - LTEXT - LFLGDB
            if (K <= 0) then
               KK = max(0, 3-NTEXT)
!%%       memcpy(&cmessc.buf[cmessi.lstrt-1L], TEXT(cmessi.itext-2L,
!%%         cmessi.lentxt+k-1), -k-kk+1L);
               BUF(LSTRT:LSTRT-K-KK)=TEXT(ITEXT-1)(LENTXT+K:LENTXT-KK)
               LSTRT = LSTRT-K-KK+1
               K = 1
            end if
            if (NTEXT > 3) then
!%%       memcpy(&cmessc.buf[cmessi.lstrt-1L], TEXT(cmessi.itext-1L,
!%%         k-1), cmessi.ntext-k-2L);
               BUF(LSTRT:LSTRT+NTEXT-K-3) = TEXT(ITEXT)(K:NTEXT-3)
               LSTRT = LSTRT + NTEXT - K - 2
            end if
         end if
         if (LFLGDB == 3) then
!  Special code to remove the "$" preceding a blank.  Only works for 1.
            do 250 L = LSTRDB + LTXTDB + max(0, KB), LSTRT
               BUF(L:L) = BUF(L+1:L+1)
  250       continue
            LFLGDB = 2
            LSTRT = LSTRT - 1
         end if
         if (C == '#') then
!                                  Output column index
!%%         sprintf(&cmessc.buf[cmessi.lstrt-1], "%*ld ",
!%%           (int)(cmessi.lbuf-cmessi.lstrt), cmessi.imag+j-1);
            write (BUF(LSTRT:LBUF), FMTC) IMAG + J - 1
            if (NTEXT /= 0) NTEXT = K
            go to 300
         end if
!                                  Set trailing blanks
!%%      if (cmessi.lstrt <= cmessi.lbuf)
!%%           for (kc=cmessi.lstrt-1; kc < cmessi.lbuf; kc++)
!%%              cmessc.buf[kc] = ' ';
         if (LSTRT <= LBUF) BUF(LSTRT:LBUF) = ' '
  300 continue
  310 return
    end subroutine messmh
!*************************************************************************

!*************************************************************************
!>
! Prints the buffer for MESS
!
! ************** Variable only used here *******************************
!
! NSCRN  Number of lines currently on CRT from messages.

    subroutine MESSPR()

      integer   NSCRN, K
      character SCRNAM*12
      save      NSCRN
!
! For comments on other variables, see the listing for MESS.
      integer   LENBUF, MEVBAS, MEVLAS
      parameter (LENBUF=250)
      parameter (MEVBAS=10)
      parameter (MEVLAS=33)
      logical          XARG, GOTFMT, XARGOK
      integer          ERRCNT, EUNIT, ICHAR0, IRC, IVAR(MEVBAS:MEVLAS), &
     &   IMAG, INC, ITEXT, IWF, IWG, KDF, KDFDEF, KDI, KDIAG, KDJ,      &
     &   KLINE, KSCRN, KSHIFT, KSPEC, KT, MAXERR, LASTI, LBUF, LENLIN,  &
     &   LENOUT, LENTRY, LENTXT, LHEAD, LINERR, LINMSG, LOCBEG, LPRINT, &
     &   LSTOP, LSTRT, LTEXT, MAXWID(2), MDAT(5), MPT, MUNIT, NCOL,     &
     &   NDIM, NFDAT, NIDAT, NMDAT, NROW, NTEXT, OUNIT, SUNIT, TABSPA
!
      character BUF*(LENBUF), DOLS*72, FMTC*7, FMTF*20, FMTG*15,        &
     &  FMTI*7, FMTIM(2)*7, FMTJ*7, FMTR*7, FMTT*15
      common /CMESSI/ SUNIT, LHEAD, KDFDEF, LINMSG, LINERR, MUNIT,      &
     &   EUNIT, KSCRN, KDIAG, MAXERR, LSTOP, LPRINT, KDF, NTEXT, NIDAT, &
     &   NFDAT, NMDAT, MDAT, TABSPA, ERRCNT, ICHAR0, IMAG, INC, IRC,    &
     &   ITEXT, IWF, IWG, KDI, KDJ, KLINE, KSHIFT, KSPEC, KT, LASTI,    &
     &   LBUF, LENLIN, LENOUT, LENTRY, LENTXT, LOCBEG, LSTRT, LTEXT,    &
     &   MAXWID, MPT, NROW, NCOL, NDIM, OUNIT, GOTFMT, XARG, XARGOK
      common /CMESSC / BUF, DOLS, FMTF, FMTG, FMTI, FMTJ, FMTT, FMTIM
      equivalence (IVAR(MEVBAS), SUNIT)
      equivalence (FMTIM(1), FMTR), (FMTIM(2), FMTC)
!
      save /CMESSI/, /CMESSC/
      data NSCRN / 0 /
!
      if (LBUF /= 0) then
   10   if (BUF(LBUF:LBUF) == ' ') then
          if (LBUF > 1) then
            LBUF = LBUF - 1
            go to 10
          end if
        end if
        if (OUNIT <= 0) then
          if (KSCRN > 0) then
            if (NSCRN >= KSCRN) then
!%%               printf( " Type 'Enter' to continue\n" );
              print '('' Type "Enter" to continue'')'
!%%               scanf( "%*[^\n]%*c" );
              read (*, *)
              NSCRN = 0
            end if
            NSCRN = NSCRN + 1
          end if
!%%      printf( "%.*s\n", (int)cmessi.lbuf, cmessc.buf);
          print '(1X, A)', BUF(1:LBUF)
          if (OUNIT == 0) go to 20
        end if
!++ CODE for ~.C. is active
        K = abs(OUNIT)
        write (K, '(A)', ERR=30) BUF(1:LBUF)
!++ CODE for .C. is inactive
!%%      fprintf(c_handle[labs(cmessi.ounit)-1], "%.*s\n",
!%%      (int)cmessi.lbuf, cmessc.buf);
!++ END
   20   LBUF = 0
      end if
      return
!++ CODE for ~.C. is active
!              See if opening fixes the error
   30 write(SCRNAM, '(A, I2.2, A)') 'MESSF_', K, '.tmp'
      open (UNIT=K, STATUS='UNKNOWN', FILE=SCRNAM)
      write (K, '(A)') BUF(1:LBUF)
      return
!++ END
    end subroutine messpr
!*************************************************************************

!*************************************************************************
!>
!  Prints FTEXT, which contains a Fortran character string, and then
!  call MESS to do the actions in MACT.  Actions in MACT can not do
!  anything other than actions that reference MACT.
!  This routine intended for use by library subroutines getting text in
!  the form of a Fortran character string.

    subroutine MESSFT(MACT, FTEXT)

      integer MACT(*)
      character FTEXT*(*)
!
      integer J, K, IDAT(1), MECONT, MEPRNT, MESUNI
!++ CODE for ~.C. is active
      character TEXT(1)*1
!++ CODE for .C. is inactive
!      character TEXT(1)*2
!++ END
      parameter (MESUNI=10, MEPRNT=21, MECONT=50)
!
      INTEGER LENBUF, MEVBAS, MEVLAS
      parameter (LENBUF=250)
      parameter (MEVBAS=10)
      parameter (MEVLAS=33)
      logical          XARG, GOTFMT, XARGOK
      integer          ERRCNT, EUNIT, ICHAR0, IRC, IVAR(MEVBAS:MEVLAS), &
     &   IMAG, INC, ITEXT, IWF, IWG, KDF, KDFDEF, KDI, KDIAG, KDJ,      &
     &   KLINE, KSCRN, KSHIFT, KSPEC, KT, MAXERR, LASTI, LBUF, LENLIN,  &
     &   LENOUT, LENTRY, LENTXT, LHEAD, LINERR, LINMSG, LOCBEG, LPRINT, &
     &   LSTOP, LSTRT, LTEXT, MAXWID(2), MDAT(5), MPT, MUNIT, NCOL,     &
     &   NDIM, NFDAT, NIDAT, NMDAT, NROW, NTEXT, OUNIT, SUNIT, TABSPA
!
      character BUF*(LENBUF), DOLS*72, FMTC*7, FMTF*20, FMTG*15,        &
     &  FMTI*7, FMTIM(2)*7, FMTJ*7, FMTR*7, FMTT*15
      common /CMESSI/ SUNIT, LHEAD, KDFDEF, LINMSG, LINERR, MUNIT,      &
     &   EUNIT, KSCRN, KDIAG, MAXERR, LSTOP, LPRINT, KDF, NTEXT, NIDAT, &
     &   NFDAT, NMDAT, MDAT, TABSPA, ERRCNT, ICHAR0, IMAG, INC, IRC,    &
     &   ITEXT, IWF, IWG, KDI, KDJ, KLINE, KSHIFT, KSPEC, KT, LASTI,    &
     &   LBUF, LENLIN, LENOUT, LENTRY, LENTXT, LOCBEG, LSTRT, LTEXT,    &
     &   MAXWID, MPT, NROW, NCOL, NDIM, OUNIT, GOTFMT, XARG, XARGOK
      common /CMESSC / BUF, DOLS, FMTF, FMTG, FMTI, FMTJ, FMTT, FMTIM
      equivalence (IVAR(MEVBAS), SUNIT)
      equivalence (FMTIM(1), FMTR), (FMTIM(2), FMTC)
!
      do 10 J = 1, 100, 2
         K = abs(MACT(J))
         if ((K > MEPRNT) .or. (K < MESUNI)) go to 20
   10 continue
   20 K = MACT(J)
      MACT(J) = MECONT
      call MESS(MACT, TEXT, IDAT)
      MACT(J) = K
!%%      k = strlen(ftext);
      K = len(FTEXT)
      NTEXT = 1
      if (K /= 0) then
         if (FTEXT(1:1) == '0') then
            NTEXT = 2
            K = K - 1
            if (LBUF == 0) then
               BUF(1:1) = ' '
               LBUF = 1
            end if
         end if
         call MESSPR
         LBUF = K
!%%      memcpy(cmessc.buf, &ftext[cmessi.ntext-1], k);
         BUF(1:K) = FTEXT(NTEXT:NTEXT+K-1)
      end if
      ICHAR0 = ICHAR('0')
      if (MACT(J) /= MECONT) call mess(MACT(J), TEXT, IDAT)
      return
    end subroutine messft
!*************************************************************************

!*************************************************************************
!>
! This subroutine is intended for the use of other library routines.
! It is used to check the storage used by options for some array.
!
! INTCHK is an array that provides information on how storage has been
!   allocated.  Information is provided in groups of three words, after
!   an initial group of 4 that contain:
!    0. INTCHK(0) / 10 is the index to use for error messages, and
!       mod(INTCHK(0), 10) defines actions as follows:
!       Value  Action on Errors  Action if no errors
!         0    Prints & Stops               Returns
!         1    Prints & Returns             Returns
!         2    Prints & Stops      Prints & Returns
!         3    Prints & Returns    Prints & Returns
!        >3  Any error message will be continued, subtract 4 from the
!            value to get one of the actions above.
!    1. Contains LAST, 1 more than the index of the last location in
!       INTCHK used to store the information as described below.
!    2. The amount of storage available in the array being checked.  If
!       this is 0, it is assumed that the user would like to know what
!       size to declare for the array, and if it is -1, it is assumed
!       that a library routine is doing the check without knowing what
!       the declared size is.
!    3. The amount of storage required in this array if no options were
!       used.  (This + 1 is first loc. available for option storage.)
!   The rest should be set as follows, for k = 5, LAST, 3:
!    k-1 Option index, 0 indicates storage required by default.  (This
!        may depend on input parameters, but doesn't depend explicitly
!        on an option.)
!    k   If >0, gives the start of space used by the option.
!        If <0, gives -(amount of space required by an option), for
!               which a starting location is to be determined.
!    k+1 If preceding entry is >0, gives space required by the option.
!        Else it is assumed that the space requested is to be found,
!        and a diagnostic will be given if space can not be found.  Else
!        INTCHK(K+1) IOPT(INTCHK(k+1)) Diagnostic on successful alloc.?
!            0           ----          No
!            <0          ----          Yes
!            >0          /= -1       Yes
!            >0          == -1       No, and IOPT(INTCHK(k+1)) is set
!                                      to starting loc. of space found.
!        When this program finds the space for an option, values of
!        INTCHK(j) for j >= LAST will be used.  INTCHK(k+1) is set
!        temporarily to -(space required) and INTCHK(k-1) is reduced by
!        2000000 to flag the fact that the location index must be saved
!        after INTCHK(LAST).  INTCHK(LAST) is assumed to contain 1 if
!        the largest space requested is required, and to contain
!        -(minimal space needed) if the large amount requested is not
!        essential.
!   On exit, INTCHK(1) is set = -LAST if there was some kind of error.
!   (Unless the user has called MESS this return can't happen.)
!   INTCHK(2) is set to suggest a value for the storage to be declared.
!   The remaining locations are changed as follows:
!    k-1 Negated if there was a problem with this option, or if the
!        following location was initially <0.
!    k   Starting location used or suggested for the option.
!    k+1 Last location used or suggested for the option.
!
!        In addition if values of INTCHK(j) for j >= LAST are used,
!        INTCHK(j), for j = LAST+1, LAST+2, ..., will be set so that
!        INTCHK(j) is equal to (one of the k's above) - 1 for an option
!        that has had a starting location assigned, and INTCHK(LAST) is
!        set to the last index value for j in the above list.
!
! IOPT  This is stored to if INTCHK(k) is < 0, see above.
! ETEXT Input text of the form 'Package_name / Argument_name$E'.  Thus
!   for argument KORD of the package DIVA, this would = 'DIVA / KORD$E'.
!
! ************************** Variable definitions **********************
!
! ERRBAD Flag to use for error severity if get an error.  17 if errors
!        are to print but not stop, 57 otherwise.
! ETEXT  (Input) Used to construct error messages, see above.
! I      Induction variable for accessing INTCHK().
! INTCHK (InOut) See above.
! IOPT   (Out) If space information is being obtained results are saved
!        here.  See above.
! ISTRT  Starting index for placing options with unknown locations.
! KEX    Points to the last place in INTCHK where indices of entries
!        that have had locations determined here are stored.
! L      Temporary index.
! LAST   First free location in INTCHK.  (=INTCHK(1))
! LNEG   Index of the first INTCHK entry that is to be positioned.
! LOPT   Location in IOPT to get location of INTCHK entry that is
!        positioned.
! LTXTAx Variables setup by PMESS defining the locations in MTXTAA where
!        various error messages start.
! LTXTEE Location in MTXTAA where data in ETEXT is stored.
! LTXTZZ Number of characters available for saving ETEXT in MTXTAA.
! LWANT  -Number of locations wanted by INTCHK entry being positioned.
! MACT   Vector used to specify error printing actions, see MESS.
!        MACT(2) flags error/diagnostic level.  = 0 none; 07 is used to
!        get diagnostics only; and ERRBAD otherwise.
! MEEMES Parameter specifying that error message is to be printed.
! MEIDAT Parameter specifying location in INTCHK to start printing
!        integer in MESS.
! MEIMAT Parameter specifying an integer matrix is to be printed.
! MENTXT Parameter specifying location in INTCHK to start printing
!        text from MTXTAA in MESS.
! MERET  Parameter specifying the end of an error message.
! MESS   Routine to print error and diagnostc messages.
! METEXT Parameter specifying that text is to be printed.
! MI     Temporary index used to save in acceptable location.
! MTXTAA Used to contain error message text and instructions, see MESS.
! MTXTAx Character variables setup by PMESS and equivalenced into ETEXT
!        used to contain parts of the error messages.
! MTXTZZ As for MTXTAx, except not setup by PMESS.  Used to hold text
!        from ETEXT.
! MV     Temporary, which contains value associated with INTCHK(MI).
! N      Temporary value.
! NERBAD Array telling what to do concerning errrors.  ERRBAD is set
!        from NERBAD(mod(INTCHK(0), 10)), and the default value for
!        MACT(2) is set from NERBAD(INTCHK(0)+4).
!
!### History
!  * 1998-11-01 OPTCHK  Krogh  ERRSEV => MACT(2) for "mangle".
!  * 1996-05-13 OPTCHK  Krogh  Changes to use .C. and C%%.
!  * 1995-03-10 OPTCHK  Krogh  Added "abs(.) just below "do 140 ..."
!  * 1994-11-11 OPTCHK  Krogh  Declared all vars.
!  * 1993-05-17 OPTCHK  Krogh  Additions for Conversion to C.
!  * 1991-11-25 OPTCHK  Krogh  More comments, little clean up of code.
!  * 1991-10-09 OPTCHK  Krogh  More comments, little clean up of code.
!  * 1991-06-27 OPTCHK  Krogh  Initial Code.

    subroutine OPTCHK(INTCHK, IOPT, ETEXT)

      integer INTCHK(0:*), IOPT(*)
      character ETEXT*(*)
      integer I, ISTRT, KEX, L, LAST, LNEG, LOPT, LWANT, MI, MV, N
! Declarations for error messages.
      integer MENTXT, MEIDAT, MECONT, MERET, MEEMES, METEXT, MEIMAT,    &
     &   LTXTEE, LTXEND
      parameter (MENTXT =23)
      parameter (MEIDAT =24)
      parameter (MECONT =50)
      parameter (MERET =51)
      parameter (MEEMES =52)
      parameter (METEXT =53)
      parameter (MEIMAT =58)
      integer MACT(16), ERRBAD, NERBAD(0:7)
!
! ********* Error message text ***************
![Last 2 letters of Param. name]  [Text generating message.]
!AA OPTCHK$B
!AB "Option #" is negated if option needs attention.$N
!   "Option 0" is for space not associated with a specific option.$N
!   "First Loc." is negated if user did not set value.$N
!   Space avail. = $I; all options must have first loc. > $I$E
!AC Option #$HFirst Loc.$HLast Loc.$E
!AD From subprogram/argument: $B
!AE Space for ETEXT.$E
      integer LTXTAA,LTXTAB,LTXTAC,LTXTAD,LTXTAE
      parameter (LTXTAA=  1,LTXTAB=  9,LTXTAC=233,LTXTAD=267,LTXTAE=295)
      character MTXTAA(2) * (156)
!                          Next 4 lines not automatically generated
!%%    #define LTXTEE  137
      parameter (LTXTEE = LTXTAE - 156 - 2)
      parameter (LTXEND = 156)
!
      data MTXTAA/'OPTCHK$B"Option #" is negated if option needs attenti&
     &on.$N"Option 0" is for space not associated with a specific option&
     &.$N"First Loc." is negated if user di','d not set value.$NSpace av&
     &ail. = $I; all options must have first loc. > $I$EOption #$HFirst$&
     & Loc.$HLast Loc.$EFrom subprogram/argument: $BSpace for ETEXT.$E'/
!
!                      1 2 3      4       5  6       7      8       9
      data MACT / MEEMES,0,1,LTXTAD, MEIDAT, 2, MENTXT,LTXTAB, METEXT,  &
     &   MEIMAT,3,3,0,LTXTAC,-1, MERET /
!            10    13     14 15     16
      data NERBAD / 57, 17, 57, 17, 0, 0, 7, 7 /
!
! *************************** Start of Executable Code *****************
!
      MACT(3) = INTCHK(0) / 10
      MACT(16)=MERET
      I = INTCHK(0) - 10*MACT(3)
      if (I > 3) then
         I = I - 4
         MACT(16) = MECONT
      end if
      ERRBAD = NERBAD(I)
      MACT(2) = NERBAD(I+4)
      LAST = INTCHK(1)
      KEX = LAST
   20 LNEG = 0
      do 100 I = 5, LAST, 3
!       Loop to sort on the low indices -- Inefficient algorithm to keep
!       code short -- LAST should never be very big.
         MI = I
         MV = INTCHK(I)
         do 50 L = I+3, LAST, 3
!                                    Find mimimum from here down.
            if (INTCHK(L) < MV) then
               MI = L
               MV = INTCHK(L)
            end if
   50    continue
         if (MI /= I) then
!                                   Interchange to get low at top.
            do 70 L = -1, 1
               N = INTCHK(I+L)
               INTCHK(I+L) = INTCHK(MI+L)
               INTCHK(MI+L) = N
   70       continue
         end if
         if (MV < 0) then
!                Save loc. of last entry that needs space to be found.
            LNEG = I
         else if (LNEG == 0) then
!        Entry I and previous entries are in their correct sorted order.
            if (INTCHK(I+1) < 0) then
               if (INTCHK(I-1) < -1000000) then
                  INTCHK(I-1) = INTCHK(I-1) + 2000000
                  INTCHK(I+1) = -INTCHK(I+1)
!                            Save INTCHK index defining allocated space.
                  KEX = KEX + 1
                  INTCHK(KEX) = I - 1
               else
!                   Error -- Got request for a negative amount of space.
                  MACT(2) = ERRBAD
                  INTCHK(I-1) = -abs(INTCHK(I-1))
               end if
            end if
!                Save final location used by the option.
            INTCHK(I+1) = INTCHK(I) + INTCHK(I+1) - 1
            if (INTCHK(I) <= INTCHK(I-2)) then
!                                           Error -- options overlap.
               INTCHK(I-1) = -abs(INTCHK(I-1))
               MACT(2) = ERRBAD
            end if
         end if
  100 continue
      if (LNEG /= 0) then
!     Find spaces that need to be allocated, starting with the smallest.
         ISTRT = LNEG
         I = LNEG
  120    LWANT = INTCHK(LNEG)
         LOPT = INTCHK(LNEG+1)
         if (I == LNEG) then
!                         Make fake entry to get started.
            INTCHK(LNEG) = 1
            INTCHK(LNEG+1) = INTCHK(3)
         end if
         do 140 ISTRT = ISTRT, LAST-3, 3
            if(INTCHK(I)+abs(INTCHK(I+1))-LWANT < INTCHK(ISTRT+3))   &
     &         go to 150
            I = ISTRT + 3
  140    continue
  150    INTCHK(LNEG) = INTCHK(I) + abs(INTCHK(I+1))
         if (LOPT /= 0) then
            if (LOPT > 0) then
               if (IOPT(LOPT) == -1) then
                  IOPT(LOPT) = INTCHK(LNEG)
                  go to 160
               end if
            end if
!                     Error -- IOPT not expecting starting loc.
            INTCHK(LNEG-1) = -abs(INTCHK(LNEG-1))
            MACT(2) = ERRBAD
         end if
  160    INTCHK(LNEG+1) = LWANT
         INTCHK(LNEG-1) = INTCHK(LNEG-1) - 2000000
         if (LNEG < 8) go to 20
         I = LNEG
         LNEG = LNEG - 3
         go to 120
      end if
      if (INTCHK(LAST-1) > INTCHK(2)) then
         if (INTCHK(2) < 0) go to 180
         if (LAST /= KEX) then
            if (INTCHK(KEX) == LAST - 3) then
               if (INTCHK(LAST) <= 0) then
                  if (INTCHK(LAST-2)-INTCHK(LAST)-1 <= INTCHK(2)) then
                     INTCHK(LAST-1) = INTCHK(2)
                     go to 180
                  end if
               end if
            end if
         end if
         INTCHK(LAST-3) = -abs(INTCHK(LAST-3))
         MACT(2) = ERRBAD
      end if
  180 if (LAST /= KEX) INTCHK(LAST) = KEX
      if (MACT(2) > 0) then
  190    if (LAST /= KEX) then
            do 200 I = LAST+1, abs(KEX)
               INTCHK(INTCHK(I)+1) = -INTCHK(INTCHK(I)+1)
  200       continue
            if (KEX < 0) go to 210
            KEX = -KEX
         end if
         MACT(13) = (LAST - 4) / 3
!%%       strcpy(&mtxtaa[1][LTXTEE-1], etext);
         MTXTAA(2)(LTXTEE:LTXEND)=ETEXT(1:)
         call MESS(MACT, MTXTAA, INTCHK(1))
         if (MACT(2) > 10) INTCHK(1) = -LAST
         if (LAST /= KEX) go to 190
      end if
  210 INTCHK(2) = INTCHK(LAST-1)
      return
    end subroutine optchk
!*************************************************************************

!*************************************************************************
    end module diva_module
!*************************************************************************