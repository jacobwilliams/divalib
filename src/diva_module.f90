!*************************************************************************
!>
!  Modernized version of the DIVA Variable Order Adams Method
!  Ordinary Differential Equation Solver From the MATH77 library.
!
!  When converting between precisions, don't forget to change the value
!  of KDIM set in parameter statements in a variety of routines, and to
!  adjust comments for the data statements associated with EIBND in
!  [[DIVACR]], and B in [[DIVAHC]].
!
!### Entries
!  * [[DIVA]]    Main entry for starting the package.
!  * [[DIVAA]]   Main program inside the package, calls the other routines,
!                and does checks for output, and noise.  Called by the user
!                if reverse communication is used.
!  * [[DIVABU]]  Back ups the solution to the current base time, if a step
!                that has been started must be taken over for some reason.
!  * [[DIVACO]]  Called by user to get certain information from the common
!                blocks.
!  * [[DIVACR]]  Corrects the solution, estimates errors, and selects order.
!  * [[DIVADB]]  Subroutine to assist in debugging codes.  Called by user to
!                get a formatted list of all the variables used in the
!                integration.  Not required in usual case.
!  * [[DIVADE]]  Needed only for delay differential equations.  This is called
!                by the user from the derivative subprogram.
!  * [[DIVAG]]   Required only if the user has G-Stops, i.e. places to call
!                his output subroutine where certain functions have zeroes.
!  * [[DIVAHC]]  Compute coefficients that depend on the step size history.
!  * [[DIVAIN]]  Used to interpolate to arbitrary points.
!  * [[DIVAOP]]  Used to process user option requests.
!  * [[DIVAPR]]  Used to update the differences and to predict the solution
!                at the end of the current step.
!  * DERIVS (formal) Name of subroutine to be called for computing
!  * [[OPTCHK]]  Used in checking storage allocation.
!  * [[DMESS]]   Used to output error messages and diaganostic messages.
!                (Just [[MESS]] if no floating point is output.)
!  * [[DZERO]]   Called only if [[DIVAG]] is used.  Iterates to find zeros of
!                arbitrary (continuous) functions.
!
!### Common blocks
!  As a left over from the distant past, some variables
!  are in common so that they would be saved.
!
!  * DIVAEV  Holds variables that depend on the environment.
!  * DIVAMC  The main common block for the package.
!  * DIVASC  The secondary common block for the package.  This contains
!            variables that are required for doing interpolation and is
!            separate to simplify saving the variables that are required
!            when the solution is being dumped (saved).
!
!### Original Copyright
!
!  Copyright (c) 1996 California Institute of Technology, Pasadena, CA.
!  ALL RIGHTS RESERVED.
!  Based on Government Sponsored Research NAS7-03001.

    module diva_module

    use messages_module
    use diva_constants

    implicit none

    contains
!*************************************************************************

!*************************************************************************
!>
!  Main routine.
!
!### Notes
!
!```
!++S Default KDIM = 16
!++  Default KDIM = 20
!++  Default MAXORD = 2, MAXSTF = 1
!++  Default INTEGO, VAREQ, OUTPUT, DUMP, GSTOP, EXTRAP
!++  Default STIFF=.F., ARGM=.F., ERRSTO=.F.
!```
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

!### Common variables and local variables
! ALPHA  (DIVAMC) Array with I-th entry = (current step size) / XI(I).
!   Used in computing integration coefficients.
! B      (DIVAHC) Array used to get started on computing integration
!   coefficients.  B(K) = 1. / (K*(K+1))
! BAKMIN (DIVADE) The largest delay at the initial point.
! BETA   (DIVAMC) Array with I-th entry = product (K=1,I-1) of
!   (current (XI(K)) / XI(K) from previous step),  BETA(1)=1.  Used in
!    updating the difference tables.
! C      (DIVAIN) Array used to hold integration/interpolation coeffs.
! C0     Parameter = 0. (in DIVAA,DE,CR,A,G,HC,IN,OP,PR)
! C1     Parameter = 1. (in DIVA,A,CR,DA,HC,IN,OP)
! C10    Parameter = 10. (in DIVAA,CR,OP)
! C1000  Parameter = 1000. (in DIVACR)
! C16    Parameter = 16. (in DIVAA,OP)
! C1M3   Parameter = .001 (in DIVAA)
! C1M5   Parameter = .00001 (in DIVAA)
! C1P125 Parameter = 1.125 (in DIVAA,HC,OP)
! C1P3   Parameter = 1.3 (in DIVAA)
! C1P4   Parameter = 1.4 (in DIVACR)
! C2     Parameter = 2. (in DIVAA,DE,BU,CR,IN,OP)
! C20    Parameter = 20. (in DIVACR)
! C2P5M3 Parameter = .0025 (in DIVAA)
! C4     Parameter = 4. (in DIVACR,OP)
! C40    Parameter = 40. (in DIVACR)
! C4096  Parameter = 4096. (in DIVAA)
! C6     Parameter = 6. (in DIVAA)
! C8M3   Parameter = .008 (in DIVAA)
! CM2    Parameter = -2. (in DIVACR)
! CM8    Parameter = -8. (in DIVACR)
! CMP5   Parameter = -.5 (in DIVACR)
! CMP75  Parameter = -.75 (in DIVAOP)
! CP0625 Parameter = .0625 (in DIVAA)
! CP1    Parameter = .1 (in DIVAA,CR,DA,HC)
! CP125  Parameter = .125 (in DIVACR)
! CP25   Parameter = .25 (in DIVAA,CR,DE,OP)
! CP3    Parameter = .3 (in DIVAA,OP)
! CP4    Parameter = .4 (in DIVAA)
! CP5    Parameter = .5 (in DIVAA,CR,DA,DE,HC,OP)
! CP5625 Parameter = .5625 (in DIVAHC)
! CP625  Parameter = .625 (in DIVAOP)
! CP75   Parameter = .75 (in DIVACR,OP)
! CP8    Parameter = .8 (in DIVACR)
! CP875  Parameter = .875 (in DIVAA, OP)
! CP9    Parameter = .9 (in DIVAOP)
! CP9375 Parameter = .9375 (in DIVACR)
! CQ3125 Parameter = .03125 (in DIVACR)
! CRBQI  Parameter = .421875 (in DIVAHC)  Initial val for computing RBQ.
! CSUM   (DIVAIN) Array used to contain partial sums of the integration
!   coefficients.  This is used to corrrect for a difference table that
!   has not yet been updated.
! D      (DIVAMC) Array to be used later to store coefficients for
!   integrating stiff equations.
!   derivatives.  Not used if option 13 is set.
! DISADJ (DIVAA) Value of stepsize when discontinuity is indicated.
! DNOISE (DIVAMC) Used in determining if noise is limiting the
!   precision.  It is usually |highest difference used in correcting|
!   of the equation with the largest error estimate.
! DS     (DIVAMC) Array to be used later to store coefficients for
!   estimating errors when integrating stiff equations.
! DVC2   (DIVADB) Array used for output of variables HC to TOUT in
!   common block DIVAMC.
! E      (DIVACR) (Estimated error) / (Requested accuracy)
! EAVE   (DIVAMC) This is a weighted average of past values of EIMAX.
!   It is adjusted to account for expected changes due to step changes.
! EEPS10 (DIVAEV) = 10. * (machine epsilon).
! EEPS16 (DIVAEV) = 16. * (machine epsilon).
! EEPS2  (DIVAEV) =  2. * (machine epsilon).
! EEPT75 (DIVAEV) = (machine epsilon) ** (.75)
! EI     (DIVACR) Estimate for what E would be if step size increased.
! EIBND  (DIVACR) Array containing limits on the estimated error with
!   the stepsize increased.  This array tends to make the code a little
!   more conservative on step size increases at low order.
! EIMAX  (DIVAMC) Estimate of (error estimate / error requested) if the
!   step size should be increased.
! EIMIN  (DIVAMC) An error estimate is small enough to allow a step
!   increase if the estimate of ((error with the step size increased) /
!   (error requested)) is less than EIMIN.
! EIMINO (DIVAA) Set to C8M3 and never changed.  When step size is being
!   reduced if EIMIN <= EIMINO then the reduction factor is set to
!   CP875.  This variable could be a parameter.
! EMAX   (DIVAMC) Largest value computed for (error estimate) / (error
!   requested).
! EOVEP2 (DIVAEV) = EEPS2 * (largest floating point number).
! EPS    (DIVACR) Current absolute error tolerance.  Also used for
!   temporary storage when computing the desired value of EPS.
! ERCOEF (DIVACR) (Error coefficient from formula) / EPS
! EREP   (DIVAMC) If EMAX > EREP, a step is repeated.  Ordinarily
!   this has the value .3.  This is set < 0 if the error tolerance is
!   specified improperly, and is set to a large value if the user
!   requests complete control over the step size.  EREP is also set
!   < 0 after a user specified discontinuity.
! EROV10 (DIVAEV) = 10. / (largest floating point number).
! ETA    (DIVAIN) Array used in computing integration/interp. coeffs.
! EVC    (DIVADB) Array used for output of variables EEPS2 to EROV10 in
!   common block DIVAEV.
! EXR    (DIVAA) Set to CP1 and never changed.  If it is estimated the
!   the (error estimate) / (error requested) on the next step will be
!   >= EXR then the step size is reduced.  Could be a parameter.
! F      (formal) Array used to store derivative values, the difference
!   tables, error tolerance requests, and values used by some other
!   options. (in DIVA,A,BU,CR,DA,DB,G,IN,PR)
! FDAT  (DIVAMC) Used to store data for error messages.  (Local array in
!   DIVAIN.)
! FOPT  (formal) in DIVAOP.  Passed as place to save floating point data
!   for options.  This package passes F in for FOPT when calling DIVAOP.
! G      (DIVAMC) Integration coefficients used for predicting solution.
!   G(I, J) gives the I-th coefficient for integrating a J-th order
!   differential equation.  G(1, 1) is equal to the step size.
! GAMMA  (DIVAIN) Array used in computing integration/interp. coeffs.
! GG     (DIVAHC) Array of length = max. differential equation order
!   allowed by code - 1.  GG(K) = (HH**(K+1)) / K!
! GNEW   (formal) in DIVAG.  Current value for vector function g, whose
!   zeroes are to be found.
! GOINT  (DIVACR) Used for assigned go to used in computing integration
!   coefficients.
! GOLD   (DIVAG) Previous value for element of G whose zero search is
!   active.
! GS     (DIVAMC) Integration coefficients used in estimating errors.
! GT     (formal) in DIVAG.  Previous value of GNEW.
! HC     (DIVAMC) Ratio of (new step size) / (old step size)
! HDEC   (DIVAMC) Default value to use for HC when reducing the step
!   size.  (Values closer to 1 may be used some of the time.)
! HH     Equivalenced to G(1,1) = current step size in DIVAA,CR,DA,G,HC.
! HI     (DIVAIN) Step length from the base value of the independent
!   variable for the interpolation.
! HINC   (DIVAMC) Default value to use for HC when increasing the step
!   size.  (Values closer to 1 may be used some of the time.)
! HINCC  (DIVAMC) Actual value used for default value of HC when
!   increasing the step size.  Set to HINC after start is considered
!   complete.  During the start HINCC is set to 1.125.
! HMAX   (DIVAMC) Largest value allowed for abs(step size).  Default
!   value is a very large number.
! HMAXP9 (DIVAMC) .9 * HMAX.
! HMIN   (DIVAMC) Smallest value allowed for abs(step size).  Default
!   value is 0.
! HNEW   (DIVADE) Value of step size when iterating at initial point
!   for delay differential equations.
! I      Used for temporary storage. (DIVAA,BU,CR,DA,DE,G,IN,OP,PR)
! IA     (DIVAOP) absolute value of first integer stored for an option.
! ICF    (DIVAMC) Final index for current loop in DIVACR.  Required by
!   option 18.
! ICI    (DIVAIN) Temporary index, = 0 for interpolation, 1 or 0 for
!   differentiation, and d-1, d-2, ... 0 for integration, where d is the
!   order of the differential equation.  Index of first location
!   in C() used is ICI + an offset.
! ICS    (DIVAMC) Starting index for current loop in DIVACR.
! ID     (formal) Array use to contain integer data from common.  Values
!   are returned in locations 1 to 5 as follows.
!   1    KEMAX = Index of equation with largest error estimate
!   2    KSTEP = Current step number
!   3    NUMDT = Number of differences used for each equation
!   4            Reserved for future use
!   5            Reserved for future use
! IDAT   (DIVAMC) Used to store integer for error messages.  (Also used
!   in DIVAA for temporary storage of KORD(2).  (Local array in DIVAIN.)
! IDE    (DIVADE - formal) Array used to contain past information so
!   that delays can stretch back indefinitely.  If the first location is
!   0, then any interpolations requested must be in the range of the
!   current difference tables.  At present, only the value 0 is allowed
!   in IDE(1).  This array is intended for the support of saving long
!   past histories.  IDE(2) must contain the declared dimension of WDE.
! IDEF   (DIVADE -- formal) Flag giving indicaion of what is going on.
!   = 0  User should compute derivatives and return to the main
!        integrator.
!   = 1  Code is computing additional values in order to get past data
!        necessary for starting.  User should compute derivatives and
!        call DIVADE.
!   < 0  Indicates an error condition.  If DIVADE is called without
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
! IDT    (DIVAIN) Used as a base index into the difference table.
! IFLAG  (formal in DIVAG) Used for communication with user.
!   = 1  Continue as if DIVAG was not called.
!   = 2  Check KORD(1) as one would do at start of OUTPUT if no G-Stops
!        were present. (Exit if in DERIVS.)
!   = 3  Return to the integrator.
!   = 4  Compute G and return to DIVAG.
!   = 5  A G-Stop has been found, and NSTOP gives its index.  (If NSTOP
!        < 0, the stop was an extrapolating stop.)
!   = 6  Same as 5, but requested accuracy was not met.
!   = 7  Same as 5, but there is a probable error in computing G.
!   = 8  Fatal error of some type.  (An error message has been printed.)
! IG     (DIVAG)  IG = KORD(2) on the initial entry (0 for extrapolating
!   G-Stops, and 1 for interpolating).
! IGFLG  (DIVAMC) Used primarily in DIVAg, but also used in DIVA to keep
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
! IGSTOP (DIVAMC) IGSTOP(k) is set in DIVAg to the index of the last G
!   with a 0, where k is one for an interpolatory G-Stop, and k is two
!   for an extrapolatory G-Stop.
! IGTYPE (DIVAMC) Array with two elements as for IGSTOP, but this saves
!   a flag giving the nature of convergence to the stop.
!   = 0  All known G-stops completely processed.
!   = 4  Need to compute next value while iterating.
!   = 5  Got good convergence.
!   = 6  Got convergence, but not to desired accuracy.
!   = 7  Problem in getting convergence.
!   = 8  A fatal error of some type.
! IHI    (DIVA) Last location used by the current option.
! ILGREP (DIVAMC) Used when correction to keep track of equations that
!   are to use a certain error tolerance.
! ILGROR (DIVACR) Index of last equation in the current group of
!   equations grouped for selecting integration order.
! ILOW   (DIVA) First location used by the current option.
! INCOM  (DIVADE) Array equivalenced to LDT in the common block DIVASC.
!   Used to simplify saving information in the common block.
! INCOP  (DIVAOP) Array containing data giving the amount of space in
!   IOPT used for each of the options.
! INGS   Current index for G-stop being examined in DIVAG.
! INICAS (DIVADE) Used to track the initialization for a delay equation.
!   = 1  Very beginning.
!   = 2  Getting derivative at the very beginning.
!   = 3  Getting derivatives at points prior to the initial point.
!   = 4  Getting derivative at initial point after iteration is started.
! INTCHK (DIVA) Array passed to OPTCHK containing information on storage
!   allocation.  See comments in OPTCHK for details.
! INTEG  (DIVAIN) Number of integrations being done. (<0 for
!   differentiations and =0 for interpolation.)  Also used as counter
!   when computing integration coefficients.
!        (DIVAPR) Number of integrations being done.
! INTEGS (DIVAPR) = -1 for equations that are not stiff, 0 for those
!   that are stiff.
! INTEGZ (DIVAIN) min(INTEG, 0)
! INTERP (DIVAIN) added to the usual integration order to get the order
!   to be used when interpolating: 3-KQMAXI, if HI=0; 1, if
!   |HI| > |XI(1)| and HI * XI(1) < 0; 0, otherwise -- the usual case.
! IOP10  (DIVAMC) Number of times diagnostic output is to be given when
!   leaving DIVAcr (the corrector).
! IOP11  (DIVAMC) Gives current step number of the method.  Tells how
!   many of certain coefficients must be computed. (Has nothing to do
!   with options.) = min(max integ order + 1, KDIM).  Also set when
!   starting to flag that certain memory locations must be set to 0.
! IOP12  (DIVAMC) Points to location in F() where user supplied values
!   of HINC, HDEC, HMIN, and HMAX.  (0 if option 12 not used.)
! IOP13  (DIVAMC) If not zero, reverse communication will be used for
!   getting the values of derivatives.  Associated with option 13.
! IOP14  (DIVAMC) If not zero, reverse communication will be used in
!   place of calls to the output routine.  Associated with option 14.
! IOP15  (DIVAMC) If not zero, a return will be made to the user after
!   the initialization.  Associated with option 15.  This might be used
!   to overlay DIVA, some of the user's code, and perhaps DIVAop.
! IOP16  (DIVAMC) Points to location in KORD() where information for
!   specifying the error tolerance is specified.  See option 16.
! IOP17  (DIVAMC) Used in initialization for option 17, afterwards this
!   cell is used by KEXIT which is equivalenced to IOP17.
! IOP18  (DIVAMC) Points to location in KORD() where information for
!   specifying a grouping of equations for derivative evaluation is
!   stored.  See option 18.
! IOP19  (DIVAMC) Points to location in KORD() where information for
!   specifying a grouping of equations for integration order control
!   is stored.  See option 19.
! IOP20  (DIVAMC) Used for option 20, gives first location in F where
!   estimated errors are to be stored.  Expected to be useful in a
!   program for solving boundary value problems using multiple shooting.
! IOP21  (DIVAMC) Was used for stiff equations option (never completely
!   coded).  The optional code still uses this (don't activate it!).
!   Now used to flag the location if F where the user has stored the
!    tolerance to use in finding G-Stops.
! IOP21S (DIVAMC) Was used for stiff equations see above.
! IOP22  (DIVAMC) Set aside for possible option for stiff equations.
! IOP3   (DIVAMC) Value set by option 3.
!   =  0 Interpolate to final point. (The default)
!   =  1 Integrate to final point.
!   = -1 Extrapolate to final point.
! IOP4   (DIVAMC) Value set by option 4.  The output routine is called
!   with KORD(1) = 4, every IOP4 steps.  (Default value for IOP4 is a
!   very large number.
! IOP5   (DIVAMC) Value provided by option 5, used to specify extra
!   output points.
! IOP6   (DIVAMC) Value provided by option 6.  If nonzero, the output
!   routine is called at the end of every step.  If > 0, there are
!   IOP6 interpolating G-Stops.
! IOP7   (DIVAMC) Value provided by option 7.  If > 0, there are K7
!   extrapolating G-Stops.
! IOP8   (DIVAMC) Value provided by option 8.  If nonzero, the output
!   routine is called with KORD(1)=8 whenever the step size is changed.
! IOP9   (DIVAMC) Value provided by option 9.  Used to specify that the
!   user wishes to save the solution.
! IOPIVA (DIVA) Used to save length of IOPT vector for error messages.
! IOPST  (DIVASC) Intended for possible use in stiff equations.
! IOPT   (formal DIVA and IVAOP) Used to specify options.
! IOPTC  (DIVAOP) In DIVAOP equivalenced so that IOPTC(3) is equivalent
!   to IOP3.
! IOPTS  (DIVAOP) Array containing the current default values to be
!   stored into IOPTC.
! IORD   (DIVACR) Index of first equation in the current group of
!   equations grouped for selecting integration order.
! IOUTKO (DIVADC) Used in DIVADI to point to KORD to keep track of
!   equation grouping for diagnostic output.
! ISVCOM (DIVADE) Used to save info. in the common block DIVASC.
! ITERS  (DIVADE) Counts iterations in starting delay differential
!   equations.  Max. value for this is arbitrarily 100.
! ITOLEP (DIVAMC) Used for temporary storage, and for the index of a
!   tolerance relative to the start of tolerances.
! IVC1   (DIVADB) Array used for output of variables IOPST to NUMDT in
!   common block DIVASC.
! IVC2   (DIVADB) Array used for output of variables ICF to NY in
!   common block DIVAMC.
! IWB    (DIVADE) Current base index for saving F values in WDE when
!   starting delay differential equations.
! IY     (DIVAMC) Used for the current index to the Y() array.  (Local
!   variable in DIVAIN used in computing IYI.)  Equivalenced to
!   IZFLAG in DIVAG.
! IYI    (DIVAIN) Y(IYI) is currently being computed.
! IYN    (DIVAIN) Y(IYN) is base Y() corresponding to Y(IYI).
! IYNI   (DIVAIN) Used as base index for computing IYN as IY is for INI.
! IYO    (DIVADE) Points to first base value of Y for current
!   interpolation when getting values for a delay differential equation.
! IZFLAG (DIVAG)  Equivalenced to IY.  Set to 0 initially, and later
!   set to the value returned by *ZERO.
!    = 0  Value set on entry at start of search.
!    = 1  Compute next g again.
!    = 2  Normal terminiation.
!    = 3  Normal termination -- error criterion not satisfied.
!    = 4  Apparent discontinuity -- no zero found.
!    = 5  Couldn't find a sign change.
!    = 6  *ZERO was called with a bad value in IZFLAG.
! J      For temporary storage. (In DIVA,A,BU,CR,DA,DB,DE,HC,IN,OP,PR)
! J1     (DIVAA & DA) Used for temporary storage.
! J2     (DIVAA) Used for temporary storage.
! JL     (DIVA) Used for checking storage.
! JLGREP (DIVACR) Contents of first location of KORD (called LGROUP in
!   DIVACR) for the current error tolerance rule.
! JLGROR (DIVACR) Contents of first location of KORD for the current
!   integration order control.
! JLIM   (DIVA) Used for checking second item in KORD list for options
!   16 and 19.
! K      For temporary storage.  (In DIVA,A,BU,CR,DA,DB,DE,HC,IN,OP,PR)
! KDIM   Parameter giving the largest number of differences supported.
!        Used in all the routines.
! KEMAX  (DIVAMC) Index associated with equation giving the largest
!   value for (estimated error) / (requested error).
! KEXIT  (DIVAMC) Equivalenced to IOP17 which is not used after
!   initialization.  Defines actions when KORD2I = -7.  (Referenced in
!   (DIVAA,DA,G).)
!   =  1  Take the step over with reduced H.
!   =  2  Take the step over.
!   =  3  Do the end of step call to OUTPUT.
!   =  4  Reset TMARK, then do same as for KEXIT = 2.
!   =  5  Reset TMARK, then do same as for KEXIT = 3.
!   =  6  Give the fatal error diagnostic.
! KFERR  (DIVA)  Temporary storage in checking for option 16.
! KGO    (DIVA)  Used to tell from whence a check is being done or an
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
! KIS    (DIVAMC) Used to check if it is time to dump the solution.
!   The check involves incrementing KIS at the end of the step, and
!   dumping the solution if KIS is 0.
!   = -1  Set in DIVAcr when it is time to dump solution
!   =  0  When starting
!   =  2  After being dumped.
!   This is set to 1000 just after a user specified discontinuity, and
!   counted up from that point.
! KMARK  (DIVAMC) Identifies the type of output associated with the next
!   output point specified by TSPECS.
! KONV   (DIVADE) Counts iterations.  Test for convergence if KONV > 1.
! KORD   (formal in DIVA,A,BU,CR,DA,DB,DE,G,IN,PR) KORD(1) is used to
!   return flags for the user to test, and KORD(2) tells what routine
!   the flag is associated with.  See KORD1I and KORD2I below and the
!   write up for the program.  KORD(3) is used for communicating extra
!   information to the user in some cases.  KORD(4) to KORD(NTE+3) are
!   used for integration order for the equations, and the rest of KORD()
!   is available for user options.
! KORD1I (DIVAMC) Helps in defining the state of the integrator.
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
!   = 24  Set on error message in DIVA, along with KORD2I = -4.
!   Also used as an index into MLOC in DIVAA when an error is being
!   processsed, see MLOC below.
! KORD2I (DIVAMC) Helps in defining the state of the integrator.
!   Frequently has the same value as KORD(2).
!   = -3  Set in DIVAg, to get a derivative evaluation.
!   = -2  Set in DIVAg, to get another entry to OUTPUT.
!   = -1  Return to calling program, done, interrupt, or got an error.
!   =  1  Calling OUTPUT or returning to user for OUTPUT type action.
!   =  0  Calling DERIVS or returning to user for DERIVS type action.
!   = -4  Error message in DIVA and in DIVAop, along with KORD1I = 24.
!   = -5  Starting
!   = -6  Starting, getting the initial derivative value or derivatives
!         for the noise test.
!   = -7  Done some extrapolation, KEXIT defines the action to take.
!         Set in DIVAg to activate KEXIT action in DIVA.
!   = -8  Set when user has requested adjustment of the difference
!         tables for a discontinutiy.
! KORDI  (DIVASC) Order of differential equation being integrated.  If
!   all orders are the same, this set once at the beginning.
! KOUNT   (DIVADE) Count of number of points back from the initial point
!   when solving a delay differential equation.
! KOUNTM  (DIVADE) Largest value currrently allowed for KOUNT.
! KOUNTX  (DIVADE) Largest value allowed for KOUNTM.
! KOUTKO  Used in DIVACR to track where output is wanted.
! KPRED  (DIVAMC) Value assigned to KORD1I when getting a predicted
!   derivative.  (1 used now, 5 planned for use with stiff equations.)
! KQD    (DIVACR) = max(2, integration order)
! KQDCON (DIVAMC) Number of coefficients computed with constant step
!   size for stiff equations.
! KQICON (DIVAMC) Number of coefficients computed with constant step
!   size for nonstiff equations.
! KQL    (DIVACR) Integration order at start of (DIVACR)
! KQLORD (DIVACR) Saved value of KQL when equations are grouped for
!   controlling the integration order.
! KQMAXD (DIVASC) Maximum integration order used for stiff equations.
! KQMAXI (DIVASC) Maximum integration order used for nonstiff equations.
! KQMAXS (DIVAMC) Maximum integration order for equations that have
!   some limit on the error that can be committed.
! KQMXDS (DIVAMC) Used to save KQMAXD in case step is repeated and the
!   solution must be dumped.
! KQMXI  (DIVAIN) Maximum integration order used for integration or
!   interpolation, = KQMAXI+INTERP-1.
! KQMXS  (DIVAIN) Maximum step number, = max(KQMXI, KQMAXD).
! KQMXIL (DIVAMC) Value of KQMAXI the last time integration coefficients
!   were computed.
! KQMXIP (DIVAMC) = KQMAXI + MAXINT, for computing integration coeffs.
! KQMXIS (DIVAMC) Used to save KQMAXI in case step is repeated and the
!   solution must be dumped.
! KQN    (DIVACR) Value of integration order at end of DIVACR.
! KQQ    Used for the integration order for current equation.  (Values
!   < 0 are intended for stiff equations.)  (In DIVA,BU,DA,IN,PR)
! KSC    (DIVAMC) Number of steps that have been taken with a constant
!   step size.
! KSOUT  (DIVAMC) When KSTEP reaches this value, the output routine is
!   called with KORD(1) = 4.  The default value is a very large number.
! KSSTRT (DIVAMC) Set when ending one derivative per step to KSTEP + 2.
!   Checked later in DIVAHC to decide whether to set the step changing
!   factors to their nominal values.
! KSTEP  (DIVAMC) Number of steps taken since the start of integration.
! L      Used for temporary storage.  In DIVAIN, L is the initial value
!   of LDT, except L=1 if LDT=-1, and MAXINT >= 0.  (Used in DIVAA,BU
!   CR,DA,DB,IN,PR.)
! LAHAG  (DIVADB) Used to get proper offset into an diagnostic message.
! LAIAG  (DIVADB) Used to get proper offset into an diagnostic message.
! LDIS   (DIVAA) Count of steps since user flagged a discontinuity.
! LDT    (DIVASC) Used to keep track of state of difference table.
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
! LEX    (DIVAMC) Indicates how to get values at next output point:
!   = -1  Extrapolate
!   =  0  Interpolate (The usual case.)
!   =  1  Integrate to the output point, integration is not continued.
! LGO    (DIVAIN) Used as an an assigned go to.  Result is to add in
!   extra correction term when LDT has been set to 2.
! LGROUP (formal) This is a part of KORD passed into DIVACR.  The first
!   location is the start of the information on the grouping of
!   equations for error control.
! LINC   (DIVAMC) Used to indicate state of step size selection.
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
! LINCD  (DIVAMC) Value of smallest k for which HINCC**k >= 2.
!   (=-2 if user is specifying all step size changes.)
! LINCQ  (DIVAMC) Value of smallest k for which HINCC**k >= 4.
! LIOPT  (DIVAOP) Value of the last index in IOPT on the last call.
!   Used so DIVA can print IOPT in error messages.
! LL     (DIVACR) Temporary variable used when equations are grouped
!   for integration order control.
! LNOTM1 (DIVAIN) Logical variable = L /= -1.  If LNOTM1 is true,
!   storage in Y() is different in some way lost to antiquity.  Such
!   a case can only arise in the case of stiff equations.
! LOCF1  (DIVADB) Gives packed data needed for output of tables by the
!   message processor MESS.  See comments there under METABL for defs.
! LOCF2  (DIVADB) As for LOCF1 above.
! LOCM   (DIVAA) Parameter = 32*256, used to unpack integers stored
!   in MLOC for use in error message processing.
! LPRINT (formal, DIVADB) Defines how much printing is to be done in
!   DIVADB.  Let |LPRINT| = 10*N1 + N2     (N1,N2 digits)
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
! LSC    (DIVAMC) Indicates if starting or if noise may be present.
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
! LX     (DIVAA) Used for temporary storage in computing TMARKA().
!        ( formal DIVADE)  An integer array containing extra
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
! LX2    (DIVADE) Value of LX(5i+2), when working on the i-th delay.
! MACT   Used in the programs which call the error message program.
!   This array difines the actions to be taken by that program.  (In
!   (DIVA,A,DA,DE,G,IN,OP)
! MACT0  (DIVADB) Used to call the message program, see MACT.
! MACT?  As for MACT, in (DIVA,CR,DB)
! MACTFV (DIVADB) As for MACT0.
! MAXDIF (DIVASC) Maximum differentiations required for stiff equations.
! MAXINT (DIVASC) Maximum integrations required.  (= max. order of
!   differential equations if equations are not stiff.)
! MAXKQ  (DIVA, BU)e
! MAXKQD (DIVAMC) Largest integration order allowed for stiff equations.
! MAXKQI (DIVAMC) Largest integ. order allowed for nonstiff equations.
! ME???? Parameters defining constants used for interaction with the
!   error message program MESS.  See comments there for definitions.
!   (In DIVA,A,DA,DE,G,IN,OP)
! METHOD (DIVAMC) Defines kind of methods being used.
!   = -1  Only stiff equations are being integrated.
!   =  0  Only nonstiff equations are being integrated.
!   =  1  Both kinds of methods are required.
! MLOC   (DIVA,A,DE) Contains locations in MTEXT for error messages.  In
!   DIVAA this data is packed using MLOC??, see below.
! MLOC?? (DIVAA) Parameters constructed to aid in making up packed data
!   for processing error messages.  Low two digits give the value of
!   KORD1I to use for the error index and later processing, the next two
!   give the error severity level, and the rest point to text used for
!   the message.
! MODF2  (DIVADB) Used in constructing the same kind of packed data as
!   described for LOCF1 above.
! MULTJ  Local to DIVAOP for calls not using F.
! MTEXT  (DIVA,A,CR,IN,OP) Text for error messages.
! MTXT?? (DIVA,A,CR,DA,DB,DE,G,IN,OP) Equivalenced into MTEXT.
! N      Used for temporary storage.  (In DIVAHC,IN,PR)
! NDTF   (DIVASC) Location in F() where difference table starts.
! NE     (DIVAMC) Number of equations in the first group.  (=NTE if
!   option 18 is not used.)
! NEDDIG (DIVADB) Parameter = -MEDDIG.
! NEPTOL (DIVAMC) Used for temporary storage and to save the value of
!   ITOLEP for error messages.
! NEQ    (formal) Total number of equations being integrated.
! NG     (DIVAMC) Used in DIVAg for the number of g's in the current
!   context.
! NGSTOP (DIVAG) Dimension 2 array equivalenced to IOP6, and IOP7.  To
!   get the number of interpolating and extrapolating G-Stops.
! NGTOT  (DIVAMC) NGTOT(1) gives the number of interpolating G-Stops,
!   and NGTOT(2) gives the number of extrapolating G-Stops.
! NKDKO  (DIVASC) If this is nonzero (option 17), it gives the location
!   in KORD() where a vector defining the order of each equation is
!   specified.
! NLX    (DIVADE) Temporary index used to keep track of interpolations
!   being done to get Y() values for a delay differential equation.
! NOISEQ (DIVAMC) max(2, order of equation for which (error estimate)/
!   (error requested) is a maximum).
! NOUTKO (DIVAMC) If nonzero, gives the index in KORD where information
!   on what equations are to be included in the diagnostic output is
!   given.   See option 10.
! NSTOP  (formal) In DIVAG.  Index of the G-stop, see IFLAG.
! NTE    (DIVASC) Total number of equations being integrated = NEQ.
! NTEXT  (formao DIVADB) Character variable containing heading text.
! NTOLF  (DIVAMC) First location in F() where tolerance specifying
!   accuracy desired is stored.
! NUMDT  (DIVASC) Maximum allowed number of differences available for
!   doing an integration.
! NXTCHK (DIVA) Equivalenced to INTCHK(1), which gives the next
!   available location in INTCHK for storing data on storage allocation.
! NY     (DIVAMC) Total order of the system.
! NYNY   (DIVASC) Location in Y() where the base value for Y() is saved.
! NYNYO  (DIVADE) Equivalenced to the saved value from common of NYNY.
! OUTPUT (formal) Name of subroutine to be called for the output of
!   data or for computing G-Stops.  Not used if option 14 is set.
! OVD10  (DIVAEV) (largest floating point number) / 10.
! OVTM75 (DIVAEV) (largest floating point number) ** (-.75)
! RBQ    (DIVAMC) Array containing data for the preliminary noise test.
! RD     (formal DIVACO) Array use to contain floating point data from
!   common.  Values are returned in locations 1 to 3 as follows.
!   1    EMAX =  Max. ratio of estimated error to requested error
!   2            Reserved for future use
!   3            Reserved for future use
! REF    (DIVACR) Array of length 3 used for translating error tolerance
!   type into the factor used for exponential averaging for that type.
! RND    (DIVACR) Usually the current estimated error.  Used in deciding
!   if noise is limiting precision.
! RNOISE (DIVACR) Value used in comparison with RBQ() for preliminary
!   noise test.
! ROBND  (DIVAMC) Used to influence the selection of integration order.
!   The larger ROBND, the harder it is to increase the order and the
!   easier it is to decrease it.
! RVC2   (DIVADB) Array used for output of variables DNOISE to SNOISE in
!   common block DIVAMC.  These are variables that don't require a great
!   deal of precision.
! S      (DIVACR) Estimate of (step size) * eigenvalue of Jacobian.
! SIGMA  (DIVAMC) The k-th entry of this array contains a factor that
!   gives the amount the k-th difference is expected to increase if the
!   step size in increased.  These numbers get bigger it there is a past
!   history of increasing the step size.
! SIGMAS (DIVAA) Saved value of SIGMA(k) from the last step, where k =
!   integration order for equation with index KEMAX.
! SNOISE (DIVAMC) Value used in comparison with RBQ() on equation with
!   largest value for (error estimate) / (error requested).
! T      (formal) in DIVAIN. T(1) contains the point to be interpolated
!   to, and T(2) is used in a check that |HI| <= |T(2)|.  When used by
!   other routines in this package, TSPECS is passed in for T.
! TB      (DIVADE) Base time for current interpolation.
! TC      (DIVADE) Original value of TN when getting past Y's for a
!   delay differential equation.
! TEMP   Used for temporary storage, in DIVAHC,PR
! TEMPA  (DIVACR) Array equivalenced to (TPS1,TPS2,TPS3,TPS4).
! TEMPAO (DIVACR) Array used to accumulate values in TEMPA.
! TG     (DIVAMC) TG(1) gives the last value of TSPECS(1) for which an
!   interpolatory G-Stop has been computed.  TG(2) is defined similarly
!   for extrapolatory G-Stops.
! TGSTOP (DIVAMC) TGSTOP(1) gives the value of TSPECS(1) where the last
!   0 for an interpolatory G-Stop was found.  TGSTOP(2) is defined
!   similarly for extrapolatory G-Stops.
! TMARK  (DIVAMC) Location of the next output point.
! TMARKA (DIVAA)  Array of length 2 equivalenced to TMARK (and TMARKX).
! TMARKX (DIVAMC) Location of the next output point to be found using
!   integration or extrapolation.  This variable must follow immediately
!   after TMARK in the common block.
! TN     (DIVASC) The value of TSPECS(1) at the conclusion of the last
!   step.
! TNEQ   (DIVADB) Array of dimension 1 equivalenced to TN so that an
!   array can be passed to *MESS.
! TOL    (formal) This is a part of F passed into DIVACR.  The first
!   location is the start of the information on the tolerances for error
!   control.
! TOLD   (DIVAG) Value of TSPECS(1) on one side of a zero.
! TOLG   (DIVAMC) Tolerance to pass to dzero when locating G-Stops.
! TOUT   (DIVAMC) Location of next output point defined by value of
!   TSPECS(3).  Such output is given with KORD(1) = 2.
! TP     (DIVA,A,DA,DE,HC) Used for temporary storage.
! TP1    (DIVAA,DA,HC,IN,PR) Used for temporary storage.
! TP2    (DIVAA,DA,HC,PR) Used for temporary storage.
! TP3    (DIVAA) Used for temporary storage.
! TPD    (DIVABU) Used for temporary storage.
! TPP    (DIVACR) Used for temporary storage.  Usually same as TPS3.
! TPS1   (DIVAA,CR) Used for temporary storage.  (In DIVACR is the
!   difference of order KQQ-2)
! TPS2   (DIVAA,CR) Used for temporary storage.  (In DIVACR is the
!   difference of order KQQ-1)
! TPS3   (DIVACR) Contains the difference of order KQQ.  This is the
!   last difference used in the corrector.
! TPS4   (DIVACR) Contains the difference of order KQQ+1.
! TPS5   (DIVACR) Temporary storage.
! TPS6   (DIVACR) Temporary storage.
! TPS7   (DIVACR) Temporary storage.
! TSAVE  (DIVAG) Value of TSPECS(1) before starting the search for a 0.
! TSPECS (formal DIVA,A,DB,DE,G)
!   TSPECS(1) is the current value of the independent variable.
!   TSPECS(2) is the current value of the step size.
!   TSPECS(3) is the increment to use between output points that give
!             output with KORD(1) = 2.
!   TSPECS(4) is the "final" output point.
! V      (DIVAMC) Array used in computing integration coefficients.
! XI     (DIVASC) XI(K) = TSPECS(1) - value of TSPECS(1) K steps
!   previous.
! W      (DIVAHC) Array used in computing integration coefficients.
! WDE    (formal, DIVADE)  Array used for working storage.  This storage
!   is used to save derivative values when iterating to get started.  To
!   be safe one should allow as much space as is allowed for differences
!   in F.  In most cases the start will not require this much space
!   however.  This array is also intended for the support of saving long
!   past histories.
! Y      (formal, DIVA,A,CR,DA,DB,DE,G,IN,PR) Array containing the
!   independent variable and all derivatives up to order one less than
!   the order of the differential equation.  Also use to save these
!   values at the beginning of the current step, the base values.
! YN     (formal, in DIVAPR)  Base values of y, these follow the
!   current values of the dependent variable, y, in Y().


    subroutine DIVA(TSPECS, Y, F, KORD, NEQ, DIVAF, DIVAO, IDIMT,     &
                    IDIMY, IDIMF, IDIMK, IOPT)

      use divaev_module
      use diva_constants
      use divasc_module
      use divamc_module

      integer NEQ, IDIMT, IDIMY, IDIMF, IDIMK
      integer KORD(*), IOPT(*)
      double precision TSPECS(*), Y(*)
      double precision F(*)
      external DIVAF, DIVAO

      integer KGO, INTCHK(0:30), NXTCHK
      integer IHI, JL, J, ILOW, K, KQQ, JLIM, KFERR
!
      equivalence (INTCHK(1), NXTCHK)

      double precision, parameter :: CM1 = -1.D0

      integer IOPIVA(2)
      save IOPIVA
!
!                      Declarations for error message processing.
!
      character TEXT1(1)*10

      integer, parameter :: MENTXT = 23
      integer, parameter :: MEIDAT = 24
      integer, parameter :: MEMDA1 = 27
      integer, parameter :: MECONT = 50
      integer, parameter :: MERET  = 51
      integer, parameter :: MEEMES = 52
      integer, parameter :: METEXT = 53
      integer, parameter :: MEIVEC = 57
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

     integer, parameter :: LTXTAA =   1
     integer, parameter :: LTXTAB =   7
     integer, parameter :: LTXTAC =  71
     integer, parameter :: LTXTAD = 183
     integer, parameter :: LTXTAE = 226
     integer, parameter :: LTXTAF = 290
     integer, parameter :: LTXTAG = 400
     integer, parameter :: LTXTAH = 427
     integer, parameter :: LTXTAI = 455
     integer, parameter :: LTXTAJ = 498
     integer, parameter :: LTXTAK = 535
     integer, parameter :: LTXTAL = 573
     integer, parameter :: LTXTAM = 620
     integer, parameter :: LTXTAN = 655
     integer, parameter :: LTXTAO =   1

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

      use divaev_module
      use diva_constants
      use divasc_module
      use divamc_module

      integer KORD(*)
      double precision TSPECS(*), Y(*)
      double precision F(*)
      external DIVAF, DIVAO

      double precision, parameter :: CMP75 = -.75D0
      double precision, parameter :: C0 = 0.D0
      double precision, parameter :: C1M5 = 1.D-5
      double precision, parameter :: C1M3 = 1.D-3
      double precision, parameter :: C2P5M3 = 2.5D-3
      double precision, parameter :: C8M3 = 8.D-3
      double precision, parameter :: CP0625 = .0625D0
      double precision, parameter :: CP1 = .1D0
      double precision, parameter :: CP25 = .25D0
      double precision, parameter :: CP3 = .3D0
      double precision, parameter :: CP4 = .4D0
      double precision, parameter :: CP5 = .5D0
      double precision, parameter :: CP875 = .875D0
      double precision, parameter :: C1 = 1.D0
      double precision, parameter :: C1P125 = 1.125D0
      double precision, parameter :: C1P3 = 1.3D0
      double precision, parameter :: C2 = 2.D0
      double precision, parameter :: C6 = 6.D0
      double precision, parameter :: C10 = 10.D0
      double precision, parameter :: C16 = 16.D0
      double precision, parameter :: C4096 = 4096.D0
!
      integer LDIS, I, J, K, J1, J2, L, LX
      double precision TP, TP1, TP2, TP3, DISADJ
      double precision SIGMAS, TPS1, TPS2, EIMINO, EXR
      double precision XP, XP1

      save EIMINO, EXR, SIGMAS, TP, TP1, TP2, TP3, TPS1, TPS2, DISADJ,  &
     &   LDIS
!
!                      Declarations for error message processing.
!
      integer MACT(17), MLOC(8)

      integer, parameter :: MENTXT = 23
      integer, parameter :: MERET  = 51
      integer, parameter :: MEEMES = 52
      integer, parameter :: METEXT = 53
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

     integer, parameter :: LTXTAA =   1
     integer, parameter :: LTXTAB =   8
     integer, parameter :: LTXTAC =  40
     integer, parameter :: LTXTAD =  80
     integer, parameter :: LTXTAE = 129
     integer, parameter :: LTXTAF = 189
     integer, parameter :: LTXTAG = 237
     integer, parameter :: LTXTAH = 272
     integer, parameter :: LTXTAI = 302
     integer, parameter :: LTXTAJ = 342
     integer, parameter :: LTXTAK = 390
     integer, parameter :: LTXTAL = 455
     integer, parameter :: LTXTAM = 484
     integer, parameter :: LTXTAN = 531

     character MTXTAA(3) * (186)
!

    integer, parameter :: LOCM = 32 * 256
!                               KORD1I   Severity   Loc. message
    integer, parameter :: MLOCAC = 23 + 32 * (99 + 256 * LTXTAC)
    integer, parameter :: MLOCAD = 13 + 32 * (38 + 256 * LTXTAD)
    integer, parameter :: MLOCAE = 22 + 32 * (38 + 256 * LTXTAE)
    integer, parameter :: MLOCAF =  3 + 32 * (14 + 256 * LTXTAF)
    integer, parameter :: MLOCAG = 21 + 32 * (38 + 256 * LTXTAG)
    integer, parameter :: MLOCAH =  2 + 32 * (25 + 256 * LTXTAH)
    integer, parameter :: MLOCAI = 11 + 32 * (38 + 256 * LTXTAI)
    integer, parameter :: MLOCAJ = 12 + 32 * (38 + 256 * LTXTAJ)
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
      call DMESS(MACT, MTXTAA, IDAT, FDAT)
      MACT(K) = MENTXT
      go to 2110

    end subroutine divaa
!*************************************************************************

!*************************************************************************
!>
!  This subroutine restores the difference table to its state
!  at the beginning of the current step.  if the integration order
!  was increased, it is reduced. the common array xi is also
!  restored to its state at the beginning of the step. if the
!  stepsize is not being changed, the array v used to compute
!  integration coefficients is restored.
!
!### History
!  * 1987-12-07 DIVABU Krogh   Initial code.

    subroutine divabu(f, kord)

    use diva_constants
    use divasc_module
    use divamc_module

    integer :: kord(*)
    double precision :: f(*)

    integer :: i, l, kqq, j, k
    double precision :: tpd

    double precision, parameter :: c0 = 0.d0
    double precision, parameter :: c2 = 2.d0

! ********
! begin loop to back up difference tables
! ********
    l = ndtf - 1
    do i = 1, nte
         kqq = kord(i + 3)
!++  code for stiff is inactive
!         if (kqq) 2302,2400,2310
!c.           equation is stiff
! 2302    if (linc>=0) go to 2310
!         if (f(l+1+i)) 2306,2308,2304
!c.     order was increased, and thus must be decreased (kqq<0)
! 2304    kqq=kqq+1
!         kord(i+3) = kqq
!         go to 2308
!c.     order was decreased
! 2306    kqq=kqq-1
! 2308    kqq=max(2,-kqq)
!         go to 2350
!++  end
!     equation is not stiff
 2310    if (kqq > 2) then
            if (f(l + kqq) == c0) then
               ! order was increased, and thus must be decreased
               kqq = kqq - 1
               kord(i + 3) = kqq
            end if
         end if
         j = min(kqq, ksc)
         kqmaxi = max(kqmaxi, kqq)
         if (kqq /= 1) f(l + kqq + 1) = 0.d0
         ! back up for backward differences
         do k = 1, j
            f(l + k) = f(l + k) - f(l + k + 1)
         end do
         if (kqq > ksc) then
            ! back up for modified divided differences
            do k = j+1, kqq
               f(l + k) = (f(l+k) - f(l+k+1)) / beta(k)
            end do
         end if
 2400    f(l + kqq + 1) = f(l + kqq + 1) / beta(kqq + 1)
         l = l + numdt
    end do
! end of loop to back up difference tables

! ********
! back up xi to beginning of the step
! ********
      i = ksc + 1
      if (i - iop11 - 1) 2420, 2440, 2450
 2420 tpd = xi(1)
      ! check below needed when starting?
      if (tpd == xi(2)) go to 2450
      do k = i, iop11
         xi(k - 1) = xi(k) - tpd
      end do
 2440 xi(iop11) = c2 * xi(iop11 - 1)
      if (iop11 /= 2) xi(iop11) = xi(iop11) - xi(iop11 - 2)
 2450 kqicon = -1
      icf = ne
      ics = 1
      ldt = 1

    end subroutine divabu
!*************************************************************************

!*************************************************************************
!>
!  THIS SUBROUTINE RETURNS THE FOLLOWING DATA FROM COMMON
!
!  * ID(1) = KEMAX  =  INDEX OF EQUATION WITH LARGEST ERROR ESTIMATE
!  * ID(2) = KSTEP  =  CURRENT STEP NUMBER
!  * ID(3) = NUMDT  =  NUMBER OF DIFFERENCES USED FOR EACH EQUATION
!  * ID(4) =           RESERVED FOR FUTURE USE
!  * ID(5) =           RESERVED FOR FUTURE USE
!  * RD(1) = EMAX   =  MAX. RATIO OF ESTIMATED ERROR TO REQUESTED ERROR
!  * RD(2) =           RESERVED FOR FUTURE USE
!  * RD(3) =           RESERVED FOR FUTURE USE
!
!### History
!  * 1987-12-07 DIVACO Krogh   Initial code.

    subroutine divaco(id, rd)

    use diva_constants
    use divasc_module
    use divamc_module

    integer :: id(5)
    double precision :: rd(3)

    id(1) = kemax
    id(2) = kstep
    id(3) = numdt
    rd(1) = emax

    end subroutine divaco
!*************************************************************************

!*************************************************************************
!>
! This subroutine:
!
!  1. Corrects y for equations which are not stiff
!  2. Estimates errors
!  3. Selects integration orders
!  4. Tests if noise limits the precision
!
!### History
!  * 1988-08-25 DIVACR Krogh   Fix bug in relative error test.
!  * 1988-01-15 DIVACR Krogh   Initial code.

    subroutine DIVACR(Y, F, KORD, TOL, LGROUP)

      use divaev_module
      use diva_constants
      use divasc_module
      use divamc_module

      integer :: lgroup(*) !! vector indicating how error tolerances are to be grouped
                           !! (and possibly how integration orders are to be grouped).
      integer :: kord(*)
      double precision :: y(*) !! vector of predicted values on entry, and of corrected
                               !! values when the return is made.
      double precision :: tol(*) !! vector containing error tolerances (and possibly relative
                                 !! error factors).
      double precision :: f(*)

!    KD = VECTOR GIVING ORDERS OF THE DIFFERENTIAL EQUATIONS
!         (IF EQUATIONS HAVE DIFFERENT ORDERS).
!    KQ = VECTOR OF INTEGRATION ORDERS.

      integer L, I, KQL, KQN, KQD, JLGREP, J, K, ILGROR, ITOLOR, JLGROR,&
              IORD, KOUTKO, KQLORD, LL, LKQMAX

      double precision, parameter :: CM8 = -8.D0
      double precision, parameter :: CM2 = -2.D0
      double precision, parameter :: CMP5 = -.5D0
      double precision, parameter :: C0 = 0.D0
      double precision, parameter :: CQ3125 = .03125D0
      double precision, parameter :: CP1 = .1D0
      double precision, parameter :: CP125 = .125D0
      double precision, parameter :: CP25 = .25D0
      double precision, parameter :: CP5 = .5D0
      double precision, parameter :: CP75 = .75D0
      double precision, parameter :: CP8 = .8D0
      double precision, parameter :: CP9375 = .9375D0
      double precision, parameter :: C1 = 1.D0
      double precision, parameter :: C1P4 = 1.4D0
      double precision, parameter :: C2 = 2.D0
      double precision, parameter :: C4 = 4.D0
      double precision, parameter :: C10 = 10.D0
      double precision, parameter :: C20 = 20.D0
      double precision, parameter :: C40 = 40.D0
      double precision, parameter :: C1000 = 1000.D0

      double precision TPP, E, EI, EPS, ERCOEF, RND, RNOISE, S !HH
      double precision TP2, TPS1, TPS2, TPS3, TPS4, TPS5, TPS6, TPS7
      double precision REF(4)
      double precision EIBND(KDIM-1)
!++  Code for INTEGO is active
      double precision TEMPA(4), TEMPAO(4)
!++  End
      save KOUTKO, LKQMAX
      equivalence (TPS1,TEMPA(1)), (TPS2,TEMPA(2)), (TPS3,TEMPA(3)),    &
     &   (TPS4, TEMPA(4))
    !  equivalence (G(1, 1), HH)
      integer MACT1(2), MACT2(12)

      ! Parameters for Interface to MESS and DMESS
      integer, parameter :: MERET  =51
      integer, parameter :: METEXT =53
      integer, parameter :: METABL =55

! ********* Error message text ***************
![Last 2 letters of Param. name]  [Text generating message.]
!AA KSTEP=$(I6) T=$(E15.8) H=$(E12.5) LSC=$(I3) $C
!   EIMIN=$(E8.2) EAVE=$G KSC=$(I2) SIGMA($J)=$G $C
!   RQ=$(E11.5)$G$E
!   $
!AB I$HKQ$HLI$HE$HEI$HEPS$HF$H$H$H$H
!   HIGH ORDER PREDICTED DIFFERENCES$HRNOISE$HSTIFF$HBETA$E
      integer, parameter :: LTXTAA = 1
      integer, parameter :: LTXTAB = 1

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

      use divaev_module
      use diva_constants
      use divasc_module
      use divamc_module

!                 K - 1 + 1 / K  is equivalent to max(1, K-1)
      double precision GG(MAXORD - 1 + 1/MAXORD), B(KDIM+MAXORD), W(KDIM+MAXORD)
      integer  K, N, J

      double precision, parameter :: C0     = 0.D0
      double precision, parameter :: CP1    = .1D0
      double precision, parameter :: CRBQI  = .421875D0
      double precision, parameter :: CP5    = .5D0
      double precision, parameter :: CP5625 = .5625D0
      double precision, parameter :: C1     = 1.D0
      double precision, parameter :: C1P125 = 1.125D0

!++  Code for STIFF is inactive
!      INTEGER          GODIF
!++  End
      double precision TP1, TP2, TEMP, TP !HH
     ! equivalence (G(1, 1), HH)
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

    use diva_constants
    use divasc_module

      integer KORD(*)
      double precision T(*), Y(*)
      double precision F(*)

      integer I, ICI, IDT, INTERP, INTEG, INTEGZ, IY, IYI, IYN, IYNI, J,&
     &    K, KQMXI, KQMXS, KQQ, L, N

      double precision, parameter :: C0 = 0.D0
      double precision, parameter :: C1 = 1.D0
      double precision, parameter :: C2 = 2.D0

      double precision C(KDIM+MAXORD-1), ETA(KDIM)
      double precision GAMMA(KDIM)
      double precision TP1, HI
      double precision CSUM(KDIM+MAXORD-1)
      double precision XP1
      logical LNOTM1
!
!              Stuff for processing error messages
      integer IDAT(1)
      double precision FDAT(6)

      integer, parameter :: MENTXT =23
      integer, parameter :: MERET  =51
      integer, parameter :: MEEMES =52
      integer, parameter :: METEXT =53

      integer MACT(8)

! ********* Error message text ***************
![Last 2 letters of Param. name]  [Text generating message.]
!AA DIVAIN$B
!AB Interpolating at T(1)=$F with $B
!AC TN=$F, T(2)=$F and H=$F.  T(1) must be in [$F, $F].$E
!AD internal variable LDT = $I.  Interpolation not allowed now.$E

      integer, parameter :: LTXTAA =  1
      integer, parameter :: LTXTAB =  9
      integer, parameter :: LTXTAC = 41
      integer, parameter :: LTXTAD = 94

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
      call DMESS(MACT, MTXTAA, IDAT, FDAT)
      if (MACT(2) < 50) go to 3820
      return
    end subroutine divain
!*************************************************************************

!*************************************************************************
!>
!  SUBROUTINE TO SET UP OPTIONS FOR DIFFERENTIAL EQUATION PACKAGE DIVA
!
!### History
!  * 1987-12-07 DIVAOP Krogh   Initial code.

    subroutine DIVAOP(IOPT, FOPT)

      use divaev_module
      use diva_constants
      use divasc_module
      use divamc_module

      double precision FOPT(*)
      integer IOPT(*)

      integer IOPTS(23), INCOP(22), I, IA, J, K, LIOPT, MULTJ

      double precision, parameter :: CMP75  = -.75D0
      double precision, parameter :: C0     = 0.D0
      double precision, parameter :: CP25   = .25D0
      double precision, parameter :: CP3    = .3D0
      double precision, parameter :: CP5    = .5D0
      double precision, parameter :: CP625  = .625D0
      double precision, parameter :: CP75   = .75D0
      double precision, parameter :: CP875  = .875D0
      double precision, parameter :: CP9    = .9D0
      double precision, parameter :: C1     = 1.D0
      double precision, parameter :: C1P125 = 1.125D0
      double precision, parameter :: C2     = 2.D0
      double precision, parameter :: C4     = 4.D0
      double precision, parameter :: C10    = 10.D0
      double precision, parameter :: C16    = 16.D0

 !     equivalence (IOPTC(3), IOP3)
      save IOPTS, LIOPT
!
!                      Declarations for error message processing.
!
      integer, parameter :: MECONT = 50
      integer, parameter :: MERET  = 51
      integer, parameter :: MEEMES = 52
      integer, parameter :: MEIVEC = 57

      integer MACT(7), MACT1(5)
!
! ********* Error message text ***************
![Last 2 letters of Param. name]  [Text generating message.]
!AA DIVAOP$B
!AB Error in IOPT() specifications: IOPT =$E
!AC HMIN = $F is > HMAX = $F.$E

      integer, parameter :: LTXTAA =  1
      integer, parameter :: LTXTAB =  9
      integer, parameter :: LTXTAC = 49

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

      use divaev_module
      use diva_constants
      use divasc_module
      use divamc_module

      integer KORD(*)
      double precision Y(*), YN(*)
      double precision F(*)

      double precision, parameter :: C0 = 0.D0
!
      integer I, INTEG, INTEGS, J, K, KQQ, L, N
      double precision TEMP(KDIM)
      double precision TP1
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
!            XP = XP + dble(G(J, INTEG)) * dble(TEMP(J))
!++  END
 4650    continue
         K = INTEG + INTEGS
         do 4660 J = K, 1, -1
!++  Code for ~{p,x} is active
            XP = XP + G(1, J) * YN(IY + J)
!++  Code for {p,x} is inactive
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

      use divaev_module
      use diva_constants
      use divasc_module
      use divamc_module

      integer LPRINT, KORD(*)
      character TEXT*(*)
      character TEXT1(1)*11, TEXT2(1)*4, TEXT3(1)*5, TEXT4(1)*4
      integer J, K, L, N1, N2
      double precision TSPECS(*), Y(*)
      double precision F(*)

    ! Declarations for error message processing.
      integer, parameter :: MEDDIG = 12
      integer, parameter :: NEDDIG = -MEDDIG
      integer, parameter :: METDIG = 22
      integer, parameter :: METABS = 32
      integer, parameter :: MERET  = 51
      integer, parameter :: METEXT = 53
      integer, parameter :: METABL = 55
      integer, parameter :: MEIVEC = 57
      integer, parameter :: MEFVEC = 61
      integer, parameter :: MEFMAT = 62

      integer MACT0(3), MACT1(2), MACT2(7), MACT3(7), MACT4(8), &
     &   MACT5(11), MACT6(3), MACT7(14), MACTFV(4)

      integer, parameter :: KPI = 400301  !! wddtrr
      integer, parameter :: KPE = 1305501 !! wwddtrr
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

     integer, parameter :: LTXTAA =   1
     integer, parameter :: LTXTAB =  14
     integer, parameter :: LTXTAC =   1
     integer, parameter :: LTXTAD =   1
     integer, parameter :: LTXTAE =   1
     integer, parameter :: LTXTAF =  21
     integer, parameter :: LTXTAG =   1
     integer, parameter :: LTXTAH =  16
     integer, parameter :: LTXTAI =  22
     integer, parameter :: LTXTAJ =   1
     integer, parameter :: LTXTAK =   1
     integer, parameter :: LTXTAL =   1
     integer, parameter :: LTXTAM = 655
     integer, parameter :: LTXTAN =   1
     integer, parameter :: LTXTAO =   1
     integer, parameter :: LTXTAP =   1
     integer, parameter :: LTXTAQ =   1

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
      call DMESS(MACTFV, TEXT1, KORD, TSPECS)
      MACTFV(3) = NY
      call DMESS(MACTFV, TEXT2, KORD, Y)
      call DMESS(MACTFV, TEXT3, KORD, Y(NYNY))
      MACTFV(3) = NTE
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
      call DMESS(MACT4, MTXTAE, KORD, F(NDTF))
!
   80 if (N2 <= 1) return
! ********
! WRITE SCALARS IN COMMON
! ********
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
      call DMESS(MACT5, MTXTAH, IVC2, RVC2)
      call DMESS(MACT1, MTXTAI, IVC2, DVC1)
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
         call DMESS(MACT7, MTXTAK, IDAT, FDAT)
  100 continue
!++  Code for STIFF is inactive
!     if (MAXDIF <= 0) return
!        Need to define MACT8 and set values
!     call DMESS(MACT8, 'D$B', IDAT, D)
!     call DMESS(MACT8, 'DS$B', IDAT, DS)
!++  End
!
      call DMESS(MACT1, MTXTAL, IDAT, EVC)
      return
    end subroutine divadb
!*************************************************************************

!*************************************************************************
!>
!  SUBROUTINE TO LOCATE OUTPUT POINTS AT ZEROS OF ARBITRARY
!  FUNCTIONS  **** GSTOPS **** FOR DIFFERENTIAL EQUATION
!  INTEGRATOR DODE (OR DIVA).
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

      use divaev_module
      use diva_constants
      use divasc_module
      use divamc_module

      integer KORD(*), IFLAG, NSTOP
      double precision TSPECS(*), Y(*), TOLD, TSAVE
      double precision F(*), GNEW(*), GT(*), GOLD
      save GOLD, TOLD, TSAVE

      integer I, IG

      ! Declarations for error message processing.
      integer, parameter :: MEMDA1 = 27
      integer, parameter :: MERET  = 51
      integer, parameter :: MEEMES = 52

      integer MACT(7)
! ********* Error message text ***************
![Last 2 letters of Param. name]  [Text generating message.]
!AA DIVAG$B
!AB Call with bad values of KORD.  KORD(1)=$I, KORD(2)=$I, when $C
!   TSPECS(1)=$F and KSTEP=$M.$E

      integer, parameter :: LTXTAA = 1
      integer, parameter :: LTXTAB = 8

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
      call DMESS(MACT, MTXTAA, KORD, TSPECS)
      IFLAG = 8
      KEXIT = 6
      go to 470
    end subroutine divag
!*************************************************************************

!*************************************************************************
!>
! Subroutine to find a bounded zero
!
!### Usage
!  Usage is as follows (of course, variations are possible.)
!  Somehow one has available X1, F1, X2, and F2 such
!  that F1 = F(X1), F2 = F(X2) and F1*F2 <= 0.
!  In addition, one should assign a value to TOL.
!
!```fortran
!     MODE = 0
! XXX call DZERO(X1,F1,X2,F2,MODE,TOL)
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
!```
!
!### History
!
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
!  * Algorithmic changes, vars. added to save stmt., Sept. 1987 by Krogh
!  * Modified for portability, April 1984 by Krogh.
!  * analysis and coding by Fred T.Krogh at the Jet Propulsion
!    Laboratory, Pasadena, Calif.  April 25, 1972.

    subroutine DZERO(X1, F1, X2, F2, MODE, TOL)

      integer :: MODE !! is a parameter used for communication between this
                      !! subroutine and the user. (The user should set MODE
                      !! only to initialize it to 0 before the first call)
                      !!
                      !! * `=1` compute F(X1) and call $ZERO
                      !! * `=2` F(X1) is approximately 0, the iteration is finished
                      !!   and the error criterion is satisfied.
                      !! * `=3` same as MODE=2, except the error criterion can
                      !!   not be satisfied.
                      !! * `=4` apparently the function has a discontinuity
                      !!   between X1 and X2 -- No zero can be found
                      !! * `=5` F1*F2 was greater than zero on the first call, and an attempt
                      !!   to bound the zero on both sides have failed.
                      !! * `=6` fatal error -- $ZERO was called after mode was set >=2.
                      !!   If $ZERO is called again, the program will be stopped.
                      !!   (Unless MODE is set to 0)
                      !! * `<0` If MODE is set <0 and $ZERO is called, no action is taken
                      !!   except that print is turned on for -MODE calls to $ZERO.
                      !!   This print gives all values of X and F used in the iteration.
      double precision :: X1 !! independent variable
      double precision :: X2 !! second value of independent variable
      double precision :: F1 !! dependent variable --  initially   F1=F(X1).
                             !! When MODE=1 (or 5) the user is to compute F(X1) given X1
      double precision :: F2 !! F(X2) on the initial entry.  When MODE = 2-4, F2=F(X2) and
                             !! F1*F2 <= 0.
      double precision :: TOL !! is the error tolerance:
                              !!
                              !! * `TOL>0` Iterate until values of X1 and X2 are known
                              !!   for which abs(X1-X2) <= tol and F1*F2 <= 0.
                              !! * `TOL<0` Iterate until a value of X1 is found for which
                              !!   abs(F1) <= abs(TOL).
                              !! * `TOL=0` Iterate until the zero is determined as
                              !!   precisely as possible.  MODE = 3 is impossible
                              !!   in this case.

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

      integer KS, KTYP, J, I

      integer, parameter :: LINIT = -40

      integer KNKP(2), LCHG(30), LMODE, LNLP(2), NP
      double precision XX, XO, XL, FF, FO, FL, FFDFO
      double precision DIV, QFM, QXM, TP, TP1, XXMXO, XXMXOL
      double precision RND, XRND, SMALL, TOLX
      double precision XXMXL, XLMXB, FFMFL, FFMFB, FLMFB
      double precision DXDFFF, DXDFFO, DFDXXX, DFDXXO

      double precision, parameter :: C0 = 0.D0
      double precision, parameter :: C1 = 1.D0
      double precision, parameter :: C2 = 2.D0
      double precision, parameter :: C4 = 4.D0
      double precision, parameter :: C8 = 8.D0
      double precision, parameter :: CP125 = 0.125D0
      double precision, parameter :: CP25 = 0.25D0
      double precision, parameter :: CP75 = 0.75D0
      double precision, parameter :: CP5 = 0.5D0
      double precision, parameter :: C1P25 = 1.25D0
      double precision, parameter :: CP01 = 0.01D0
      double precision, parameter :: CP001 = 0.001D0
      double precision, parameter :: CP99 = 0.99D0
      double precision, parameter :: C1P031 = 1.03125D0
!
!                      Declarations for error message processing.
!

      double precision FDAT(4)
      integer MACT(5), MACT1(2), MLOC(4), IDAT(2)
      save DIV, FL, FLMFB, FO, KNKP, KS, KTYP, LCHG, LMODE,             &
     &   LNLP, MACT, NP, RND, SMALL, XL, XLMXB, XO, XX, XXMXOL

      integer, parameter :: MERET  = 51
      integer, parameter :: MEEMES = 52
      integer, parameter :: METEXT = 53

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

      integer, parameter :: LTXTAA =   1
      integer, parameter :: LTXTAB =   8
      integer, parameter :: LTXTAC =  63
      integer, parameter :: LTXTAD = 112
      integer, parameter :: LTXTAE = 169
      integer, parameter :: LTXTAF =   1
      integer, parameter :: LTXTAG =   1

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
      integer, parameter :: MENTXT = 23
      integer, parameter :: MEIDAT = 24
      integer, parameter :: MECONT = 50
      integer, parameter :: MERET  = 51
      integer, parameter :: MEEMES = 52
      integer, parameter :: METEXT = 53
      integer, parameter :: MEIMAT = 58

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

      integer, parameter :: LTXTAA =   1
      integer, parameter :: LTXTAB =   9
      integer, parameter :: LTXTAC = 233
      integer, parameter :: LTXTAD = 267
      integer, parameter :: LTXTAE = 295

      character MTXTAA(2) * (156)
!                          Next 4 lines not automatically generated
      integer, parameter :: LTXTEE = LTXTAE - 156 - 2
      integer, parameter :: LTXEND = 156
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