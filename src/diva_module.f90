subroutine diva(tspecs, y, f, kord, neq, divaf, divao, idimt, &
   idimy, idimf, idimk, iopt)
! Copyright (c) 1996 California Institute of Technology, Pasadena, CA.
! ALL RIGHTS RESERVED.
! Based on Government Sponsored Research NAS7-03001.
!>> 2015-03-15 DIVA  Krogh  Removed extra call divabu after noise test
!>> 2015-03-15 DIVA  Krogh  Forced restart needs more reduction in h.
!>> 2010-02-20 DIVA  Krogh  Fixed calling DIVAOP with array other than F.
!>> 2009-11-03 DIVA  Krogh  Added option 11, more variables initialized.
!>> 2009-10-30 DIVA  Krogh  Gave KSSTRT and ROBND initial values.
!>> 2009-10-30 DIVA  Krogh  Fixed reference to undefined location in F.
!>> 2009-10-21 DIVA  Krogh  Got rid of NaN in diag. print when LSC=3.
!>> 2009-10-15 DIVA  Krogh  A few changes on how noise is handled.
!>> 2002-11-12 DIVA  Krogh  Fixed problem integrating to final output pt
!>> 2002-08-29 DIVA  Krogh  Added test for invalid HMIN/HMAX.
!>> 2002-07-26 DIVA  Krogh  Added KOUTKO to fully support Option 10.
!>> 2002-05-14 DIVA  Krogh  Fix starting prob. for Option 18.
!>> 2002-05-13 DIVA  Krogh  Put exponent letter in  numbers missing them
!>> 2002-05-12 DIVA  Krogh  Added error message for bad option 5 usage.
!>> 2001-09-07 DIVA  Krogh  Changes to allow user tol on G-Stops.
!>> 2001-05-25 DIVA  Krogh  Minor change for making .f90 version.
!>> 2001-05-18 DIVA  Krogh  Less computing with no error test
!>> 2001-05-17 DIVA  Krogh  Fixed so with no error test can't start dump
!>> 2001-04-24 DIVA  Krogh  Inserted comments from ivacom.
!>> 2000-12-01 DIVA  Krogh  Removed (some of) unused C1, MAXSTF, METEXT.
!>> 1999-12-28 DIVA  Krogh  Saved S in DIVACR for output consistency.
!>> 1999-08-19 DIVA  Krogh  Removed superfluous test above label 3520.
!>> 1997-04-22 DIVA  Krogh  Got rid of assigned go to's. F=0 if diag.
!>> 1996-08-26 DIVA  Krogh  Initialize F to 0 if dumping solution.
!>> 1996-08-23 DIVA  Krogh  Print TN not TSPECS(1) in error messages.
!>> 1996-05-30 DIVA  Krogh  Changed DERIVS/OUTPUT to  DIVAF/DIVAO.
!>> 1996-04-27 DIVA  Krogh  Changes to use .C. and C%%.
!>> 1996-03-30 DIVA  Krogh  Added external statement.
!>> 1996-03-25 DIVA  Krogh  Introduced TEXT1 to comply with F77.
!>> 1996-02-27 DIVA  Krogh  Fixed so DUMP not affected by ignored eqs.
!>> 1995-12-18 DIVA  Krogh  Fixed so no solution dump on 0 length integ.
!>> 1995-11-09 DIVA  Krogh  Fixed so char. data at col. 72 is not ' '.
!>> 1995-06-19 DIVA  Krogh  Fixed prob. with discon. just after restart.
!>> 1995-05-09 DIVA  Krogh  Fixed G-Stop/discontinuity code interaction
!>> 1995-04-26 DIVA  Krogh  Use KQMAXS instead of KQMAXI when LDIS>1000.
!>> 1995-04-26 DIVA  Krogh  Keep current KQL on discontinutiy.
!>> 1994-12-16 DIVA  Krogh  Fixed option 12 with K12 < 0.
!>> 1994-11-11 DIVA  Krogh  Declared all vars.
!>> 1994-11-02 DIVA  Krogh  Changes to use M77CON
!>> 1994-09-08 DIVA  Krogh  Added CHGTYP code.
!>> 1994-07-11 DIVA  Krogh  Fix to get same state with/without var. eqs.
!>> 1994-03-07 DIVA  Krogh  Allow larger order in single precision.
!>> 1994-01-14 DIVA  Krogh  Minor change to allow changing TFINAL.
!>> 1993-04-27 DIVA  Krogh  Additions for Conversion to C.
!>> 1993-04-12 DIVA  Krogh  Converted to use slightly altered MESS.
!>> 1993-04-12 DIVA  Krogh  Fixed LSC so sol. saved when HMAX is small.
!>> 1992-10-13 DIVA  Krogh  Fixed G-Stop/discontinuity code interaction.
!>> 1992-09-21 DIVA  Krogh  Fixed bug in discontinuity code.
!>> 1992-09-09 DIVA  Krogh  Fixed bug - Var. Eqs. with discontinuities.
!>> 1992-08-07 DIVA  Krogh  Storage map printed only if option 10 .ne. 0
!>> 1992-07-16 DIVA  Krogh  Restored correct discontinuity code.
!>> 1992-06-16 DIVA  Krogh  Eliminate reuse of storage for option 12.
!>> 1992-04-08 DIVA  Krogh  Removed unused labels, 1020, 2120.
!>> 1992-03-30 DIVA  Krogh  Fixed bug in DIVAOP error message.
!>> 1992-03-12 DIVA  Krogh  Simplified DIVABU, more digits in B's.
!>> 1992-01-16 DIVA  Krogh  Fixed minor bug in error messages.
!>> 1991-12-03 DIVA  Krogh  Major change for improved error checks.
!>> 1991-06-17 DIVA  Krogh  Fixed bug in checking storage allocation.
!>> 1991-04-11 DIVA  Krogh  Fixed minor bug re. option 12 in DIVAOP.
!>> 1991-03-28 DIVA  Krogh  Removed check at label 650 for KORD2I<0.
!>> 1991-02-08 DIVA  Krogh  Changed some floats to generics
!>> 1990-11-08 DIVA  Krogh  Fixed bug on TSPECS on discon.
!>> 1990-09-14 DIVA  Krogh  Fixed bug when discon. and sol. save.
!>> 1990-09-13 DIVA  Krogh  Increased dimension of BETA by 1.
!>> 1990-09-13 DIVA  Krogh  Added one more poss. on rel. error test.
!>> 1990-09-11 DIVA  Krogh  Recent change messed up getting dump output.
!>> 1990-06-05 DIVA  Krogh  Fixed bug in noise test, comments in IVACOM.
!>> 1990-05-08 DIVA  Krogh  Fixed new bug when TMARK hit in DIVAG.
!>> 1990-04-17 DIVA  Krogh  Fixed minor problem in DIVAIN error msg.
!>> 1990-04-10 DIVA  Krogh  Fixed interaction between discon. & dump.
!>> 1990-03-23 DIVA  Krogh  Fixed bug on option "-2", see 1989-12-07.
!>> 1990-03-20 DIVA  Krogh  Fixed rarely occuring loop.
!>> 1990-01-29 DIVA  Krogh  Removed unneeded labels.
!>> 1989-12-14 DIVA  Krogh  Saved common block DIVAEV.
!>> 1989-12-07 DIVA  Krogh  Added option "2" to DIVAOP.
!>> 1989-11-09 DIVA  Krogh  Made GG a save var. in DIVAHC
!>> 1989-08-21 DIVA  Krogh  Fix out of bounds ref. to V in DIVABU
!>> 1989-07-26 DIVA  Krogh  Fix bug in initial dim. check
!>> 1989-07-21 DIVA  Krogh  Code for integrating discontinuities
!>> 1987-12-07 DIVA  Krogh  Initial code.
!
!--D replaces "?": ?IVA,?IVAA,?IVABU,?IVACO,?IVACR,?IVAEV,?IVAF,?IVAHC,
!-- & ?IVAG,?IVAIN,?IVAMC,?IVAO,?IVAOP,?IVAPR,?IVASC,?IVACE,?IVAIE,
!-- & ?IVAPE,?MESS
!
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
!   reduced if EIMIN .le. EIMINO then the reduction factor is set to
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
!   .ge. EXR then the step size is reduced.  Could be a parameter.
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
!   of LDT, except L=1 if LDT=-1, and MAXINT .ge. 0.  (Used in *IVAA,BU
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
! LINCD  (*IVAMC) Value of smallest k for which HINCC**k .ge. 2.
!   (=-2 if user is specifying all step size changes.)
! LINCQ  (*IVAMC) Value of smallest k for which HINCC**k .ge. 4.
! LIOPT  (*IVAOP) Value of the last index in IOPT on the last call.
!   Used so *IVA can print IOPT in error messages.
! LL     (*IVACR) Temporary variable used when equations are grouped
!   for integration order control.
! LNOTM1 (*IVAIN) Logical variable = L .ne. -1.  If LNOTM1 is true,
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
!   to, and T(2) is used in a check that |HI| .le. |T(2)|.  When used by
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
!
!++S Default KDIM = 16
!++  Default KDIM = 20
!++  Default MAXORD = 2, MAXSTF = 1
!++  Default INTEGO, VAREQ, OUTPUT, DUMP, GSTOP, EXTRAP
!++  Default STIFF=.F., ARGM=.F., ERRSTO=.F.
!
integer neq, idimt, idimy, idimf, idimk
integer kord(*), iopt(*)
!--D Next line special: P=>D, X=>Q
double precision tspecs(*), y(*)
double precision f(*)
external divaf, divao
!
! *********************** Internal Variables ***************************
!
! Comments for variables used in this package can be found in the file
!   IVACOM.
!
! *********************** Type Declarations ****************************
!
integer kdim, maxord, maxstf
!++ Substitute for KDIM, MAXORD, MAXSTF below
parameter (kdim = 20, maxord = 2, maxstf = 1)
!--D Next line special: P=>D, X=>Q
double precision tn
double precision xi(kdim)
!
!--D Next line special: P=>D, X=>Q
double precision tg(2), tgstop(2), tmark, tmarkx, tout, tolg
double precision alpha(kdim), beta(kdim+1)
double precision  d(maxstf+maxord,maxord), g(kdim,maxord)
double precision v(kdim+maxord)
double precision hc, hdec, hinc, hincc, hmax, hmaxp9, hmin
double precision fdat(11)
!
double precision ds(maxstf+maxord, maxord), gs(kdim)
double precision sigma(kdim), rbq(kdim), dnoise
double precision eave, eimax, eimin, emax, erep, robnd, snoise
!
!.    SPECIFICATION OF ENVIRONMENTAL CONSTANTS.
double precision eeps10, eeps16, erov10, eeps2
double precision eept75, eovep2, ovtm75, ovd10
common / divaev / eeps2, eept75, eovep2, ovtm75, ovd10, eeps10, &
   eeps16, erov10
save / divaev /
integer iopst, kordi, kqmaxd, kqmaxi, ldt, maxdif, maxint, nkdko, &
   nte, nyny, ndtf, numdt
common / divasc / tn, xi, iopst, kordi, kqmaxd, kqmaxi, ldt, &
   maxdif, maxint, nkdko, nte, nyny, ndtf, numdt
!
integer icf,ics,igflg,igtype(2),igstop(2),ilgrep,ings,iop3,iop4, &
   iop5,iop6,iop7,iop8,iop9,iop10,iop11,iop12,iop13,iop14,iop15, &
   iop16,iop17,iop18,iop19,iop20,iop21,iop22,iop21s,itolep,iy, &
   kemax,kis,kmark,kord1i,kord2i,kpred,kqdcon,kqicon,kqmaxs, &
   kqmxds,kqmxil,kqmxip,kqmxis,ksc,ksout,ksstrt,kstep,lex,linc, &
   lincd,lincq,lsc,maxkqd,maxkqi,method,ne,neptol,ng,ngtot, &
   noiseq,noutko,ntolf,ny,idat(6)
common /divamc/ tg,tgstop,tmark,tmarkx,tout,tolg,hc,hdec,hinc, &
   hincc,hmax,hmaxp9,hmin,alpha,beta,d,g,v,ds,gs,sigma,rbq,dnoise, &
   eave,eimax,eimin,emax,erep,robnd,snoise,fdat,icf,ics,igflg, &
   igtype,igstop,ilgrep,ings,iop3,iop4,iop5,iop6,iop7,iop8,iop9, &
   iop10,iop11,iop12,iop13,iop14,iop15,iop16,iop17,iop18,iop19, &
   iop20,iop21,iop22,iop21s,itolep,iy,kemax,kis,kmark,kord1i, &
   kord2i,kpred,kqdcon,kqicon,kqmaxs,kqmxds,kqmxil,kqmxip,kqmxis, &
   ksc,ksout,ksstrt,kstep,lex,linc,lincd,lincq,lsc,maxkqd,maxkqi, &
   method,ne,neptol,ng,ngtot,noiseq,noutko,ntolf,ny,idat
save / divamc / , / divasc /
!
integer kgo, intchk(0:30), nxtchk
integer ihi, jl, j, ilow, k, kqq, jlim, kferr
!
equivalence (intchk(1), nxtchk)
double precision cm1
parameter (cm1 = (-1.d0))
integer iopiva(2)
save iopiva
!
!                      Declarations for error message processing.
!
character text1(1)*10
integer  mentxt,meidat,memda1,mecont,meret,meemes,metext,meivec
parameter (mentxt =23)
parameter (meidat =24)
parameter (memda1 =27)
parameter (mecont =50)
parameter (meret  =51)
parameter (meemes =52)
parameter (metext =53)
parameter (meivec =57)
!
integer mact(16), mloc(12), mact1(4)
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
!AN Option 5 argument must be .le. 0 or .gt. 4.$E
!   $
!AO KORD values for this option starting at KORD($M) are:$E
integer ltxtaa,ltxtab,ltxtac,ltxtad,ltxtae,ltxtaf,ltxtag,ltxtah, &
 ltxtai,ltxtaj,ltxtak,ltxtal,ltxtam,ltxtan,ltxtao
parameter (ltxtaa=  1,ltxtab=  7,ltxtac= 71,ltxtad=183,ltxtae=226, &
 ltxtaf=290,ltxtag=400,ltxtah=427,ltxtai=455,ltxtaj=498, &
 ltxtak=535,ltxtal=573,ltxtam=620,ltxtan=655,ltxtao=  1)
character mtxtaa(3) * (233)
character mtxtab(1) * (55)
data mtxtaa/'DIVA$BThe interval [1, 10**6], bounds the allowed val &
ues for nte=$i.$EFor option $i, the interval [$i, $i], bounds the$ &
 allowed values for the integration order which is set to $i.$EOpt &
ion 16 must be used for error control.$ef($i) = ','$f, but it must &
 be -1.0 when skipping the error check.$EFor option $i, the interv &
al [$i, $i] bounds the allowed values for kord($i)=$i, which is us &
ed to specify an $Boutput type for printing.$Eoutput group for pri &
nting.$Eequation gro','up for variational equations.$Eorder for a$ &
 differential equation.$Eequation group for diagnostic print.$Eequ &
ation group for integration order control.$Eequation group for err &
or control.$EOption 5 argument must be .le. 0 or .gt. 4.$e'/
data mtxtab/'KORD values for this option starting at KORD($M) are: &
$e'/

! End of automatically generated error message code.
!
!        for KGO =     1      2      3      4      5      6      7
data mloc / ltxtai,ltxtak,ltxtal,ltxtam,ltxtaj,ltxtag,ltxtah, &
            ltxtab,ltxtac,ltxtad,ltxtae, ltxtan /
!           KGO        8      9     10     11      12
!
!                      1  2  3 4       5  6       7       8       9 10
data mact / meemes,38,24,0, mentxt, 0, metext, mecont, memda1,0, &
  metext, meidat,0, meivec,0, mecont /
!           11      12 13     14 15     16
data mact1 / metext, meivec, 0, meret /
data iopiva(1) / 1111 /
data text1 / 'IOPT()= $B' /
!
! ************** START OF EXECUTABLE CODE ******************
!
!     **** TEST IF CONTINUING AN INTEGRATION
if (kord(1) .ne. 0) go to 330
!     **** INITIALIZE VARIOUS SCALARS
kstep = 0
kqmxis = 0
kord2i = -5
kord(2) = -1
nte = neq
ne = nte
tolg = 0.d0
!     **** SET UP OPTIONS
if (iopt(1) .ne. 0) call divaop(iopt, f)
call divaop(iopiva, f)
if (iopt(1) .eq. 0) iopiva(2) = 1
!
if ((ne .le. 0) .or. (ne .gt. 1000000)) then
   idat(1) = ne
   kgo = 8
   go to 650
end if
!                         Set up diagnostic print on storage allocation.
intchk(0) = 245
if (iop10 .ne. 0) intchk(0) = 247
!
!     **** CHECK TSPECS STORAGE ALLOCATION
intchk(2) = idimt
intchk(3) = 4
nxtchk = 4
if (iop5 .ne. 0) then
   intchk(4) = 5
   intchk(5) = 5
   if (iop5 .gt. 0) then
      intchk(6) = iop5 - 4
      if (iop5 .lt. 5) then
         kgo = 12
         go to 600
      end if
   else
      ihi = -iop5
      jl = 4
      do 15 ihi = ihi, idimk-3, 3
         j = abs(kord(ihi))
         if (j .eq. 0) go to 20
         if (abs(kord(ihi + 2)) .gt. 1) then
            idat(2) = -1
            idat(3) = 1
            kgo = 6
            go to 600
         end if
         if ((j .le. jl) .or. (j .gt. kord(ihi+1))) then
            kgo = 7
            idat(2) = jl + 1
            idat(3) = kord(ihi+1)
            go to 610
         end if
         jl = kord(ihi+1)
15       continue
      if (kord(ihi) .ne. 0) ihi = ihi + 3
20       intchk(6) = jl - 4
   end if
   nxtchk = 7
end if
25 call optchk(intchk, iopt, 'DIVA / TSPECS$E')
if (nxtchk .lt. 0) kord2i = -4
!
!     **** CHECK KORD STORAGE ALLOCATION
intchk(2) = idimk
intchk(3) = ne + 3
nxtchk = 4
if (iop5 .lt. 0) then
   intchk(4) = 5
   intchk(5) = -iop5
   intchk(6) = ihi + iop5
   nxtchk = 7
end if
!
!++  Code for VAREQ is active
if (iop18 .ne. 0) then
   ne = abs(kord(iop18))
   intchk(nxtchk) = 18
   ilow = iop18
   kgo = 1
!.       **** CHECK OPTION FOR VALID INPUT
   go to 430
end if
!++  End
30 continue
if (nkdko .ne. 0) then
!                        **** STORAGE ALLOCATED FOR ODE ORDERS
   intchk(nxtchk) = 17
   intchk(nxtchk+1) = nkdko
   intchk(nxtchk+2) = nte
   nxtchk = nxtchk + 3
end if
!++  Code for STIFF is inactive
!      IF (IOPST .ne. 0) then
!         INTCHK(NXTCHK) = 17
!         INTCHK(NXTCHK+1) = IOPST
!         INTCHK(NXTCHK+2) = NTE
!         NXTCHK = NXTCHK + 3
!      end if
!++  End
!
! **** SET INITIAL INTEGRATION ORDERS, TEST ODE ORDERS ****
!
maxint = 0
maxdif = 0
ny = 0
do 80 k = 1, nte
   if (nkdko .ne. 0) kordi = kord(nkdko + k - 1)
   ny = ny + abs(kordi)
!++  Code for STIFF is inactive
!      IF (IOPST .EQ. 0) GO TO 60
!c.    **** CHECK FOR POTENTIAL STIFF EQUATION
!      JS = abs(KORD(IOPST+K-1)) - 1
!      IF ( JS ) 52,60,54
!c.    **** EQUATION IS NOT ACTIVE
!   52 KQQ = 0
!      GO TO 56
!c.    **** EQUATION USES IMPLICIT METHOD
!   54 KQQ = -1
!      IF (JS .GT. abs(KORDI)) then
!        Set up an error message.
!      end if
!      MAXINT = max(MAXINT, abs(KORDI) - JS)
!   56 IF (KORDI .GE. 0) GO TO 70
!      KORDI = -1 - KORDI
!      JS = JS - 1
!      MAXDIF = max(MAXDIF, JS, 1)
!      GO TO 70
!++  End
!     **** EQUATION IS TO USE AN EXPLICIT METHOD
60    kqq = 1
   maxint = max(maxint, kordi)
70    if ((kordi .gt. maxord) .or. (kordi .le. 0)) then
!                    Set up error message.  KORDI is out of range.
      idat(1) = 17
      idat(2) = 1
      idat(3) = maxord
      if (nkdko .ne. 0) then
         kgo = 5
         ilow = nkdko
         ihi = nkdko + k - 1
         go to 640
      else
         kgo = 9
         idat(4) = kordi
         go to 650
      end if
   end if
   kord(k + 3) = kqq
80 continue
!     **** SET FLAGS WHICH DEPEND ON METHOD USED
!++  Code for STIFF is inactive
!      METHOD = 1
!      IF (MAXINT .GT. 0) IF (MAXDIF) 85,90,85
!      METHOD = -1
!   85 CONTINUE
!      KPRED = 5
!      GO TO 100
!++  End
90 method = 0
kpred = 1
100 continue
!
! ******* CHECK KORD FOR DIAGNOSTIC OUTPUT CONTROL *********
!
!++  Code for OUTPUT is active
if (iop10 .gt. 0) then
   if (noutko .ne. 0) then
      intchk(nxtchk) = 10
      ilow = noutko
!.    **** Check option for valid input
      kgo = 2
      go to 430
   end if
end if
!++  End
110 continue
!
! ********** CHECK KORD FOR INTEGRATION ORDER CONTROL ******
!
!++  Code for INTEGO is active
if (iop19 .ne. 0) then
!.           **** Check option for valid input
   intchk(nxtchk) = 19
   ilow = iop19
   jlim = -30
   kgo = 3
   go to 430
end if
!++  End
120 continue
!
! ********** CHECK SET UP FOR ERROR TOLERANCES *************
!
intchk(nxtchk) = 16
ilow = iop16
jlim = -5
kgo = 4
if (iop16 .ne. 0) go to 430
!.                      **** IN CURRENT CODE, IOP16=0 IS AN ERROR
kgo = 10
go to 650
150 continue
!     **** CHECK KORD STORAGE ALLOCATION
call optchk(intchk, iopt, 'DIVA / KORD$E')
if (nxtchk .lt. 0) kord2i = -4
!
!     ******** DONE CHECKING KORD STORAGE ALLOCATION *******
!
!     **** CHECK  Y  STORAGE ALLOCATION
intchk(2) = idimy
intchk(3) = ny + ny
nxtchk = 4
nyny = ny + 1
call optchk(intchk, iopt, 'DIVA / Y$E')
if (nxtchk .lt. 0) kord2i = -4
!
!     **** CHECK  F  STORAGE ALLOCATION
intchk(2) = idimf
intchk(3) = nte
nxtchk = 4
if (iop16 .ne. 0) then
!                                Error tolerance info.
   intchk(4) = 16
   intchk(5) = ntolf
   intchk(6) = ihi - iop16 + 1
   nxtchk = 7
end if
if (iop12 .gt. 0) then
   intchk(nxtchk) = 12
   intchk(nxtchk+1) = iop12
   intchk(nxtchk+2) = 4
   nxtchk = nxtchk + 3
end if
if (iop21 .gt. 0) then
   intchk(nxtchk) = 21
   intchk(nxtchk+1) = iop21
   intchk(nxtchk+2) = 1
   nxtchk = nxtchk + 3
end if
!
!++  Code for ERRSTO is inactive
!      IF (IOP20 .ne. 0) then
!c.                                Space for saving error estimates
!         INTCHK(NXTCHK) = 20
!         INTCHK(NXTCHK) = IOP20
!         INTCHK(NXTCHK) = NTE
!         NXTCHK = NXTCHK + 3
!      end if
!++  Code for STIFF is inactive
!      if (IOP21 .gt. 0) then
!c.                               Info. for stiff equations
!         INTCHK(NXTCHK) = 21
!         INTCHK(NXTCHK+1) = IOP21
!         INTCHK(NXTCHK+2) = IOP21S
!         NXTCHK = NXTCHK + 3
!      end if
!      MAXKQD = min(MAXKQI, 6)
!++  End
!                          Set aside space for the difference tables.
intchk(nxtchk) = 0
intchk(nxtchk+1) = -kdim * nte
intchk(nxtchk+2) = 0
nxtchk = nxtchk + 3
intchk(nxtchk) = -5 * nte
call optchk(intchk, iopt, 'DIVA / F$E')
if (nxtchk .lt. 0) then
   kord2i = -4
else if (kord2i .ne. -4) then
   do 290 k = nxtchk+1, intchk(nxtchk)
      if (intchk(intchk(k)) .eq. 0) then
         ndtf = intchk(intchk(k)+1)
         numdt = min(kdim, (intchk(intchk(k)+2)-ndtf+1) / nte)
         maxkqi = numdt - 1
      else
!         Take a quick return if needed space was not specified by user.
         kord2i = -4
      end if
290    continue
end if
if (iop9 + abs(iop10) + iop11 .ne. 0) then
! Insure user doesn't get in trouble with F not iniitalized.
   do 300 k = ndtf, ndtf + nte*numdt - 1
      f(k) = 0.d0
300    continue
end if
320 continue
if ((kord2i .eq. -4) .or. (iop10 .ne. 0)) then
   mact1(3) = iopiva(2)
   call mess(mact1, text1, iopt)
   kord1i = 24
   kord(1) = 24
end if
tmark = tspecs(1)
tmarkx = tspecs(4) + tspecs(2)
!
!     **** DONE WITH INITIALIZATION AND CHECKING INPUTS
if (iop13 + iop14 + iop15 .ne. 0) return
330 call divaa(tspecs, y, f, kord, divaf, divao)
return
!
! ************ LOOP TO CHECK OPTION SPECIFICATIONS *********
!
430 jl = 0
do 560 ihi = ilow, idimk
   j = kord(ihi)
   go to (460, 480, 490, 490), kgo
!     **** CHECK ON VARIATIONAL EQUATIONS
460    continue
!++  Code for VAREQ is active
   if (j - nte) 470, 565, 620
470    if (j .eq. 0) go to 560
   if (j .le. jl) go to 620
!++  End
!     **** Check on diagnostic output option
480    continue
!++  Code for OUTPUT is active
!.    **** CHECK IF DONE
   if (j .ge. nte) go to 565
   if (j .le. jl) go to 620
   go to 550
!++  End
490    continue
!     **** Check integration order control (KGO=3) and
!     **** error tolerance equation grouping (KGO=4).
   if (j - nte) 500, 565, 620
500    if (j) 510, 530, 540
510    if ((jl .le. 0) .and. (ihi .ne. ilow)) go to 620
   if (j .lt. jlim) then
!                         Output an error message.
      idat(2) = jlim
      idat(3) = 0
      go to 630
   end if
520    jl = -jl
   go to 560
530    if (kgo .eq. 3) go to 520
   kferr = ntolf + ihi - ilow
   if (f(kferr) .eq. cm1) go to 510
!                         Set up error message, TOL must be -1.
      idat(1) = kferr
      kgo = 11
      go to 650
540    if (abs(jl) .ge. abs(j)) go to 620
550    jl = j
560    continue
565 nxtchk = nxtchk + 3
intchk(nxtchk-2) = ilow
intchk(nxtchk-1) = ihi - ilow + 1
go to (30, 110, 120, 150), kgo
!
!     **** AN ERROR HAS BEEN MADE
!                  Error in setting up TSPECS for extra output
600 ihi = ihi + 2
610 ilow = -iop5
go to 630
!                  Error in KORD indices
620 idat(2) = abs(jl) + 1
idat(3) = nte
!                  Set up for print of message about KORD
630 idat(1) = intchk(nxtchk)
640 idat(4) = ihi
idat(5) = kord(ihi)
!
! ***************** Process Errors *************************************
!
650 kord2i = -4
mact(4) = ltxtaf
if (kgo .ge. 8) mact(4) = -1
mact(6) = mloc(kgo)
!--D Next line special: P=>S, X=>D
call dmess(mact, mtxtaa, idat, fdat)
if (kgo .lt. 8) then
   mact(10) = ilow
   mact(13) = ilow
   mact(15) = -min(ihi+2, idimk)
   call mess(mact(9), mtxtab, kord)
   if (kgo .le. 4) go to 565
end if
!              5   6   7    8    9   10   11  12
go to (100, 25, 25, 320, 100, 150, 660, 25), kgo - 4
660 kgo = 4
go to 565
end
!   End of DIVA

subroutine divaa(tspecs, y, f, kord, divaf, divao)
!>> 1989-02-24 DIVAA  Krogh   Big error with BETA(2)=1+epsilon -- looped
!>> 1988-07-27 DIVAA  Krogh   Fixed to allow restart on a restart.
!>> 1988-03-07 DIVAA  Krogh   Initial code.
!
!  MAIN SUBROUTINE FOR VARIABLE ORDER INTEGRATION OF ORDINARY
!  DIFFERENTIAL EQUATIONS
!
integer kord(*)
!--D Next line special: P=>D, X=>Q
double precision tspecs(*), y(*)
double precision f(*)
external divaf, divao
!
integer kdim, maxord, maxstf
!++ Substitute for KDIM, MAXORD, MAXSTF below
parameter (kdim = 20, maxord = 2, maxstf = 1)
!--D Next line special: P=>D, X=>Q
double precision tn
double precision xi(kdim)
!
!--D Next line special: P=>D, X=>Q
double precision tg(2), tgstop(2), tmark, tmarkx, tout, tolg
double precision alpha(kdim), beta(kdim+1)
double precision  d(maxstf+maxord,maxord), g(kdim,maxord)
double precision v(kdim+maxord)
double precision hc, hdec, hinc, hincc, hmax, hmaxp9, hmin
double precision fdat(11)
!
double precision ds(maxstf+maxord, maxord), gs(kdim)
double precision sigma(kdim), rbq(kdim), dnoise
double precision eave, eimax, eimin, emax, erep, robnd, snoise
!
!.    SPECIFICATION OF ENVIRONMENTAL CONSTANTS.
double precision eeps10, eeps16, erov10, eeps2
double precision eept75, eovep2, ovtm75, ovd10
common / divaev / eeps2, eept75, eovep2, ovtm75, ovd10, eeps10, &
   eeps16, erov10
save / divaev /
integer iopst, kordi, kqmaxd, kqmaxi, ldt, maxdif, maxint, nkdko, &
   nte, nyny, ndtf, numdt
common / divasc / tn, xi, iopst, kordi, kqmaxd, kqmaxi, ldt, &
   maxdif, maxint, nkdko, nte, nyny, ndtf, numdt
!
integer icf,ics,igflg,igtype(2),igstop(2),ilgrep,ings,iop3,iop4, &
   iop5,iop6,iop7,iop8,iop9,iop10,iop11,iop12,iop13,iop14,iop15, &
   iop16,iop17,iop18,iop19,iop20,iop21,iop22,iop21s,itolep,iy, &
   kemax,kis,kmark,kord1i,kord2i,kpred,kqdcon,kqicon,kqmaxs, &
   kqmxds,kqmxil,kqmxip,kqmxis,ksc,ksout,ksstrt,kstep,lex,linc, &
   lincd,lincq,lsc,maxkqd,maxkqi,method,ne,neptol,ng,ngtot, &
   noiseq,noutko,ntolf,ny,idat(6)
common /divamc/ tg,tgstop,tmark,tmarkx,tout,tolg,hc,hdec,hinc, &
   hincc,hmax,hmaxp9,hmin,alpha,beta,d,g,v,ds,gs,sigma,rbq,dnoise, &
   eave,eimax,eimin,emax,erep,robnd,snoise,fdat,icf,ics,igflg, &
   igtype,igstop,ilgrep,ings,iop3,iop4,iop5,iop6,iop7,iop8,iop9, &
   iop10,iop11,iop12,iop13,iop14,iop15,iop16,iop17,iop18,iop19, &
   iop20,iop21,iop22,iop21s,itolep,iy,kemax,kis,kmark,kord1i, &
   kord2i,kpred,kqdcon,kqicon,kqmaxs,kqmxds,kqmxil,kqmxip,kqmxis, &
   ksc,ksout,ksstrt,kstep,lex,linc,lincd,lincq,lsc,maxkqd,maxkqi, &
   method,ne,neptol,ng,ngtot,noiseq,noutko,ntolf,ny,idat
save / divamc / , / divasc /
!
double precision cmp75, c0, c1m5, c1m3, c2p5m3, c8m3, cp0625, cp1
double precision cp25, cp3, cp4, cp5, cp875, c1, c1p125, c1p3, c2
double precision c6, c10, c16, c4096
parameter (cmp75 = -.75d0)
parameter (c0 = 0.d0)
parameter (c1m5 = 1.d-5)
parameter (c1m3 = 1.d-3)
parameter (c2p5m3 = 2.5d-3)
parameter (c8m3 = 8.d-3)
parameter (cp0625 = .0625d0)
parameter (cp1 = .1d0)
parameter (cp25 = .25d0)
parameter (cp3 = .3d0)
parameter (cp4 = .4d0)
parameter (cp5 = .5d0)
parameter (cp875 = .875d0)
parameter (c1 = 1.d0)
parameter (c1p125 = 1.125d0)
parameter (c1p3 = 1.3d0)
parameter (c2 = 2.d0)
parameter (c6 = 6.d0)
parameter (c10 = 10.d0)
parameter (c16 = 16.d0)
parameter (c4096 = 4096.d0)
!
integer ldis, kexit, i, j, k, j1, j2, l, lx
double precision tp, tp1, tp2, tp3, hh, disadj
double precision sigmas, tps1, tps2, eimino, exr
!--D Next line special: P=>D, X=>Q
double precision  tmarka(2), xp, xp1
equivalence (g(1, 1), hh), (tmarka(1), tmark)
equivalence (kexit, iop17)
save eimino, exr, sigmas, tp, tp1, tp2, tp3, tps1, tps2, disadj, &
   ldis
!
!                      Declarations for error message processing.
!
integer mact(17), mloc(8), mentxt, meret, meemes, metext
parameter (mentxt =23)
parameter (meret  =51)
parameter (meemes =52)
parameter (metext =53)
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
integer ltxtaa,ltxtab,ltxtac,ltxtad,ltxtae,ltxtaf,ltxtag,ltxtah, &
 ltxtai,ltxtaj,ltxtak,ltxtal,ltxtam,ltxtan
parameter (ltxtaa=  1,ltxtab=  8,ltxtac= 40,ltxtad= 80,ltxtae=129, &
 ltxtaf=189,ltxtag=237,ltxtah=272,ltxtai=302,ltxtaj=342, &
 ltxtak=390,ltxtal=455,ltxtam=484,ltxtan=531)
character mtxtaa(3) * (186)
!
integer locm, mlocac, mlocad, mlocae, mlocaf, mlocag, mlocah, &
   mlocai, mlocaj

parameter (locm = 32 * 256)
!                     KORD1I   Severity   Loc. message
parameter (mlocac = 23 + 32 * (99 + 256 * ltxtac))
parameter (mlocad = 13 + 32 * (38 + 256 * ltxtad))
parameter (mlocae = 22 + 32 * (38 + 256 * ltxtae))
parameter (mlocaf =  3 + 32 * (14 + 256 * ltxtaf))
parameter (mlocag = 21 + 32 * (38 + 256 * ltxtag))
parameter (mlocah =  2 + 32 * (25 + 256 * ltxtah))
parameter (mlocai = 11 + 32 * (38 + 256 * ltxtai))
parameter (mlocaj = 12 + 32 * (38 + 256 * ltxtaj))
!
data mtxtaa/'DIVAA$BAt: TN=$F, KSTEP=$I, with H=$F$EA previously r &
eported error was fatal.$EPrint points not properly ordered: tspec &
($i)=$f$EAn error tolerance of 0 requires setting special flags. $ &
 ','$BStep size reduced too fast, doing a restart.  $bh is so smal &
l that tn + h = tn.  $BError tolerance too small.  $BStep size at$ &
 end of start < hmin=$f, $BError estimates require a steps','ize < &
 hmin=$f, $b(Estimated Error) / (Requested Error) for equation $i$ &
 is $f.  $BTolerance $i is f($i) = $f.$BTolerance $i is f($i) * f( &
$i) = $f * $f = $f.$b  Replacing f($i) with $f.$e'/
data mloc / mlocac, mlocad, mlocae, mlocaf, mlocag, mlocah, &
   mlocai, mlocaj /
!
!                      1 2 3 4       5 6       7       8 9      10
data mact / meemes,0,0,0, mentxt,0, metext, mentxt,0, metext, &
  mentxt, 0, metext, mentxt, ltxtan, metext, meret /
!           11 12      13      14      15      16     17
!
data exr, eimino / cp1, c8m3 /
data ldis / 0 /
!
! ************** START OF EXECUTABLE CODE ******************************
!
660 if (kord2i) 670, 1380, 1840
670 if (kord2i .eq. -1) go to 2140
if (kord2i .eq. -5) go to 720
if (kord2i+8) 1840, 710, 1840
!     **** SPECIAL OUTPUT CASE (EXTRAPOLATION OR GSTOP)
680 go to (1220, 1190, 830, 690, 690, 2190), kexit
!     **** RESET TMARK BEFORE GOING WHERE DIRECTED BY KEXIT
690 kexit = kexit - 2
go to 1890
!     **** USER HAS REQUESTED A RESTART
700 kord2i = -5
igflg = 0
if ((kord(2) .ge. 2) .and. (lsc .lt. 3)) kord2i = -8
if ((kord1i .le. 3) .or. (kord1i .eq. 5)) go to 1890
if (kord2i .eq. -5) go to 720
!                Set up for a discontinuity
710 xp = tspecs(1)
if (kord(2) .eq. 3) then
!                Adjust for discontinuity in Y
   j = 0
   do 712 i = 1, nte
      if (nkdko .ne. 0) kordi = kord(nkdko + i -1)
      k = 1
      j = j + kordi
      xp1 = y(j)
711       y(nyny + j - k) = y(nyny + j - k) + xp1
      if (k .lt. kordi) then
         xp1 = y(j-k) + (tn - xp) * xp1 / dble(k)
         k = k + 1
         go to 711
      end if
712    continue
else
   igflg = -3
end if
disadj = hh
xp1 = (xp - tn) / xi(1)
if (xp1 .lt. cp25) then
   k = 1
   if (xp1 .lt. cmp75) k = 2
   tspecs(1) = tn - xi(k)
   tspecs(2) = 2.d0*(tn - tspecs(1))
   call divain(tspecs(1), y, f, kord)
   do 713 j = 1, ny
      y(nyny + j - 1) = y(j)
713    continue
!          Move difference tables back one step
   tn = tspecs(1)
714    call divabu(f, kord)
   if (k .eq. 2) then
      ksc = max(ksc - 1, 1)
      do 715 k = max(1, ksc), iop11-1
         beta(k+1) = beta(k) * (xi(k) / (xi(k+1) - xi(1)))
715       continue
      k = 0
      go to 714
   end if
end if
!      Take step to discontinuity.
hh = (xp - tn)
erep = -abs(erep)
kis = 1000
ldis = 1
lsc = 0
linc = 0
hincc = c1p125
ksstrt = kstep + 2
go to 1115
! ********
! INITIALIZE FOR STARTING AN INTEGRATION
! ********
720 hh = tspecs(2)
lsc = 8
linc = -5
lincd = 64
ldt = -4
kqmaxi = 0
eimin = cp3
eave = c10
xi(1) = c0
igflg = 0
kis = 0
ksstrt = 0
robnd = 0.d0


! GO COMPUTE INITIAL DERIVATIVES
730 kord2i = -6
!++  Code for VAREQ is active
ics = 0
icf = ne
!++  End
740 kord1i = kpred
go to 1360
!   RETURN AFTER COMPUTING INITIAL (OR NOISE TEST) DERIVATIVES
750 continue
!++  Code for VAREQ is active
if (iop18 .eq. 0) go to 790
!.    **** SPECIAL LOGIC TO COMPUTE VARIATIONAL DERIVATIVES
760 if (kord1i .eq. 3) go to 770
kord1i = 3
kord(3) = 0
770 if (icf .eq. nte) go to 790
if (kord2i .eq. -6) then
   if (ics .eq. 1) go to 790
end if
ics = icf + 1
780 kord(3) = kord(3) + 1
icf = iop18 + kord(3)
icf = abs(kord(icf))
if (icf) 1360, 780, 1360
!++  End
790 ics = 1
icf = ne
!++  Code for VAREQ is active
if (kord2i .eq. 0) if (erep) 2220, 2220, 1430
!++  End
if (linc + 5) 1490, 810, 1490
! END OF SPECIAL CODE FOR INITIALIZATION AT THE START
! ********
! UPDATE VARIABLES TO PREPARE FOR NEXT STEP
! ********
800 ldt = 0
eimin = c2p5m3+eimin*(c6*eave+eimax)/((eimin+c1)*(c6*eimax+eave))
810 do 820 j = 1, ny
820    y(nyny + j - 1) = y(j)
tn = tspecs(1)
! ********
! TEST FOR VARIOUS TYPES OF OUTPUT
! ********
!     TEST FOR END OF STEP OUTPUT (OR IF DUMP OUTPUT TO BE TESTED FOR)
if (iop6 .eq. 0) go to 840
!   SET UP FOR END OF STEP OUTPUT
kord1i = 6
go to 1810
!     SET UP AFTER OTHER COMPUTATIONS HAVE BEEN MADE (TSPECS(1).NE.TN)
830 kord1i = 6
go to 860
!     TEST FOR AN OUTPUT POINT
840 if (hh * (tmark - tn)) 1760, 1760, 850
!     TEST FOR TOO MANY STEPS OUTPUT
850 continue
if (ksout .gt. kstep) go to 890
kord1i = 4
860 if (tspecs(1) .eq. tn) go to 1810
!     GO INTERPOLATE VALUES AT END OF LAST STEP
870 tspecs(1) = tn
go to 1780
!     CONTINUE AFTER TOO MANY STEPS OUTPUT
880 continue
ksout = kstep + iop4
890 continue
!++  Code for DUMP is active
if (iop9 .eq. 0) go to 920
!++  Code for DUMP & STIFF is inactive
!      KQMXDS=KQMAXD
!++  Code for DUMP is active
kqmxis = kqmaxi
kis = kis + 1
if (kis .ne. 0) go to 920
!.   TIME TO DUMP THE SOLUTION
900 kord1i = 9
!.    SET TO UPDATE THE DIFFERENCE TABLE
if (ldt .eq. 1) go to 1810
! Note that doing this update can lead to very small differences in the
! results because of round off differences.
ldt = -2
go to 1780
!++  End
! RETURN AFTER DUMPING THE SOLUTION
910 continue
!++  Code for DUMP is active
kis = 2
!.    TEST IF SOLUTION DUMP DUE TO RESTART, END, OR
!.    DROP IN INTEG. ORDER
if (linc .lt. 0) if (linc + 6) 1860, 1750, 1180
!++  End
! END OF TESTING FOR VARIOUS TYPES OF OUTPUT
!   TEST IF STEPSIZE MAY BE INCREASED OR IF TESTS SHOULD BE MADE
!   FOR DECREASING THE STEPSIZE (OR IF STARTING OR USER SELECTING H)
920 kstep = kstep + 1
if (linc) 930, 980, 940
930 continue
!     **** ONLY POSSIBLE VALUES AT THIS POINT ARE
!          LINC = -2 OR -5
if (linc + 5) 1120, 1110, 1120
! ********
! ERROR ESTIMATES INDICATE STEPSIZE CAN BE INCREASED
! ********
940 hc = hincc
if (linc .gt. 1) hc = hc ** linc
hh = hc * hh
!     TEST IF NEW STEPSIZE IS TOO BIG
if (abs(hh) .gt. hmax) if (hmax) 970, 970, 960
950 eave = eimax
robnd = cp3 + hincc
go to 1110
!     NEW STEPSIZE IS TOO BIG
960 if (abs(xi(1)) .ge. hmaxp9) go to 970
hh = sign(hmax, hh)
go to 950
!     RESTORE THE OLD STEPSIZE
970 hh = xi(1)
linc = 0
go to 1150
! END OF CODE FOR CASE WHEN ERROR ESTIMATES INDICATE STEPSIZE INCREASE
! ********
! TEST IF ESTIMATED ERRORS INDICATE STEPSIZE SHOULD BE DECREASED
! ********
980 robnd = c1p3
if (eimax .le. eave) go to 990
eave = eave + cp4 * (eimax - eave)
if ((eimax * emax) - c1m3) 1000, 1010, 1010
990 eave = eimax
if ((eimax * emax) .ge. eimin) go to 1010
1000 robnd = cp3 + (sigma(kqmaxs)/sigma(kqmaxs-1))
go to 1180
!     TEST IF STEPSIZE SHOULD BE REDUCED
1010 if (emax * eimax .lt. exr * eave) go to 1180
! ********
! ERROR ESTIMATES INDICATE STEPSIZE SHOULD BE REDUCED
! ********
hc = hdec
if (eimin .le. eimino) go to 1030
eimin = eimino
hc = cp875
1030 hh = hc * xi(1)
if (lsc - 1) 1040, 1080, 1090
1040 if (abs(hh) .ge. hmin) go to 1090
if (abs(cp875*xi(1)) .le. hmin) if (linc) 1050, 970, 970
hh = sign(hmin, hh)
go to 1090
!     STEPSIZE IS TOO SMALL TO BE REDUCED
!     SET UP ERROR INDICATORS AND PREPARE FOR RETURN TO USER
1050 kord1i = 8
go to 2240
!     PROCEED WITH CURRENT STEPSIZE DESPITE ERROR BEING TOO BIG
1070 hh = xi(1)
tspecs(2) = hh
emax = c0
linc = 0
if (kord1i - 2) 1420, 1430, 1430
!     SET LSC TO END STARTING PHASE
1080 lsc = 2
!     CHECK IF REPEATING A STEP
1090 if (linc .ne. -1) go to 1110
!   WHEN REPEATING A STEP, BACK UP THE DIFFERENCES AND STEPSIZE INFO.
1100 call divabu(f, kord)
!   TEST IF NOISE TEST (LINC = -7) OR IF H IS NOT
!     BEING CHANGED (LINC = -4)
if (linc + 4) 1780, 1180, 1110
! ********
! STEPSIZE IS TO BE CHANGED
! ********
1110 continue
! MODIFY STEPSIZE TO REDUCE ROUNDOFF ERROR IN ACCUMULATING INDEP. VAR.
tp = c2 * abs(tn) + c4096 * abs(hh)
tp = (tp + abs(hh)) - tp
if (tp .ne. c0) hh = sign(tp, hh)
!     TEST IF NEW STEPSIZE SELECTED ACTUALLY GIVES A CHANGE
if (hh .eq. tspecs(2)) go to 1140
1115 tspecs(2) = hh
if (iop8 .eq. 0) go to 1140
!     SETUP TO TELL USER ABOUT STEPSIZE CHANGE (OR TO CHANGE STEPSIZE)
1120 kord1i = 8
go to 860
!     RETURN AFTER TELLING USER ABOUT STEPSIZE CHANGE
1130 hh = tspecs(2)
1140 if (hh .ne. xi(1)) kqicon = -1
1150 hc = min(eovep2, abs(hh)) / eeps2
! ********
! PREPARE FOR BEGINNING A NEW STEP
! ********
if (linc .gt. 0) then
   linc = min(linc, lincq) + lincq
   go to 1190
end if
1180 linc = lincd
1190 if (hc .gt. abs(tn)) go to 1200
!     **** GIVE SINGULARITY DIAGNOSTIC
kord1i = 5
go to 2240
1200 tspecs(1) = tn + hh
if (lex .eq. 0) go to 1250
if (hh * (tspecs(1) - tmarkx) .lt. c0) go to 1250
tspecs(1) = tmarkx
hh = tmarkx - tn
linc = 64
if (lex .gt. 0) go to 1240
if ((lsc .lt. 4) .and. (hh / xi(1) .lt. cp3)) go to 1230
1220 hh = cp875 * hh
go to 1110
!     **** GIVE OUTPUT AT CURRENT TMARK (WITH EXTRAPOLATION)
1230 kord1i = -kmark
go to 1770
!     **** INTEGRATE TO TMARKX
1240 kqicon = -1
!   TEST IF SUBROUTINE FOR COMPUTING INTEGRATION COEFF. SHOULD BE CALLED
1250 continue
!++  Code for STIFF is inactive
!      IF ((KQMAXI .LT.KQICON) .OR. (KQMAXD.LT.KQDCON)) GO TO 1320
!++  Code for ~STIFF is active
if (kqmaxi .lt. kqicon) go to 1320
!++  End
!   GO COMPUTE COEFFICIENTS REQUIRED FOR THE INTEGRATION
!     TEST IF STARTING
if (lsc .lt. 7) go to 1310
1260 kqmaxi = 2
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
1270 kqmaxd = 0
call divahc
j = ndtf
do 1300 i = 1, nte
!++  Code for STIFF is inactive
!         if (KORD(I + 3) .le. 0) go to 1290
!++  End
   kord(i + 3) = 1
!     INITIALIZE THE DIFFERENCE TABLE
   if (ldt .eq. -4) f(j) = f(i)
   f(j + 1) = c0
   f(j + 2) = c0
1290    continue
   j = j + numdt
1300    continue
if (lsc .eq. 5) go to 1340
lsc = 7
ldt = 1
go to 1330
!   INTEGRATION IS NOT BEING STARTED
1310 k = kord(kemax + 3)
sigmas = sigma(k)
call divahc
!     **** ADJUST EAVE
tps1 = beta(k)
if (tps1 .gt. c1) tps1 = cp5 * tps1 + cp5
eave = eave * tps1 * (sigma(k) / sigmas)
!     TEST BELOW USED TO GET SAME RESULTS WITH/WITHOUT EXTRA EQUATIONS
if (k .gt. kqicon) lsc = max(lsc, -3)
! END OF SPECIAL LOGIC FOR CASE WHEN INTEG. COEFF. ROUTINE IS CALLED
1320 continue
! ********
! PREDICT Y
! ********
1330 continue
!++  Code for ~ARGM is active
call divapr(y, y(nyny), f, kord)
!++  Code for ARGM is inactive
!      CALL DIVAPE
!++  End
!     GO GET PREDICTED DERIVATIVES
1340 kord1i = kpred
! ********
! CALL DIVAF  (OR RETURN)
! ********
1350 kord2i = 0
1360 kord(1) = kord1i
kord(2) = 0
if (iop13 .ne. 0) return
call divaf(tspecs(1), y, f, kord(1))
!     TEST FOR SPECIAL USER RETURN
1380 if (kord(1) .lt. 0) go to 2130
!     TEST FOR SPECIAL CASE
if (kord2i .ne. 0) go to 660
! ********
! TRANSFER CONTROL TO PROPER PLACE AFTER COMPUTING DERIVATIVES
! ********
if (kord1i - 2) 1400, 800, 1390
1390 continue
!++  Code for VAREQ is active
if (ics - icf) 1410, 1410, 760
!++  End
! ********
! PREPARE FOR CORRECTING, AND CORRECT Y
! ********
1400 itolep = 0
ilgrep = 0
iy = 1
eimax = c1m5
emax = c0
kqmaxi = 2
kqmaxs = 2
!++  Code for STIFF is inactive
!      IF (METHOD) 1404,1410,1406
! 1404 KQMAXI=0
! 1406 KQMAXD=2
!++  End
1410 continue
!++  Code for ~ARGM is active
call divacr(y, f, kord, f(ntolf), kord(iop16))
!++  Code for ARGM is inactive
!      CALL DIVACE
!++  End
!     TEST IF ESTIMATED ERROR IS TOO BIG (OR IF DIAGNOSTIC CALLED FOR)
if (emax .gt. erep) if (erep) 2210, 2210, 1670
1420 continue
!++  Code for VAREQ is active
if (iop18 .ne. 0) go to 760
!++  End
1430 kord1i = 2
!     TEST IF NOISE APPEARS TO LIMIT PRECISION
if (emax .lt. c0) go to 1470
!++  Code for ~STIFF is active
if (lsc) 1450, 1360, 1610
!++  Code for STIFF is inactive
!      IF (LSC) 1450,1460,1610
!++  End
!     SET LSC=0 IF NOISE NO LONGER APPEARS TO LIMIT PRECISION
!     OR IF THE END OF THE STARTING PHASE HAS BEEN REACHED
1450 lsc = 0
1460 if (method) 800, 1350, 1350
! ********
! NOISE APPEARS TO BE LIMITING PRECISION
! ********
1470 continue
if (lsc .le. 0) lsc = max(lsc - 1, -kqmaxs)
if (lsc .eq. -1) go to 1460
if (abs(emax) .lt. exr) go to 1590
linc = -7
tps2 = (c1 + beta(noiseq - 1)) ** noiseq
if (snoise .lt. eeps10 * tps2) go to 1550
tp = sign(eept75 * abs(tn) + ovtm75, hh)
if (abs(tp) .gt. abs(hh)) go to 1550
tspecs(1) = tn + tp
kord1i = 0
!     **** GO TO BACK UP THE DIFFERENCES AND GET F(TSPECS(1))
go to 1100
!     **** SOLUTION HAS BEEN INTERPOLATED AND F COMPUTED
1490 continue
kord1i = 0
linc = linc - 1
if (linc + 9) 1510, 1520, 1500
1500 tspecs(1) = tn + (tp + tp)
tp1 = f(kemax)
tp2 = f(ndtf + numdt * kemax - numdt)
go to 1780
!     **** COMPUTE 2-ND DIFFERENCE AT CLOSE SPACED T VALUES
1510 tp2 = tp3
1520 tp3 = f(kemax)
tps1 = abs((tp3 - tp1) - (tp1 - tp2))
if ((c16 * tps1 * tps2) .ge. dnoise) if (linc + 9) 1550, 870, 1550
1530 continue
tps2 = cp25 * snoise / rbq(noiseq)
do 1540 k = 2, numdt
   tps1 = tps1 + tps1
   rbq(k) = max(tps1, tps2 * rbq(k))
1540    continue
linc = 0

!FTK Next two lines added 2009-10-15
if (abs(emax) .lt. erep) go to 1460
!FTK  LINC = -1  And then on 2015-03-14 commented out this line

hh = cp875 * hh
go to 1040
!     **** SET UP TO GIVE NOISE DIAGNOSTIC
1550 kord1i = 6
go to 2240
!     **** AFTER GIVING NOISE DIAGNOSTIC
1560 kord1i = 2
if (kord(2) .ge. 0) then
  tps1 = eeps10
!FTK Next line added 2009-10-15
  if (tps1 .lt. .49d0 * rbq(2)) go to 1530
end if
!     **** SET NEW VALUE FOR OFFENDING TOL
f(ntolf + itolep - 1) = fdat(7)
if (linc + 7) 1180, 1570, 1180
1570 linc = 0
1580 if (lsc) 1460, 1460, 1610
!     **** CHANGE HINCC AND ADJUST SIGMA( )
1590 if (lsc .ne. -4) go to 1580
if (hincc .eq. c1p125) go to 1580
tps1 = c1p125 / hincc
tps2 = 1.0d0
do 1600 k = 2, iop11
   tps2 = tps2 * tps1
   sigma(k) = sigma(k) * tps2
1600    continue
eave = eave * tps1 ** (1-kord(kemax+3))
lincd = 6
lincq = 12
hincc = c1p125
go to 1460
!   END OF CODE FOR CASE WHEN NOISE APPEARS TO LIMIT PRECISION
! ********
! SPECIAL LOGIC FOR STARTING THE INTEGRATION
! ********
1610 if (lsc .eq. 1) go to 800
lsc = lsc - 1
if (lsc - 2) 1620, 1640, 1650
1620 if (eimax .le. (cp0625*eave*(sigma(kqmaxs)/sigmas)*(beta(kqmaxs+ &
   1))**2)) go to 800
1630 ksstrt = kstep + 2
!   TEST IF STEPSIZE IS TOO SMALL BEFORE ENDING STARTING PHASE
if (abs(hh) .ge. hmin) go to 1450
!     GIVE DIAGNOSTIC FOR STEPSIZE TOO SMALL AT END OF START
kord1i = 7
go to 2240
!   SET LSC TO DO ONE DERIVATIVE EVAL. PER STEP
1640 lsc = 1
go to 800
!     TEST IF FIRST TIME THROUGH THE FIRST STEP
1650 if (lsc .eq. 6) go to 1340
!     END STARTING PHASE IF CONVERGENCE OF CORRECTOR ITERATES TOO SLOW
if (ldt .eq. -5) go to 1660
lsc = min(lsc, 4)
go to 800
1660 ldt = 0
if (lsc - 4) 1260, 1630, 1260
! END OF SPECIAL LOGIC FOR STARTING THE INTEGRATION
! ********
! ESTIMATED ERROR IS TOO BIG
! ********
1670 if (beta(2) - c1) 1690, 1730, 1680
1680 hc = c1 / beta(2)
if (beta(2) .ge. c1p125) go to 1740
1690 if (beta(2) .gt. cp1) go to 1730
!   REQUIRED STEPSIZE REDUCTION IS TOO RAPID -- GIVE A DIAGNOSTIC
kord1i = 4
go to 2240
!
!     TEST KORD(2) AFTER ABOVE DIAGNOSTIC OR A DISCONTINUITY DIAGNOSTIC
1700 continue
if (kord(2) .eq. 0) go to 1730
!  TEST IF SOLUTION MUST BE DUMPED BEFORE A RESTART
1710 linc = -1
!++  Code for DUMP is active
if (iop9 .eq. 0) go to 1750
if (kis .eq. 2) go to 1750
linc = -6
!.    GO DUMP SOLUTION BEFORE REPEATING THE STEP
1720 kqmaxi = kqmxis
!++  Code for DUMP & STIFF is inactive
!      KQMAXD=KQMXDS
!++  Code for DUMP is active
call divabu(f, kord)
go to 900
!++  End
!   SET UP TO REPEAT THE STEP
1730 hc = cp5
1740 linc = -1
if (lsc .le. 3) go to 1030
!   RESTART THE INTEGRATION IF ERROR IS TOO BIG ON FIRST OR SECOND STEP
! LOOP TO SELECT A NEW INITIAL STEPSIZE
1750 lsc = 7
1755 hh = hh * cp5
emax = emax * cp25
if (emax .ge. cp3) go to 1755
go to 1090
!   END OF SELECTING A NEW INITIAL STEPSIZE
! END OF LOGIC FOR CASE WHEN ESTIMATED ERROR IS TOO BIG
! ********
! INTEGRATION HAS REACHED AN OUTPUT POINT
! ********
1760 if (kmark .eq. 0) go to 1920
kord1i = min(kmark, 5)
kord(3) = kmark
if (tspecs(1) .eq. tmark) go to 1790
1770 tspecs(1) = tmark
1780 call divain(tspecs(1), y, f, kord)
1790 continue
if (kord1i) 1800, 730, 1810
!   OUTPUT POINT IS OBTAINED BY EXTRAPOLATION
1800 continue
!++  Code for EXTRAP is active
kord1i = -kord1i
kord2i = -7
kexit = 4
!.  TEST IF GSTOP-S ARE PRESENT
!++  Code for EXTRAP &  GSTOP is active
if (ngtot .eq. 0) go to 1820
igflg = 4
kexit = 2
kord1i = 7
if (iop7) 740, 1820, 740
!++  End
! ********
! CALL DIVAO  (OR RETURN)
! ********
1810 kord2i = 1
1820 kord(1) = kord1i
kord(2) = 1
if (iop14 .ne. 0) return
call divao(tspecs(1), y, f, kord(1))
!     TEST FOR SPECIAL USER RETURN OR FOR A RESTART
!++  Code for ~DUMP is inactive
! 1840 IF (KORD(1)) 2130,700,1880
!++  Code for DUMP is active
1840 if (kord(1) .gt. 0) go to 1880
1850 if (iop9 .eq. 0) go to 1870
!.    **** GO DUMP THE SOLUTION
linc = -7
itolep = kord(1)
idat(1) = kord(2)
neptol = kord1i
if (lsc .ne. 8) go to 900
1860 linc = min(0, lincd)
kord1i = neptol
kord(1) = itolep
kord(2) = idat(1)
1870 if (kord(1)) 2130, 700, 2100
!++  End
1880 if (kord2i .lt. 0) go to (2140, 1810, 1350, 2110, 720, 750, 680, &
   710), -kord2i
if (kord2i .eq. 0) go to 1380
! ********
! TRANSFER CONTROL TO PROPER PLACE AFTER OUTPUT
! ********
1890 if (kord1i - 5) 1910, 1930, 1900
1900 if (kord1i - 8) 840, 1130, 910
1910 if (kord1i - 3) 1920, 1930, 880
!   GET NEW TOUT
1920 tout = tspecs(1) + tspecs(3)
! GET NEW TMARK (NEXT INDEP. VAR. OUTPUT POINT)
1930 xp = tmark
k = kmark
tmark = tout
kmark = 2
lex = 0
if (iop5) 1940, 1980, 1970
1940 i = -iop5
1950 i = i + 3
j1 = kord(i - 3)
if (j1) 1950, 1980, 1960
1960 j2 = kord(i - 2)
l = kord(i - 1)
go to 1990
1970 j1 = 5
j2 = iop5
l = 0
if (j2 .ge. j1) go to 1990
1980 j1 = 4
j2 = 4
l = iop3
!
!     **** LOOP TO SET NEW TMARK (AND TMARKX)
1990 do 2060 j = j1, j2
!        **** TEST IF EXTRAPOLATION NOT POSSIBLE
   if (l .eq. 0) go to 2010
   lx = 2
   if (lex) 2020, 2030, 2020
2000    lex = l
2010    lx = 1
2020    if (hh * (tspecs(j) - tmarka(lx))) 2030, 2060, 2060
2030    if (j .eq. 4) go to 2050
   if (hh * (tspecs(j) - xp)) 2060, 2040, 2050
2040    if ((k .ge. j) .or. (k .eq. 3)) go to 2060
2050    tmarka(lx) = tspecs(j)
   if (lx .eq. 2) go to 2000
   kmark = j
2060    continue
if (iop5 .lt. 0) go to 1950
if (j1 .ne. 4) go to 1980
if (kmark .eq. 4) kmark = 3
!     **** TEST IF NEW TMARK IS ACCEPTABLE
if (hh * (xp - tmark)) 2070, 2080, 2090
2070 if (kord2i - 1) 670, 840, 670
2080 if (k .ne. kmark) go to 2070
!++  Code for DUMP is active
if (kord1i .eq. 3) go to 1850
!++  Code for ~DUMP is inactive
!      IF (KORD1I .EQ. 3) GO TO 2100
!++  End
2090 if (kord1i .eq. 13) go to 2190
! SETUP TO INDICATE ERROR IN SPECIFICATION OF OUTPUT POINTS
kord1i = 2
idat(2) = kmark
if (kmark .le. 3) idat(2) = kmark + 1
fdat(3) = tspecs(idat(2))
go to 2240
!     SET KORD1I=1 TO INDICATE THAT END OF INTEGRATION HAS BEEN REACHED
2100 kord1i = 1
! ********
! RETURN TO USER
! ********
2110 kord2i = -1
kord(1) = kord1i
2130 kord(2) = -1
return
! ********
! TRANSFER CONTROL TO PROPER PLACE AFTER RETURN TO USER
! ********
2140 if (kord1i - 2) 2150, 1560, 2160
2150 kord2i = 1
go to 1930
2160 if (kord1i - 4) 1700, 2200, 2170
2170 if (kord1i - 13) 2180, 1930, 2190
2180 if (abs(hh) .ge. hmin) if (kord1i - 11) 1030, 1450, 1030
if (kord(2) .eq. 0) if (kord1i - 11) 1070, 800, 1070
!   ERROR MESSAGES HAVE BEEN IGNORED -- COMPUTATION CAN NOT CONTINUE
2190 kord1i = 1
go to 2240
!
!        AFTER A DISCONTINUITY RETURN
2200 linc = -4
if (kord(2)) 1710, 1730, 1100
! ********
! PROBLEM ENCOUNTERED WHEN CORRECTING
! ********
2210 if (ldis .eq. 0) go to 2230
!           Extra checks when had a user specified discontinuity.
!++  Code for VAREQ is active
if (iop18 .ne. 0) go to 760
!++  End
2220 kord1i = 2
ldis = ldis + 1
tp = disadj / hh
if (kis .ge. 1000) then
   if (ldis .eq. 2) then
      if (kqmaxs .le. 3) then
         ldis = 0
         erep = abs(erep)
         tspecs(2) = hh*min(min(tp, tp**2), &
            (cp25 * exr/emax)**.333333333d0)
         go to 720
      end if
      linc = -5
      if (iop9 .eq. 0) kis = 1001
      go to 800
   end if
   if (iop9 .eq. 0) kis = kis + 1
   if (kqmaxs .le. ldis + 2) kis = ldis + 1
   linc = min(linc, ldis-2)
end if
if (ldis .gt. 2*kqmaxs) then
   erep = abs(erep)
   ldis = 0
   if (emax .gt. erep) go to 1670
   go to 1430
end if
if (tp .ge. hincc**(linc+2)) then
   if ((ldis .ne. 3) .and. (tp .gt. dble(kqmaxs))) lsc = 1
   eimin = cp5
   eave = eave * tp**8
end if
if (lsc .eq. 2) go to 1630
if (emax .gt. exr) go to 1730
go to 1430
!
2230 erep = abs(erep)
!++  Code for DUMP is active
if (linc .lt. -3) go to 1720
!++  End
!     BAD TOL
kord1i = 3
! ********
! ERROR PROCESSING
! ********
2240 fdat(1) = tn
fdat(2) = hh
idat(1) = kstep
itolep = max(neptol, -neptol - 1)
j = 3
if (kord1i .ge. 7) then
   j = 4
   fdat(3) = hmin
end if
if (kord1i .le. 3) then
   if (kord1i .lt. 3) then
      k = 8
   else
      mact(9) = ltxtal
      fdat(3) = c0
      idat(2) = itolep
      idat(3) = itolep + ntolf - 1
      k = 11
   end if
else
   mact(9) = ltxtak
   fdat(j) = emax
   idat(2) = kemax
   idat(3) = itolep
   idat(4) = itolep + ntolf - 1
   fdat(j+1) = f(idat(4))
   k = 14
   if (kord1i .eq. 6) then
      k = 17
      idat(5) = idat(4)
      fdat(7) = 32.d0 * abs(emax) * fdat(j+1)
      fdat(j+2) = fdat(7)
   end if
   mact(12) = ltxtal
   if (neptol .lt. 0) then
      mact(12) = ltxtam
      idat(6) = idat(4)
      idat(5) = idat(4) + 1
      fdat(j+2) = f(idat(5))
      fdat(j+3) = fdat(j+1) * fdat(j+2)
   end if
end if
! Set the location for the first part of the message that varies, set
! the error severity, and the index number, print the error and
! return or stop.
l = mloc(kord1i)
mact(6) = l / locm
mact(2) = (l - mact(6) * locm) / 32
kord1i = mod(l, 32)
mact(3) = kord1i
mact(k) = meret
!--D Next line special: P=>S, X=>D
call dmess(mact, mtxtaa, idat, fdat)
mact(k) = mentxt
go to 2110
!
end
!   End of DIVAA

subroutine divabu(f, kord)
!>> 1987-12-07 DIVABU Krogh   Initial code.
!
! THIS SUBROUTINE RESTORES THE DIFFERENCE TABLE TO ITS STATE
! AT THE BEGINNING OF THE CURRENT STEP.  IF THE INTEGRATION ORDER
! WAS INCREASED, IT IS REDUCED. THE COMMON ARRAY XI IS ALSO
! RESTORED TO ITS STATE AT THE BEGINNING OF THE STEP. IF THE
! STEPSIZE IS NOT BEING CHANGED, THE ARRAY V USED TO COMPUTE
! INTEGRATION COEFFICIENTS IS RESTORED.
!
integer kord(*)
double precision f(*)
!
integer kdim, maxord, maxstf
!++ Substitute for KDIM, MAXORD, MAXSTF below
parameter (kdim = 20, maxord = 2, maxstf = 1)
!--D Next line special: P=>D, X=>Q
double precision tn
double precision xi(kdim)
!
!--D Next line special: P=>D, X=>Q
double precision tg(2), tgstop(2), tmark, tmarkx, tout, tolg
double precision alpha(kdim), beta(kdim+1)
double precision  d(maxstf+maxord,maxord), g(kdim,maxord)
double precision v(kdim+maxord)
double precision hc, hdec, hinc, hincc, hmax, hmaxp9, hmin
double precision fdat(11)
!
double precision ds(maxstf+maxord, maxord), gs(kdim)
double precision sigma(kdim), rbq(kdim), dnoise
double precision eave, eimax, eimin, emax, erep, robnd, snoise
!
integer iopst, kordi, kqmaxd, kqmaxi, ldt, maxdif, maxint, nkdko, &
   nte, nyny, ndtf, numdt
common / divasc / tn, xi, iopst, kordi, kqmaxd, kqmaxi, ldt, &
   maxdif, maxint, nkdko, nte, nyny, ndtf, numdt
!
integer icf,ics,igflg,igtype(2),igstop(2),ilgrep,ings,iop3,iop4, &
   iop5,iop6,iop7,iop8,iop9,iop10,iop11,iop12,iop13,iop14,iop15, &
   iop16,iop17,iop18,iop19,iop20,iop21,iop22,iop21s,itolep,iy, &
   kemax,kis,kmark,kord1i,kord2i,kpred,kqdcon,kqicon,kqmaxs, &
   kqmxds,kqmxil,kqmxip,kqmxis,ksc,ksout,ksstrt,kstep,lex,linc, &
   lincd,lincq,lsc,maxkqd,maxkqi,method,ne,neptol,ng,ngtot, &
   noiseq,noutko,ntolf,ny,idat(6)
common /divamc/ tg,tgstop,tmark,tmarkx,tout,tolg,hc,hdec,hinc, &
   hincc,hmax,hmaxp9,hmin,alpha,beta,d,g,v,ds,gs,sigma,rbq,dnoise, &
   eave,eimax,eimin,emax,erep,robnd,snoise,fdat,icf,ics,igflg, &
   igtype,igstop,ilgrep,ings,iop3,iop4,iop5,iop6,iop7,iop8,iop9, &
   iop10,iop11,iop12,iop13,iop14,iop15,iop16,iop17,iop18,iop19, &
   iop20,iop21,iop22,iop21s,itolep,iy,kemax,kis,kmark,kord1i, &
   kord2i,kpred,kqdcon,kqicon,kqmaxs,kqmxds,kqmxil,kqmxip,kqmxis, &
   ksc,ksout,ksstrt,kstep,lex,linc,lincd,lincq,lsc,maxkqd,maxkqi, &
   method,ne,neptol,ng,ngtot,noiseq,noutko,ntolf,ny,idat
save / divamc / , / divasc /
!
integer i, l, kqq, j, k
double precision tpd, c0, c2
parameter (c0 = 0.d0)
parameter (c2 = 2.d0)
! ********* START OF EXECUTABLE CODE **********
!
! ********
! BEGIN LOOP TO BACK UP DIFFERENCE TABLES
! ********
l = ndtf - 1
do 2410 i = 1, nte
   kqq = kord(i + 3)
!++  Code for STIFF is inactive
!         IF (KQQ) 2302,2400,2310
!c.           EQUATION IS STIFF
! 2302    IF (LINC.GE.0) GO TO 2310
!         IF (F(L+1+I)) 2306,2308,2304
!c.     ORDER WAS INCREASED, AND THUS MUST BE DECREASED (KQQ.LT.0)
! 2304    KQQ=KQQ+1
!         KORD(I+3) = KQQ
!         GO TO 2308
!c.     ORDER WAS DECREASED
! 2306    KQQ=KQQ-1
! 2308    KQQ=max(2,-KQQ)
!         GO TO 2350
!++  End
!     EQUATION IS NOT STIFF
2310    if (kqq .gt. 2) then
      if (f(l + kqq) .eq. c0) then
!                 ORDER WAS INCREASED, AND THUS MUST BE DECREASED
         kqq = kqq - 1
         kord(i + 3) = kqq
      end if
   end if
   j = min(kqq, ksc)
   kqmaxi = max(kqmaxi, kqq)
   if (kqq .ne. 1) f(l + kqq + 1) = 0.d0
!           BACK UP FOR BACKWARD DIFFERENCES
   do 2360 k = 1, j
      f(l + k) = f(l + k) - f(l + k + 1)
2360    continue
   if (kqq .gt. ksc) then
!           BACK UP FOR MODIFIED DIVIDED DIFFERENCES
      do 2390 k = j+1, kqq
         f(l + k) = (f(l+k) - f(l+k+1)) / beta(k)
2390       continue
   end if
2400    f(l + kqq + 1) = f(l + kqq + 1) / beta(kqq + 1)
   l = l + numdt
2410 continue
! END OF LOOP TO BACK UP DIFFERENCE TABLES
! ********
! BACK UP XI TO BEGINNING OF THE STEP
! ********
i = ksc + 1
if (i - iop11 - 1) 2420, 2440, 2450
2420 tpd = xi(1)
!                Check below needed when starting?
if (tpd .eq. xi(2)) go to 2450
do 2430 k = i, iop11
2430    xi(k - 1) = xi(k) - tpd
2440 xi(iop11) = c2 * xi(iop11 - 1)
if (iop11 .ne. 2) xi(iop11) = xi(iop11) - xi(iop11 - 2)
2450 kqicon = -1
icf = ne
ics = 1
ldt = 1
return
end
!   End of DIVABU

subroutine divaco(id, rd)
!>> 1987-12-07 DIVACO Krogh   Initial code.
!
! THIS SUBROUTINE RETURNS THE FOLLOWING DATA FROM COMMON
! ID(1) = KEMAX  =  INDEX OF EQUATION WITH LARGEST ERROR ESTIMATE
! ID(2) = KSTEP  =  CURRENT STEP NUMBER
! ID(3) = NUMDT  =  NUMBER OF DIFFERENCES USED FOR EACH EQUATION
! ID(4) =           RESERVED FOR FUTURE USE
! ID(5) =           RESERVED FOR FUTURE USE
! RD(1) = EMAX   =  MAX. RATIO OF ESTIMATED ERROR TO REQUESTED ERROR
! RD(2) =           RESERVED FOR FUTURE USE
! RD(3) =           RESERVED FOR FUTURE USE
!
integer id(5)
double precision rd(3)
!
integer kdim, maxord, maxstf
!++ Substitute for KDIM, MAXORD, MAXSTF below
parameter (kdim = 20, maxord = 2, maxstf = 1)
!--D Next line special: P=>D, X=>Q
double precision tn
double precision xi(kdim)
!
!--D Next line special: P=>D, X=>Q
double precision tg(2), tgstop(2), tmark, tmarkx, tout, tolg
double precision alpha(kdim), beta(kdim+1)
double precision  d(maxstf+maxord,maxord), g(kdim,maxord)
double precision v(kdim+maxord)
double precision hc, hdec, hinc, hincc, hmax, hmaxp9, hmin
double precision fdat(11)
!
double precision ds(maxstf+maxord, maxord), gs(kdim)
double precision sigma(kdim), rbq(kdim), dnoise
double precision eave, eimax, eimin, emax, erep, robnd, snoise
!
integer iopst, kordi, kqmaxd, kqmaxi, ldt, maxdif, maxint, nkdko, &
   nte, nyny, ndtf, numdt
common / divasc / tn, xi, iopst, kordi, kqmaxd, kqmaxi, ldt, &
   maxdif, maxint, nkdko, nte, nyny, ndtf, numdt
!
integer icf,ics,igflg,igtype(2),igstop(2),ilgrep,ings,iop3,iop4, &
   iop5,iop6,iop7,iop8,iop9,iop10,iop11,iop12,iop13,iop14,iop15, &
   iop16,iop17,iop18,iop19,iop20,iop21,iop22,iop21s,itolep,iy, &
   kemax,kis,kmark,kord1i,kord2i,kpred,kqdcon,kqicon,kqmaxs, &
   kqmxds,kqmxil,kqmxip,kqmxis,ksc,ksout,ksstrt,kstep,lex,linc, &
   lincd,lincq,lsc,maxkqd,maxkqi,method,ne,neptol,ng,ngtot, &
   noiseq,noutko,ntolf,ny,idat(6)
common /divamc/ tg,tgstop,tmark,tmarkx,tout,tolg,hc,hdec,hinc, &
   hincc,hmax,hmaxp9,hmin,alpha,beta,d,g,v,ds,gs,sigma,rbq,dnoise, &
   eave,eimax,eimin,emax,erep,robnd,snoise,fdat,icf,ics,igflg, &
   igtype,igstop,ilgrep,ings,iop3,iop4,iop5,iop6,iop7,iop8,iop9, &
   iop10,iop11,iop12,iop13,iop14,iop15,iop16,iop17,iop18,iop19, &
   iop20,iop21,iop22,iop21s,itolep,iy,kemax,kis,kmark,kord1i, &
   kord2i,kpred,kqdcon,kqicon,kqmaxs,kqmxds,kqmxil,kqmxip,kqmxis, &
   ksc,ksout,ksstrt,kstep,lex,linc,lincd,lincq,lsc,maxkqd,maxkqi, &
   method,ne,neptol,ng,ngtot,noiseq,noutko,ntolf,ny,idat
save / divamc / , / divasc /
!
id(1) = kemax
id(2) = kstep
id(3) = numdt
rd(1) = emax
return
end
!   End of DIVACO

subroutine divacr(y, f, kord, tol, lgroup)
!>> 1988-08-25 DIVACR Krogh   Fix bug in relative error test.
!>> 1988-01-15 DIVACR Krogh   Initial code.
!
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
integer lgroup(*), kord(*)
!--D Next line special: P=>D, X=>Q
double precision y(*)
double precision tol(*), f(*)
!
integer kdim, maxord, maxstf
!++ Substitute for KDIM, MAXORD, MAXSTF below
parameter (kdim = 20, maxord = 2, maxstf = 1)
!--D Next line special: P=>D, X=>Q
double precision tn
double precision xi(kdim)
!
!--D Next line special: P=>D, X=>Q
double precision tg(2), tgstop(2), tmark, tmarkx, tout, tolg
double precision alpha(kdim), beta(kdim+1)
double precision  d(maxstf+maxord,maxord), g(kdim,maxord)
double precision v(kdim+maxord)
double precision hc, hdec, hinc, hincc, hmax, hmaxp9, hmin
double precision fdat(11)
!
double precision ds(maxstf+maxord, maxord), gs(kdim)
double precision sigma(kdim), rbq(kdim), dnoise
double precision eave, eimax, eimin, emax, erep, robnd, snoise
!
integer iopst, kordi, kqmaxd, kqmaxi, ldt, maxdif, maxint, nkdko, &
   nte, nyny, ndtf, numdt
common / divasc / tn, xi, iopst, kordi, kqmaxd, kqmaxi, ldt, &
   maxdif, maxint, nkdko, nte, nyny, ndtf, numdt
!
integer icf,ics,igflg,igtype(2),igstop(2),ilgrep,ings,iop3,iop4, &
   iop5,iop6,iop7,iop8,iop9,iop10,iop11,iop12,iop13,iop14,iop15, &
   iop16,iop17,iop18,iop19,iop20,iop21,iop22,iop21s,itolep,iy, &
   kemax,kis,kmark,kord1i,kord2i,kpred,kqdcon,kqicon,kqmaxs, &
   kqmxds,kqmxil,kqmxip,kqmxis,ksc,ksout,ksstrt,kstep,lex,linc, &
   lincd,lincq,lsc,maxkqd,maxkqi,method,ne,neptol,ng,ngtot, &
   noiseq,noutko,ntolf,ny,idat(6)
common /divamc/ tg,tgstop,tmark,tmarkx,tout,tolg,hc,hdec,hinc, &
   hincc,hmax,hmaxp9,hmin,alpha,beta,d,g,v,ds,gs,sigma,rbq,dnoise, &
   eave,eimax,eimin,emax,erep,robnd,snoise,fdat,icf,ics,igflg, &
   igtype,igstop,ilgrep,ings,iop3,iop4,iop5,iop6,iop7,iop8,iop9, &
   iop10,iop11,iop12,iop13,iop14,iop15,iop16,iop17,iop18,iop19, &
   iop20,iop21,iop22,iop21s,itolep,iy,kemax,kis,kmark,kord1i, &
   kord2i,kpred,kqdcon,kqicon,kqmaxs,kqmxds,kqmxil,kqmxip,kqmxis, &
   ksc,ksout,ksstrt,kstep,lex,linc,lincd,lincq,lsc,maxkqd,maxkqi, &
   method,ne,neptol,ng,ngtot,noiseq,noutko,ntolf,ny,idat
save / divamc / , / divasc /
!.    SPECIFICATION OF ENVIRONMENTAL CONSTANTS.
double precision eeps10, eeps16, erov10, eeps2, eept75, eovep2
double precision ovtm75, ovd10
common / divaev / eeps2, eept75, eovep2, ovtm75, ovd10, eeps10, &
   eeps16, erov10
save / divaev /
!
integer l, i, kql, kqn, kqd, jlgrep, j, k, ilgror, itolor, jlgror, &
   iord, koutko, kqlord, ll, lkqmax
double precision cm8, cm2, cmp5, c0, cq3125, cp1, cp125, cp25, cp5
double precision cp75, cp8, cp9375, c1, c1p4, c2, c4, c10, c20
double precision c1000, c40
parameter (cm8 = -8.d0)
parameter (cm2 = -2.d0)
parameter (cmp5 = -.5d0)
parameter (c0 = 0.d0)
parameter (cq3125 = .03125d0)
parameter (cp1 = .1d0)
parameter (cp125 = .125d0)
parameter (cp25 = .25d0)
parameter (cp5 = .5d0)
parameter (cp75 = .75d0)
parameter (cp8 = .8d0)
parameter (cp9375 = .9375d0)
parameter (c1 = 1.d0)
parameter (c1p4 = 1.4d0)
parameter (c2 = 2.d0)
parameter (c4 = 4.d0)
parameter (c10 = 10.d0)
parameter (c20 = 20.d0)
parameter (c40 = 40.d0)
parameter (c1000 = 1000.d0)
double precision tpp, hh, e, ei, eps, ercoef, rnd, rnoise, s
double precision tp2, tps1, tps2, tps3, tps4, tps5, tps6, tps7
double precision ref(4)
double precision eibnd(kdim-1)
!++  Code for INTEGO is active
double precision tempa(4), tempao(4)
!++  End
save koutko, lkqmax
equivalence (tps1,tempa(1)), (tps2,tempa(2)), (tps3,tempa(3)), &
   (tps4, tempa(4))
equivalence (g(1, 1), hh)
integer mact1(2), mact2(12)
!             Parameters for Interface to MESS and DMESS
integer meret, metext, metabl
parameter (meret  =51)
parameter (metext =53)
parameter (metabl =55)
! ********* Error message text ***************
![Last 2 letters of Param. name]  [Text generating message.]
!AA KSTEP=$(I6) T=$(E15.8) H=$(E12.5) LSC=$(I3) $C
!   EIMIN=$(E8.2) EAVE=$G KSC=$(I2) SIGMA($J)=$G $C
!   RQ=$(E11.5)$G$E
!   $
!AB I$HKQ$HLI$HE$HEI$HEPS$HF$H$H$H$H
!   HIGH ORDER PREDICTED DIFFERENCES$HRNOISE$HSTIFF$HBETA$E
integer ltxtaa,ltxtab
parameter (ltxtaa=  1,ltxtab=  1)
character mtxtaa(1) * (104)
character mtxtab(1) * (88)
data mtxtaa/'KSTEP=$(I6) T=$(E15.8) H=$(E12.5) LSC=$(I3) EIMIN=$(E &
8.2) eave=$g ksc=$(i2) sigma($j)=$g rq=$(e11.5)$g$e'/
data mtxtab/'I$HKQ$HLI$HE$HEI$HEPS$HF$H$H$H$HHIGH ORDER PREDICTED$ &
 differences$hrnoise$hstiff$hbeta$e'/
!
data mact1 / metext, meret /
! (rr=repeat, t=3/5 for I/E format)  wwddtrr  wwddtrr  wwddtrr
data mact2 / metabl, 1, 0, 14, 0400201, 0300202, 0801503, &
    1507501, 1103504, 0901501, 1002501, 1205501 /
!         wwddtrr  wwddtrr  wwddtrr  wwddtrr  wwddtrr
!          End of stuff for interface to message processor
!
data ref(1), ref(2), ref(3), ref(4) / c1, cp9375, cp75, cp5 /
!++ Save data by elements if ~.C.
!++ Of next 20 lines, only the first KDIM-1 are active
data eibnd(1) / .1d0 /
data eibnd(2) / .1d0 /
data eibnd(3) / .14d0 /
data eibnd(4) / .19d0 /
data eibnd(5) / .26d0 /
data eibnd(6) / .36d0 /
data eibnd(7) / .50d0 /
data eibnd(8) / .69d0 /
data eibnd(9) / .94d0 /
data eibnd(10) / c1 /
data eibnd(11) / c1 /
data eibnd(12) / c1 /
data eibnd(13) / c1 /
data eibnd(14) / c1 /
data eibnd(15) / c1 /
data eibnd(16) / c1 /
data eibnd(17) / c1 /
data eibnd(18) / c1 /
data eibnd(19) / c1 /
!     data EIBND(20) / C1 /
!
!++  Code for ARGM is inactive
!      RETURN
!      ENTRY DIVACE
!++  End
! ********
! START OF CODE
! ********
l = ndtf - 1
if (ics .ne. 1) l = l + (ics - 1) * numdt
do 3340 i = ics, icf
   if (nkdko .ne. 0) kordi = kord(nkdko + i - 1)
   iy = iy + abs(kordi)
   kql = kord(i + 3)
   kqn = abs(kql)
   kqd = max(2, kqn)
! ********
! OBTAIN ERROR TOLERANCE SPECIFIED BY THE USER
! ********
   if (i .le. ilgrep) if (kql) 2600, 3310, 2610
   itolep = abs(itolep) + 1
   eps = tol(itolep)
   ilgrep = lgroup(itolep)
!   TEST IF SIMPLE ABSOLUTE ERROR TEST IS BEING USED
   if (ilgrep .gt. 0) go to 2580
   jlgrep = ilgrep
!     GET OLD RELATIVE ERROR FACTOR
   tps6 = tol(itolep + 1)
   ilgrep = lgroup(itolep + 1)
   itolep = -itolep - 1
!
   if (jlgrep + 1) 2540, 2570, 2510
!   NO CHECK ON THE ERROR ESTIMATE IS TO BE MADE
2510    if (eps + c1) 2520, 2590, 2520
!   ERROR TOLERANCE IS SPECIFIED IMPROPERLY
2520    kemax = i
   neptol = itolep
   linc = -3
   erep = -abs(erep)
   return
!   COMPUTE NEW RELATIVE ERROR FACTOR
2540    continue
   tps1 = c0
   do 2550 j = i, ilgrep
      tps1 = tps1 + abs(f(j))
2550       continue
   tps1 = abs(hh) * tps1 / dble(ilgrep - i + 1)
   if (lsc .le. 2) go to 2560
!     ON FIRST 3 STEPS INCREASE TPS6 WHEN COMPUTING REL. ERROR FACTOR
   tps6 = max(c4 * tps1, tps6)
!     ON 1-ST TIME THROUGH THE FIRST STEP, REL. ERR. FAC. IS NOT STORED
   if (lsc .eq. 7) go to 2570
2560    continue
   tps6 = max(tps1, tps6)
!   STORE NEW RELATIVE ERROR FACTOR
   tol(-itolep) = tps6 * ref(-jlgrep - 1)
!   COMPUTE ABSOLUTE ERROR TOLERANCE
2570    eps = eps * tps6
2580    if (eps .le. c0) go to 2520
2590    if (kql) 2600, 3330, 2610
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
!      IF (KQD.EQ.2) TPS1=Y(IY-1)
!      E=ABS(TPS3)+ABS(TPS4)
!      EI=E+ABS(TPS2)
!      RND=EI
!      IF (KORDI.GE.0) GO TO 2604
!c.    EQUATION IS IMPLICIT
!      JSI=JSI-1
!      IF (JSI.NE.0) GO TO 2604
!      IF (KORDI.EQ.-1) GO TO 2602
!      ERCOEF=GS(KQN+1)
!      GO TO 2606
! 2602 ERCOEF=.5D0*DS(KQD,1)
!      JSI=1
!      GO TO 2606
!c.    END OF SPECIAL CODE FOR IMPLICIT EQUATIONS
! 2604 ERCOEF = DS(KQD,JSI)
! 2606 ERCOEF = ABS(ERCOEF) / EPS
!      IF (LSC.LE.2)  GO TO 2710
!      IF (LSC-5) 2650,2710,2710
!c.  END OF CODE FOR STIFF EQUATIONS
!++  End
!
! EQUATION IS NOT STIFF
2610    tpp = f(i) - f(l + 1)
   tps3 = tpp
   tps4 = tpp - f(l + kqd + 1)
   tps2 = tpp + f(l + kqd)
   tps1 = tpp + f(l + kqd - 1)
   e = abs(tps3) + abs(tps4)
   rnd = e
   ei = e + abs(tps2)
   ercoef = abs(gs(kqn + 1)) / eps
   if (kql .ge. 4) go to 2710
!   TEST IF STARTING OR IF INTEGRATION ORDER IS ONE
   if (lsc .le. 2) if (kql - 2) 2660, 2710, 2710
! ********
! LOGIC ASSOCIATED WITH STARTING THE INTEGRATION
! ********
   tps4 = c0
   if (lsc - 4) 2650, 2640, 2620
! FIRST STEP
2620    e = e * cq3125
   tps3 = c0
   f(l + 4) = c0
   s = c0
!   TEST IF FIRST TIME THROUGH THE FIRST STEP
   if (lsc .eq. 7) go to 2690
!   COMPUTE S=ESTIMATE OF H * EIGENVALUE OF JACOBIAN = 2*(F(A)-F(B))/
!   (F(B)-F(C)) WHERE F(A)=CURRENT F(I), AND F(B) AND F(C) PRECEDING
!   VALUES OR ESTIMATES OF F(I)
   tpp = f(i) - f(l + 5)
   tps4 = tpp
   e = c2 * abs(tps4)
   if (s .ne. c0) s = (tps4 + tps4) / s
   if (s + cp125) 2630, 2700, 2700
!     SET LDT=-5  TO INDICATE POSSIBLE PROBLEMS DUE TO INSTABILITY
2630    ldt = -5
   go to 2690
!   ADJUST CORRECTION MADE ON SECOND STEP
2640    tpp = cp8 * tpp
!   ADJUST ESTIMATED ERRORS ON SECOND AND THIRD STEPS
2650    e = abs(tps3)
   rnd = c4 * e
   go to 2710
! END OF SPECIAL LOGIC FOR STARTING
! ********
! INTEGRATION ORDER =1 IS TREATED AS A SPECIAL CASE
! ********
2660    tpp = tpp + f(l + 2)
   if (beta(2) .ge. c1p4) ei = ei * c1000
!   ESTIMATE NEW VALUE FOR S
   s = f(l + 4)
   if (s .eq. c0) go to 2680
   s = max(cm8, c2 * beta(2) * (tps1 - tps2 - f(l + 5)) / s)
   if (s .ge. cmp5) go to 2670
!   MODIFY TPP (TO GET BETTER STABILITY CHARACTERISTICS)
   tpp = tpp * max(cp25, (cm2 - c2 * s) / (s * s))
2670    tps4 = tps4 * abs(s)
2680    e = cp25 * (e + abs(tps4))
   ei = ei + abs(tps4 * s)
!     STORE INFORMATION REQUIRED TO ESTIMATE S ON NEXT STEP
2690    f(l + 4) = tpp
2700    f(l + 5) = f(i)
! END OF SPECIAL CODE FOR INTEGRATION ORDER =1
! ********
! CODE FOR NOISE TEST AND GETTING ERROR ESTIMATE
! ********
2710    e = e * ercoef
   rnoise = c0
   if (eps .lt. c0) go to 2810
   tps5 = abs(f(l + 2)) + abs(f(i))
   if (tps5 .eq. c0) go to 2760
2720    rnoise = rnd / tps5
   if (rnoise .gt. rbq(kqd)) if (rnoise - c1) 2760, 2750, 2750
!   NOISE IS APPARENTLY SLOWING CONVERGENCE OF THE DIFFERENCES
!     REDUCE EI
   ei = rnd
   tps5 = abs(eeps2 * y(iy - 1)) / eps
   if (tps5 .lt. abs(e)) if (lsc) 2730, 2730, 2760
   e = tps5
   rnoise = c0
2730    e = -abs(e)
   if (eimin .gt. cp1) ei = (c10 * eimin) * ei
!     COMPUTE REDUCTION TO BE MADE IN EI
   if (rnoise .gt. (c20 * rbq(kqd))) go to 2760
   k = -6 - lsc
2740    if (k .le. 0) go to 2760
!     REDUCE EI WHEN NOISE APPARENTLY LIMITS PRECISION
   k = k - 1
   ei = cp5 * ei
   if (ei .gt. eimin) go to 2740
   go to 2760
2750    tps4 = 1.1d0 * rnd
   tps3 = rnd
2760    continue
!   TEST FOR STIFFNESS GOES HERE WHEN IMPLEMENTED
! *       INGREDIENTS OF TEST MAY INCLUDE --
! *       RNOISE, WHETHER (ABS(TPS4).GT.ABS(TPS3)),
! *       WHETHER EMAX IS INCREASING, RESULT OF TEST ON
! *       PREVIOUS STEPS, ETC.
!
! ********
! COMPUTE ERROR ESTIMATES AND INFORMATION FOR SELECTING THE STEPSIZE
! ********
   if (e .ge. abs(emax)) go to 2770
   if (-e .le. abs(emax)) go to 2780
   snoise = rnoise
   dnoise = rnd
   noiseq = kqd
!   STORE PARAMETERS ASSOCIATED WITH LARGEST VALUE OF E
2770    emax = e
   kemax = i
   neptol = itolep
!   DETERMINE HOW MUCH STEPSIZE CAN BE INCREASED
2780    ei = ei * ercoef * sigma(kqd)
   eimax = max(eimax, ei)
   if (linc .le. 0) go to 2810
   k = 0
2790    if (ei .ge. min(eimin, eibnd(kqn))) go to 2800
   k = k + 1
   if (k .eq. linc) go to 2810
   ei = ei * sigma(kqd)
   go to 2790
2800    linc = k
! END OF COMPUTING ERROR ESTIMATES
2810    continue
!++  Code for ERRSTO is inactive
!      IF (IOP20 .EQ. 0) GO TO 780
!c.********
!c.STORE ERROR ESTIMATE (OPTIONAL)
!c.********
!      F(IOP20+I-1)=TPS3*GS(KQN+1)
!c.END OF STORING ERROR ESTIMATE
!++  Code for INTEGO | ERRSTO is active
   if (iop19 .eq. 0) go to 3090
!.********
!.EQUATIONS ARE GROUPED TO USE SAME INTEGRATION METHOD (OPTIONAL)
!.********
!++  Code for INTEGO is active
   if (i .gt. 1) if (i - ilgror) 2900, 2900, 2830
   itolor = iop19
2830    jlgror = kord(itolor)
   itolor = itolor + 1
   if (jlgror .gt. 0) go to 2870
   ilgror = kord(itolor)
   itolor = itolor + 1
   if (jlgror + 1) 2840, 2850, 2890
2840    if (jlgror .lt. -2) if (kqd + jlgror) 2850, 2880, 2880
!.INITIALIZE FOR ACCUMULATING VARIABLES USED IN ORDER SELECTION
2850    iord = i
   kqlord = kql
   do 2860 k = 1, 4
2860       tempao(k) = abs(tempa(k))
   go to 2930
!.ORDERS IN CURRENT GROUP CAN BE DIFFERENT
2870    ilgror = jlgror
   go to 3090
!.ORDER IS NOT GOING TO BE CHANGED
2880    jlgror = 0
2890    if (kql) 3240, 3270, 3270
!.TAKE ACTION FOR EQUATION WHICH IS NOT THE FIRST IN THE GROUP
2900    if (jlgror) 2910, 2890, 3090
!.ACCUMULATE VARIABLES USED IN ORDER SELECTION
2910    do 2920 k = 1, 4
2920       tempao(k) = tempao(k) + abs(tempa(k))
!.    TEST IF THIS IS LAST EQUATION IN THE GROUP
2930    if (i .ne. ilgror) if (kql) 3310, 3290, 3290
!.SET UP TO GO SELECT INTEGRATION ORDER
   kql = 0
   do 2940 k = 1, 4
2940       tempa(k) = tempao(k)
   go to 3090
!.INTEGRATION ORDER HAS BEEN SELECTED
!++  Code for INTEGO | STIFF is active
2950    continue
!++  Code for INTEGO is active
   kql = kqlord
   if (kqn - abs(kql)) 2960, 2980, 3020
!.  TEST IF ORDER CAN BE DECREASED
2960    if (jlgror .ge. -2) if (kql) 3010, 3040, 3040
!.    INTEGRATION ORDER WAS SELECTED OUTSIDE PERMITTED RANGE
2970    kqn = abs(kql)
!.    INTEGRATION ORDER IS NOT GOING TO BE CHANGED
2980    if ((kql .ne. 1) .or. (lsc .gt. 0)) if (kql) 3030, 3040, 3040
!.    SET  4-TH ENTRY IN DIFFERENCE TABLES SO THAT STANDARD ADAMS
!.    METHOD IS USED WHEN KQL=1
2990    do 3000 k = iord, i
3000       f(ndtf + k*numdt - numdt + 3) = c0
   go to 3270
!.  ORDER FOR STIFF EQUATION WAS REDUCED
3010    continue
!++  Code for INTEGO & STIFF is inactive
!      IF (KQN.LT.JSI) GO TO 990
!      TPP=-C1
!      GO TO 1090
!c.  TEST IF ORDER CAN BE INCREASED
!++  Code for INTEGO is active
3020    if (jlgror .eq. -2) go to 2970
!++  Code for INTEGO & STIFF is inactive
!      IF (KQL.GE.0) GO TO 1140
!      IF ((JSI.NE.0).AND.(KQN.GT.(MAXKQD+JSI))) GO TO 990
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
3040    ll = ndtf + numdt * iord - numdt
   do 3080 j = iord, i
      kord(j + 3) = kqn
      if (kqn - kql) 3050, 3070, 3060
3050       f(ll + kqd - 1) = f(ll + kqd - 1) + (f(j) - f(ll))
      go to 3080
3060       f(ll + kqn) = f(ll + kqd)
3070       f(ll + kqd) = c0
3080       ll = ll + numdt
   if (kqn - 1) 3270, 2990, 3270
!++  End
!.********
!.SELECT INTEGRATION ORDER
!.********
3090    if (lsc .le. 0) go to 3120
!. SPECIAL ORDER SELECTION WHEN STARTING
   if (lsc - 3) 3110, 3210, 3100
3100    if (lsc .eq. 5) if (s + .125d0) 3160, 3130, 3130
   if (lsc - 6) 3130, 3130, 3210
3110    if (c40 * min(abs(tps4), abs(tps3)) .gt. abs(tps2)) then
      if (eps .ne. -c1) lsc = 2
   end if
   if (abs(tps4) .lt. abs(tps3)) if (c4 * abs(tps4) - abs(tps2)) &
      3130, 3130, 3210
!.  CHECK IF ORDER CAN BE INCREASED OR SHOULD BE DECREASED
3120    tps5 = robnd * abs(tps4)
   tps6 = robnd * (tps5 + abs(tps3))
   tps7 = abs(tps1) + abs(tps2)
   if (tps5 .ge. abs(tps3)) go to 3140
   if (tps6 .ge. tps7) go to 3210
3130    if (kqn .ge. maxkqi) go to 3210
!.    INCREASE THE INTEGRATION ORDER
   kqn = kqn + 1
!++  Code for INTEGO | STIFF is active
   if (kql) 3230, 2950, 3250
!++  Code for ~(INTEGO | STIFF) is inactive
!      GO TO 3250
!++  End
!.  CHECK IF ORDER SHOULD BE DECREASED
3140    if (tps6 .lt. tps7) go to 3210
   if (tps5 .lt. abs(tps3 - tps4)) go to 3210
   if ((tps3.eq.tps4) .and. (lsc.le.0)) go to 3210
   if (kqn - 2) 3210, 3160, 3180
3160    kqn = 1
!++  Code for INTEGO | STIFF is active
   if (kql) 3220, 2950, 3170
!++  End
!.    WHEN ORDER IS REDUCED TO 1 WITH ADAMS METHOD SET F(L+4)=0
3170    f(l + 4) = c0
   go to 3260
!.    DECREASE THE INTEGRATION ORDER
3180    kqn = kqn - 1
!++  Code for INTEGO | STIFF is active
   if (kql) 3220, 2950, 3200
!++  End
3200    f(l+kqd) = f(l+kqd) + tpp
   go to 3260
!   NO CHANGE IN INTEGRATION ORDER IS BEING MADE
3210    continue
!++  Code for INTEGO is active
   if (kql) 3240, 2950, 3270
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
!      IF (KQN.LT.JSI) GO TO 3236
!      F(L+1)=-C1
!      GO TO 3233
!c.    ORDER WAS INCREASED
!++  Code for INTEGO |  STIFF  is active
3230    continue
!++  Code for STIFF is inactive
!      IF ((JSI.NE.0).AND.(KQN.GT.(MAXKQD+JSI))) GO TO 3236
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
! 3245 IF (JSI.NE.0) KQMAXD=max(KQN,KQMAXD)
!      IF (JS.LT.abs(KORDI)) KQMAXI=max(KQN,KQMAXI)
!      GO TO 3290
!++  End
! EQUATION IS NOT STIFF
!     ORDER INCREASED
3250    f(l + kqn + 1) = -f(l + kqd + 1)
   if (lsc .gt. 0) f(l + kqn + 1) = f(l + 1) - f(i)
!     ORDER CHANGED
3260    kord(i + 3) = kqn
3270    kqmaxi = max(kqn, kqmaxi)
   if (eps .gt. c0) kqmaxs = max(kqn, kqmaxs)
   f(l + kqd + 1) = c0
3290    continue
   if (kqn .gt. kis) go to 3310
!.********
!.DETERMINE IF TIME TO STORE SOLUTION (OPTIONAL)
!.********
   if (kis .ge. 1000) then
      tp2 = max(1.5d0, dble(kqn) * c2 ** (1001 - kis)) * abs(tps4)
3295       if (tp2 .gt. abs(f(l+kqn))) then
         if (kqn .le. kql) then
            kqn = kqn - 1
            if (kqn .gt. 1) go to 3295
            kqn = 1
         end if
      end if
      kord(i+3) = kqn
      if (i .eq. 1) lkqmax = 0
      lkqmax = max(kqn, lkqmax)
      kqmaxi = lkqmax
      if (kis .eq. 1000) then
         if (i .eq. kemax) emax = dble(8 + kqn**2) * abs(emax)
         go to 3325
      end if
!++  Code for DUMP is active
   else if ((e .ne. c0) .and. (eps .gt. c0)) then
      if (iop9 .gt. 0) if((abs(e)*dble(kis-kqn+2)**(kqn+1))-1.d-2) &
         3310, 3310, 3300
3300       kis = -1
!++  End
   end if
3310    continue
! ********
! CORRECT
! ********
   do 3320 k = 1, kordi
!++  Code for ~{p,x} is active
      y(iy - k) = y(iy - k) + g(kql + 1, k) * tpp
!++  Code for {p,x} is inactive
!c--D Next line special: P=>D, X=>Q
!            Y(IY - K) = Y(IY - K) + dble(G(KQL + 1, K)) * dble(TPP)
!++  END
3320    continue
! END OF CORRECTING
3325 continue
!++  Code for OUTPUT is active
if (iop10 .gt. 0) then
   if (i .eq. 1) then
      idat(1) = kstep
      idat(2) = lsc
      idat(3) = ksc
      idat(4) = iop11
      fdat(1) = tn
      fdat(2) = hh
      fdat(3) = eimin
      fdat(4) = eave
      fdat(5) = sigma(iop11)
      fdat(6) = robnd
      mact2(3) = nte
!--D Next line special: P=>S, X=>D
      call dmess(mact1, mtxtaa, idat, fdat)
      koutko = noutko
   end if
   if (koutko .ne. 0) then
      if (kord(koutko) .gt. 0) then
         if (i .lt. kord(koutko)) go to 3328
         koutko = koutko + 1
      else
         if (i .ge. abs(kord(koutko))) koutko = koutko + 1
      end if
   end if
   idat(1) = i
   idat(2) = kql
   idat(3) = linc
   fdat(1) = e
   fdat(2) = ei
   fdat(3) = eps
   fdat(4) = f(i)
   fdat(5) = tps1
   fdat(6) = tps2
   fdat(7) = tps3
   fdat(8) = tps4
   fdat(9) = rnoise
   fdat(10) = 0.d0
   if (kql .eq. 1) fdat(10) = s
   fdat(11) = beta(kqd)
!--D Next line special: P=>S, X=>D
   call dmess(mact2, mtxtab, idat, fdat)
3328    if (i .eq. nte) iop10 = iop10 - 1
end if
!++  End
3330    l = l + numdt
3340    continue
return
end
!   End of DIVACR

subroutine divahc
!>> 1988-05-20 DIVAHC Krogh   Initial code.
!
! SUBROUTINE TO COMPUTE COEFFICIENTS REQUIRED FOR INTEGRATING
! ORDINARY DIFFERENTIAL EQUATIONS
!
integer kdim, maxord, maxstf
!++ Substitute for KDIM, MAXORD, MAXSTF below
parameter (kdim = 20, maxord = 2, maxstf = 1)
!--D Next line special: P=>D, X=>Q
double precision tn
double precision xi(kdim)
!
!--D Next line special: P=>D, X=>Q
double precision tg(2), tgstop(2), tmark, tmarkx, tout, tolg
double precision alpha(kdim), beta(kdim+1)
double precision  d(maxstf+maxord,maxord), g(kdim,maxord)
double precision v(kdim+maxord)
double precision hc, hdec, hinc, hincc, hmax, hmaxp9, hmin
double precision fdat(11)
!
double precision ds(maxstf+maxord, maxord), gs(kdim)
double precision sigma(kdim), rbq(kdim), dnoise
double precision eave, eimax, eimin, emax, erep, robnd, snoise
!
integer iopst, kordi, kqmaxd, kqmaxi, ldt, maxdif, maxint, nkdko, &
   nte, nyny, ndtf, numdt
common / divasc / tn, xi, iopst, kordi, kqmaxd, kqmaxi, ldt, &
   maxdif, maxint, nkdko, nte, nyny, ndtf, numdt
!
integer icf,ics,igflg,igtype(2),igstop(2),ilgrep,ings,iop3,iop4, &
   iop5,iop6,iop7,iop8,iop9,iop10,iop11,iop12,iop13,iop14,iop15, &
   iop16,iop17,iop18,iop19,iop20,iop21,iop22,iop21s,itolep,iy, &
   kemax,kis,kmark,kord1i,kord2i,kpred,kqdcon,kqicon,kqmaxs, &
   kqmxds,kqmxil,kqmxip,kqmxis,ksc,ksout,ksstrt,kstep,lex,linc, &
   lincd,lincq,lsc,maxkqd,maxkqi,method,ne,neptol,ng,ngtot, &
   noiseq,noutko,ntolf,ny,idat(6)
common /divamc/ tg,tgstop,tmark,tmarkx,tout,tolg,hc,hdec,hinc, &
   hincc,hmax,hmaxp9,hmin,alpha,beta,d,g,v,ds,gs,sigma,rbq,dnoise, &
   eave,eimax,eimin,emax,erep,robnd,snoise,fdat,icf,ics,igflg, &
   igtype,igstop,ilgrep,ings,iop3,iop4,iop5,iop6,iop7,iop8,iop9, &
   iop10,iop11,iop12,iop13,iop14,iop15,iop16,iop17,iop18,iop19, &
   iop20,iop21,iop22,iop21s,itolep,iy,kemax,kis,kmark,kord1i, &
   kord2i,kpred,kqdcon,kqicon,kqmaxs,kqmxds,kqmxil,kqmxip,kqmxis, &
   ksc,ksout,ksstrt,kstep,lex,linc,lincd,lincq,lsc,maxkqd,maxkqi, &
   method,ne,neptol,ng,ngtot,noiseq,noutko,ntolf,ny,idat
save / divamc / , / divasc /
!.    SPECIFICATION OF ENVIRONMENTAL CONSTANTS.
double precision eeps10, eeps16, erov10, eeps2, eept75, eovep2
double precision ovtm75, ovd10
common / divaev / eeps2, eept75, eovep2, ovtm75, ovd10, eeps10, &
   eeps16, erov10
save / divaev /
!                 K - 1 + 1 / K  is equivalent to max(1, K-1)
double precision gg(maxord - 1 + 1/maxord), b(kdim+maxord), &
   w(kdim+maxord)
integer  k, n, j
double precision c0, cp1, crbqi, cp5, cp5625, c1, c1p125
parameter (c0 = 0.d0)
parameter (cp1 = .1d0)
parameter (crbqi = .421875d0)
parameter (cp5 = .5d0)
parameter (cp5625 = .5625d0)
parameter (c1 = 1.d0)
parameter (c1p125 = 1.125d0)
!++  Code for STIFF is inactive
!      INTEGER          GODIF
!++  End
double precision tp1, tp2, hh, temp, tp
equivalence (g(1, 1), hh)
!
save gg, w
!
!  B(K)= 1/(K*(K+1))
!++ Save data by elements if ~.C.
!++ Of next 23 lines, only the first KDIM+MAXORD are active
data b(1)  / 5.000000000000000000000000000000000000000d-1 /
data b(2)  / 1.666666666666666666666666666666666666667d-1 /
data b(3)  / 8.333333333333333333333333333333333333333d-2 /
data b(4)  / 5.000000000000000000000000000000000000000d-2 /
data b(5)  / 3.333333333333333333333333333333333333333d-2 /
data b(6)  / 2.380952380952380952380952380952380952381d-2 /
data b(7)  / 1.785714285714285714285714285714285714286d-2 /
data b(8)  / 1.388888888888888888888888888888888888889d-2 /
data b(9)  / 1.111111111111111111111111111111111111111d-2 /
data b(10) / 9.090909090909090909090909090909090909091d-3 /
data b(11) / 7.575757575757575757575757575757575757576d-3 /
data b(12) / 6.410256410256410256410256410256410256410d-3 /
data b(13) / 5.494505494505494505494505494505494505495d-3 /
data b(14) / 4.761904761904761904761904761904761904762d-3 /
data b(15) / 4.166666666666666666666666666666666666667d-3 /
data b(16) / 3.676470588235294117647058823529411764706d-3 /
data b(17) / 3.267973856209150326797385620915032679739d-3 /
data b(18) / 2.923976608187134502923976608187134502924d-3 /
data b(19) / 2.631578947368421052631578947368421052632d-3 /
data b(20) / 2.380952380952380952380952380952380952381d-3 /
data b(21) / 2.164502164502164502164502164502164502165d-3 /
data b(22) / 1.976284584980237154150197628458498023715d-3 /
!     data B(23) / 1.811594202898550724637681159420289855072D-3 /
!
! ********
! START OF CODE
! ********
!     SET STEP NUMBER OF METHOD
!++  Code for STIFF is inactive
!      IOP11 = MIN(max(KQMAXI,KQMAXD) + 1), KDIM)
!++  Code for ~STIFF is active
iop11 = min(kqmaxi + 1, kdim)
!++  End
!     TEST IF STEPSIZE WAS CHANGED
if (kqicon .ge. 0) go to 3510
! ********
! STEPSIZE JUST CHANGED
! ********
!     SET CONSTANTS DEPENDING ON NEW STEPSIZE
kqmxil = kqmaxi
tp1 = hh
gg(1) = tp1 * tp1
g(1, 2) = gg(1) * cp5
if (maxint .le. 2) go to 3450
do 3440 k = 3, maxint
   gg(k - 1) = g(1, k - 1) * tp1
   g(1, k) = gg(k - 1) / dble(k)
3440    continue
!     SET CONSTANTS INDICATING STEP CHANGE
3450 kqicon = 0
!++  Code for STIFF is inactive
!      KQDCON=0
!++  End
kqmxip = 1
ksc = 1
if (lsc .lt. 7) go to 3490
!     SPECIAL SET-UP OF CONSTANTS ON THE VERY FIRST STEP
hincc = c1p125
lincd = 6
lincq = 12
if (hinc .gt. c0) go to 3460
lincd = -2
linc = -2
robnd = c1
3460 sigma(1) = 1.0d0
beta(1) = c1
do 3470 n = 1, iop11
!++  Code for STIFF is inactive
!      D(1,N)=C0
!++  End
   xi(n) = tp1
   alpha(n) = c1
   beta(n + 1) = c1
   sigma(n + 1) = dble(n + 1) * sigma(n) * hincc
3470    continue
temp = eeps16
rbq(1) = c1
rbq(2) = cp1
tp = crbqi
!     **** IN THE LOOP BELOW RBQ(K) IS COMPUTED TO BE
!          APPROXIMATELY (3/4 ** ((K-1) ** 2 - 1) / 10
!          .5625 = (3/4) ** 2    TP = (3/4) ** (2*K -3)
do 3480 k = 3, kdim
   temp = temp + temp
   rbq(k) = max(temp, rbq(k - 1) * tp)
3480    tp = tp * cp5625
go to 3560
!     SET-UP AFTER THE FIRST STEP
3490 tp2 = xi(1)
xi(1) = tp1
beta(2) = tp1 / tp2
k = 2
if (hincc .eq. hinc) go to 3540
if ((lsc .ne. 0) .or. ((kstep-ksstrt-kqmaxs) .lt. 10)) go to 3540
hincc = c1
lincd = 0
3500 lincd = lincd + 1
hincc = hincc * hinc
if (hincc .lt. 2.d0) go to 3500
linc = (linc * (lincd + lincd)) / lincq
lincq = lincd + lincd
hincc = hinc
go to 3540
! END OF LOGIC FOR CASE WHEN STEPSIZE JUST CHANGED
!     TEST IF MAXIMUM INTEGRATION ORDER DID NOT INCREASE
3510 if (kqmaxi .gt. kqmxil) then
! ********
! INTEGRATION ORDER WAS INCREASED -- GET NEW V'S
! ********
   kqmxil = kqmaxi
   kqmxip = kqmxil + maxint
   k = kqmxip
   v(k) = b(k)
   if (kqicon .eq. 1) go to 3530
!     if (KQICON .eq. K) KQICON = KQICON - 1 --- Removed 1999-08-19
   do 3520 n = 2, kqicon
      k = k - 1
3520       v(k) = v(k) - alpha(n) * v(k + 1)
! END OF GETTING NEW V'S
else
   iop11 = max(iop11, kqmxil+1)
end if
3530 if (iop11 .le. ksc) go to 3560
! ********
! COMPUTE PARAMETERS WHICH ARE STILL CHANGING AS A RESULT OF
! A CHANGE IN THE STEPSIZE
! ********
tp2 = xi(ksc)
!     UPDATE CONSTANT STEP COUNTER
ksc = ksc + 1
k = ksc
beta(k) = c1
3540 continue
temp = hincc
!
!   LOOP TO COMPUTE NEW VALUES OF PARAMETERS
do 3550 n = k, iop11
   tp1 = tp2 + hh
   tp2 = xi(n)
   xi(n) = tp1
   alpha(n) = hh / tp1
   beta(n + 1) = beta(n) * (tp1 / tp2)
   temp = max(temp, dble(n) * (alpha(n) * hincc))
   sigma(n) = sigma(n - 1) * temp
3550    continue
if (iop11 .ne. kdim) xi(iop11 + 1) = tp2 + hh
! END OF CODE FOR COMPUTING PARAMETERS WHICH ARE STILL CHANGING
!
3560 if (kqicon .ge. kqmxip) go to 3690
! ********
! COMPUTE INTEGRATION COEFFICIENTS WHICH ARE STILL CHANGING
! ********
kqmxil = max(kqmaxi, kqmxil)
kqmxip = kqmxil + maxint
j = kqmxip - kqicon
n = kqicon + 1
kqicon = n
if (n .ne. 1) go to 3580
! INITIALIZE V AND W
do 3570 k = 1, j
   v(k) = b(k)
3570    w(k) = v(k)
go to 3600
! UPDATE V AND INITIALIZE W
3580 if (n .eq. kdim) go to 3690
do 3590 k = 1, j
   v(k) = v(k) - alpha(n) * v(k + 1)
3590    w(k) = v(k)
! SET TRANSFER FOR LOOP BELOW DEPENDING ON VALUE OF MAXINT
3600 continue
go to 3660
!
3640 j = j - 1
! INNER LOOP FOR COMPUTING INTEGRATION COEFFICIENTS
do 3650 k = 1, j
3650    w(k) = w(k) - alpha(n) * w(k + 1)
!     STORE INTEGRATION COEFFICIENTS
3660 g(n + 1, 1) = hh * w(1)
gs(n + 1) = g(n + 1, 1) - g(n, 1)
!++  Code for MAXORD >= 2 is active
if (maxint .ge. 2) then
   g(n + 1, 2) = gg(1) * w(2)
!++  Code for MAXORD >= 3 is inactive
!        if (MAXINT .gt. 2) then
!           DO 3665 K=3,MAXINT
!3665          G(N+1,K)=GG(K-1)*W(K)
!        end if
!++  Code for MAXORD >= 2 is active
end if
!++  End
n = n + 1
if (n .le. kqmxil) go to 3640
! END OF COMPUTING INTEGRATION COEFFICIENTS
!
3690 continue
!++  Code for STIFF is inactive
!      IF (KQDCON.GT.KQMAXD) GO TO 4662
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
end
!   End of DIVAHC

subroutine divain(t, y, f, kord)
!>> 1988-01-14 DIVAIN Krogh   Initial code.
!
!  SUBROUTINE TO DO INTERPOLATION FOR VARIABLE ORDER INTEG. ROUTINE
!
integer kord(*)
!--D Next line special: P=>D, X=>Q
double precision t(*), y(*)
double precision f(*)
integer kdim, maxord
!++ Substitute for KDIM, MAXORD below
parameter (kdim = 20, maxord = 2)
!--D Next line special: P=>D, X=>Q
double precision tn
double precision xi(kdim)
!
integer iopst, kordi, kqmaxd, kqmaxi, ldt, maxdif, maxint, nkdko, &
   nte, nyny, ndtf, numdt
common / divasc / tn, xi, iopst, kordi, kqmaxd, kqmaxi, ldt, &
   maxdif, maxint, nkdko, nte, nyny, ndtf, numdt
save / divasc /
integer i, ici, idt, interp, integ, integz, iy, iyi, iyn, iyni, j, &
    k, kqmxi, kqmxs, kqq, l, n
double precision c0, c1, c2
parameter (c0 = 0.d0)
parameter (c1 = 1.d0)
parameter (c2 = 2.d0)
double precision c(kdim+maxord-1), eta(kdim)
double precision gamma(kdim)
double precision tp1, hi
double precision csum(kdim+maxord-1)
!--D Next line special: P=>D, X=>Q
double precision xp1
logical lnotm1
!
!              Stuff for processing error messages
integer idat(1)
double precision fdat(6)
integer mentxt, meret, meemes, metext
parameter (mentxt =23)
parameter (meret  =51)
parameter (meemes =52)
parameter (metext =53)
integer mact(8)
! ********* Error message text ***************
![Last 2 letters of Param. name]  [Text generating message.]
!AA DIVAIN$B
!AB Interpolating at T(1)=$F with $B
!AC TN=$F, T(2)=$F and H=$F.  T(1) must be in [$F, $F].$E
!AD internal variable LDT = $I.  Interpolation not allowed now.$E
integer ltxtaa,ltxtab,ltxtac,ltxtad
parameter (ltxtaa=  1,ltxtab=  9,ltxtac= 41,ltxtad= 94)
character mtxtaa(1) * (154)
data mtxtaa/'DIVAIN$BInterpolating at T(1)=$F with $BTN=$F, T(2)=$ &
f and h=$f.  t(1) must be in [$f, $f].$Einternal variable ldt = $i &
.  Interpolation not allowed now.$e'/
!
!                      1 2 3 4       5 6       7      8
data mact / meemes,0,0,0, mentxt,0, metext, meret /
!
!++  Code for ARGM is inactive
!      ENTRY DIVAIE
!++  End
! ********
! START OF CODE -- CHECK ON STATE OF DIFFERENCE TABLE
! ********
l = ldt
if (l) 3710, 3730, 3780
3710 if (l + 2) 4170, 3730, 3720
3720 if (maxint .ge. 0) l = 1
go to 3840
! ********
! UPDATE DIFFERENCE TABLE TO START OF NEXT STEP
! ********
3730 k = ndtf
do 3770 i = 1, nte
   kqq = kord(i + 3)
   if (kqq .le. 0) go to 3760
! EQUATION IS NOT STIFF
   tp1 = f(i) - f(k)
! LOOP TO DO UPDATING
   n = k + max(abs(kqq), 2)
   do 3750 j = k, n
3750       f(j) = f(j) + tp1
3760    continue
   k = k + numdt
3770    continue
ldt = 1
if (l .ne. 0) return
! END OF UPDATING DIFFERENCE TABLE
! ********
! INITIALIZE FOR COMPUTATION OF COEFFICIENTS
! ********
3780 interp = 0
hi = t(1) - tn
gamma(1) = hi / xi(1)
if (gamma(1)) 3790, 3800, 3810
3790 if (gamma(1) .ge. -c1) go to 3820
interp = 1
if (abs(hi) - abs(t(2))) 3820, 3820, 4180
3800 interp = 2 - kqmaxi
go to 3820
3810 if (gamma(1) .gt. c2) if (ldt - 2) 4180, 3820, 3820
3820 kqmxi = kqmaxi + interp - 1
!++  Code for STIFF is inactive
!      KQMXS=max(KQMXI,KQMAXD)
!++  Code for ~STIFF is active
kqmxs = kqmxi
!++  End
do 3830 n = 2, kqmxs
3830    gamma(n) = (hi + xi(n-1)) / xi(n)
3840 lnotm1 = l .ne. -1
integ = maxint
if (integ .le. 0) if (integ + maxdif) 4160, 3950, 3950
! ********
! COMPUTE INTEGRATION COEFFICIENTS
! ********
!     INITIAL SET-UP
!         COMPUTE INITIAL C VALUES
do 3850 n = 1, integ
   c(n) = hi / dble(n)
3850 continue
i = integ + 1
integ = integ + kqmxi
do 3860 n = i, integ
   c(n) = c(n - 1) * (dble(n - maxint) / dble(n))
3860    continue
!         COMPUTE ETA'S
do 3870 n = 1, kqmxi
3870    eta(n) = hi / xi(n)
!         COMPUTE C(K)'S TO CORRESPOND TO G(K-MAXINT+1,MAXINT),
!         K=MAXINT, MAXINT+1,..., MAXINT+KQMXI-1
i = integ
3880 j = integ
integ = j - 1
if (integ .le. maxint) go to 3900
do 3890 n = j, i
3890    c(n) = eta(n - integ) * c(n) + c(n - 1)
go to 3880
3900 do 3910 n = j, i
3910    c(n) = eta(n - integ) * c(n)
!         END OF COMPUTING  G(---,MAXINT)
integz = 0
go to 3940
!         COMPUTE C(K)-S TO CORRESPOND TO G(K-INTEG+1,INTEG),
!         K=INTEG+1,INTEG+2,..., INTEG+KQMXI
3920 do 3930 n = 1, kqmxi
3930    c(integ+n) = gamma(n)*c(integ+n-1) - eta(n)*c(integ+n)
3940 ici = integ - 1
go to 4020
! END OF COMPUTING INTEGRATION COEFFICIENTS
! ********
! COMPUTE COEFFICIENTS FOR INTERPOLATION
! ********
3950 c(1) = c1
ici = 0
do 3960 n = 1, kqmxs
3960    c(n + 1) = gamma(n) * c(n)
if (integ + 1) 3970, 3990, 4010
! END OF COMPUTING INTERPOLATION COEFFICIENTS
!
!     SET-UP TO COMPUTE DIFFERENTIATION COEFFICIENTS REQUIRED
!     IN ORDER TO GET COEFFICIENTS ACTUALLY USED
3970 integ = 0
ici = 1
3980 integ = integ - 1
if (integ .eq. maxint) ici = 0
! ********
! COMPUTE DIFFERENTIATION COEFFICIENTS
! ********
3990 interp = max(interp, 0)
tp1 = dble(-integ)
c(1) = tp1 * c(1) / xi(-integ)
j = kqmaxd + integ
do 4000 n = 1, j
4000    c(n + 1) = (tp1*c(n)) / xi(n - integ) + gamma(n - integ) * c(n)
!     C(N) NOW CORRESPONDS TO THE DIFFERENTIAL COEFFICIENT
!          D(N-INTEG,-INTEG)
4010 integz = integ
if (ici .ne. 0) go to 3980
! END OF COMPUTING DIFFERENTIATION COEFFICIENTS
! ********
! BEGINNING OF LOOP TO DO
!         INTEGRATION       (INTEG.GT.0)
!         INTERPOLATION     (INTEG.EQ.0)
!         DIFFERENTIATION   (INTEG.LT.0)
! TO THE POINT INDICATED BY T.
! ********
!     SET UP INITIAL INDICES
4020 if (nyny .lt. 0) then
   iy = -nyny
   iyni = nyny + ici + 1
   if (ldt .eq. 2) then
      csum(ici+1) = c(ici+1)
      do 4025 j = ici+2, integ+kqmxi
      csum(j) = csum(j-1) + c(j)
4025       continue
   end if
else
   iy = 1
   iyni = nyny + ici - 1
end if
idt = ndtf - integz
do 4140 i = 1, nte
   if (nkdko .ne. 0) kordi = kord(nkdko + i - 1)
   iy = iy + abs(kordi)
   kqq = kord(i + 3)
!         GET INDEX OF HIGHEST ORDER DIFFERENCE TO BE USED
   k = max(abs(kqq) + interp, 2)
   iyi = -integ
   if (kqq) 4030, 4130, 4040
! EQUATION IS STIFF
4030    continue
!++  Code for STIFF is inactive
!      JS=abs(KORD(NJSKO+I-1))-1
!      IYI=IYI-JS
!      IF(LNOTM1) IF (IYI) 4034,4032,4130
!      IF (KORDI.LT.0) IYI=IYI+1
!      IYI=IYI+MAXINT-abs(KORDI)
!      IF (IYI) 4034,4130,4130
!c.      IF EQUATION IS IMPLICIT DO NOT COMPUTE AN F
! 4032 IF (KORDI.LT.0) GO TO 4130
!c.      TEST IF INTEG TOO BIG FOR THIS EQUATION
! 4034 IF (abs(KORDI).LT.-IYI) GO TO 4130
!      IYI=IYI+IY
!      IYN=IYI+IYNI
!c. COMPUTE INNER PRODUCT FOR STIFF EQUATIONS
!      IF (INTEGZ.EQ.0) GO TO ???
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
!      IF (INTEG.EQ.0) GO TO 4120
!      TP1=TP1 + C(ICI+1)*Y(IYN+1)
!++  End
   go to 4100
! END OF SPECIAL CODE FOR STIFF EQUATIONS
!
! EQUATION IS NOT STIFF
4040    if (lnotm1) if (iyi) 4050, 4060, 4130
   iyi = iyi + maxint - kordi
   if (iyi .ge. 0) go to 4130
!       TEST IF INTEG TOO BIG FOR THIS EQUATION
4050    if (kordi .lt. -iyi) go to 4130
4060    iyi = iyi + iy
   iyn = iyi + iyni
!  COMPUTE INNER PRODUCT FOR EQUATION WHICH IS NOT STIFF
   xp1 = c0
   if (ldt .eq. 2) then
      if (kqq .ne. kqmaxi) xp1 = csum(k+integz+ici) * &
         f(idt+integz+numdt-1)
   end if
   do 4070 j = k + integz + ici, ici + 1, -1
      xp1 = xp1 + c(j) * f(idt - ici - 1 + j)
4070       continue
   if (integ) 4080, 4090, 4100
! STORE FINAL RESULT IN Y WHEN DIFFERENTIATING
4080    continue
   y(iyi) = xp1
   go to 4130
! STORE INTERPOLATED VALUE IN F (OR STIFF DIFFERENTIATION)
4090    f(i) = xp1
   go to 4130
! PICK UP EXTRA STUFF TO ADD TO INNER PRODUCT WHEN INTEGRATING
4100    k = ici
   if (k .eq. 0) go to 4120
4110    continue
   xp1 = c(k) * (xp1 + y(iyn))
   iyn = iyn - 1
   k = k - 1
   if (k .ne. 0) go to 4110
! STORE FINAL RESULT IN Y WHEN INTEGRATING (OR STIFF INTERPOLATION)
4120    y(iyi) = xp1 + y(iyn)
4130    continue
   idt = idt + numdt
4140    continue
!
integ = integ - 1
if (integ .ge. -maxdif) if (integ) 3990, 3950, 3920
4160 return
! ********
! ERROR PROCESSING
! ********
4170 mact(2) = 68
mact(3) = 11
mact(6) = ltxtad
idat(1) = ldt
go to 4190
4180 mact(2) = 28
mact(3) = 1
mact(6) = ltxtac
fdat(2) = tn
fdat(3) = t(2)
fdat(4) = xi(1)
fdat(5) = tn - t(2)
fdat(6) = tn + c2 * xi(1)
if (xi(1) .lt. 0) then
   fdat(5) = fdat(6)
   fdat(6) = tn - t(2)
end if
4190 fdat(1) = t(1)
!--D Next line special: P=>S, X=>D
call dmess(mact, mtxtaa, idat, fdat)
if (mact(2) .lt. 50) go to 3820
return
end
!   End of DIVAIN

subroutine divaop(iopt, fopt)
!>> 1987-12-07 DIVAOP Krogh   Initial code.
!
!  SUBROUTINE TO SET UP OPTIONS FOR DIFFERENTIAL EQUATION  PACKAGE -IVA
double precision fopt(*)
integer iopt(*)
!
integer kdim, maxord, maxstf
!++ Substitute for KDIM, MAXORD, MAXSTF below
parameter (kdim = 20, maxord = 2, maxstf = 1)
!--D Next line special: P=>D, X=>Q
double precision tn
double precision xi(kdim)
!
!--D Next line special: P=>D, X=>Q
double precision tg(2), tgstop(2), tmark, tmarkx, tout, tolg
double precision alpha(kdim), beta(kdim+1)
double precision  d(maxstf+maxord,maxord), g(kdim,maxord)
double precision v(kdim+maxord)
double precision hc, hdec, hinc, hincc, hmax, hmaxp9, hmin
double precision fdat(11)
!
double precision ds(maxstf+maxord, maxord), gs(kdim)
double precision sigma(kdim), rbq(kdim), dnoise
double precision eave, eimax, eimin, emax, erep, robnd, snoise
!
integer iopst, kordi, kqmaxd, kqmaxi, ldt, maxdif, maxint, nkdko, &
   nte, nyny, ndtf, numdt
common / divasc / tn, xi, iopst, kordi, kqmaxd, kqmaxi, ldt, &
   maxdif, maxint, nkdko, nte, nyny, ndtf, numdt
!
integer icf,ics,igflg,igtype(2),igstop(2),ilgrep,ings,iop3,iop4, &
   iop5,iop6,iop7,iop8,iop9,iop10,iop11,iop12,iop13,iop14,iop15, &
   iop16,iop17,iop18,iop19,iop20,iop21,iop22,iop21s,itolep,iy, &
   kemax,kis,kmark,kord1i,kord2i,kpred,kqdcon,kqicon,kqmaxs, &
   kqmxds,kqmxil,kqmxip,kqmxis,ksc,ksout,ksstrt,kstep,lex,linc, &
   lincd,lincq,lsc,maxkqd,maxkqi,method,ne,neptol,ng,ngtot, &
   noiseq,noutko,ntolf,ny,idat(6)
common /divamc/ tg,tgstop,tmark,tmarkx,tout,tolg,hc,hdec,hinc, &
   hincc,hmax,hmaxp9,hmin,alpha,beta,d,g,v,ds,gs,sigma,rbq,dnoise, &
   eave,eimax,eimin,emax,erep,robnd,snoise,fdat,icf,ics,igflg, &
   igtype,igstop,ilgrep,ings,iop3,iop4,iop5,iop6,iop7,iop8,iop9, &
   iop10,iop11,iop12,iop13,iop14,iop15,iop16,iop17,iop18,iop19, &
   iop20,iop21,iop22,iop21s,itolep,iy,kemax,kis,kmark,kord1i, &
   kord2i,kpred,kqdcon,kqicon,kqmaxs,kqmxds,kqmxil,kqmxip,kqmxis, &
   ksc,ksout,ksstrt,kstep,lex,linc,lincd,lincq,lsc,maxkqd,maxkqi, &
   method,ne,neptol,ng,ngtot,noiseq,noutko,ntolf,ny,idat
save / divamc / , / divasc /
!.    SPECIFICATION OF ENVIRONMENTAL CONSTANTS.
double precision eeps10, eeps16, erov10, eeps2, eept75, eovep2
double precision ovtm75, ovd10
common / divaev / eeps2, eept75, eovep2, ovtm75, ovd10, eeps10, &
   eeps16, erov10
save / divaev /
!
integer iopts(23), incop(22), ioptc(23), i, ia, j, k, liopt, multj
double precision cmp75, c0, cp25, cp3, cp5, cp625, cp75, cp875, &
   cp9, c1, c1p125, c2, c4, c10, c16
parameter (cmp75 = (-.75d0))
parameter (c0 = 0.d0)
parameter (cp25 = .25d0)
parameter (cp3 = .3d0)
parameter (cp5 = .5d0)
parameter (cp625 = .625d0)
parameter (cp75 = .75d0)
parameter (cp875 = .875d0)
parameter (cp9 = .9d0)
parameter (c1 = 1.d0)
parameter (c1p125 = 1.125d0)
parameter (c2 = 2.d0)
parameter (c4 = 4.d0)
parameter (c10 = 10.d0)
parameter (c16 = 16.d0)
external d1mach
double precision d1mach
equivalence (ioptc(3), iop3)
save iopts, liopt
!
!                      Declarations for error message processing.
!
integer mecont, meret, meemes, meivec
parameter (mecont =50)
parameter (meret  =51)
parameter (meemes =52)
parameter (meivec =57)
integer mact(7), mact1(5)
!
! ********* Error message text ***************
![Last 2 letters of Param. name]  [Text generating message.]
!AA DIVAOP$B
!AB Error in IOPT() specifications: IOPT =$E
!AC HMIN = $F is > HMAX = $F.$E
integer ltxtaa,ltxtab,ltxtac
parameter (ltxtaa= 1,ltxtab= 9,ltxtac=49)
character mtxtaa(1) * (76)
data mtxtaa/'DIVAOP$BError in IOPT() specifications: IOPT =$EHMIN$ &
 = $f is > hmax = $f.$e'/
! **** End of text generated by pmess
!                      1   2  3   4       5  6      7
data mact / meemes, 88, 24, 0, meivec, 0, meret /
data mact1 / meemes, 28, 24, ltxtac, meret /
!
!                      IOP4       IOP17
data iopts / 3*0, 500000, 12*0, 1, 6*0 /
!
!                  1  2    3  8  9 10 11 12  13  16   17 21 22
data incop / 1, 3, 5*2, 1, 2, 3, 1, 2, 3*1, 3, 4*2, 2, 2 /
!
! ********* START OF EXECUTABLE CODE ***********************
!
multj = 1
k = 1
4200 i = iopt(k)
ia = abs(i)
! 1 and 6 lines below need 21 changed if more options are added.
if (ia .le. 21) if (i) 4220, 4520, 4280
if (ia .ne. 1111) go to 4490
if (i .lt. 0) then
  multj = -1
  k = k + 1
  go to 4200
end if
iopt(2) = liopt
!
!     ****  INITIALIZE FOR STARTING A NEW INTEGRATION
do 4210 j = 3, 23
4210    ioptc(j) = iopts(j)
ksout = iopts(4)
kmark = 1 - iopts(1)
kordi = iopts(17)
nkdko = max(-kordi, 0)
iopst = iopts(22)
go to 4260
!
!     **** SET A NOMINAL VALUE
4220 iopts(ia) = 0
if (ia .eq. 12) go to 4420
if (ia - 2) 4400, 4240, 4230
4230 if (ia .eq. 4) iopts(4) = 500000
if (ia .eq. 21) tolg = 0.d0
go to 4390
!
!     **** SET ALL OPTIONS TO THEIR NOMINAL VALUES
4240 ia = 1
iopts(1) = 0
do 4250 j = 3, 22
   iopts(j) = 0
   ioptc(j) = 0
4250 continue
iopts(4) = 500000
ioptc(4) = iopts(4)
iopts(17) = 1
tolg = 0.d0
4260 ngtot = iopts(7) + max(iopts(6), 0)
if (iopts(12) .eq. 0) go to 4420
4270 return
!
!     **** SET SPECIFIED OPTION
4280 j = iopt(k + 1)
if (incop(ia) - 2) 4290, 4330, 4300
!     **** OPTION INVOLVES NO EXTRA PARAMETERS
4290 iopts(ia) = 1
if (ia - 2) 4400, 4400, 4390
!     **** TAKE CARE OF SECOND EXTRA PARAMETER
4300 if (ia .ne. 10) go to 4310
noutko = iopt(k + 2)
if (noutko) 4500, 4350, 4350
4310 if (ia .ne. 16) go to 4320
ntolf = iopt(k + 2)
if (ntolf) 4500, 4500, 4350
4320 if (j .eq. 3) then
  if (kmark .ne. 3) then
     if (xi(1)*(fopt(iopt(k+2)) - tmark) .ge. c0) go to 4400
  end if
end if
tmark = fopt(iopt(k+2))
kmark = j
go to 4400
!     **** TAKE CARE OF FIRST EXTRA PARAMETER
4330 continue
if (ia .eq. 12) go to 4410
if (ia .eq. 4) ksout = j
if (ia .eq. 21) tolg = fopt(j)
4350 iopts(ia) = j * multj
if (abs(ia - 7) .gt. 1) go to 4360
!     **** SET SPECIAL PARAMETERS FOR GSTOP-S
igflg = 0
ngtot = iopts(7) + max(iopts(6), 0)
!     **** TEST FOR ERROR
if (j .gt. 500) go to 4500
4360 if (j .gt. 0) go to 4390
if ((ia .eq. 5) .or. (ia .eq. 17)) go to 4390
if (j + 1) 4500, 4370, 4380
4370 if (ia .eq. 7) go to 4500
4380 if ((ia .eq. 4) .or. (ia .eq. 11) .or. (ia .ge. 16)) go to 4500
!     **** STORE SAVED VALUE IN COMMON
4390 ioptc(ia) = iopts(ia)
!
!     **** INCREMENT K TO GET NEXT OPTION
4400 k = k + incop(ia)
go to 4200
!
! ******* SET UP INFORMATION FOR CHANGING STEPSIZE *********
!
!     **** TEST IF VALUES ARE ALREADY SET
4410 if (iopts(12) .ne. 0) go to 4430
!     **** SET NOMINAL VALUES FOR VARIABLES ONLY SET ONCE
4420 erep = cp3
!     **** SET NOMINAL VALUES FOR STEPSIZE CONTROL AND ENV. CONSTANTS
eeps2 = d1mach(4)
eeps16 = c16 * eeps2
eeps10 = cp625 * eeps16
eept75 = eeps2 ** cp75
eeps2 = eeps2 + eeps2
ovd10 = d1mach(2)
erov10 = c10 / ovd10
eovep2 = ovd10 * eeps2
ovtm75 = ovd10 ** cmp75
ovd10 = ovd10 / c10
hinc = c2
hdec = cp5
hmin = erov10
hmax = ovd10
if (i .ne. 12) go to 4470
4430 iopts(12) = j * multj
if (j) 4450, 4470, 4460
!     **** SET UP TO GIVE USER COMPLETE STEPSIZE CONTROL
4450 erep = c1 / erov10
hinc = -c2
iop8 = 1
!## Recent code 12/16/94
lincd = -2
linc = -2
robnd = c1
!## Recent code 9/6/2001
tolg = 0.d0
!## End of recent code
go to 4480
!     **** SET USER VALUES FOR STEPSIZE CONTROL
4460 if (fopt(j) .ne. c0) hinc = max(c1p125, min(fopt(j),c4))
if (fopt(j + 1) .ne. c0) hdec = min(cp875, max(fopt(j + 1), cp25))
if (fopt(j + 2) .ne. c0) hmin = fopt(j + 2)
if (fopt(j + 3) .ne. c0) hmax = fopt(j + 3)
if ((hmin .gt. hmax) .or. (hmax .le. 0.d0)) then
   call dmess(mact1, mtxtaa, iopt, fopt(j+2))
   kord2i = -4
end if
4470 kqicon = -1
4480 hmaxp9 = hmax * cp9
if (i - 1111) 4400, 4270, 4400
!
! ***************** ERROR  IN  IOPT ************************
!
4490 ia = 1
4500 mact(6) = k + incop(ia) - 1
call mess(mact, mtxtaa, iopt)
kord1i = 24
kord2i = -4
! Usual return with no error is here.
4520 liopt = k
return
end
!   End of DIVAOP

subroutine divapr(y, yn, f, kord)
!>> 1988-01-13 DIVAPR Krogh   Initial code.
!
! THIS SUBROUTINE
!   1. UPDATES THE DIFFERENCE TABLE FROM THE PREVIOUS STEP (IF NOT
!      DONE ALREADY).
!   2. PREDICTS WHAT THE VALUES OF THE DEPENDENT VARIABLES, Y, AND
!      THE DIFFERENCE TABLE, DT, WILL BE AT THE END OF THE CURRENT STEP.
!
!   Y = VECTOR OF PREDICTED VALUES COMPUTED BY THIS SUBROUTINE.
!   YN= VECTOR OF VALUES OF Y COMPUTED ON THE LAST STEP.
!   F = VECTOR OF DERIVATIVE VALUES.
!   DT= ARRAY CONTAINING DIFFERENCE TABLES.
!   KD= VECTOR GIVING ORDERS OF THE DIFFERENTIAL EQUATIONS (IF
!       EQUATIONS HAVE DIFFERENT ORDERS).
!   KQ= VECTOR OF INTEGRATION ORDERS.
!
integer kord(*)
!--D Next line special: P=>D, X=>Q
double precision y(*), yn(*)
double precision f(*)
!
integer kdim, maxord, maxstf
!++ Substitute for KDIM, MAXORD, MAXSTF below
parameter (kdim = 20, maxord = 2, maxstf = 1)
!--D Next line special: P=>D, X=>Q
double precision tn
double precision xi(kdim)
!
!--D Next line special: P=>D, X=>Q
double precision tg(2), tgstop(2), tmark, tmarkx, tout, tolg
double precision alpha(kdim), beta(kdim+1)
double precision  d(maxstf+maxord,maxord), g(kdim,maxord)
double precision v(kdim+maxord)
double precision hc, hdec, hinc, hincc, hmax, hmaxp9, hmin
double precision fdat(11)
!
double precision ds(maxstf+maxord, maxord), gs(kdim)
double precision sigma(kdim), rbq(kdim), dnoise
double precision eave, eimax, eimin, emax, erep, robnd, snoise
!
!.    SPECIFICATION OF ENVIRONMENTAL CONSTANTS.
double precision eeps10, eeps16, erov10, eeps2
double precision eept75, eovep2, ovtm75, ovd10
common / divaev / eeps2, eept75, eovep2, ovtm75, ovd10, eeps10, &
   eeps16, erov10
save / divaev /
integer iopst, kordi, kqmaxd, kqmaxi, ldt, maxdif, maxint, nkdko, &
   nte, nyny, ndtf, numdt
common / divasc / tn, xi, iopst, kordi, kqmaxd, kqmaxi, ldt, &
   maxdif, maxint, nkdko, nte, nyny, ndtf, numdt
!
integer icf,ics,igflg,igtype(2),igstop(2),ilgrep,ings,iop3,iop4, &
   iop5,iop6,iop7,iop8,iop9,iop10,iop11,iop12,iop13,iop14,iop15, &
   iop16,iop17,iop18,iop19,iop20,iop21,iop22,iop21s,itolep,iy, &
   kemax,kis,kmark,kord1i,kord2i,kpred,kqdcon,kqicon,kqmaxs, &
   kqmxds,kqmxil,kqmxip,kqmxis,ksc,ksout,ksstrt,kstep,lex,linc, &
   lincd,lincq,lsc,maxkqd,maxkqi,method,ne,neptol,ng,ngtot, &
   noiseq,noutko,ntolf,ny,idat(6)
common /divamc/ tg,tgstop,tmark,tmarkx,tout,tolg,hc,hdec,hinc, &
   hincc,hmax,hmaxp9,hmin,alpha,beta,d,g,v,ds,gs,sigma,rbq,dnoise, &
   eave,eimax,eimin,emax,erep,robnd,snoise,fdat,icf,ics,igflg, &
   igtype,igstop,ilgrep,ings,iop3,iop4,iop5,iop6,iop7,iop8,iop9, &
   iop10,iop11,iop12,iop13,iop14,iop15,iop16,iop17,iop18,iop19, &
   iop20,iop21,iop22,iop21s,itolep,iy,kemax,kis,kmark,kord1i, &
   kord2i,kpred,kqdcon,kqicon,kqmaxs,kqmxds,kqmxil,kqmxip,kqmxis, &
   ksc,ksout,ksstrt,kstep,lex,linc,lincd,lincq,lsc,maxkqd,maxkqi, &
   method,ne,neptol,ng,ngtot,noiseq,noutko,ntolf,ny,idat
save / divamc / , / divasc /
double precision c0
parameter (c0 = 0.d0)
!
integer i, integ, integs, j, k, kqq, l, n
double precision temp(kdim)
double precision tp1
!--D Next line special: P=>D, X=>Q
double precision xp
data integs / -1 /
!
!++  Code for ARGM is inactive
!      RETURN
!      ENTRY DIVAPE
!++  End
! ********
! START OF CODE
! ********
iy = 0
l = ndtf - 1
do 4680 i = 1, nte
   integ = kordi
   if (nkdko .ne. 0) integ = kord(nkdko + i - 1)
   kqq = kord(i + 3)
   k = max(abs(kqq), 2)
   if (kqq) 4530, 4520, 4540
4520    iy = iy + abs(integ)
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
4540    n = kqq
   if (ldt .ne. 0) if (k - ksc) 4570, 4570, 4550
!     DIFFERENCE TABLE HAS NOT BEEN UPDATED
   tp1 = f(i) - f(l + 1)
   if (k - ksc) 4610, 4610, 4590
! END OF SET-UP FOR EQUATIONS WHICH ARE NOT STIFF
! ********
! GET PREDICTED DIFFERENCES FROM UPDATED DIFFERENCE TABLE
! ********
4550    f(l + k + 1) = f(l + k + 1) * beta(k + 1)
   temp(k) = f(l + k) * beta(k)
   f(l + k) = temp(k)
! LOOP FOR MODIFIED DIVIDED DIFFERENCES
4560    k = k - 1
   if (k .le. ksc) go to 4580
   temp(k) = f(l + k) * beta(k)
   f(l + k) = temp(k) + f(l + k + 1)
   go to 4560
! CODE FOR BACKWARD DIFFERENCES
4570    f(l + k + 1) = f(l + k + 1)
   temp(k) = f(l + k)
   k = k - 1
!
4580    temp(k) = f(l + k)
   f(l + k) = temp(k) + f(l + k + 1)
   k = k - 1
   if (k .ne. 0) go to 4580
   go to 4630
! ********
! UPDATE DIFFERENCE TABLE AND GET PREDICTED DIFFERENCES
! ********
! CODE FOR MODIFIED DIVIDED DIFFERENCES
4590    f(l + k + 1) = (f(l+k+1) + tp1) * beta(k + 1)
   temp(k) = (f(l + k) + tp1) * beta(k)
   f(l + k) = temp(k)
4600    k = k - 1
   if (k .le. ksc) go to 4620
   temp(k) = (f(l + k) + tp1) * beta(k)
   f(l + k) = temp(k) + f(l + k + 1)
   go to 4600
! CODE FOR BACKWARD DIFFERENCES
4610    f(l + k + 1) = (f(l+k+1) + tp1)
   temp(k) = f(l + k) + tp1
   f(l + k) = temp(k)
   k = k - 1
!
4620    temp(k) = f(l + k) + tp1
   f(l + k) = temp(k) + f(l + k + 1)
   k = k - 1
   if (k .ne. 0) go to 4620
! ********
! COMPUTE Y-S OBTAINED USING INTEGRATION
! ********
!     TEST IF NEXT Y TO BE OBTAINED BY INTERPOLATION
4630    continue
!++  Code for STIFF is inactive
!      IF (INTEG.EQ.0) GO TO 4662
!++  End
   iy = iy + 1
!     FORM INNER PRODUCT
   xp = c0
   do 4650 j = integs + n + 1, integs + 2, -1
!++  Code for ~{p,x} is active
      xp = xp + g(j, integ) * temp(j)
!++  Code for {p,x} is inactive
!c--D Next line special: P=>D, X=>Q
!            XP = XP + dble(G(J, INTEG)) * dble(TEMP(J))
!++  END
4650    continue
   k = integ + integs
   do 4660 j = k, 1, -1
!++  Code for ~{p,x} is active
      xp = xp + g(1, j) * yn(iy + j)
!++  Code for {p,x} is inactive
!c--D Next line special: P=>D, X=>Q
!            XP = XP + dble(G(1, J)) * dble(YN(IY + J))
!++  END
4660    continue
   y(iy) = yn(iy) + xp
   integ = integ - 1
   if (k) 4670, 4670, 4630
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
!      IF (KQQ.EQ.1) Y(IY)=YN(IY)
! 4663 INTEG=INTEG+1
!      IF (INTEG.EQ.JS) IF (IMPLIC) 4680,4680,4664
!c.    COMPUTE INTEG-TH DERIVATIVE
!      XP = C0
! 4664 DO 4666 J = KQQ+1, INTEG+1, -1
!         XP = XP + D(J, INTEG) * TEMP(J)
! 4666 CONTINUE
!      IF (INTEG.EQ.JS) GO TO 4667
!      IY=IY+1
!      Y(IY)=XP
!      GO TO 4663
!c.STORE PREDICTED VALUE FOR F
! 4667 CONTINUE
!      F(L+NUMDT)=XP
!++  End
4670    l = l + numdt
4680    continue
ldt = -3
return
end
subroutine divadb(lprint, tspecs, y, f, kord, text)
! Copyright (c) 1996 California Institute of Technology, Pasadena, CA.
! ALL RIGHTS RESERVED.
! Based on Government Sponsored Research NAS7-03001.
!>> 2009-11-04 DIVADB Krogh Included TOLG, initilized the unitialized.
!>> 2000-12-01 DIVADB Krogh Removed unused parameter METXTF.
!>> 1996-07-02 DIVADB Krogh Transpose flag for matrix output in C.
!>> 1996-03-25 DIVADB Krogh Introduced TEXT1-TEXT4 to comply with F77.
!>> 1996-01-19 DIVADB Krogh Changed NTEXT to TEXT to agree with doc.
!>> 1995-04-26 DIVADB Krogh Fixed print of V & G's for high order eqs.
!>> 1994-11-11 DIVADB Krogh Declared all vars.
!>> 1994-10-20 DIVADB Krogh Changes to use M77CON
!>> 1994-09-12 DIVADB Krogh Added CHGTYP code.
!>> 1994-03-07 DIVADB Krogh Allow larger order in single precision.
!>> 1993-05-03 DIVADB Krogh Additions for Conversion to C.
!>> 1993-04-14 DIVADB Krogh Changes for new MESS usage.
!>> 1992-04-08 DIVADB Krogh Unused labels 10 and 60 removed.
!>> 1992-03-10 DIVADB Krogh Fixed value for KDIM in single p. version.
!>> 1992-02-17 DIVADB Krogh Made tabs depend on # digits output.
!>> 1991-11-04 DIVADB Krogh Switched to use MESS, DMESS
!>> 1990-03-08 DIVADB Krogh Unused stiff vars. set to 0.
!>> 1989-07-21 DIVADB Krogh Code for integrating discontinuities
!>> 1988-06-07 DIVADB Krogh Dim. of IVC2 and DVC2 upped by 1 (old bug)
!>> 1987-12-07 DIVADB Krogh Initial code.
!
!--D replaces "?": ?IVADB, ?IVAEV, ?IVAMC, ?IVASC, ?MESS
!
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
integer lprint, kord(*)
character text*(*)
character text1(1)*11, text2(1)*4, text3(1)*5, text4(1)*4
integer ivc1(12), ivc2(65), j, k, l, n1, n2
double precision  dvc2(7), rvc2(8), evc(8)
!--D Next line special: P=>D, X=>Q
double precision tspecs(*), y(*), tneq(1), dvc1(7)
double precision f(*)
!
!++S Default KDIM = 16
!++  Default KDIM = 20
!++  Default MAXORD = 2, MAXSTF = 1
!++  Default STIFF=.F.
integer kdim, maxord, maxstf
!++ Substitute for KDIM, MAXORD, MAXSTF below
parameter (kdim = 20, maxord = 2, maxstf = 1)
!--D Next line special: P=>D, X=>Q
double precision tn
double precision xi(kdim)
!
!--D Next line special: P=>D, X=>Q
double precision tg(2), tgstop(2), tmark, tmarkx, tout, tolg
double precision alpha(kdim), beta(kdim+1)
double precision  d(maxstf+maxord,maxord), g(kdim,maxord)
double precision v(kdim+maxord)
double precision hc, hdec, hinc, hincc, hmax, hmaxp9, hmin
double precision fdat(11)
!
double precision ds(maxstf+maxord, maxord), gs(kdim)
double precision sigma(kdim), rbq(kdim), dnoise
double precision eave, eimax, eimin, emax, erep, robnd, snoise
!
!.    SPECIFICATION OF ENVIRONMENTAL CONSTANTS.
double precision eeps10, eeps16, erov10, eeps2
double precision eept75, eovep2, ovtm75, ovd10
common / divaev / eeps2, eept75, eovep2, ovtm75, ovd10, eeps10, &
   eeps16, erov10
save / divaev /
integer iopst, kordi, kqmaxd, kqmaxi, ldt, maxdif, maxint, nkdko, &
   nte, nyny, ndtf, numdt
common / divasc / tn, xi, iopst, kordi, kqmaxd, kqmaxi, ldt, &
   maxdif, maxint, nkdko, nte, nyny, ndtf, numdt
!
integer icf,ics,igflg,igtype(2),igstop(2),ilgrep,ings,iop3,iop4, &
   iop5,iop6,iop7,iop8,iop9,iop10,iop11,iop12,iop13,iop14,iop15, &
   iop16,iop17,iop18,iop19,iop20,iop21,iop22,iop21s,itolep,iy, &
   kemax,kis,kmark,kord1i,kord2i,kpred,kqdcon,kqicon,kqmaxs, &
   kqmxds,kqmxil,kqmxip,kqmxis,ksc,ksout,ksstrt,kstep,lex,linc, &
   lincd,lincq,lsc,maxkqd,maxkqi,method,ne,neptol,ng,ngtot, &
   noiseq,noutko,ntolf,ny,idat(6)
common /divamc/ tg,tgstop,tmark,tmarkx,tout,tolg,hc,hdec,hinc, &
   hincc,hmax,hmaxp9,hmin,alpha,beta,d,g,v,ds,gs,sigma,rbq,dnoise, &
   eave,eimax,eimin,emax,erep,robnd,snoise,fdat,icf,ics,igflg, &
   igtype,igstop,ilgrep,ings,iop3,iop4,iop5,iop6,iop7,iop8,iop9, &
   iop10,iop11,iop12,iop13,iop14,iop15,iop16,iop17,iop18,iop19, &
   iop20,iop21,iop22,iop21s,itolep,iy,kemax,kis,kmark,kord1i, &
   kord2i,kpred,kqdcon,kqicon,kqmaxs,kqmxds,kqmxil,kqmxip,kqmxis, &
   ksc,ksout,ksstrt,kstep,lex,linc,lincd,lincq,lsc,maxkqd,maxkqi, &
   method,ne,neptol,ng,ngtot,noiseq,noutko,ntolf,ny,idat
save / divamc / , / divasc /
equivalence (ivc1(1), iopst), (ivc2(1), icf)
equivalence (tneq, tn)
equivalence (rvc2(1), dnoise), (dvc1(1), tg), (dvc2(1), hc), &
   (evc(1), eeps2)
!
!                      Declarations for error message processing.
integer meddig, neddig, metdig, metabs, meret, metext, &
   metabl, meivec, mefvec, mefmat
parameter (meddig =12)
parameter (neddig =-meddig)
parameter (metdig =22)
parameter (metabs =32)
parameter (meret  =51)
parameter (metext =53)
parameter (metabl =55)
parameter (meivec =57)
parameter (mefvec =61)
parameter (mefmat =62)
integer mact0(3), mact1(2), mact2(7), mact3(7), mact4(8), &
   mact5(11), mact6(3), mact7(14), mactfv(4)
integer kpi, kpe
!                      wddtrr        wwddtrr
parameter (kpi = 400301, kpe = 1305501)
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
integer ltxtaa,ltxtab,ltxtac,ltxtad,ltxtae,ltxtaf,ltxtag,ltxtah, &
 ltxtai,ltxtaj,ltxtak,ltxtal,ltxtam,ltxtan,ltxtao,ltxtap,ltxtaq
parameter (ltxtaa=  1,ltxtab= 14,ltxtac=  1,ltxtad=  1,ltxtae=  1, &
 ltxtaf= 21,ltxtag=  1,ltxtah= 16,ltxtai= 22,ltxtaj=  1, &
 ltxtak=  1,ltxtal=  1,ltxtam=655,ltxtan=  1,ltxtao=  1, &
 ltxtap=  1,ltxtaq=  1)
character mtxtaa(1) * (26)
character mtxtab(1) * (13)
character mtxtac(1) * (13)
character mtxtad(1) * (34)
character mtxtae(1) * (28)
character mtxtaf(1) * (9)
character mtxtag(1) * (120)
character mtxtah(3) * (242)
character mtxtai(1) * (80)
character mtxtaj(1) * (68)
character mtxtak(1) * (84)
character mtxtal(1) * (88)
data mtxtaa/'$NKORD:    $BInt. Ord.: $B'/
data mtxtab/'D.E. Ord.: $B'/
data mtxtac/'Meth.Type: $B'/
data mtxtad/'Tolerance Groups: $BTolerances: $B'/
data mtxtae/'$NDifferences$BEq. $#Ord. $#'/
data mtxtaf/'$NTN=$F$E'/
data mtxtag/'$NIOPST=$I$TKORDI=$I$TKQMAXD=$I$TKQMAXI=$I$TLDT=$I$TM &
axdif=$i$tmaxint=$i$tnkdko=$i$tnte=$i$tnyny=$i$tndtf=$i$tnumdt=$i$ &
e'/
data mtxtah/'$NICF=$I$TICS=$I$TIGFLG=$I$TIGTYPE(1)=$I$TIGTYPE(2)=$ &
i$tigstop(1)=$i$tigstop(2)=$i$tilgrep=$i$tings=$i$tiop3=$i$tiop4=$ &
i$tiop5=$i$tiop6=$i$tiop7=$i$tiop8=$i$tiop9=$i$tiop10=$i$tiop11=$i &
$tiop12=$i$tiop13=$i$tiop14=$i$tiop15=$i$tiop16=$i$tiop17','=$i$ti &
op18=$i$tiop19=$i$tiop20=$i$tiop21=$i$tiop22=$i$tiop21s=$i$titolep &
=$i$tiy=$i$tkemax=$i$tkis=$i$tkmark=$i$tkord1i=$i$tkord2i=$i$tkpre &
d=$i$tkqdcon=$i$tkqicon=$i$tkqmaxs=$i$tkqmxds=$i$tkqmxil=$i$tkqmxi &
p=$i$tkqmxis=$i$tksc=$i$tksout=$i$tkss','trt=$i$tkstep=$i$tlex=$i$ &
tlinc=$i$tlincd=$i$tlincq=$i$tlsc=$i$tmaxkqd=$i$tmaxkqi=$i$tmethod &
=$i$tne=$i$tneptol=$i$tng=$i$tngtot=$i$tnoiseq=$i$tnoutko=$i$tntol &
f=$i$tny=$i$e$ndnoise=$f$teave=$f$teimax=$f$teimin=$f$temax=$f$ter &
ep=$f$trobnd=$f$e  '/
data mtxtai/'$NTG(1)=$F$TTG(2)=$F$TTGSTOP(1)=$F$TTGSTOP(2)=$F$TTMA &
rk=$f$ttmarkx=$f$ttout=$f$e'/
data mtxtaj/'HC=$F$THDEC=$F$THINC=$F$THINCC=$F$THMAX=$F$THMAXP9=$F &
$thmin=$f$t$n$e'/
data mtxtak/'K$HXI(K)$HBETA(K)$HALPHA(K)$HG(K,1)$HRBQ(K)$HSIGMA(K) &
$hgs(k)$hv(k)$hg(k,2..maxint)$e'/
data mtxtal/'$NEEPS2=$F$TEEPT75=$F$TEOVEP2=$F$TOVTM75=$F$TOVD10=$F &
$teeps10=$f$teeps16=$f$terov10=$f$e'/
!
data mact0 / metabs, 10, meret /
data mact1 / metext, meret /
!                       1       2  3       4       5  6      7
data mact2 / metext, meivec, 3, metext, meivec, 0, meret /
!                       1       2  3       4       5  6      7
data mact3 / metext, meivec, 0, metext, mefvec, 0, meret /
!                       1       2  3  4  5       6       7
data mact4 / metext, mefmat, 0, 0, 0, ltxtai, ltxtah, meret /
!                       1   2       3       4   5       6  7       8
data mact5 / metabs, 12, metext, metabs, 18, metdig, 5, metext, &
             metabs, 0, meret /
!                       9 10     11
data mact6 / neddig, 0, meret /
!                         2 3 4   5 6 7 8 9  10  11  12 13 14
data mact7 / metabl,0,0,0,kpi,0,0,0,0,kpe,kpe,kpe, 0, 0 /
!                        1       2  3      4
data mactfv / metext, mefvec, 3, meret /
!
data text1 / '$NTSPECS:$B' /
data text2 / 'Y:$B' /
data text3 / 'YN:$B' /
data text4 / 'F:$B' /
!
! ********
! START OF CODE -- PRINT TEXT AND SET INDEX FOR F
! ********
!    Getting variables that are not yet assigned some values.
!++  Code for ~STIFF is active
kqdcon = 0
kqmxds = 0
maxkqd = 0
!++  End
gs(1) = 1.d0
if (iop6 .eq. 0) then
  igtype(1) = 0
  igstop(1) = 0
  tg(1) = 0.d0
  tgstop(1) = 0.d0
end if
if (iop7 .eq. 0) then
  igtype(2) = 0
  igstop(2) = 0
  tg(2) = 0.d0
  tgstop(2) = 0.d0
end if
if (iop6 + iop7 .eq. 0) then
  ings = 0
  ng = 0
end if
if (iop10 .eq. 0) noutko = 0
j = 0
call messft(mact0, text)
!
n1 = lprint / 10
n2 = lprint - 10 * n1
if (n1 .le. 1) go to 80
! ********
! PRINT ALL EXTERNAL VARIABLES EXCEPT FOR THE DIFFERENCES
! ********
mactfv(3) = max(iop5, 4)
!--D Next line special: P=>D, X=>Q
call dmess(mactfv, text1, kord, tspecs)
mactfv(3) = ny
!--D Next line special: P=>D, X=>Q
call dmess(mactfv, text2, kord, y)
!--D Next line special: P=>D, X=>Q
call dmess(mactfv, text3, kord, y(nyny))
mactfv(3) = nte
!--D Next line special: P=>S, X=>D
call dmess(mactfv, text4, kord, f)
mact2(6) = nte
call mess(mact2, mtxtaa, kord)
if (nkdko .gt. 0) call mess(mact2(4), mtxtab, kord(nkdko))
if (iopst .gt. 0) call mess(mact2(4), mtxtac, kord(iopst))
! WRITE TOL
k = iop16
70 if (kord(k) .lt. 0) k = k + 1
k = k + 1
if (kord(k - 1) .lt. nte) go to 70
mact3(3) = k - iop16
mact3(6) = mact3(3)
!--D Next line special: P=>S, X=>D
call dmess(mact3, mtxtad, kord(iop16), f(ntolf))
if (n1 .eq. 2) go to 80
! ********
! WRITE THE DIFFERENCE TABLES
! ********
k = numdt
if (n1 .eq. 3) k = kqmaxs
mact4(3) = numdt
mact4(4) = -k
mact4(5) = nte
!--D Next line special: P=>S, X=>D
call dmess(mact4, mtxtae, kord, f(ndtf))
!
80 if (n2 .le. 1) return
! ********
! WRITE SCALARS IN COMMON
! ********
!--D Next line special: P=>D, X=>Q
call dmess(mact1, mtxtaf, kord, tneq)
!
! ===== COMMON 1  -- INTEGER
!
call mess(mact1, mtxtag, ivc1)
if (n2 .eq. 2) return
call mess(mact6, mtxtaa, idat)
mact5(10) = mact6(2) + 14
!
! ===== COMMON 2  -- INTEGER AND FLOATING POINT
!
!--D Next line special: P=>S, X=>D
call dmess(mact5, mtxtah, ivc2, rvc2)
!--D Next line special: P=>D, X=>Q
call dmess(mact1, mtxtai, ivc2, dvc1)
!--D Next line special: P=>S, X=>D
call dmess(mact1, mtxtaj, ivc2, dvc2)
if (n2 .eq. 3) return
!         wddtrr              wddtrr
j = 101000 * mact6(2) + 800501
mact7(2) = 1
mact7(3) = kqmaxs
mact7(4) = 8
do 90 k = 6, 9
   mact7(k) = j
90 continue
if (n2 .gt. 0) then
   mact7(4) = 8 + maxint
   mact7(13) = j
   l = min(maxint, 4)
   mact7(14) = j + l - 2
end if
do 100 k = 1, mact7(3)
   fdat(1) = xi(k)
   fdat(2) = beta(k)
   fdat(3) = alpha(k)
   fdat(4) = g(k, 1)
   fdat(5) = rbq(k)
   fdat(6) = sigma(k)
   fdat(7) = gs(k)
   if (n2 .ge. 4) then
      fdat(8) = v(k)
      do 95 j = 2, l
         fdat(7+j) = g(k, j)
95       continue
   end if
!--D Next line special: P=>S, X=>D
   call dmess(mact7, mtxtak, idat, fdat)
100 continue
!++  Code for STIFF is inactive
!     if (MAXDIF .le. 0) return
!        Need to define MACT8 and set values
!c--D Next line special: P=>S, X=>D
!     call DMESS(MACT8, 'D$B', IDAT, D)
!c--D Next line special: P=>S, X=>D
!     call DMESS(MACT8, 'DS$B', IDAT, DS)
!++  End
!
!--D Next line special: P=>S, X=>D
call dmess(mact1, mtxtal, idat, evc)
return
end
subroutine divag(tspecs, y, f, kord, iflag, nstop, gnew, gt)
! Copyright (c) 1996 California Institute of Technology, Pasadena, CA.
! ALL RIGHTS RESERVED.
! Based on Government Sponsored Research NAS7-03001.
!>> 2001-09-07 DIVAG  Krogh  Changes to allow user tol on G-Stops.
!>> 1995-06-20 DIVAG  Krogh  Fixed problem introduced with last change.
!>> 1995-05-09 DIVAG  Krogh  Fixed G-Stop/discontinuity code interaction
!>> 1994-11-11 DIVAG  Krogh  Declared all vars.
!>> 1994-10-20 DIVAG  Krogh  Changes to use M77CON
!>> 1994-09-12 DIVAG  Krogh  Added CHGTYP code.
!>> 1994-08-17 DIVAG  Krogh  Modified internal comment.
!>> 1994-03-07 DIVAG  Krogh  Allow larger order in single precision.
!>> 1993-04-27 DIVAG  Krogh  Additions for Conversion to C.
!>> 1992-10-12 DIVAG  Krogh  Fixed G-Stop/discontinuity code interaction
!>> 1992-09-17 DIVAG  Krogh  Slight change in check for sign change.
!>> 1992-04-08 DIVAG  Krogh  Unused labels 140,150,230, and 250 removed.
!>> 1992-03-10 DIVAG  Krogh  Fixed value for KDIM in single p. version.
!>> 1990-01-29 DIVAG  Krogh  Added arg to call to DERMN.
!>> 1988-03-04 DIVAG  Krogh  Initial code.
!
!--D replaces "?": ?IVAG,?IVABU,?IVAEV,?IVAIN,?IVAMC,?IVASC,?MESS,?ZERO
!
!     -XXXG(TSPECS, Y, F, KORD, IFLAG, NSTOP, GNEW, GT)
!
!  SUBROUTINE TO LOCATE OUTPUT POINTS AT ZEROS OF ARBITRARY
!  FUNCTIONS  **** GSTOPS **** FOR DIFFERENTIAL EQUATION
!  INTEGRATOR -ODE (OR -IVA).
!
integer kord(*), iflag, nstop
!--D Next line special: P=>D, X=>Q
double precision tspecs(*), y(*), told, tsave
double precision f(*), gnew(*), gt(*), gold
save gold, told, tsave
!
!++SP Default KDIM = 16
!++  Default KDIM = 20
!++  Default MAXORD = 2, MAXSTF = 1
integer kdim, maxord, maxstf
!++ Substitute for KDIM, MAXORD, MAXSTF below
parameter (kdim = 20, maxord = 2, maxstf = 1)
!--D Next line special: P=>D, X=>Q
double precision tn
double precision xi(kdim)
!
!--D Next line special: P=>D, X=>Q
double precision tg(2), tgstop(2), tmark, tmarkx, tout, tolg
double precision alpha(kdim), beta(kdim+1)
double precision  d(maxstf+maxord,maxord), g(kdim,maxord)
double precision v(kdim+maxord)
double precision hc, hdec, hinc, hincc, hmax, hmaxp9, hmin
double precision fdat(11)
!
double precision ds(maxstf+maxord, maxord), gs(kdim)
double precision sigma(kdim), rbq(kdim), dnoise
double precision eave, eimax, eimin, emax, erep, robnd, snoise
!
!.    SPECIFICATION OF ENVIRONMENTAL CONSTANTS.
double precision eeps10, eeps16, erov10, eeps2
double precision eept75, eovep2, ovtm75, ovd10
common / divaev / eeps2, eept75, eovep2, ovtm75, ovd10, eeps10, &
   eeps16, erov10
save / divaev /
integer iopst, kordi, kqmaxd, kqmaxi, ldt, maxdif, maxint, nkdko, &
   nte, nyny, ndtf, numdt
common / divasc / tn, xi, iopst, kordi, kqmaxd, kqmaxi, ldt, &
   maxdif, maxint, nkdko, nte, nyny, ndtf, numdt
!
integer icf,ics,igflg,igtype(2),igstop(2),ilgrep,ings,iop3,iop4, &
   iop5,iop6,iop7,iop8,iop9,iop10,iop11,iop12,iop13,iop14,iop15, &
   iop16,iop17,iop18,iop19,iop20,iop21,iop22,iop21s,itolep,iy, &
   kemax,kis,kmark,kord1i,kord2i,kpred,kqdcon,kqicon,kqmaxs, &
   kqmxds,kqmxil,kqmxip,kqmxis,ksc,ksout,ksstrt,kstep,lex,linc, &
   lincd,lincq,lsc,maxkqd,maxkqi,method,ne,neptol,ng,ngtot, &
   noiseq,noutko,ntolf,ny,idat(6)
common /divamc/ tg,tgstop,tmark,tmarkx,tout,tolg,hc,hdec,hinc, &
   hincc,hmax,hmaxp9,hmin,alpha,beta,d,g,v,ds,gs,sigma,rbq,dnoise, &
   eave,eimax,eimin,emax,erep,robnd,snoise,fdat,icf,ics,igflg, &
   igtype,igstop,ilgrep,ings,iop3,iop4,iop5,iop6,iop7,iop8,iop9, &
   iop10,iop11,iop12,iop13,iop14,iop15,iop16,iop17,iop18,iop19, &
   iop20,iop21,iop22,iop21s,itolep,iy,kemax,kis,kmark,kord1i, &
   kord2i,kpred,kqdcon,kqicon,kqmaxs,kqmxds,kqmxil,kqmxip,kqmxis, &
   ksc,ksout,ksstrt,kstep,lex,linc,lincd,lincq,lsc,maxkqd,maxkqi, &
   method,ne,neptol,ng,ngtot,noiseq,noutko,ntolf,ny,idat
save / divamc / , / divasc /
!
integer i, ig, igflgs, izflag, kexit, ngstop(2)
double precision hh
equivalence (g(1,1), hh), (ngstop(1), iop6), (kexit, iop17), &
   (izflag, iy), (igflgs, itolep)
!
!                      Declarations for error message processing.
!
integer memda1, meret, meemes
parameter (memda1 =27)
parameter (meret  =51)
parameter (meemes =52)
integer mact(7)
! ********* Error message text ***************
![Last 2 letters of Param. name]  [Text generating message.]
!AA DIVAG$B
!AB Call with bad values of KORD.  KORD(1)=$I, KORD(2)=$I, when $C
!   TSPECS(1)=$F and KSTEP=$M.$E
integer ltxtaa,ltxtab
parameter (ltxtaa= 1,ltxtab= 8)
character mtxtaa(1) * (95)
data mtxtaa/'DIVAG$BCall with bad values of KORD.  KORD(1)=$I, KOR &
d(2)=$i, when tspecs(1)=$f and kstep=$m.$e'/
data mact / memda1, 0, meemes, 68, 24, 0, meret /
!
! ****************** START OF EXECUTABLE CODE **************
!
iflag = 1
ig = kord(2)
if ((ig .ne. 0) .and. (ig .ne. 1)) go to 500
if (igflg - 3) 10, 80, 70
10 if (igflg - 1) 20, 210, 60
!
! ******************** INITIAL POINT ***********************
!
20 if (igflg .eq. -3) then
   igflg = 5
   return
end if
igtype(ig + 1) = 0
if ((ngstop(ig + 1) .le. 0) .or. (ig + igflg .eq. -1)) go to 30
igflg = ig - 2
go to 40
30 igflg = 5
40 ng = ngstop(2 - ig)
do 50 i = 1, ng
50    gt(i) = gnew(i)
go to 480
!
!     **** USER HAS BEEN TOLD THAT A GSTOP WAS FOUND
!          TEST IF CALLED FROM OUTPUT WHEN CALL SHOULD
!          BE FROM DERIVS
60 if ((igtype(1) .eq. 0) .and. (ig .ne. 0)) go to 420
!     **** PROTECT AGAINST NOISEY G NEAR THE ZERO
if (gnew(ings) .eq. 0.d0) gnew(ings) = gt(ings)
!
! ****** TEST FOR CHANGE IN THE SIGN OF A G ****************
!
70 ng = ngstop(2 - ig)
ings = 0
80 ings = ings + 1
if (ings .gt. ng) if (igflg - 4) 400, 380, 480
if (gnew(ings)) 90, 100, 110
90 if (gt(ings)) 120, 120, 130
100 if (gt(ings)) 130, 80, 130
110 if (gt(ings)) 130, 120, 120
120 gt(ings) = gnew(ings)
go to 80
!
! ********* A SIGN CHANGE HAS BEEN FOUND *******************
!
130 nstop = ings
if (ig .eq. 0) nstop = -ings
if (igflg .ne. 5) go to 200
!     **** USUAL CASE -- TEST IF OUTPUT POINT PRECEDES THE
!          SIGN CHANGE, OR IF PREDICTING, CORRECTING, OR
!          NOT THROUGH THE FIRST STEP.
!     **** TEST IF AN INTERPOLATED G WAS WHAT CHANGED SIGN
if (ig .ne. 0) go to 180
!     **** BACK UP DIFFERENCES AND OTHER STUFF TO BEGINNING
!          OF THE STEP
call divabu(f, kord)
!     **** TEST IF CORRECTING
if (kord1i .eq. 2) go to 170
!     **** TEST IF THROUGH THE FIRST STEP
if (lsc .lt. 4) go to 180
!     **** IF FIRST DERIVATIVE EVALUATION OF THE FIRST
!          STEP, FIND THE GSTOP, AND USE IT TO GET A NEW
!          INITIAL STEPSIZE
if (lsc .eq. 7) go to 200
!     **** SET NEW STEPSIZE AFTER SIGN CHANGE WHEN STARTING
160 hh = tspecs(1) - tn
!     **** SET KEXIT TO TRY NEW STEPSIZE
170 kexit = 1
tspecs(1) = tn
go to 460
!     **** SET KEXIT FOR USUAL CASE
180 kexit = ig + 2
!     **** TEST IF SIGN CHANGE IN G PRECEDES NEXT OUTPUT PT.
if (hh * (tspecs(1) - tmark)) 200, 200, 190
!     **** SET UP TO EVALUATE G AT OUTPUT POINT
190 igflg = 4
tspecs(1) = tmark
nstop = 0
go to 240
!
! ***************** FIND THE ZERO OF G *********************
!
!     **** INITIALIZE ZERO FINDER
200 told = tg(2 - ig)
gold = gt(ings)
tsave = tspecs(1)
igflgs = igflg
igflg = 1
izflag = 0
go to 220
!     **** TEST IF ZERO ALREADY FOUND
210 if (izflag - 1) 350, 220, 310
220 continue
call dzero(tspecs(1), gnew(ings), told, gold, izflag, tolg)
!     **** TEST FOR CONVERGENCE
if (izflag .ne. 1) go to 260
!     **** INTERPOLATE NEW Y, AND GO COMPUTE G AGAIN
240 call divain(tspecs(1), y, f, kord)
iflag = 4
kord2i = ig - 3
return
!     **** CONVERGENCE -- CHOOSE TOLD TO GIVE A CHANGE
!          IN SIGN
260 if (gnew(ings) .eq. 0.d0) go to 290
if (tspecs(1) - told) 270, 300, 280
270 if (hh) 290, 300, 300
280 if (hh) 300, 300, 290
290 told = tspecs(1)
!     **** CHECK IF SIGN CHANGE DUE TO NOISE
300 tspecs(1) = told + xi(1)
go to 240
310 tspecs(1) = told
if (gnew(ings)) 320, 340, 330
320 if (gt(ings)) 340, 360, 360
330 if (gt(ings)) 360, 360, 340
!     **** ZERO WAS EVIDENTLY DUE TO NOISE
340 tspecs(1) = tsave
izflag = 0
go to 370
350 igflg = igflgs
!     SET KORD2I TO INITIAL VALUE TO AVOID LOOP
kord2i = ig
go to 80
!     **** SAVE INFORMATION ABOUT THE STOP
360 igflg = 3
tgstop(2 - ig) = tspecs(1)
igtype(2 - ig) = izflag + 3
igstop(2 - ig) = nstop
370 nstop = 0
go to 240
!
! ************** AFTER SEARCH FOR A SIGN CHANGE ************
!
!     **** NO SIGN CHANGE AT A T OUTPUT POINT
!     TEST IF CALLED FROM OUTPUT
380 if (ig .ne. 0) go to 390
!     SET UP FOR CALL TO OUTPUT
kord1i = 7
kord2i = -2
iflag = 3
go to 480
!     **** ADJUST KEXIT AND SET UP TO GIVE OUTPUT
390 kexit = kexit + 2
kord1i = min(kmark, 5)
kord(3) = kmark
kord(1) = kord1i
igflg = 5
iflag = 2
go to 470
!     **** TEST IF USER HAS BEEN TOLD OF GSTOP
400 if (igflg .eq. 2) go to 450
!     **** A GSTOP HAS BEEN FOUND
!     TEST IF STARTING
if (lsc .eq. 7) go to 160
iflag = igtype(2 - ig)
nstop = igstop(2 - ig)
ings = abs(nstop)
if (ings .eq. 0) go to 410
gt(ings) = -gt(ings)
if (ig .eq. 0) go to 430
igflg = 2
!     If interpolated GSTOP was found set to check again in case of
!     multiple stops at exactly the same point.
if (igtype(1) .ne. 0) go to 440
!     **** TELL USER OF AN EXTRAPOLATED GSTOP
410 iflag = igtype(2)
nstop = igstop(2)
ings = abs(nstop)
!     **** SET SO DERIVS IS CALLED WITH KORD(1) = KPRED
420 kord1i = kpred
kord2i = -3
return
!     **** AN EXTRAPOLATED GSTOP WAS FOUND, SET UP TO CHECK
!          INTERPOLATED STOPS (IF ANY)
430 ng = ngstop(1)
ings = 0
iflag = 3
nstop = igstop(2)
!     **** SET SO OUTPUT IS CALLED WITH KORD(1) = 7
440 kord1i = 7
kord2i = -2
go to 490
!     **** CHECK THAT AN EXTRAPOLATED G-STOP IS NOT MISSED
450 if ((ig .eq. 0) .or. (igtype(2) .eq. 0)) go to 460
!     SET TO CHECK FOR INTERPOLATED G-S.
tg(1) = tspecs(1)
igtype(1) = 0
tspecs(1) = tgstop(2)
ings = 0
igflg = 3
go to 240
!     **** SET SO INTEGRATOR GOES TO PLACE DIRECTED BY KEXIT
460 nstop = 0
igflg = 5
iflag = 3
470 kord2i = -7
!     **** STORE INFO. ON LAST G COMPUTED
480 igtype(2 - ig) = 0
490 tg(2 - ig) = tspecs(1)
return
!
! ********************** ERROR PROCESSING ******************
!
500 mact(2) = kstep
!--D Next line special: P=>S, X=>D
call dmess(mact, mtxtaa, kord, tspecs)
iflag = 8
kexit = 6
go to 470
end
subroutine dmess (mact, text, idat, fdat)
! Copyright (c) 1996 California Institute of Technology, Pasadena, CA.
! ALL RIGHTS RESERVED.
! Based on Government Sponsored Research NAS7-03001.
!>> 2009-09-27 DMESS Krogh  Same as below, in another place.
!>> 2009-07-23 DMESS Krogh  Changed ,1x to :1x in write to FMTF.
!>> 2008-06-13 DMESS Krogh  Changed -0's to 0.
!>> 2007-09-08 DMESS Krogh  Fixed definitions of MEVLAS.
!>> 2006-10-08 DMESS Krogh  Another try, see 2005-05-26
!>> 2006-10-08 DMESS Krogh  Fixed minor problem in matrix/vector output.
!>> 2006-10-01 DMESS Krogh  Print NaN's and infity (at least for g77).
!>> 2006-07-01 DMESS Krogh  messxc => dmessxc (and not static) (for C)
!>> 2006-04-07 DMESS Krogh  Major rewrite of code for F.Pt. format.
!>> 2006-04-04 DMESS Krogh  Fixes in C code for vectors & matrices.
!>> 2006-04-02 DMESS Krogh  Added code for output of sparse vector.
!>> 2005-07-10 DMESS Krogh  Small adjustment for last correction.
!>> 2005-05-26 DMESS Krogh  Fixed "*****" output in boundary case.
!>> 2002-05-16 DMESS Krogh  Added way for user to get error count.
!>> 2002-03-27 DMESS Krogh  Fixed crash when number is -INF.
!>> 2001-06-08 DMESS Krogh  Eliminated Hollerith in formats.
!>> 2001-05-25 DMESS Krogh  Added a couple of commas in formats.
!>> 1997-06-17 DMESS Krogh  In C code made messxc, static.
!>> 1996-07-12 DMESS Krogh  Changes to use .C. and C%%.
!>> 1996-03-30 DMESS Krogh  Added external statement.
!>> 1994-10-20 DMESS Krogh  Changes to use M77CON
!>> 1994-09-21 DMESS Krogh  Added CHGTYP code.
!>> 1994-09-08 DMESS Krogh  Added new matrix/vector capabilities.
!>> 1994-08-17 DMESS Krogh  Removed duplicate save statement.
!>> 1994-04-19 DMESS Krogh  Removed blank line from DMESS.
!>> 1993-05-14 DMESS Krogh  Changed TEXT to array of character strings.
!>> 1993-04-14 DMESS Krogh  Fixes for conversion to C. (C%% comments.)
!>> 1992-07-12 DMESS Krogh  Fixed so negative KDFDEF works.
!>> 1992-05-27 DMESS Krogh  Initialized LDFDEF in a data statement.
!>> 1992-05-14 DMESS Krogh  Put common blocks in save statement.
!>> 1992-04-28 DMESS Krogh  Corrected minor error in floating pt. format
!>> 1992-02-28 DMESS Krogh  Initial Code.
!
!--D replaces "?": ?MESS,?MESSXC
!
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
external         d1mach
integer          mact(*), idat(*)
double precision fdat(*)
character        text(*)*(*)
character        fmtsp*29
integer          icol, id, j, k, ksma, kbig, kexe, ldfdef, neg
double precision fbig, fout, fsma, d1mach
save ldfdef, fmtsp
save /cmessi/, /cmessc/
!++ CODE for .C. is inactive
!      integer  kciwid, kccwid, kcrwid, lbeg, lend, lfprec, lgprec
!      common /MESSCC/ kciwid,kccwid,kcrwid,lbeg,lend,lfprec,lgprec
!++ END
!
! ************************** Data from common block ********************
!
! For comments on these variables, see the listing for MESS.
!
integer   lenbuf, mevbas, mevlas
parameter (lenbuf=250)
parameter (mevbas=10)
parameter (mevlas=33)
logical          xarg, gotfmt, xargok
integer          errcnt, eunit, ichar0, irc, ivar(mevbas:mevlas), &
   imag, inc, itext, iwf, iwg, kdf, kdfdef, kdi, kdiag, kdj, &
   kline, kscrn, kshift, kspec, kt, maxerr, lasti, lbuf, lenlin, &
   lenout, lentry, lentxt, lhead, linerr, linmsg, locbeg, lprint, &
   lstop, lstrt, ltext, maxwid(2), mdat(5), mpt, munit, ncol, &
   ndim, nfdat, nidat, nmdat, nrow, ntext, ounit, sunit, tabspa
!
character buf*(lenbuf), dols*72, fmtc*7, fmtf*20, fmtg*15, &
  fmti*7, fmtim(2)*7, fmtj*7, fmtr*7, fmtt*15
common /cmessi/ sunit, lhead, kdfdef, linmsg, linerr, munit, &
   eunit, kscrn, kdiag, maxerr, lstop, lprint, kdf, ntext, nidat, &
   nfdat, nmdat, mdat, tabspa, errcnt, ichar0, imag, inc, irc, &
   itext, iwf, iwg, kdi, kdj, kline, kshift, kspec, kt, lasti, &
   lbuf, lenlin, lenout, lentry, lentxt, locbeg, lstrt, ltext, &
   maxwid, mpt, nrow, ncol, ndim, ounit, gotfmt, xarg, xargok
common /cmessc / buf, dols, fmtf, fmtg, fmti, fmtj, fmtt, fmtim
equivalence (ivar(mevbas), sunit)
equivalence (fmtim(1), fmtr), (fmtim(2), fmtc)
!
data ldfdef / 0 /
!

! ************************* Start of Executable Code *******************
!
xargok = .true.
if (ldfdef .eq. 0) then
   ldfdef = 1 - int(log10(d1mach(3)))
end if
kdfdef = ldfdef
kdf = kdfdef
10 call mess (mact, text, idat)
!             4    5    6    7    8    9   10   11
go to (20, 100, 100, 200, 300, 400, 100, 500), lentry-3
xargok = .false.
ldfdef = kdfdef
return
!                                      Print from FDAT
20 j = lbuf + 1
fout = fdat(nfdat)
nfdat = nfdat + 1
if (kspec .ge. 8) then
   lbuf = lbuf + iwg
!%% messcc.lend = cmessi.lbuf;
!%% cmessc.buf[messcc.lend] = ' ';
!%% if ((j > 1) && (cmessc.buf[j-2] >= '0') &&
!%%    (cmessc.buf[j-2] <= '9')) cmessc.buf[j++ - 1] = ' ';
!%% sprintf(&cmessc.buf[j-1], cmessc.fmtg, cmessi.iwg,
!%%    messcc.lgprec, fout);
!%% if (cmessc.buf[messcc.lend] != 0) {messcc.lbeg=j; dmessxc(kexe);}
   write (buf(j:lbuf), fmtg) fout
   go to 10
end if
if (fout .le. 0.d0) then
  if (fout .eq. 0.d0) then
    fdat(nfdat-1) = 0.d0
    fout = 0.d0
  else
    neg = 1
  end if
else if (fout .gt. 0.d0) then
  neg = 0
else
!               Must be a Nan
  neg = 0
  fbig = 1.0
  fsma = 1.0
  iwf = 2
  go to 40
end if

fbig = abs(fout)
fsma = fbig
iwf = 2
!                                      Get the format.
40 continue
if (kdf .eq. 0) kdf = kdfdef
kexe = 0
if (fbig .ne. 0.d0) then
  if (fsma .eq. 0.d0) fsma = 1.d0
  fbig = fbig * (1.d0 + .5d0 * .1d0**abs(kdf))
  iwf = iwf + neg
  if (kdf .lt. 0) then
    ksma = 0
  else
    ksma = -log10(fsma)
    if (fsma .lt. 1.d0) ksma = ksma + 1
  end if
  kbig = log10(fbig)
  if (fbig .lt. 1.d0) then
    kbig = kbig - 1
    if (fbig .gt. 1.d0 - .1d0**abs(kdf-1)) kbig = 0
  end if
!         This is to get ininities out (at least with g77)
  if ((kbig .lt. -1000) .or. (kbig .gt. 1000)) kbig = 8
  if ((ksma .lt. -1000) .or. (ksma .gt. 1000)) ksma = 8
  if (max(kbig, 0) + max(ksma,0) .ge. 4) then
!               Want to use an "E" format
    kexe = 3 + max(0, int(log10(dble(max(kbig,abs(ksma))+1.d-5))))
    if (kdf .lt. 0) then
      id = -kdf
    else
      id = kdf - 1
    end if
    iwf = iwf + id + kexe
!++ CODE for ~.C. is active
    if (lentry .eq. 10) iwf = iwf - 1
    write (fmtf, '(''(1P,99(E'',I2,''.'',I2,''E'',I1,'':1X))'')') &
      iwf, id, kexe - 2
!++ CODE for .C. is inactive
!c WATCOM C and others (?) man need an extra 1 here??
!%%    strcpy(cmessc.fmtf, "%*.*E ");
!      lfprec = id
!++ END
    go to 60
  end if
else
  ksma = 1
  kbig = -1
end if
!               Want to use an "F" format
if (kdf .lt. 0) then
  id = -kdf
else
  id = kdf + ksma - 1
end if
!++ CODE for ~.C. is active
iwf = iwf + id + max(kbig, -1)
write (fmtf, '(''(0P,99(F'',I2,''.'',I2,'':1X))'')') iwf,id
!++ CODE for .C. is inactive
!      IWF = IWF + ID + max(KBIG, 0)
!%%    strcpy(cmessc.fmtf, "%*.*f ");
!      lfprec = id
!++ END
60 if (lentry .ne. 4) then
  iwf = iwf + 1
  if (lentry .ne. 10) go to 10
!               Take care of setup for sparse vector
  imag = 0
  do 70 j = locbeg, lasti
    imag = max(abs(imag), idat(j))
70   continue
  call messfi
!  Format forms:     12345678901234567890   123456789012345678  1234567
!                    (1P,99(Edd.ddEd:1X))   (0P,99(Fxx.xx:1X))  (99Idd)
!++ CODE for ~.C. is active
  if (fmtf(8:8) .eq. 'F') then
    fmtsp= &
      '(99(' // fmti(4:6) // ''') '',0P,' // fmtf(8:18)
  else
    fmtsp= &
     '(99(' // fmti(4:6) // ''')'' ,1P,' // fmtf(8:20)
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
  iwf = iwf + kdi + 1
  go to 10
end if
!
lbuf = lbuf + iwf
!%% messcc.lend = cmessi.lbuf;
!%% cmessc.buf[messcc.lend] = ' ';
!%% if ((j > 1) && (cmessc.buf[j-2] >= '0') &&
!%%    (cmessc.buf[j-2] <= '9')) cmessc.buf[j++ - 1] = ' ';
!%% sprintf(&cmessc.buf[j-1], cmessc.fmtf, cmessi.iwf,
!%%    messcc.lfprec, fout);
!%% if (cmessc.buf[messcc.lend] != 0) {messcc.lbeg=j; dmessxc(kexe);}
write (buf(j:lbuf),fmtf) fout
go to 10
!                                     Get format for a vector or matrix
100 icol = 1
if (fdat(locbeg) .lt. 0.d0) then
  neg = 1
else if (fdat(locbeg) .ge. 0.d0) then
  neg = 0
else
!               Must be a Nan
  neg = 0
  fbig = 1.0
  fsma = 1.0
  go to 110
end if

fbig = abs(fdat(locbeg))
fsma = fbig
110 do 120 j = locbeg, lasti, inc
  if (fdat(j) .le. 0.d0) then
    if (fdat(j) .eq. 0.d0) then
      fdat(j) = 0.d0
    else
      neg = 1
    end if
  end if
  fbig = max(fbig, abs(fdat(j)))
  if (fsma .eq. 0.d0) then
    fsma = abs(fdat(j))
  else if (fdat(j) .ne. 0.d0) then
    fsma = min(fsma, abs(fdat(j)))
  end if
120 continue
if (ncol .ne. 0) then
   icol = icol + 1
   locbeg = locbeg + ndim
   lasti = lasti + ndim
   if (icol .le. ncol) go to 110
end if
iwf = 2
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
write(buf(lstrt:lbuf),fmtf)(fdat(k),k=mpt,mpt+inc*(kline-1),inc)


!      print '(/A/)', BUF(1:LBUF)

mpt = mpt + kline * inc
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
write (buf(lstrt:lbuf), fmtf) (fdat(k), k = mpt, lasti, ndim)
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
write (buf(lstrt:lbuf), fmtf) (fdat(k), k = mpt, mpt+kline-1)
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
write (buf(1:lbuf), fmtsp) (idat(k),fdat(k),k=mpt,mpt+kline-1)
mpt = mpt + kline
go to 10

end

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
subroutine dzero(x1, f1, x2, f2, mode, tol)
! Copyright (c) 1996 California Institute of Technology, Pasadena, CA.
! ALL RIGHTS RESERVED.
! Based on Government Sponsored Research NAS7-03001.
!>> 2010-04-14 DZERO  Krogh  No discontinuity message if |F1-F2| small.
!>> 2010-04-12 DZERO  Krogh  Fixed KNKP to get discontinuity diagnostic.
!>> 2010-02-20 DZERO  Krogh  $G => $F for print out of iterations
!>> 2008-03-01 DZERO  Krogh  Minor change in diagnostic print.
!>> 2000-12-01 DZERO  Krogh  Removed unused variable C1P01.
!>> 1998-11-01 DZERO  Krogh  Set so errors stop less easily.
!>> 1998-11-01 DZERO  Krogh  For "mangle", INDIC replaced with MACT(3).
!>> 1996-03-30 DZERO  Krogh  Added external statement.
!>> 1995-11-09 DZERO  Krogh  Fixed so char. data at col. 72 is not ' '.
!>> 1994-11-11 DZERO  Krogh  Declared all vars.
!>> 1994-10-20 DZERO  Krogh  Changes to use M77CON
!>> 1994-09-08 DZERO  Krogh  Added CHGTYP code.
!>> 1993-04-27 DZERO  Krogh  Additions for Conversion to C.
!>> 1993-04-13 DZERO  Krogh  Minor change for new MESS.
!>> 1992-04-08 DZERO  Krogh  Unused label 400 removed.
!>> 1992-01-09 DZERO  Krogh  Moved calc. of XXMXO up (for error msg.)
!>> 1991-11-26 DZERO  Krogh  Converted to new error processor.
!>> 1988-08-14 DZERO  Krogh  Labels runumbered.
!>> 1988-03-07 DZERO  Krogh  Initial code.
!
!--D replaces "?": ?ZERO, ?MESS
!
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
!        F1*F2 .le. 0.
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
!      =6  fatal error -- $ZERO was called after mode was set .ge.2.
!          If $ZERO is called again, the program will be stopped.
!          (Unless MODE is set to 0)
!      <0  If MODE is set <0 and $ZERO is called, no action is taken
!          except that print is turned on for -MODE calls to $ZERO.
!          This print gives all values of X and F used in the iteration.
!  TOL    is the error tolerance
!     TOL.GT.0  Iterate until values of X1 and X2 are known
!              for which abs(X1-X2) .le. tol and F1*F2 .le. 0.
!     TOL.LT.0  Iterate until a value of X1 is found for which
!              abs(F1) .le. abs(TOL).
!     TOL  = 0  Iterate until the zero is determined as
!              precisely as possible.  MODE = 3 is impossible
!              in this case.
!
! Parameters in the calling sequence have the following types
!
integer mode
double precision x1, x2, f1, f2, tol
!
! Usage is as follows (of course, variations are possible.)
!         Somehow one has available X1, F1, X2, and F2 such
!         that F1 = F(X1), F2 = F(X2) and F1*F2 .le. 0.
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
! J      This is 1 if FF .le. 0., and is 2 if FF > 0.
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
external d1mach
integer linit, ks, ktyp, j, i
parameter (linit = -40)
integer knkp(2), lchg(30), lmode, lnlp(2), np
double precision xx, xo, xl, ff, fo, fl, ffdfo
double precision div, qfm, qxm, tp, tp1, xxmxo, xxmxol
double precision rnd, xrnd, small, tolx
double precision xxmxl, xlmxb, ffmfl, ffmfb, flmfb
double precision dxdfff, dxdffo, dfdxxx, dfdxxo
double precision c0, c1, c2, c4, cp125, cp25, cp5, cp75, c1p25
double precision c8, cp01, cp99, cp001, c1p031
double precision d1mach
!
parameter (c0 = 0.d0, c1 = 1.d0, c2 = 2.d0, c4 = 4.d0)
parameter (c8 = 8.d0)
parameter (cp125 = 0.125d0, cp25 = 0.25d0, cp75 = 0.75d0)
parameter (cp5 = 0.5d0)
parameter (c1p25 = 1.25d0)
parameter (cp01 = 0.01d0)
parameter (cp001 = 0.001d0)
parameter (cp99 = 0.99d0)
parameter (c1p031 = 1.03125d0)
!
!                      Declarations for error message processing.
!
integer meret, meemes, metext
double precision fdat(4)
integer mact(5), mact1(2), mloc(4), idat(2)
save div, fl, flmfb, fo, knkp, ks, ktyp, lchg, lmode, &
   lnlp, mact, np, rnd, small, xl, xlmxb, xo, xx, xxmxol
parameter (meret  =51)
parameter (meemes =52)
parameter (metext =53)
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
integer ltxtaa,ltxtab,ltxtac,ltxtad,ltxtae,ltxtaf,ltxtag
parameter (ltxtaa=  1,ltxtab=  8,ltxtac= 63,ltxtad=112,ltxtae=169, &
 ltxtaf=  1,ltxtag=  1)
character mtxtaa(1) * (193)
character mtxtab(1) * (46)
character mtxtac(1) * (25)
data mtxtaa/'DZERO$BBest bound for zero is [$F, $F], but tolerance &
 is $f.$EApparent discontinuity in function near x = $f.$ECan not$ &
 find a sign change: x1=$f, x2=$f, f1=$f, f2=$f$ECalled with mode$ &
 = $i.$e'/
data mtxtab/'In DZERO -- X1=$F F1=$F KTYP=$I DIV=$F KS=$I$E'/
data mtxtac/'            X2=$F F2=$F$E'/
! **** End of automatically generated text
!                      1  2  3  4      5
data mact / meemes, 0, 0, 0, meret /
data mact1 / metext, meret /
data mloc / ltxtab, ltxtac, ltxtad, ltxtae /
!
data rnd / c0 /
data ks, ktyp, lmode, div / 0, 0, 2, c0 /
data lchg / &
   1, 2, 0, 0, 0, &
   1, 1, 4, 5, 0, &
   1, 2, 3, 3, 0, &
   1, 4, 4, 3, 0, &
   1, 4, 5, 5, 0, &
   1, 5, 5, 5, 0 /
data np / 0 /

!
! INITIALIZE
!
if (mode .lt. 0) then
   np = -1 - mode
   return
end if
if (np .gt. 0) then
   np = np - 1
   fdat(1) = x1
   fdat(2) = f1
   fdat(3) = div
   idat(1) = ktyp
   idat(2) = ks
   call dmess(mact1, mtxtab, idat, fdat)
   if (mode .ne. 0) if (lmode - 1) 70, 80, 450
   fdat(1) = x2
   fdat(2) = f2
   call dmess(mact1, mtxtac, idat, fdat)
else if (mode .ne. 0) then
   if (lmode - 1) 70, 80, 450
end if
!
if (rnd .eq. c0) then
   rnd = d1mach(4)
   small = cp5 / (rnd * d1mach(2))
end if
xl = x2
fl = f2
30 tp = c1
mode = 1
mact(3) = 2
xxmxol = c0
knkp(1) = 0
knkp(2) = 0
lnlp(1) = 0
lnlp(2) = 0
ks = -1
xx = x1
ff = f1
if (fl) 40, 75, 50
40 if (ff) 60, 230, 100
50 if (ff) 100, 230, 60
60 lmode = 0
!             Take care of points on same side of zero.
70 ff = f1
xx = x1
tp = ff / fl
if (tp .lt. c0) go to 30
lmode = lmode - 1
if (lmode .lt. linit) then
   mact(3) = 5
   fdat(1) = xx
   fdat(2) = xl
   fdat(3) = ff
   fdat(4) = fl
   go to 250
end if
if (tp .gt. c1) then
   ff = fl
   xx = xl
   fl = f1
   xl = x1
end if
if (abs(ff) .ge. c8 * abs(fl-ff)) then
   tp = c8
else
   tp = max(-cp25*dble(lmode), ff / (fl - ff))
end if
fl = ff
xo = xl
xl = xx
if (xx .eq. xo) xo = c1p031 * xx + sign(cp001, xx)
xx = xx + tp * (xx - xo)
x1 = xx
mode = 1
return
!
75 x1 = xl
f1 = fl
go to 250
! END OF INITIALIZATION
!
!
! ENTRY AFTER COMPUTING F FOR THE LAST ITERATE
80 ff = f1
tp = ff / fl
if (tp) 90, 230, 110
90 tp = ff / fo
ks = 0
100 fo = fl
xo = xl
go to 120
110 ks = ks + 1
120 j = 1
if (ff .gt. c0) j = 2
if (tp - c1) 150, 140, 130
130 i = 4
if (tp .gt. c4) i = 5
go to 160
140 i = 3
go to 160
150 i = 2
if (tp .lt. cp01) i = 1
if (tp .lt. cp99) go to 170
160 knkp(j) = knkp(j) + 1
go to 180
170 knkp(j) = 0
180 xxmxo = xx - xo
lnlp(j) = lchg(5*lnlp(j) + i)
if (lnlp(j) .ge. 4) then
   if (lnlp(3 - j) .ge. 4) go to 210
end if
! XXMXO GIVES THE LENGTH OF THE INTERVAL INSIDE WHICH
! THE ZERO IS KNOWN TO LIE.
if (c2 * abs(xxmxo) .lt. abs(xxmxol)) then
   knkp(j) = max(0, knkp(1) - 3)
end if
xxmxol = xxmxo
xrnd = rnd * (abs(xx) + abs(xo) + small)
!
! TEST FOR CONVERGENCE
if (tol) 190, 200, 200
190 continue
if (abs(ff) .le. abs(tol)) go to 220
200 continue
tolx = max(tol, xrnd)
if (abs(xxmxo) .gt. tolx) go to 310
!
! CONVERGENCE -- PREPARE FOR FINAL EXIT
210 if ((abs(xxmxo) .gt. tol) .and. (tol .ne. c0)) then
   mact(3) = 3
   fdat(3) = tol
   if (xxmxo .gt. 0) then
      fdat(2) = xx
      fdat(1) = xo
   else
      fdat(1) = xx
      fdat(2) = xo
   end if
end if
! SET FINAL VALUES FOR X1,F1,X2,AND F2
220 continue
if (abs(ff) .le. abs(fo)) go to 240
f1 = fo
x1 = xo
230 fo = ff
xo = xx
240 x2 = xo
f2 = fo
! TEST FOR DISCONTINUITY
if ((knkp(1) .gt. 5) .or. (knkp(2) .gt. 5)) then
  if (abs(f1 - f2) .gt. rnd * max(x1, 1.d0)) then
    mact(3) = 4
    fdat(1) = xx
  end if
end if
250 mode = mact(3)
if (mact(3) - 2) 420, 420, 430
! END OF CODE FOR FINAL EXIT
!
! F NOT DECREASING (OR THE FIRST ITERATE)
! PREPARE TO DIVIDE THE INTERVAL
260 tp = c1
if (ks) 370, 280, 270
270 if (ktyp .eq. 0) go to 290
280 div = c2
290 continue
div = max(div, ffdfo)
! KTYP=0 IF AND ONLY IF THE INTERVAL WAS DIVIDED (USING DIV)
! ON THE LAST ITERATION
if (ktyp .eq. 0) div = div * (c1p25 / (c1p25 - tp))
! DIVIDE THE INTERVAL AS SPECIFIED BY DIV
300 tp1 = -xxmxo * (div/(div+c1))
ktyp = 0
go to 410
!
310 continue
xxmxl = xx - xl
ffmfl = ff - fl
ffdfo = abs(ff / fo)
tolx = cp5 * tolx
if (tp .ge. c1) go to 260
! DIVIDE THE INTERVAL IF F HAS HAD THE SAME SIGN FOR
! FOUR OR MORE TIMES IN SUCCESSION
if (ks - 4) 320, 340, 290
320 continue
if (flmfb .eq. c0) go to 340
! BEGINNING OF CODE TO DETERMINE IF INVERSE QUADRATIC
! INTERPOLATION IS TO BE USED.
ffmfb = ffmfl + flmfb
if (ffmfb .eq. c0) go to 330
qfm = c1 - (ffmfl / flmfb) * (xlmxb / xxmxl)
qxm = c1 - (xxmxl / xlmxb) * (flmfb / ffmfl)
dxdfff = c1 + (ffmfl / ffmfb) * qfm
dxdffo = dxdfff + c2 * ((fo - ff) / ffmfb) * qfm
tp1 = xxmxl + xlmxb
dfdxxx = c1 + (xxmxl / tp1) * qxm
dfdxxo = dfdxxx + c2 * ((xo - xx) / tp1) * qxm
tp1 = dxdfff * dfdxxx
if ((tp1 .le. cp25) .or. (tp1 .ge. c4)) go to 330
tp1 = dxdffo * dfdxxo
if ((tp1 .gt. cp25) .and. (tp1 .lt. c4)) go to 380
!
! DERIVATIVES DO NOT MATCH WELL ENOUGH
330 continue
if (ks .eq. 0) if (ffdfo - c1) 350, 370, 360
340 continue
if ((ktyp .eq. 0) .and. (tp .ge. cp75)) go to 290
continue
tp = c1 - tp
if (tp .le. ffdfo) go to 280
ffdfo = ffdfo / tp
div = cp125
go to 290
350 continue
div = cp5 * max(max(cp25, ffdfo), tp / (c1p25 - min(tp, c1)))
go to 300
360 continue
div = min(c4, cp5 * ffdfo)
go to 300
! INTERPOLATE WITH SECANT METHOD
370 tp1 = -xxmxl
go to 390
!
! DERIVATIVES MATCH UP PRETTY WELL.
380 continue
! INTERPOLATE USING THE INVERSE QUADRATIC
tp1 = xxmxl * (qfm * (fl / ffmfb) - c1)
390 tp1 = (ff/ffmfl) * tp1
ktyp = 1
!
! EXIT TO GET F(X)
410 continue
fl = ff
flmfb = ffmfl
xlmxb = xxmxl
xl = xx
! COMPUTE X1, INSURING THAT IT IS NOT TOO CLOSE TO THE
! ENDS OF THE INTERVAL
xx = min(max(xl + tp1, min(xl, xo) + tolx), max(xl, xo) - tolx)
x1 = xx
420 lmode = mode
return
!
430 mact(2) = 11*mact(3)  - 19
440 mact(4) = mloc(mact(3)-2)
call dmess(mact, mtxtaa, idat, fdat)
go to 420
!
! A CALL TO THE SUBROUTINE HAS BEEN MADE WITH MODE.NE.1
450 idat(1) = mode
mact(3) = 6
mode = 6
if (lmode .ne. 6) go to 430
mact(2) = 99
go to 440
end
!++ CODE for .C. is inactive
!%% static FILE *c_handle[2], *scratch_file;
!%% static char *c_fname[2]={"MESSF-xx", "MESSF-xx"};
!%% char *ctmp;
!++ END
subroutine mess(mact, text, idat)
! Copyright (c) 1996 California Institute of Technology, Pasadena, CA.
! ALL RIGHTS RESERVED.
! Based on Government Sponsored Research NAS7-03001.
!>> 2010-02-22 MESS  Krogh  Moved NSKIP=0 to start of code.
!>> 2009-10-30 MESS  Krogh  Defined DSCRN.
!>> 2009-02-28 MESS  Krogh  Added FMTT = ' ' for NAG compiler.
!>> 2009-02-28 MESS  Krogh  Fixed "f" format for C code.
!>> 2007-09-08 MESS  Krogh  Fixed definitions of MEVLAS.
!>> 2006-07-27 MESS  Krogh  Fixed boundary case in printing long text.
!>> 2006-03-20 MESS  Krogh  Added code for output of sparse vector.
!>> 2005-04-07 MESS  Krogh  Declared LFLGDB integer in MESSMH.
!>> 2004-12-15 MESS  Krogh  Added " - 1" at end of line on label 410.
!>> 2002-05-17 MESS  Krogh  Added way for user to get error count.
!>> 2001-12-28 MESS  Krogh  Added NSKIP for more flexible output values.
!>> 2000-12-30 MESS  Krogh  Fixed some types/casts in C code.
!>> 1997-12-12 MESS  Krogh  Prefixed 0P edit descriptor to F format.
!>> 1996-07-11 MESS  Krogh  Transpose matrix output for C.
!>> 1996-06-27 MESS  Krogh  fprintf(stdout, => printf( & memset now used
!>> 1996-06-18 MESS  Krogh  "Saved" NTEXTR.
!>> 1996-05-15 MESS  Krogh  Changes to use .C. and C%%.
!>> 1996-03-30 MESS  Krogh  Added external statement.
!>> 1996-01-24 MESS  Krogh  Fixed minor bug introduced with "$ " stuff.
!>> 1996-01-23 MESS  Krogh  Minor changes for C conversion.
!>> 1995-11-10 MESS  Krogh  Add code to change "$ " to " " in headings.
!>> 1995-08-11 MESS  Krogh  Made code default to not using UMESS.
!>> 1995-01-20 MESS  Krogh  Fixed unusual case in matrix output.
!>> 1994-12-15 MESS  Krogh  Removed block data for Cray T3D.
!>> 1994-11-11 MESS  Krogh  Declared all vars.
!>> 1994-09-14 MESS  Krogh  Fixed to get 1 more "$" in C output.
!>> 1994-09-08 MESS  Krogh  Added new matrix/vector capabilities.
!>> 1994-08-22 MESS  Krogh  Fix for conversion to C for new converter.
!>> 1994-07-05 MESS  Krogh  Fixed bug, KDI and FMTI could be inconsist.
!>> 1994-05-20 MESS  Krogh  Changes to MESSFT so line 1 can go to file.
!>> 1994-05-20 MESS  Krogh  Changes to setting output unit.
!>> 1994-05-09 MESS  Krogh  Integer vectors had overflow & space probs.
!>> 1993-05-19 MESS  Krogh  Changed TEXT to array of character strings.
!>> 1993-04-14 MESS  Krogh  Fixes for conversion to C. (C%% comments.)
!>> 1993-03-10 MESS  Krogh  Broke into smaller pieces.
!>> 1992-12-02 MESS  Krogh  Added save statement to block data subpr.
!>> 1992-07-13 MESS  Krogh  Add checks in heading set up.
!>> 1992-07-12 MESS  Krogh  Fixed so $$ prints a single $ in TEXT.
!>> 1992-07-12 MESS  Krogh  Set out of bound inputs to limit values.
!>> 1992-07-12 MESS  Krogh  Fixed so output works to alternate files.
!>> 1992-07-12 MESS  Krogh  Added integer declarations for parameters.
!>> 1992-06-24 MESS  Krogh  More blanks allowed on break of long lines.
!>> 1992-06-10 MESS  Krogh  Minor fix to vector output.
!>> 1992-05-27 MESS  Krogh  Fixed bug on line width setting.
!>> 1992-05-14 MESS  Krogh  Put common blocks in save statement.
!>> 1992-05-11 MESS  Krogh  Added label to assigned go to & a comment.
!>> 1992-04-08 MESS  Krogh  Unused labels 60, 220 and 320 removed.
!>> 1992-03-20 MESS  Krogh  Changed status on open to SCRATCH.
!>> 1992-03-10 MESS  Krogh  1 line below label 690 changed max to min.
!>> 1992-02-05 MESS  Krogh  Fixed bugs in printing matrix labels.
!>> 1992-01-29 MESS  Krogh  Added UMESS and multiple print option.
!>> 1991-12-09 MESS  Krogh  Fine tuning of vector output.
!>> 1991-10-10 MESS  Krogh  Insure no stop if stop level = 9.
!>> 1991-06-26 MESS  Krogh  Initial Code.
! Processes Messages -- Actions are controlled by MACT().
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
! MESUNI=10  (0 .le. K10 .le. 99) Set the unit to use for a scratch
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
! MEHEAD=11  (0 .le. K11 .le. 1) Defines the print that surrounds an
!            error message.  K11=0 gives nothing, and 1 gives the first
!            4 characters in TEXT repeated 18 times.  If this is not
!            used, one gets 72 $'s.  (To get a blank line use 1 with
!            TEXT = '    '.)
! MEDDIG=12  (-50 .le. K12 .le. 50) Set default digits to print for
!            floating point.  If K12 > 0 then K12 significant digits
!            will be printed, if K12 < 0, then -K12 digits will be
!            printed after the decimal point, and if K12 = 0, the
!            default will be used, which is the full machine precision.
!            Setting or getting this value will only work properly if
!            the action is taken by calling SMESS or DMESS as
!            appropriate.
! MEMLIN=13  (39 .le. K13 .le. 500) Set message line length to K13.
!            (Default is 128.)
! MEELIN=14  (39 .le. K14 .le. 500) Set error message line length to
!            K14. (Default is 79)
! MEMUNI=15  (-99 .le. K15 .le. 99) Messages go to unit K15.  If K15 = 0
!            (default), 'print' is used.  If K15 < 0, messages go to
!            both 'print' and to unit abs(K15).  If a write can not be
!            done to unit abs(K15), this unit will be opened with file
!            name MESS_Fxx.tmp, where xx is the value of abs(K15).
! MEEUNI=16  (-99 .le. K16 .le. 99) As for MEMUNI, except for Error
!            Messages.
! MESCRN=17  (0 .le. K17 .le. 100000000) Set number of lines to print to
!            standard output before pausing for "go" from user.  Default
!            is 0, which never stops.
! MEDIAG=18  (0 .le. K18 .le. 1000000000) Set the diagnostic level
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
! MEMAXE=19  (0 .le. K19 .le. 1000000000) Set the maximum error value.
!            When retrieving this value, it is the maximum value seen
!            for 10000*s + 1000*p + i, where s, p, and i are the stop
!            and print levels, and the index on the last error message
!            processed, respectively.  See MEEMES below.
! MESTOP=20  (0 .le. K20 .le. 8) Set the stop level for error messages.
!            If an error message has a stop index > min(K20, 8), the
!            program is stopped after processing the message.  The
!            default value is K20=3.
! MEPRNT=21  (0 .le. K21 .le. 8) Set the print level for error messages.
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
! METDIG=22  (-50 .le. K22 .le. 50) As for MEDDIG, except the value here
!            is temporary, lasting until the return, or next use of this
!            action.  If 0, the internal value for K12 is used instead.
! MENTXT=23  (1 .le. K23 .le. 10000000) Set value of NTEXT to K23.
! MEIDAT=24  (1 .le. K24 .le. 1000000000) Set value of NIDAT to K24.
! MEFDAT=25  (1 .le. K25 .le. 1000000000) Set value of NFDAT to K25.
! MEMDAT=26  (1 .le. K26 .le. 5) Set value of NMDAT to K26.
! MEMDA1=27  (K27) set MDAT(1) to K27.  See description of NMDAT above.
! MEMDA2=28  (K28) set MDAT(2) to K28.
! MEMDA3=29  (K29) set MDAT(3) to K29.
! MEMDA4=30  (K30) set MDAT(4) to K30.
! MEMDA5=31  (K31) set MDAT(5) to K31.
! METABS=32  (1 .le. K32 .le. 100) set spacing for tabs to K32.
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
!            p the print level, and should have 10 > p .ge. s .ge. 0.
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
!            0 .le. rr,dd,ww < 100, and 0 .le. t < 10.
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
! LPRINT For error messages with a print level .le. LPRINT nothing is
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
integer lnmsg, lnerr
parameter (lnmsg=128)
parameter (lnerr=79)
!
! ************** Parameters Defining Actions (See Above) ***************
!
integer   mesuni, mehead, meddig, memlin, meelin, memuni, meeuni, &
  mescrn, mediag, memaxe, mestop, meprnt, metdig, mentxt, meidat, &
  mefdat, memdat, memda1, memda2, memda3, memda4, memda5, metabs, &
  meerrs, mecont, meret , meemes, metext, metabl, meres3, meivci, &
  meivec, meimat, mejvci, mejvec, mejmat, mefvci, mefvec, mefmat, &
  megvci, megvec, megmat, memaxi, megbas, mevbas, mevlas, mefspv
! Parameters for changing the environment.
parameter (mesuni=10,mehead=11,meddig=12,memlin=13,meelin=14, &
 memuni=15,meeuni=16,mescrn=17,mediag=18,memaxe=19,mestop=20, &
 meprnt=21,metdig=22,mentxt=23,meidat=24,mefdat=25,memdat=26, &
 memda1=27,memda2=28,memda3=29,memda4=30,memda5=31,metabs=32, &
 meerrs=33)
! Parameters for actions.
parameter (mecont=50, meret=51,meemes=52,metext=53,mefspv=54, &
 metabl=55,meres3=56,meivec=57,meimat=58,mejvec=59,mejmat=60, &
 mefvec=61,mefmat=62,megvec=63,megmat=64,meivci=65,mejvci=66, &
 mefvci=67,megvci=68)
! Parameter derived from those above.
parameter (memaxi=68,megbas=49,mevbas=10,mevlas=33)
!
! ************************** Variable Declarations *********************
!
external messgs
integer    mact(*), idat(*)
character  text(*)*(*)
!
integer    i, icol, incm(mecont:meimat), inerr, iout, irow, irow1, &
    itextr, itxtsv, j, jj, k, k1, k2, kdilab, kk, knt, kolwid, kp, &
    ks, lasknt, lbuf1, lbuf2, lenbuf, linstr, m, &
    mbndhi(mevbas:mevlas), mbndlo(mevbas:mevlas), mtext(2), &
    mtextc, mtextr, nline, nroco(2), nskip, ntextr, ntxtsv
integer messgs
logical   getw, first
character ermsg*63, ermsg1*27
character sc, c
parameter (sc='$')
save  first, i, icol, inerr, irow, irow1, itxtsv, kdilab, knt, &
   lasknt, m, mtext, nline, nskip, ntextr, ntxtsv
save /cmessi/, /cmessc/
equivalence (mtext(1), mtextr), (mtext(2), mtextc)
!
! ************************** Data from common block ********************
!
parameter (lenbuf=250)
logical          xarg, gotfmt, xargok
integer          errcnt, eunit, ichar0, irc, ivar(mevbas:mevlas), &
   imag, inc, itext, iwf, iwg, kdf, kdfdef, kdi, kdiag, kdj, &
   kline, kscrn, kshift, kspec, kt, maxerr, lasti, lbuf, lenlin, &
   lenout, lentry, lentxt, lhead, linerr, linmsg, locbeg, lprint, &
   lstop, lstrt, ltext, maxwid(2), mdat(5), mpt, munit, ncol, &
   ndim, nfdat, nidat, nmdat, nrow, ntext, ounit, sunit, tabspa
!
character buf*(lenbuf), dols*72, fmtc*7, fmtf*20, fmtg*15, &
  fmti*7, fmtim(2)*7, fmtj*7, fmtr*7, fmtt*15
common /cmessi/ sunit, lhead, kdfdef, linmsg, linerr, munit, &
   eunit, kscrn, kdiag, maxerr, lstop, lprint, kdf, ntext, nidat, &
   nfdat, nmdat, mdat, tabspa, errcnt, ichar0, imag, inc, irc, &
   itext, iwf, iwg, kdi, kdj, kline, kshift, kspec, kt, lasti, &
   lbuf, lenlin, lenout, lentry, lentxt, locbeg, lstrt, ltext, &
   maxwid, mpt, nrow, ncol, ndim, ounit, gotfmt, xarg, xargok
common /cmessc / buf, dols, fmtf, fmtg, fmti, fmtj, fmtt, fmtim
equivalence (ivar(mevbas), sunit)
equivalence (fmtim(1), fmtr), (fmtim(2), fmtc)
!
equivalence (nroco, nrow)
equivalence (maxwid(1), linstr), (maxwid(2), kolwid)
! ************************** End of stuff from common block ************
!
data inerr, first / 0, .true. /
data ermsg / &
' reports error: Stop level = x, Print level = y, Error index = '/
data ermsg1 / ': Print level = y, Index = ' /
!                 50  51, 52  53  54 55 56 57 58
data incm /  1,  1,  4,  1,  2, 0, 0, 2, 6 /
data mbndlo /  0, 0, -50,  39,  39, -99, -99,         0, &
       0,          0, 0, 0, -50,        1,          1, &
       1, 1, -1000000000, -1000000000, -1000000000, -1000000000, &
       -1000000000,   1,  0 /
data mbndhi / 99, 1,  50, 500, 500,  99,  99, 100000000, &
  1000000000, 1000000000, 8, 8,  50, 10000000, 1000000000, &
  1000000000, 5, 1000000000, 1000000000, 1000000000, 1000000000, &
  1000000000, 100, 1000000000 /
!
! ************************* Start of Executable Code *******************
!
!
nskip = 0
if (first) then
   first = .false.
! Initialize common block
   sunit = -1
   lhead = 1
   linmsg = lnmsg
   linerr = lnerr
   munit = 0
   eunit = 0
   kscrn = 0
   maxerr = 0
   tabspa = 6
   lstop = 3
   lprint = 3
   errcnt = 0
   ichar0 = ichar('0')
   kdi = 1
   kdj = 6
   lenlin = lnmsg
   lentry = 1
   ounit = 0
!++ CODE for ~.C. is active
   dols(1:40) = '$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$'
   dols(41:72) ='$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$'
   fmti = '(99I01)'
   fmtj = '(99I06)'
   fmtg = '(1P,99Exx.xx)  '
!++ CODE for .C. is inactive
!%%    memset(cmessc.dols,'$',72);
!      FMTI = '%*d'
!      FMTJ = '%*d\0'
!      FMTG = '%*.*E\0'
!++ END
else
!               1  2  3   4    5    6    7    8   9   10   11
   go to (5,10,20,850,1160,1620,1130,1530,960,1210,1220), lentry
end if
!                             First entry for a message
5 lbuf = 0
!                             Usual continuation entry
10 i = 1
ntext = 1
itext = 1
lentxt = len(text(1))
nidat = 1
nfdat = 1
nmdat = 1
go to 120
!                     Continuation entry when have a non printing error
! Skip all actions -- Inside non-printing error message.
20 i = 1
30 k = mact(i)
if (k .le. meret) then
   if (k .eq. meret) go to 120
   if (k .eq. mecont) return
   if (k .le. -megbas) go to 180
   i = i + 2
else
   if (k .gt. meimat) then
      if (k .gt. memaxi) go to 180
      k = meivec + mod(k - meivec, 2)
   end if
   i = i + incm(k)
end if
go to 30
!
! Print BUF
40 call messpr
!                             Usual place to end an action request.
100 i = i + incm(m)
!                             Pick up the next action request
120 m = mact(i)
if (m .gt. megbas) go to 140
i = i + 2
if (abs(m) .gt. mevlas) go to 180
if (m .gt. 0) then
   ivar(m) = mact(i-1)
   if (ivar(m) .lt. mbndlo(m)) then
      ivar(m) = mbndlo(m)
   else if (ivar(m) .gt. mbndhi(m)) then
      ivar(m) = mbndhi(m)
   end if
!            MEHEAD, MEDDIG, MEMLIN, MEELIN, MEMUNI, MEEUNI
   go to (122,    124,    126,    126,    128,    128), m - mesuni
   if (m .ne. mentxt) go to 120
   itext = (ntext-1) / lentxt
   ntext = ntext - lentxt*itext
   itext = itext + 1
   go to 120
122    if (lhead .ne. 0) then
   end if
   go to 120
124    kdf = kdfdef
   go to 120
126    lenlin = linmsg
   go to 120
128    if (ivar(m) .ne. 0) then
!%%          k = labs(cmessi.ounit);
!%%          c_fname[m-15][6] = k / 10 + '0';
!%%          c_fname[m-15][7] = k % 10 + '0';
!%%          if (strcmp(&c_fname[16-m][6], &c_fname[m-15][6]))
!%%             c_handle[m-15] = fopen(c_fname[m-15],"w");
!%%          else
!%%             c_handle[m-15] = c_handle[16-m];
      k = abs(ivar(m))
   end if
   ounit = munit
   go to 120
end if
if (m .eq. -mesuni) then
!%%      if (cmessi.sunit == -1L) {
!%%          scratch_file = tmpfile();
!%%          cmessi.sunit = 1L;}
   if (sunit .le. 0) sunit = messgs()
end if
!
mact(i-1) = ivar(-m)
go to 120
!  ME ..    CONT  RET EMES ETXT  FSPV TABL
140 go to (170, 200, 310, 400, 1200, 910, 180), m-megbas
if (m .le. megvci) go to 1000
go to 180
!
! Action MECONT -- Continue message on next entry
170 lentry = 2
return
!
! Some kind of error in message specification.
180 continue
!++ CODE for ~.C. is active
buf(1:57) = &
   'Actions in MESS terminated due to error in usage of MESS.'
!++ CODE for .C. is inactive
!%%   memcpy(cmessc.buf,
!%%   "Actions in MESS terminated due to error in usage of MESS.",57);
!++ END
lbuf = 57
!
! Action MERET -- Finish a message.
200 lentry = 1
j = inerr
inerr = 0
if (j .ge. 2) inerr = -2
if (j .gt. 0) go to 330
!                       Finish print before exit.
call messpr
return
!
! Action MEEMES -- Start an error message
310 lentry = 3
errcnt = errcnt + 1
!++  Code for UMESS is inactive
!      call UMESS(TEXT, MACT(I+1), IVAR)
!++  End
imag = max( 0, min(999, mact(i+2)))
k = mact(i+1)
maxerr = max(maxerr, 1000*k + imag)
ks = k / 10
kp = k - 10 * ks
if (ks .le. min(lstop, 8)) then
   if (kp .le. lprint) then
      inerr = -1
      go to 20
   end if
   inerr = 1
else
   inerr = 2
end if
ounit = eunit
lenlin = linerr
!                        Output a blank line.
buf(1:1) = ' '
lbuf = 1
330 call messpr
!                        Put out line of $'s
if (lhead .ne. 0) then
   lbuf = min(len(dols), lenlin)
!++ CODE for ~.C. is active
   buf(1:lbuf) = dols(1:lbuf)
   if (inerr.lt.0) buf(5:37)=' Fatal error -- Program stopped. '
!++ CODE for .C. is inactive
!%%      memcpy(cmessc.buf, cmessc.dols, cmessi.lbuf);
!%%      if (inerr < 0L)
!%%      memcpy(&cmessc.buf[4]," Fatal error -- Program stopped. ",34);
!++ END
   call messpr
end if
if (inerr .le. 0) then
!                                 Just finished an error message
   if (inerr .ne. 0) stop
   ounit = munit
   lenlin = linmsg
   return
end if
!                     Just starting an error message get program name
ntextr = 0
go to 410
!                     Got the program name in BUF.
370 lbuf = min(lbuf, 40)
if (ks .eq. 0) then
   ermsg1(17:17) = char(kp + ichar0)
!%%       memcpy(&cmessc.buf[cmessi.lbuf], ermsg1, strlen(ermsg1));
   buf(lbuf+1:lbuf+len(ermsg1)) = ermsg1
   lbuf = lbuf + len(ermsg1)
else
   ermsg(30:30) = char(ks + ichar0)
   ermsg(47:47) = char(kp + ichar0)
!%%       memcpy(&cmessc.buf[cmessi.lbuf], ermsg, strlen(ermsg));
   buf(lbuf+1:lbuf+len(ermsg)) = ermsg
   lbuf = lbuf + len(ermsg)
end if
lstrt = lbuf + 1
call messfi
lbuf = lbuf + kdi
!%%   sprintf(&cmessc.buf[cmessi.lstrt-1L], "%*ld",
!%%           (int)messcc.kciwid, cmessi.imag);
write (buf(lstrt:lbuf), fmti) imag
!          Finish up the start error message action.
if (mact(i+3) .lt. 0) go to 40
if (mact(i+3) .ne. 0) then
   itext = (mact(i+3)-1) / lentxt
   ntext = mact(i+3) - lentxt*itext
   itext = itext + 1
end if
kspec = 13
go to 480
!                  Take care of any left over print from error header
390 if (lbuf .ne. 0) call messpr
!
! Action METEXT -- Print string from TEXT
400 lentry = 4
ntextr = ntext
itextr = itext
!                  Continue with print from TEXT
! K     take at most K-1 chars., but if 0 take max number
! K1    is last loc. used from TEXT if LENTXT is BIG.
! NEXT  is first character location in TEXT(ITEXT)
! K2    is last character location in TEXT(ITEXT)
! LSTRT is first character position in BUF
! LBUF  is last used character position in BUF

410 lstrt = lbuf + 1
k2 = min(lentxt, ntext + (lenbuf - lstrt))
!%%       if ((ctmp=memchr(TEXT(cmessi.itext-1L,cmessi.ntext-1), SC,
!%%          k2 - cmessi.ntext + 1)) == NULL)
!%%             k = 0;
!%%       else
!%%             k = ctmp - TEXT(cmessi.itext-1L,cmessi.ntext-1) + 1;
k = index(text(itext)(ntext:k2), sc)
if (k .eq. 0) then
! Want to take all that we can.
   lbuf = lstrt + k2 - ntext
!%%       memcpy(&cmessc.buf[cmessi.lstrt-1L], TEXT(cmessi.itext-1L,
!%%         cmessi.ntext-1), k2 - cmessi.ntext + 1L);
   buf(lstrt:lbuf) = text(itext)(ntext:k2)
   if (k2 .eq. lentxt) then
     itext = itext + 1
     ntext = 1
     if (lbuf .le. lenlin) go to 410
   else
     ntext = k2 + 1
   end if
   kspec = 12
   if (itext - itextr .lt. 4000) go to 480
   kspec = 2
   go to 430
end if
lbuf = lbuf + k - 1
!%%   if (k >= 2) memcpy(&cmessc.buf[cmessi.lstrt-1],
!%%     TEXT(cmessi.itext-1L, cmessi.ntext-1), k - 1L);
if (k .ge. 2) buf(lstrt:lbuf) = text(itext)(ntext:ntext+k-2)
!        Jump to location below if get $ after computing an NSKIP.
415 continue
ntext = ntext + k + 1
if (ntext .gt. lentxt) then
   itext = itext + 1
   if (ntext .eq. lentxt + 1) then
      c = text(itext-1)(lentxt:lentxt)
      ntext = 1
   else
      c = text(itext)(1:1)
      ntext = 2
   end if
else
   c = text(itext)(ntext-1:ntext-1)
end if
if (c .eq. ' ') then
!                Special code to take care of " " following "$".
   ntext = ntext - 1
   if (ntext .eq. 0) then
      itext = itext - 1
      ntext = lentxt
   end if
   go to 410
end if
if (ntextr .eq. 0) then
   if (lentry .eq. 3) go to 370
   go to 1510
end if
kspec = index('BERNIMFJG(T', c)
430 if (lbuf .gt. lenlin) go to 480
!              1   2   3   4   5   6   7   8   9  10  11  12, 13
!              B   E   R   N   I   M   F   J   G   (   T done end err
go to (455,480,450,460,700,680,900,700,900,600,690,410,390), kspec
!               No match  -- Check for setting NSKIP
if (((c .ge. '0') .and. (c .le. '9')) .or. (c .eq. '-')) then
   nskip = 0
   k1 = 1
   if (c .ne. '-') go to 436
   k1 = -1
433    c = text(itext)(ntext:ntext)
   ntext = ntext + 1
   if (ntext .ge. lentxt) then
      itext = itext + 1
      ntext = 1
   end if
436    if ((c .ge. '0') .and. (c .le. '9')) then
      nskip = 10 * nskip + k1 * (ichar(c) - ichar0)
      go to 433
   end if
   if (c .eq. '$') then
      k = 0
      go to 415
   end if
end if
!
! Continue with the text.
440 lbuf = lbuf + 1
buf(lbuf:lbuf) = c
go to 410
!                        Reset NTEXT for $R
450 ntext = ntextr
itext = itextr
!                             Done with METEXT action.
455 nmdat = 1
go to 100
!           At this point want to output all in BUF
460 do 470 lbuf = lbuf, 1, -1
   if (buf(lbuf:lbuf) .ne. ' ') go to 480
470 continue
480 lbuf2 = lbuf
if (lbuf2 .eq. 0) then
   lbuf = 1
   buf(1:1) = ' '
else if (lbuf .gt. lenlin) then
   do 485 k = lenlin+1, lenlin/3, -1
      if (buf(k:k) .eq. ' ') then
         lbuf = k - 1
         go to 490
      end if
485    continue
   lbuf = lenlin
end if
490 lbuf1 = lbuf
call messpr
if (lbuf1 .ge. lbuf2) then
!                       The entire buffer has been printed.
   if (kspec .le. 2) go to 455
   if (kspec .ne. 4) go to 430
   go to 410
end if
!                       Remove trailing blanks
do 510 lbuf1 = lbuf1+1, lbuf2
   if (buf(lbuf1:lbuf1) .ne. ' ') go to 520
510 continue
!                       Shift the contents of the buffer.
520 lbuf = lbuf2-lbuf1+1
lstrt = 1
530 if (lbuf .ge. lbuf1) then
!                              Take care of overlap.
   k = 2*lbuf1 - lstrt
!%%memcpy(&cmessc.buf[cmessi.lstrt-1],&cmessc.buf[lbuf1-1],k-lbuf1);
   buf(lstrt:lbuf1-1) = buf(lbuf1:k-1)
   lstrt = lbuf1
   lbuf1 = k
   go to 530
end if
!%% if (cmessi.lbuf>=cmessi.lstrt) memcpy(&cmessc.buf[cmessi.lstrt-1],
!%%       &cmessc.buf[lbuf1-1L], lbuf2-lbuf1+1);
if (lbuf .ge. lstrt) buf(lstrt:lbuf) = buf(lbuf1:lbuf2)
go to 430
!
! Get information on user format
600 kspec = 8
!              I,   i,   F,   f,   E,   e,   G,   g
go to (604, 604, 601, 601, 602, 602, 602, 602), &
   index('IiFfEeGg',text(itext)(ntext:ntext))
go to 180
601 continue
!++ CODE for ~.C. is active
fmtg='(0P,99F'
!++ CODE for .C. is inactive
!%%   strcpy(cmessc.fmtg, "%*.*f\0");
!%%   messcc.lgprec = 0;
!++ END
go to 603
602 continue
!++ CODE for ~.C. is active
fmtg='(1P,99'//text(itext)(ntext:ntext)
!++ CODE for .C. is inactive
!%%   strcpy(cmessc.fmtg, "%*.*E\0");
!      FMTG(5:5) = TEXT(ITEXT)(NTEXT:NTEXT)
!%%   messcc.lgprec = 0;
!++ END
603 kspec = 9
604 imag = 0
getw = .true.
k = ntext
fmtt = ' '
606 continue
   ntext = ntext + 1
   if (ntext .gt. lentxt) then
      itext = itext + 1
      ntext = 1
   end if
!++ CODE for ~.C. is active
   fmtt(ntext-k:ntext-k) = text(itext)(ntext:ntext)
!++ END
   jj = ichar(text(itext)(ntext:ntext)) - ichar0
   if (getw) then
      if ((jj .ge. 0) .and. (jj .le. 9)) then
         imag = 10*imag + jj
      else
         if (text(itext)(ntext:ntext) .eq. ')')  go to 610
         if (text(itext)(ntext:ntext) .ne. '.')  go to 180
         getw = .false.
      end if
   else
      if (text(itext)(ntext:ntext) .eq. ')') go to 610
      if ((jj .lt. 0) .or. (jj .gt. 9)) go to 180
!++ CODE for .C. is inactive
!%%         messcc.lgprec = 10*messcc.lgprec + jj;
!++ END
   end if
go to 606
!
610 ntext = ntext + 1
if (ntext .gt. lentxt) then
   itext = itext + 1
   ntext = 1
end if
!++ CODE for ~.C. is active
if (kspec .eq. 8) then
   kdj = imag
   fmtj(5:7) = fmtt
else
   iwg = imag
   fmtg(8:15) = fmtt
end if
!++ CODE for .C. is inactive
!%%   if (cmessi.kspec == 8)
!%%       cmessi.kdj = cmessi.imag;
!%%   else
!%%       cmessi.iwg = cmessi.imag;
!++ END
if (text(itext)(ntext:ntext) .eq. sc) go to 410
if (kspec .eq. 8) go to 700
if (xargok) return
go to 440
!
!                         Print from MDAT
680 iout = mdat(nmdat)
if (nmdat .ge. 6) then
   mdat(nmdat) = mdat(nmdat) + 1
else
   nmdat = nmdat + 1
end if
go to 720
!
!                         Process a tab
690 lstrt = lbuf + 1
lbuf = min(lbuf + tabspa - mod(lbuf, tabspa), lenlin+1)
!%%  for (kc=cmessi.lstrt-1; kc<cmessi.lbuf; kc++) cmessc.buf[kc]=' ';
buf(lstrt:lbuf) = ' '
go to 850
!                         Print from IDAT
700 nidat = nidat + nskip
nskip = 0
iout = idat(nidat)
nidat = nidat + 1
720 lstrt = lbuf + 1
imag = iout
if (kspec .ge. 8) then
   lbuf = lbuf + kdj
!%%   sprintf(&cmessc.buf[cmessi.lstrt-1],"%*ld",(int)cmessi.kdj, iout);
write (buf(lstrt:lbuf), fmtj) iout
   go to 850
end if
!
!                Get format for integer output.
call messfi
lbuf = lbuf + kdi
!%% sprintf(&cmessc.buf[cmessi.lstrt-1],"%*ld",(int)messcc.kciwid,iout);
write (buf(lstrt:lbuf), fmti) iout
!                         Entry here to check line after numeric output.
850 if (lbuf .le. lenlin) go to 410
kspec = 12
go to 480
!
!                          Take care of output for extra argument.
900 if (xargok) return
go to 180
!
! Action METABL -- Start a table
910 gotfmt = mact(i+1) .ne. 1
if (.not. gotfmt) then
   irow = 0
   kolwid = 0
end if
lentry = 9
if (lbuf .ne. 0) call messpr
920 continue
!%%   memset(cmessc.buf,' ',LENBUF);
buf = ' '
nrow = 1
ncol = mact(i+3)
icol = i + 3
940 icol = icol + 1
jj = mact(icol)
kline = mod(jj, 100)
lenout = jj / 100000
ncol = ncol - max(kline, 1)
if (gotfmt) then
!                                 Print the data
   lstrt = lbuf + 1
   lbuf = min(lbuf + kline * lenout, lenbuf)
   jj = jj / 100
   kk = mod(jj, 10)
!              Text,   I   I',   F    E    G
   go to (948, 941, 941, 943, 945, 944), kk
   go to 180
!                             Integer output
941    continue
!++ CODE for ~.C. is active
   kdi = lenout
   fmti(5:5) = char(lenout / 10 + ichar0)
   fmti(6:6) = char(mod(lenout, 10) + ichar0)
!++ END
   if (kk .eq. 3) then
!%%         sprintf(&cmessc.buf[cmessi.lstrt-1], "%*ld",
!%%            (int)cmessi.lenout, mact[i]);
      write (buf(lstrt:lbuf), fmti) mact(i+1)
      go to 960
   end if
!                            Regular integer output
   nidat = nidat + nskip
   nskip = 0
!++ CODE for ~.C. is active
   write (buf(lstrt:lbuf), fmti) (idat(k), k = nidat, &
      nidat+kline-1)
   nidat = nidat + kline
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
   fmtf = '(0P,99F  .  )'
!++ END
   go to 946
944    continue
!++ CODE for ~.C. is active
   fmtf = '(1P,99G  .  )'
!++ END
   go to 946
945    continue
!++ CODE for ~.C. is active
   fmtf = '(1P,99E  .  )'
!++ END
946    jj = mod(jj/10, 100)
!++ CODE for ~.C. is active
   fmtf(8:8) = char(ichar0 + lenout / 10)
   fmtf(9:9) = char(ichar0 + mod(lenout, 10))
   fmtf(11:11) = char(ichar0 + jj / 10)
   fmtf(12:12) = char(ichar0 + mod(jj, 10))
!++ CODE for .C. is inactive
!%%      strcpy(cmessc.fmtf, "%*.*E\0");
!        IWF = LENOUT
!        lfprec = JJ
!++ END
   if (.not. xargok) go to 180
   mpt = nfdat
   nfdat = nfdat + kline
   return
!                           Text output
948    k1 = ntext + lbuf - lstrt
!%%    memcpy(&cmessc.buf[cmessi.lstrt-1], TEXT(cmessi.itext-1,
!%%       cmessi.ntext -1), k1 - cmessi.ntext);
   buf(lstrt:lbuf) = text(itext)(ntext:k1-1)
   ntext = k1
else
!                                 Print the heading
   kt = 1
   call messmh(text)
   if (kt .lt. 0) go to 180
end if
960 if ((lbuf .le. mdat(nrow)) .and. (ncol .gt. 0)) go to 940
if (nrow .eq. 1) then
   jj = lbuf
   lbuf = mdat(1)
   call messpr
   lbuf = jj
else
   if (irow .eq. 0) then
      if (nrow .eq. 2) then
!++ CODE for ~.C. is active
         if (sunit .le. 0) sunit = messgs()
         rewind(sunit)
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
   write(sunit) buf(5:mdat(nrow))
end if
if (lbuf .gt. mdat(nrow)) then
!%%  memcpy(&cmessc.buf[4], &cmessc.buf[cmessi.mdat[cmessi.nrow-1]],
!%%     cmessi.lbuf - cmessi.mdat[cmessi.nrow-1]);
   buf(5:lbuf - mdat(nrow) + 4) = buf(mdat(nrow)+1:lbuf)
   lbuf = lbuf - mdat(nrow) + 4
   nrow = nrow + 1
   if (.not. gotfmt) then
      if (nrow .gt. 5) go to 180
      mdat(nrow) = lbuf
   end if
   if (ncol .eq. 0) go to 960
   go to 940
end if
lbuf = 0
if (.not. gotfmt) then
   gotfmt = .true.
   irow = irow - 1
   go to 920
end if
mact(i+1) = mact(i+1) + 1
if (mact(i+1) .le. mact(i+2)) go to 999
mact(i+1) = 1
if (nrow .eq. 1) go to 999
!%%    fputc(EOF, scratch_file);
endfile sunit
kk = 1
994 kk = kk + 1
if (kk .gt. nrow) go to 999
!%%   rewind(scratch_file);
rewind(sunit)
irow = -1
k = kk
995 lbuf = 5
irow = irow + 1
if (irow .ne. 0) then
!%%      sprintf(cmessc.buf, "%4ld",  irow%10000);
   write(buf(1:4), '(I4)') mod(irow, 10000)
else
!%%    memset(cmessc.buf,' ',4);
   buf(1:4) = ' '
end if
do 996 j = 2, k
   if (j .eq. k) lbuf = mdat(kk)
!%%       if (fread(&cmessc.buf[4], cmessi.lbuf-4, 1,
!%%         scratch_file) == 0) goto L_994;
   read(sunit, end = 994) buf(5:lbuf)
996 continue
k = nrow
call messpr
go to 995
999 lentry = 1
return
!
!                          Get started with vector or matrix output
1000 inc = 1
locbeg = nidat
if (m .gt. megmat) then
! Have a user set increment between entries of a vector.
  m = meivec + 2 * (m - meivci)
  i = i + 1
  inc = mact(i)
end if
xarg = m .gt. mejmat
if (xarg) then
   m = m - 4
   locbeg = nfdat
   if (.not. xargok) go to 40
end if
gotfmt = m .gt. meimat
if (gotfmt) m = m - 2
locbeg = locbeg + nskip
nskip = 0
mpt = locbeg
if (m .eq. meimat) go to 1300
!                           Take care of setup for vector output
knt = 0
lasknt = mact(i+1)
if (lasknt .le. 0) then
   lasknt = -lasknt
   knt = locbeg - 1
   if (lasknt .le. knt) go to 40
end if
imag = lasknt
lasti = locbeg + inc * (lasknt - 1 - knt)
ncol = 0
!                          Get format for label output.
call messfi
!%%   messcc.kcrwid = messcc.kciwid;
fmtr = fmti
kdilab = kdi+1
linstr = 2*kdilab+2
if (xarg) then
   if (.not. gotfmt) go to 1150
   iwf = iwg
   fmtf = fmtg
!++ CODE for .C. is inactive
!%%      cmessi.iwf = cmessi.iwg;
!%%      messcc.lfprec = messcc.lgprec;
!++ END
   go to 1160
end if
call messfd(idat)
!                                          After integer format
lenout = kdi
nidat = lasti + 1
!                                          Common code continues here
1080 nline = (lenlin - linstr + 1) / lenout
if (lbuf .eq. 0) go to 1090
k = max(linstr, lbuf+1)
if (((k-linstr)/lenout + (lenlin-k+1)/lenout) .lt. nline) k = k + &
   lenout - mod(k-linstr, lenout)
kline = (lenlin - k + 1) / lenout
if (kline .lt. min(lasknt-knt, nline/2)) go to 1085
linstr = k - lenout * ((k - linstr) / lenout)
if (kline .ge. lasknt-knt)  then
   kline = lasknt - knt
   k = lbuf + 1
end if
knt = knt + kline
!%%    for (kc=cmessi.lbuf; kc < k; kc++) cmessc.buf[kc] = ' ';
buf(lbuf+1:k) = ' '
lbuf = k
go to 1110
1085 call messpr
1090 continue
!++ CODE for ~.C. is active
buf = ' '
write (buf(1:kdilab), fmtr) knt+1
!++ CODE for .C. is inactive
!%%   memset(cmessc.buf,' ',LENBUF);
!%%   sprintf(cmessc.buf, "%*ld", (int)messcc.kcrwid, knt+1);
!++ END
buf(kdilab:kdilab) = '-'
kline = min(nline, lasknt - knt)
knt = knt + kline
!%%    sprintf(&cmessc.buf[kdilab], "%*ld", (int)messcc.kcrwid, knt);
write (buf(kdilab+1:2*kdilab), fmtr) knt
!%%    cmessc.buf[kdilab*2L-1] = ':';
!%%    for (kc=kdilab*2L; kc < *linstr-1; kc++) cmessc.buf[kc] = ' ';
buf(2*kdilab:linstr-1) = ':'
lbuf = linstr
1110 lstrt = lbuf
lbuf = lbuf + lenout * kline - 1
if (xarg) return
!                                    Integer output
!++ CODE for ~.C. is active
write (buf(lstrt:lbuf), fmti) (idat(k), k = mpt, &
    mpt+inc*(kline-1), inc)
!++ CODE for .C. is inactive
!%%   for (k=cmessi.mpt; k<=cmessi.mpt+cmessi.kline-1; k++)
!%%  sprintf(&cmessc.buf[cmessi.lstrt+messcc.kciwid*(k-cmessi.mpt)-1],
!%%      "%*ld", (int)messcc.kciwid, idat[cmessi.inc*k-1]);
!++ END
mpt = mpt + kline * inc
!
!                                     Entry here after vector output.
1130 if (mpt .le. lasti) go to 1085
go to 40
!                                          Get other format
1150 lentry = 5
return
!                                          After other format
1160 lenout = iwf
lentry = 7
nfdat = lasti + 1
go to 1080


!                         Sparse vector output.
1200 xarg = .true.
if (.not. xargok) go to 40
gotfmt = .false.
mpt = 1
locbeg = 1
inc = 1
lasknt = mact(i+1)
lasti = lasknt
lentry = 10
return

!                Entry after getting format for sparse data output.
1210 lenout = iwf
lentry = 11
nline = lenlin / iwf

1220 call messpr
kline = min(lasknt - mpt + 1, nline)
if (kline .le. 0) go to 40
lbuf = lenout * kline
return

!
!                           Take care of setup for matrix output
1300 continue
ndim = mact(i+1)
if (ndim .le. 0) then
   if (ndim .eq. 0) go to 40
   inc = -ndim
   ndim = 1
end if
icol = 1
irow1 = 1
nrow = mact(i+2)
if (nrow .le. 0) then
   if (nrow .eq. 0) go to 40
   irow1 = -nrow / 100000
   nrow = -nrow - 99999 * irow1 - 1
end if
ncol = mact(i+3)
if (ncol .le. 0) then
   if (ncol .eq. 0) go to 40
   icol = -ncol / 100000
   ncol = -ncol - 99999 * irow1 - 1
end if
ntxtsv = ntext
itxtsv = itext
irc = 1
!                        Compute widths for row and column labels
1320 maxwid(irc) = 0
mtext(irc) = mact(i+irc+3)
imag = nroco(irc)
kline = imag
1330 ntext = mtext(irc)
if (ntext .ge. 0) then
   if (ntext .eq. 0) then
      ltext = 5
   else
!                        Go get row/column widths
      kt = 2
      call messmh(text)
      if (kt .lt. 0) then
         mtext(irc) = 0
         go to 1330
      end if
   end if
   call messfi
   maxwid(irc) = max(maxwid(irc), ltext + kdi+1)
!%%      if (cmessi.irc == 1)
!%%         messcc.kcrwid = cmessi.kdi;
!%%      else
!%%         messcc.kccwid = cmessi.kdi;
   fmtim(irc) = fmti
end if
irc = irc + 1
if (irc .eq. 2) go to 1320
!                 Widths for Row and column titles have been computed.
kshift = 1
lasti = locbeg + inc * (nrow - irow1)
if (xarg) then
   if (.not. gotfmt) go to 1610
!++ CODE for ~.C. is active
   iwf = iwg
   fmtf = fmtg
!++ CODE for .C. is inactive
!%%      cmessi.iwf = cmessi.iwg;
!%%      messcc.lfprec = messcc.lgprec;
!++ END
   go to 1620
end if
call messfd(idat)
!
If (kdi .ge. kolwid) then
   lenout = kdi
else
   kshift = (kolwid - kdi + 2) /2
   lenout = kolwid
!++ CODE for ~.C. is active
   kdi = kolwid
   fmti(5:5) = char(ichar0 + kolwid / 10)
   fmti(6:6) = char(ichar0 + mod(kolwid, 10))
!++ CODE for .C. is inactive
!%%  messcc.kciwid = *kolwid;
!++ END
end if
nidat = nidat + ndim*ncol
!                              Continue with commmon code
1390 nline = (lenlin - linstr) / lenout
if (lbuf .le. linstr) go to 1420
1400 call messpr
1420 irow = irow1
kline = min(nline, ncol-icol+1)
!                       Output column labels (if any)
if (mtextc .lt. 0) go to 1480
ntext = mtextc
imag = icol
kt = 3
call messmh(text)
if (kt .lt. 0) go to 180
!                       Return from output of column labels.
mtextc = ntext
1480 icol = icol + kline
1490 call messpr
!
!                      Output row labels (if any)
if (mtextr .lt. 0) go to 1520
if (mtextr .eq. 0) then
!%%       memcpy(&cmessc.buf[cmessi.lbuf],"Row ", 4);
   buf(lbuf+1:lbuf+4) = 'Row '
   lbuf = lbuf + 4
   go to 1515
end if
ntext = mtextr
itext = (ntext-1) / lentxt
ntext = ntext - itext * lentxt
itext = itext + 1

!                     Go get text for row label
ntextr = 0
go to 410
!                     Return from getting text for row label
1510 if (c .ne. '#') then
   mtextr = ntext + lentxt * (itext-1)
!%%    for (kc=cmessi.lbuf; kc < *linstr; kc++) cmessc.buf[kc] = ' ';
   buf(lbuf+1:linstr) = ' '
   go to 1520
end if
1515 continue
!%%   sprintf(&cmessc.buf[cmessi.lbuf],"%*ld",(int)messcc.kcrwid,irow);
!%%    for (kc=cmessi.lbuf+messcc.kcrwid;
!%%       kc < *linstr; kc++) cmessc.buf[kc] = ' ';
write (buf(lbuf+1:linstr), fmtr) irow
1520 lstrt = linstr + 1
lbuf = linstr + lenout*kline
lasti = mpt + ndim * kline - 1
if (xarg) return
!                                    Integer output
!%% for (k=cmessi.mpt; k<=cmessi.lasti; k+=cmessi.ndim)
!%%  sprintf(&cmessc.buf[cmessi.lstrt + messcc.kciwid*(k-cmessi.mpt)/
!%%     cmessi.ndim - 1], "%*ld", (int)messcc.kciwid, idat[k-1]);
   write (buf(lstrt:lbuf), fmti) (idat(k), k=mpt,lasti,ndim)
!
!                                     Entry here after matrix output.
1530 mpt = mpt + inc
irow = irow + 1
!
if (irow .le. nrow) go to 1490
if (icol .gt. ncol) then
   ntext = ntxtsv
   itext = itxtsv
   go to 40
end if
mpt = ndim*(icol-1) + 1
mtextr = mact(i+4)
call messpr
lbuf = 1
buf(1:1) = ' '
go to 1400
!                                Need to get format for matrix print.
1610 lentry = 6
return
!                                Entry after got format for matrix print
1620 If (iwf .ge. kolwid) then
   lenout = iwf
else
   kshift = (kolwid - iwf + 2) /2
   lenout = kolwid
!%%      cmessi.iwf = *kolwid;
!%%      strcpy(cmessc.fmtf, "%*.*E\0");
   write (fmtf(7:8), '(I2)') kolwid
end if
nfdat = nfdat + ndim*ncol
lentry = 8
go to 1390
end

subroutine messfd(idat)
! Get the format for data to be printed in vectors and arrays.
!
! ************** Variable only used here *******************************
!
! K      Temporary index.
! J      Temporary index.
! IDAT   Input array to MESS
! IMAX   Used when computing largest integer in array.
! IMIN   Used when computing smallest integer in array.
!
integer j, k, idat(*), imax, imin
!
! For comments on other variables, see the listing for MESS.
integer   lenbuf, mevbas, mevlas
parameter (lenbuf=250)
parameter (mevbas=10)
parameter (mevlas=33)
logical          xarg, gotfmt, xargok
integer          errcnt, eunit, ichar0, irc, ivar(mevbas:mevlas), &
   imag, inc, itext, iwf, iwg, kdf, kdfdef, kdi, kdiag, kdj, &
   kline, kscrn, kshift, kspec, kt, maxerr, lasti, lbuf, lenlin, &
   lenout, lentry, lentxt, lhead, linerr, linmsg, locbeg, lprint, &
   lstop, lstrt, ltext, maxwid(2), mdat(5), mpt, munit, ncol, &
   ndim, nfdat, nidat, nmdat, nrow, ntext, ounit, sunit, tabspa
!
character buf*(lenbuf), dols*72, fmtc*7, fmtf*20, fmtg*15, &
  fmti*7, fmtim(2)*7, fmtj*7, fmtr*7, fmtt*15
common /cmessi/ sunit, lhead, kdfdef, linmsg, linerr, munit, &
   eunit, kscrn, kdiag, maxerr, lstop, lprint, kdf, ntext, nidat, &
   nfdat, nmdat, mdat, tabspa, errcnt, ichar0, imag, inc, irc, &
   itext, iwf, iwg, kdi, kdj, kline, kshift, kspec, kt, lasti, &
   lbuf, lenlin, lenout, lentry, lentxt, locbeg, lstrt, ltext, &
   maxwid, mpt, nrow, ncol, ndim, ounit, gotfmt, xarg, xargok
common /cmessc / buf, dols, fmtf, fmtg, fmti, fmtj, fmtt, fmtim
equivalence (ivar(mevbas), sunit)
equivalence (fmtim(1), fmtr), (fmtim(2), fmtc)
!
save /cmessi/, /cmessc/
!
if (gotfmt) then
   kdi = kdj
!%%      messcc.kciwid = cmessi.kdj;
   fmti = fmtj
   return
end if
k = 1
imax = 1
imin = 0
10 do 20 j = locbeg, lasti, inc
   imax = max(imax, idat(j))
   imin = min(imin, idat(j))
20 continue
if (ncol .ne. 0) then
   k = k + 1
   locbeg = locbeg + ndim
   lasti = lasti + ndim
   if (k .le. ncol) go to 10
end if
imag = imax
if ((imag/10) + imin .lt. 0) imag = imin
kdi = -kdi
call messfi
return
end


subroutine messfi
! Get the format for the integer IMAG.
!
! ************** Variable only used here *******************************
!
! I, K, KD are used in determining number of characters needed to
!          represent IMAG.
!
integer i, k, kd
!
! For comments on other variables, see the listing for MESS.
integer   lenbuf, mevbas, mevlas
parameter (lenbuf=250)
parameter (mevbas=10)
parameter (mevlas=33)
logical          xarg, gotfmt, xargok
integer          errcnt, eunit, ichar0, irc, ivar(mevbas:mevlas), &
   imag, inc, itext, iwf, iwg, kdf, kdfdef, kdi, kdiag, kdj, &
   kline, kscrn, kshift, kspec, kt, maxerr, lasti, lbuf, lenlin, &
   lenout, lentry, lentxt, lhead, linerr, linmsg, locbeg, lprint, &
   lstop, lstrt, ltext, maxwid(2), mdat(5), mpt, munit, ncol, &
   ndim, nfdat, nidat, nmdat, nrow, ntext, ounit, sunit, tabspa
!
character buf*(lenbuf), dols*72, fmtc*7, fmtf*20, fmtg*15, &
  fmti*7, fmtim(2)*7, fmtj*7, fmtr*7, fmtt*15
common /cmessi/ sunit, lhead, kdfdef, linmsg, linerr, munit, &
   eunit, kscrn, kdiag, maxerr, lstop, lprint, kdf, ntext, nidat, &
   nfdat, nmdat, mdat, tabspa, errcnt, ichar0, imag, inc, irc, &
   itext, iwf, iwg, kdi, kdj, kline, kshift, kspec, kt, lasti, &
   lbuf, lenlin, lenout, lentry, lentxt, locbeg, lstrt, ltext, &
   maxwid, mpt, nrow, ncol, ndim, ounit, gotfmt, xarg, xargok
common /cmessc / buf, dols, fmtf, fmtg, fmti, fmtj, fmtt, fmtim
equivalence (ivar(mevbas), sunit)
equivalence (fmtim(1), fmtr), (fmtim(2), fmtc)
!
save /cmessi/, /cmessc/
!
kd = 1
if (kdi .lt. 0) then
!              KDI < 0 to flag need for extra space -- avoids overflows
   kdi = -kdi
   kd = 2
end if
k = 1
if (imag .lt. 0) then
   imag = -imag
   kd = kd + 1
end if
i = imag / 10
if (i .ne. 0) then
10    k = 10 * k
   kd = kd + 1
   if (i .ge. k) go to 10
end if
if (kd .ne. kdi) then
   kdi = kd
!++ CODE for ~.C. is active
   fmti(5:5) = char(ichar0 + kdi / 10)
   fmti(6:6) = char(ichar0 + mod(kdi, 10))
!++ CODE for .C. is inactive
!%%      messcc.kciwid = cmessi.kdi;
!++ END
end if
return
end

!++ CODE for ~.C. is active
integer function messgs()
!                                 Get a scratch unit assigned.
integer j
!
messgs = 31
10 messgs = messgs - 1
if (messgs .eq. 0) stop 'Could not assign scratch unit in MESS.'
open (messgs, status='SCRATCH', access='SEQUENTIAL', &
    form='UNFORMATTED', iostat=j)
if (j .ne. 0) go to 10
return
end
!++ END

subroutine messmh(text)
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
!
integer j, k, kb, kk, l, lflgdb, lstrdb, ltxtdb
character*(*)  text(*)
character sc, c
parameter (sc='$')
! For comments on other variables, see the listing for MESS.
integer   kolwid, linstr, lenbuf, mevbas, mevlas
parameter (lenbuf=250)
parameter (mevbas=10)
parameter (mevlas=33)
logical          xarg, gotfmt, xargok
integer          errcnt, eunit, ichar0, irc, ivar(mevbas:mevlas), &
   imag, inc, itext, iwf, iwg, kdf, kdfdef, kdi, kdiag, kdj, &
   kline, kscrn, kshift, kspec, kt, maxerr, lasti, lbuf, lenlin, &
   lenout, lentry, lentxt, lhead, linerr, linmsg, locbeg, lprint, &
   lstop, lstrt, ltext, maxwid(2), mdat(5), mpt, munit, ncol, &
   ndim, nfdat, nidat, nmdat, nrow, ntext, ounit, sunit, tabspa
!
character buf*(lenbuf), dols*72, fmtc*7, fmtf*20, fmtg*15, &
  fmti*7, fmtim(2)*7, fmtj*7, fmtr*7, fmtt*15
common /cmessi/ sunit, lhead, kdfdef, linmsg, linerr, munit, &
   eunit, kscrn, kdiag, maxerr, lstop, lprint, kdf, ntext, nidat, &
   nfdat, nmdat, mdat, tabspa, errcnt, ichar0, imag, inc, irc, &
   itext, iwf, iwg, kdi, kdj, kline, kshift, kspec, kt, lasti, &
   lbuf, lenlin, lenout, lentry, lentxt, locbeg, lstrt, ltext, &
   maxwid, mpt, nrow, ncol, ndim, ounit, gotfmt, xarg, xargok
common /cmessc / buf, dols, fmtf, fmtg, fmti, fmtj, fmtt, fmtim
equivalence (ivar(mevbas), sunit)
equivalence (fmtim(1), fmtr), (fmtim(2), fmtc)
!
equivalence (maxwid(1), linstr), (maxwid(2), kolwid)

save /cmessi/, /cmessc/
!++ CODE for .C. is inactive
!%%      long int kc;
!++ END
save lflgdb
data lflgdb / 2 /
!
if (ntext .ne. 0) then
   itext = (ntext-1) / lentxt
   ntext = ntext - itext * lentxt
   itext = itext + 1
end if
do 300 j = 1, max(1,kline)
   if (ntext .eq. 0) then
      k = kolwid
      go to 210
   end if
   lflgdb = 2
   ltext = 0
110    continue
!%%       ctmp=memchr(TEXT(cmessi.itext-1L,cmessi.ntext-1), SC,
!%%          cmessi.lentxt - cmessi.ntext + 1);
!%%       if (ctmp == NULL)
!%%             l = 0;
!%%       else
!%%             l = ctmp - TEXT(cmessi.itext-1L,cmessi.ntext-1) + 1;
   l = index(text(itext)(ntext:lentxt), sc)
   if (l .eq. 0) then
      ltext = ltext + lentxt - ntext + 1
      if (ltext .lt. 80) then
         itext = itext + 1
         ntext = 1
         go to 110
      end if
      ltext = 0
      if (kt .eq. 3) go to 310
      go to 160
   end if
   ntext = ntext + l + 1
   ltext = l + ltext - 1
   if (ntext .gt. lentxt) then
      itext = itext + 1
      if (ntext .eq. lentxt + 1) then
         c = text(itext-1)(lentxt:lentxt)
         ntext = 1
      else
         c = text(itext)(1:1)
         ntext = 2
      end if
   else
      c = text(itext)(ntext-1:ntext-1)
   end if
   if (c .eq. 'H') go to (180, 190, 200), kt
   if (c .eq. 'E') go to (180, 310, 200), kt
   if (c .eq. '#') go to (140, 310, 200), kt
   if (c .eq. ' ') then
!  Special code to set for removing the "$" preceding a blank.
      lstrdb = lstrt
      lflgdb = 3
      ltxtdb = ltext
      ltext = ltext + 1
      go to 110
   end if
   if (kt .ne. 1) go to 160
140    ltext = ltext + 2
   go to 110
160    kt = -kt
   go to 310
!
180    kolwid = kolwid + lenout
   if (ltext .eq. 0) go to 300
   kb = kolwid-ltext
   if (kb .lt. 0) stop &
   'Stopped in MESS -- Column width too small in a heading.'
   if (xarg)  kb = 1 + kb/2
   lstrt = lbuf + kb + 1
   lbuf = lbuf + kolwid
   if (lbuf .le. lenlin) mdat(nrow) = lbuf
   kolwid = 0
   go to 220
!
!                                  Set up column widths
190    maxwid(irc) = max(maxwid(irc), ltext)
   go to 300
!
!                                  Output matrix column
200    k = kolwid
   if (c .ne. '#') k = ltext
210    kb = lenout - kolwid
   if (j .eq. 1) then
!                        Special setup for the first column.
      if (xarg) kb = (kb + 1) / 2
      kb = kb + kshift + linstr - lbuf
   end if
   kb = kb + kolwid - k
   lstrt = lbuf + kb + 1
   lbuf = lstrt + k - 1
!                                  Set initial blanks
220    continue
!%%      if (kb > 0) for (kc=cmessi.lstrt-kb-1; kc<cmessi.lstrt-1; kc++)
!%%         cmessc.buf[kc] = ' ';
   if (kb .gt. 0) buf(lstrt-kb:lstrt-1) = ' '
!                                  Move characters
   if (ntext .eq. 0) then
!%%       memcpy(&cmessc.buf[cmessi.lstrt-1],"Col ", 4);
      buf(lstrt:lstrt+3) = 'Col '
      c = '#'
      lstrt = lstrt+4
   else
      k = ntext - ltext - lflgdb
      if (k .le. 0) then
         kk = max(0, 3-ntext)
!%%       memcpy(&cmessc.buf[cmessi.lstrt-1L], TEXT(cmessi.itext-2L,
!%%         cmessi.lentxt+k-1), -k-kk+1L);
         buf(lstrt:lstrt-k-kk)=text(itext-1)(lentxt+k:lentxt-kk)
         lstrt = lstrt-k-kk+1
         k = 1
      end if
      if (ntext .gt. 3) then
!%%       memcpy(&cmessc.buf[cmessi.lstrt-1L], TEXT(cmessi.itext-1L,
!%%         k-1), cmessi.ntext-k-2L);
         buf(lstrt:lstrt+ntext-k-3) = text(itext)(k:ntext-3)
         lstrt = lstrt + ntext - k - 2
      end if
   end if
   if (lflgdb .eq. 3) then
!  Special code to remove the "$" preceding a blank.  Only works for 1.
      do 250 l = lstrdb + ltxtdb + max(0, kb), lstrt
         buf(l:l) = buf(l+1:l+1)
250       continue
      lflgdb = 2
      lstrt = lstrt - 1
   end if
   if (c .eq. '#') then
!                                  Output column index
!%%         sprintf(&cmessc.buf[cmessi.lstrt-1], "%*ld ",
!%%           (int)(cmessi.lbuf-cmessi.lstrt), cmessi.imag+j-1);
      write (buf(lstrt:lbuf), fmtc) imag + j - 1
      if (ntext .ne. 0) ntext = k
      go to 300
   end if
!                                  Set trailing blanks
!%%      if (cmessi.lstrt <= cmessi.lbuf)
!%%           for (kc=cmessi.lstrt-1; kc < cmessi.lbuf; kc++)
!%%              cmessc.buf[kc] = ' ';
   if (lstrt .le. lbuf) buf(lstrt:lbuf) = ' '
300 continue
310 return
end

subroutine messpr
! Prints the buffer for MESS
!
! ************** Variable only used here *******************************
!
! NSCRN  Number of lines currently on CRT from messages.
!
integer   nscrn, k
character scrnam*12
save      nscrn
!
! For comments on other variables, see the listing for MESS.
integer   lenbuf, mevbas, mevlas
parameter (lenbuf=250)
parameter (mevbas=10)
parameter (mevlas=33)
logical          xarg, gotfmt, xargok
integer          errcnt, eunit, ichar0, irc, ivar(mevbas:mevlas), &
   imag, inc, itext, iwf, iwg, kdf, kdfdef, kdi, kdiag, kdj, &
   kline, kscrn, kshift, kspec, kt, maxerr, lasti, lbuf, lenlin, &
   lenout, lentry, lentxt, lhead, linerr, linmsg, locbeg, lprint, &
   lstop, lstrt, ltext, maxwid(2), mdat(5), mpt, munit, ncol, &
   ndim, nfdat, nidat, nmdat, nrow, ntext, ounit, sunit, tabspa
!
character buf*(lenbuf), dols*72, fmtc*7, fmtf*20, fmtg*15, &
  fmti*7, fmtim(2)*7, fmtj*7, fmtr*7, fmtt*15
common /cmessi/ sunit, lhead, kdfdef, linmsg, linerr, munit, &
   eunit, kscrn, kdiag, maxerr, lstop, lprint, kdf, ntext, nidat, &
   nfdat, nmdat, mdat, tabspa, errcnt, ichar0, imag, inc, irc, &
   itext, iwf, iwg, kdi, kdj, kline, kshift, kspec, kt, lasti, &
   lbuf, lenlin, lenout, lentry, lentxt, locbeg, lstrt, ltext, &
   maxwid, mpt, nrow, ncol, ndim, ounit, gotfmt, xarg, xargok
common /cmessc / buf, dols, fmtf, fmtg, fmti, fmtj, fmtt, fmtim
equivalence (ivar(mevbas), sunit)
equivalence (fmtim(1), fmtr), (fmtim(2), fmtc)
!
save /cmessi/, /cmessc/
data nscrn / 0 /
!
if (lbuf .ne. 0) then
10   if (buf(lbuf:lbuf) .eq. ' ') then
    if (lbuf .gt. 1) then
      lbuf = lbuf - 1
      go to 10
    end if
  end if
  if (ounit .le. 0) then
    if (kscrn .gt. 0) then
      if (nscrn .ge. kscrn) then
!%%               printf( " Type 'Enter' to continue\n" );
        print '('' Type "Enter" to continue'')'
!%%               scanf( "%*[^\n]%*c" );
        read (*, *)
        nscrn = 0
      end if
      nscrn = nscrn + 1
    end if
!%%      printf( "%.*s\n", (int)cmessi.lbuf, cmessc.buf);
    print '(1X, A)', buf(1:lbuf)
    if (ounit .eq. 0) go to 20
  end if
!++ CODE for ~.C. is active
  k = abs(ounit)
  write (k, '(A)', err=30) buf(1:lbuf)
!++ CODE for .C. is inactive
!%%      fprintf(c_handle[labs(cmessi.ounit)-1], "%.*s\n",
!%%      (int)cmessi.lbuf, cmessc.buf);
!++ END
20   lbuf = 0
end if
return
!++ CODE for ~.C. is active
!              See if opening fixes the error
30 write(scrnam, '(A, I2.2, A)') 'MESSF_', k, '.tmp'
open (unit=k, status='UNKNOWN', file=scrnam)
write (k, '(A)') buf(1:lbuf)
return
!++ END
end

subroutine messft(mact, ftext)
!  Prints FTEXT, which contains a Fortran character string, and then
!  call MESS to do the actions in MACT.  Actions in MACT can not do
!  anything other than actions that reference MACT.
!  This routine intended for use by library subroutines getting text in
!  the form of a Fortran character string.
!
integer mact(*)
character ftext*(*)
!
integer j, k, idat(1), mecont, meprnt, mesuni
!++ CODE for ~.C. is active
character text(1)*1
!++ CODE for .C. is inactive
!      character TEXT(1)*2
!++ END
parameter (mesuni=10, meprnt=21, mecont=50)
!
integer lenbuf, mevbas, mevlas
parameter (lenbuf=250)
parameter (mevbas=10)
parameter (mevlas=33)
logical          xarg, gotfmt, xargok
integer          errcnt, eunit, ichar0, irc, ivar(mevbas:mevlas), &
   imag, inc, itext, iwf, iwg, kdf, kdfdef, kdi, kdiag, kdj, &
   kline, kscrn, kshift, kspec, kt, maxerr, lasti, lbuf, lenlin, &
   lenout, lentry, lentxt, lhead, linerr, linmsg, locbeg, lprint, &
   lstop, lstrt, ltext, maxwid(2), mdat(5), mpt, munit, ncol, &
   ndim, nfdat, nidat, nmdat, nrow, ntext, ounit, sunit, tabspa
!
character buf*(lenbuf), dols*72, fmtc*7, fmtf*20, fmtg*15, &
  fmti*7, fmtim(2)*7, fmtj*7, fmtr*7, fmtt*15
common /cmessi/ sunit, lhead, kdfdef, linmsg, linerr, munit, &
   eunit, kscrn, kdiag, maxerr, lstop, lprint, kdf, ntext, nidat, &
   nfdat, nmdat, mdat, tabspa, errcnt, ichar0, imag, inc, irc, &
   itext, iwf, iwg, kdi, kdj, kline, kshift, kspec, kt, lasti, &
   lbuf, lenlin, lenout, lentry, lentxt, locbeg, lstrt, ltext, &
   maxwid, mpt, nrow, ncol, ndim, ounit, gotfmt, xarg, xargok
common /cmessc / buf, dols, fmtf, fmtg, fmti, fmtj, fmtt, fmtim
equivalence (ivar(mevbas), sunit)
equivalence (fmtim(1), fmtr), (fmtim(2), fmtc)
!
do 10 j = 1, 100, 2
   k = abs(mact(j))
   if ((k .gt. meprnt) .or. (k .lt. mesuni)) go to 20
10 continue
20 k = mact(j)
mact(j) = mecont
call mess(mact, text, idat)
mact(j) = k
!%%      k = strlen(ftext);
k = len(ftext)
ntext = 1
if (k .ne. 0) then
   if (ftext(1:1) .eq. '0') then
      ntext = 2
      k = k - 1
      if (lbuf .eq. 0) then
         buf(1:1) = ' '
         lbuf = 1
      end if
   end if
   call messpr
   lbuf = k
!%%      memcpy(cmessc.buf, &ftext[cmessi.ntext-1], k);
   buf(1:k) = ftext(ntext:ntext+k-1)
end if
ichar0 = ichar('0')
if (mact(j) .ne. mecont) call mess(mact(j), text, idat)
return
end
subroutine optchk(intchk, iopt, etext)
! Copyright (c) 1996 California Institute of Technology, Pasadena, CA.
! ALL RIGHTS RESERVED.
! Based on Government Sponsored Research NAS7-03001.
!>> 1998-11-01 OPTCHK  Krogh  ERRSEV => MACT(2) for "mangle".
!>> 1996-05-13 OPTCHK  Krogh  Changes to use .C. and C%%.
!>> 1995-03-10 OPTCHK  Krogh  Added "abs(.) just below "do 140 ..."
!>> 1994-11-11 OPTCHK  Krogh  Declared all vars.
!>> 1993-05-17 OPTCHK  Krogh  Additions for Conversion to C.
!>> 1991-11-25 OPTCHK  Krogh  More comments, little clean up of code.
!>> 1991-10-09 OPTCHK  Krogh  More comments, little clean up of code.
!>> 1991-06-27 OPTCHK  Krogh  Initial Code.
!
! OPTCHK -- Fred T. Krogh, Jet Propulsion Lab., Pasadena, CA.
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
!            >0          .ne. -1       Yes
!            >0          .eq. -1       No, and IOPT(INTCHK(k+1)) is set
!                                      to starting loc. of space found.
!        When this program finds the space for an option, values of
!        INTCHK(j) for j .ge. LAST will be used.  INTCHK(k+1) is set
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
!        In addition if values of INTCHK(j) for j .ge. LAST are used,
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
! ************************** Variable Declarations *********************
!
integer intchk(0:*), iopt(*)
character etext*(*)
integer i, istrt, kex, l, last, lneg, lopt, lwant, mi, mv, n
! Declarations for error messages.
integer mentxt, meidat, mecont, meret, meemes, metext, meimat, &
   ltxtee, ltxend
parameter (mentxt =23)
parameter (meidat =24)
parameter (mecont =50)
parameter (meret =51)
parameter (meemes =52)
parameter (metext =53)
parameter (meimat =58)
integer mact(16), errbad, nerbad(0:7)
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
integer ltxtaa,ltxtab,ltxtac,ltxtad,ltxtae
parameter (ltxtaa=  1,ltxtab=  9,ltxtac=233,ltxtad=267,ltxtae=295)
character mtxtaa(2) * (156)
!                          Next 4 lines not automatically generated
!%%    #define LTXTEE  137
parameter (ltxtee = ltxtae - 156 - 2)
parameter (ltxend = 156)
!
data mtxtaa/'OPTCHK$B"Option #" is negated if option needs attenti &
on.$n"Option 0" is for space not associated with a specific option &
.$n"First Loc." is negated if user di','d not set value.$NSpace av &
ail. = $i; all options must have first loc. > $i$EOption #$HFirst$ &
 Loc.$HLast Loc.$EFrom subprogram/argument: $BSpace for etext.$e'/
!
!                      1 2 3      4       5  6       7      8       9
data mact / meemes,0,1,ltxtad, meidat, 2, mentxt,ltxtab, metext, &
   meimat,3,3,0,ltxtac,-1, meret /
!            10    13     14 15     16
data nerbad / 57, 17, 57, 17, 0, 0, 7, 7 /
!
! *************************** Start of Executable Code *****************
!
mact(3) = intchk(0) / 10
mact(16)=meret
i = intchk(0) - 10*mact(3)
if (i .gt. 3) then
   i = i - 4
   mact(16) = mecont
end if
errbad = nerbad(i)
mact(2) = nerbad(i+4)
last = intchk(1)
kex = last
20 lneg = 0
do 100 i = 5, last, 3
!       Loop to sort on the low indices -- Inefficient algorithm to keep
!       code short -- LAST should never be very big.
   mi = i
   mv = intchk(i)
   do 50 l = i+3, last, 3
!                                    Find mimimum from here down.
      if (intchk(l) .lt. mv) then
         mi = l
         mv = intchk(l)
      end if
50    continue
   if (mi .ne. i) then
!                                   Interchange to get low at top.
      do 70 l = -1, 1
         n = intchk(i+l)
         intchk(i+l) = intchk(mi+l)
         intchk(mi+l) = n
70       continue
   end if
   if (mv .lt. 0) then
!                Save loc. of last entry that needs space to be found.
      lneg = i
   else if (lneg .eq. 0) then
!        Entry I and previous entries are in their correct sorted order.
      if (intchk(i+1) .lt. 0) then
         if (intchk(i-1) .lt. -1000000) then
            intchk(i-1) = intchk(i-1) + 2000000
            intchk(i+1) = -intchk(i+1)
!                            Save INTCHK index defining allocated space.
            kex = kex + 1
            intchk(kex) = i - 1
         else
!                   Error -- Got request for a negative amount of space.
            mact(2) = errbad
            intchk(i-1) = -abs(intchk(i-1))
         end if
      end if
!                Save final location used by the option.
      intchk(i+1) = intchk(i) + intchk(i+1) - 1
      if (intchk(i) .le. intchk(i-2)) then
!                                           Error -- options overlap.
         intchk(i-1) = -abs(intchk(i-1))
         mact(2) = errbad
      end if
   end if
100 continue
if (lneg .ne. 0) then
!     Find spaces that need to be allocated, starting with the smallest.
   istrt = lneg
   i = lneg
120    lwant = intchk(lneg)
   lopt = intchk(lneg+1)
   if (i .eq. lneg) then
!                         Make fake entry to get started.
      intchk(lneg) = 1
      intchk(lneg+1) = intchk(3)
   end if
   do 140 istrt = istrt, last-3, 3
      if(intchk(i)+abs(intchk(i+1))-lwant .lt. intchk(istrt+3)) &
         go to 150
      i = istrt + 3
140    continue
150    intchk(lneg) = intchk(i) + abs(intchk(i+1))
   if (lopt .ne. 0) then
      if (lopt .gt. 0) then
         if (iopt(lopt) .eq. -1) then
            iopt(lopt) = intchk(lneg)
            go to 160
         end if
      end if
!                     Error -- IOPT not expecting starting loc.
      intchk(lneg-1) = -abs(intchk(lneg-1))
      mact(2) = errbad
   end if
160    intchk(lneg+1) = lwant
   intchk(lneg-1) = intchk(lneg-1) - 2000000
   if (lneg .lt. 8) go to 20
   i = lneg
   lneg = lneg - 3
   go to 120
end if
if (intchk(last-1) .gt. intchk(2)) then
   if (intchk(2) .lt. 0) go to 180
   if (last .ne. kex) then
      if (intchk(kex) .eq. last - 3) then
         if (intchk(last) .le. 0) then
            if (intchk(last-2)-intchk(last)-1 .le. intchk(2)) then
               intchk(last-1) = intchk(2)
               go to 180
            end if
         end if
      end if
   end if
   intchk(last-3) = -abs(intchk(last-3))
   mact(2) = errbad
end if
180 if (last .ne. kex) intchk(last) = kex
if (mact(2) .gt. 0) then
190    if (last .ne. kex) then
      do 200 i = last+1, abs(kex)
         intchk(intchk(i)+1) = -intchk(intchk(i)+1)
200       continue
      if (kex .lt. 0) go to 210
      kex = -kex
   end if
   mact(13) = (last - 4) / 3
!%%       strcpy(&mtxtaa[1][LTXTEE-1], etext);
   mtxtaa(2)(ltxtee:ltxend)=etext(1:)
   call mess(mact, mtxtaa, intchk(1))
   if (mact(2) .gt. 10) intchk(1) = -last
   if (last .ne. kex) go to 190
end if
210 intchk(2) = intchk(last-1)
return
end
