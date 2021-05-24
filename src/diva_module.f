      subroutine DIVA(TSPECS, Y, F, KORD, NEQ, DIVAF, DIVAO, IDIMT,
     1   IDIMY, IDIMF, IDIMK, IOPT)
c Copyright (c) 1996 California Institute of Technology, Pasadena, CA.
c ALL RIGHTS RESERVED.
c Based on Government Sponsored Research NAS7-03001.
c>> 2015-03-15 DIVA  Krogh  Removed extra call divabu after noise test
c>> 2015-03-15 DIVA  Krogh  Forced restart needs more reduction in h.
c>> 2010-02-20 DIVA  Krogh  Fixed calling DIVAOP with array other than F.
c>> 2009-11-03 DIVA  Krogh  Added option 11, more variables initialized.
c>> 2009-10-30 DIVA  Krogh  Gave KSSTRT and ROBND initial values.
c>> 2009-10-30 DIVA  Krogh  Fixed reference to undefined location in F.
c>> 2009-10-21 DIVA  Krogh  Got rid of NaN in diag. print when LSC=3.
c>> 2009-10-15 DIVA  Krogh  A few changes on how noise is handled.
c>> 2002-11-12 DIVA  Krogh  Fixed problem integrating to final output pt
c>> 2002-08-29 DIVA  Krogh  Added test for invalid HMIN/HMAX.
c>> 2002-07-26 DIVA  Krogh  Added KOUTKO to fully support Option 10.
c>> 2002-05-14 DIVA  Krogh  Fix starting prob. for Option 18.
c>> 2002-05-13 DIVA  Krogh  Put exponent letter in  numbers missing them
c>> 2002-05-12 DIVA  Krogh  Added error message for bad option 5 usage.
c>> 2001-09-07 DIVA  Krogh  Changes to allow user tol on G-Stops.
C>> 2001-05-25 DIVA  Krogh  Minor change for making .f90 version.
c>> 2001-05-18 DIVA  Krogh  Less computing with no error test
c>> 2001-05-17 DIVA  Krogh  Fixed so with no error test can't start dump
c>> 2001-04-24 DIVA  Krogh  Inserted comments from ivacom.
c>> 2000-12-01 DIVA  Krogh  Removed (some of) unused C1, MAXSTF, METEXT.
c>> 1999-12-28 DIVA  Krogh  Saved S in DIVACR for output consistency.
c>> 1999-08-19 DIVA  Krogh  Removed superfluous test above label 3520.
c>> 1997-04-22 DIVA  Krogh  Got rid of assigned go to's. F=0 if diag.
c>> 1996-08-26 DIVA  Krogh  Initialize F to 0 if dumping solution.
c>> 1996-08-23 DIVA  Krogh  Print TN not TSPECS(1) in error messages.
c>> 1996-05-30 DIVA  Krogh  Changed DERIVS/OUTPUT to  DIVAF/DIVAO.
c>> 1996-04-27 DIVA  Krogh  Changes to use .C. and C%%.
c>> 1996-03-30 DIVA  Krogh  Added external statement.
C>> 1996-03-25 DIVA  Krogh  Introduced TEXT1 to comply with F77.
c>> 1996-02-27 DIVA  Krogh  Fixed so DUMP not affected by ignored eqs.
c>> 1995-12-18 DIVA  Krogh  Fixed so no solution dump on 0 length integ.
c>> 1995-11-09 DIVA  Krogh  Fixed so char. data at col. 72 is not ' '.
c>> 1995-06-19 DIVA  Krogh  Fixed prob. with discon. just after restart.
c>> 1995-05-09 DIVA  Krogh  Fixed G-Stop/discontinuity code interaction
C>> 1995-04-26 DIVA  Krogh  Use KQMAXS instead of KQMAXI when LDIS>1000.
C>> 1995-04-26 DIVA  Krogh  Keep current KQL on discontinutiy.
C>> 1994-12-16 DIVA  Krogh  Fixed option 12 with K12 < 0.
C>> 1994-11-11 DIVA  Krogh  Declared all vars.
c>> 1994-11-02 DIVA  Krogh  Changes to use M77CON
c>> 1994-09-08 DIVA  Krogh  Added CHGTYP code.
c>> 1994-07-11 DIVA  Krogh  Fix to get same state with/without var. eqs.
c>> 1994-03-07 DIVA  Krogh  Allow larger order in single precision.
c>> 1994-01-14 DIVA  Krogh  Minor change to allow changing TFINAL.
c>> 1993-04-27 DIVA  Krogh  Additions for Conversion to C.
c>> 1993-04-12 DIVA  Krogh  Converted to use slightly altered MESS.
c>> 1993-04-12 DIVA  Krogh  Fixed LSC so sol. saved when HMAX is small.
c>> 1992-10-13 DIVA  Krogh  Fixed G-Stop/discontinuity code interaction.
c>> 1992-09-21 DIVA  Krogh  Fixed bug in discontinuity code.
c>> 1992-09-09 DIVA  Krogh  Fixed bug - Var. Eqs. with discontinuities.
c>> 1992-08-07 DIVA  Krogh  Storage map printed only if option 10 .ne. 0
c>> 1992-07-16 DIVA  Krogh  Restored correct discontinuity code.
c>> 1992-06-16 DIVA  Krogh  Eliminate reuse of storage for option 12.
c>> 1992-04-08 DIVA  Krogh  Removed unused labels, 1020, 2120.
c>> 1992-03-30 DIVA  Krogh  Fixed bug in DIVAOP error message.
c>> 1992-03-12 DIVA  Krogh  Simplified DIVABU, more digits in B's.
c>> 1992-01-16 DIVA  Krogh  Fixed minor bug in error messages.
c>> 1991-12-03 DIVA  Krogh  Major change for improved error checks.
c>> 1991-06-17 DIVA  Krogh  Fixed bug in checking storage allocation.
c>> 1991-04-11 DIVA  Krogh  Fixed minor bug re. option 12 in DIVAOP.
c>> 1991-03-28 DIVA  Krogh  Removed check at label 650 for KORD2I<0.
c>> 1991-02-08 DIVA  Krogh  Changed some floats to generics
c>> 1990-11-08 DIVA  Krogh  Fixed bug on TSPECS on discon.
c>> 1990-09-14 DIVA  Krogh  Fixed bug when discon. and sol. save.
c>> 1990-09-13 DIVA  Krogh  Increased dimension of BETA by 1.
c>> 1990-09-13 DIVA  Krogh  Added one more poss. on rel. error test.
c>> 1990-09-11 DIVA  Krogh  Recent change messed up getting dump output.
c>> 1990-06-05 DIVA  Krogh  Fixed bug in noise test, comments in IVACOM.
c>> 1990-05-08 DIVA  Krogh  Fixed new bug when TMARK hit in DIVAG.
c>> 1990-04-17 DIVA  Krogh  Fixed minor problem in DIVAIN error msg.
c>> 1990-04-10 DIVA  Krogh  Fixed interaction between discon. & dump.
c>> 1990-03-23 DIVA  Krogh  Fixed bug on option "-2", see 1989-12-07.
c>> 1990-03-20 DIVA  Krogh  Fixed rarely occuring loop.
c>> 1990-01-29 DIVA  Krogh  Removed unneeded labels.
c>> 1989-12-14 DIVA  Krogh  Saved common block DIVAEV.
c>> 1989-12-07 DIVA  Krogh  Added option "2" to DIVAOP.
c>> 1989-11-09 DIVA  Krogh  Made GG a save var. in DIVAHC
c>> 1989-08-21 DIVA  Krogh  Fix out of bounds ref. to V in DIVABU
c>> 1989-07-26 DIVA  Krogh  Fix bug in initial dim. check
c>> 1989-07-21 DIVA  Krogh  Code for integrating discontinuities
c>> 1987-12-07 DIVA  Krogh  Initial code.
c
c--D replaces "?": ?IVA,?IVAA,?IVABU,?IVACO,?IVACR,?IVAEV,?IVAF,?IVAHC,
c-- & ?IVAG,?IVAIN,?IVAMC,?IVAO,?IVAOP,?IVAPR,?IVASC,?IVACE,?IVAIE,
c-- & ?IVAPE,?MESS
c
c
c Note a "*" at the start of a name is used to indicate "D" for the
c double precision version and "S" for the single precision version.
c
c When converting between precisions, don't forget to change the value
c of KDIM set in parameter statements in a variety of routines, and to
c adjust comments for the data statements associated with EIBND in
c *IVACR, and B in *IVAHC.
c
c Entries
c  *IVA    Main entry for starting the package.
c  *IVAA   Main program inside the package, calls the other routines,
c          and does checks for output, and noise.  Called by the user
c          if reverse communication is used.
c  *IVABU  Back ups the solution to the current base time, if a step
c          that has been started must be taken over for some reason.
c  *IVACO  Called by user to get certain information from the common
c          blocks.
c  *IVACR  Corrects the solution, estimates errors, and selects order.
c  *IVADB  Subroutine to assist in debugging codes.  Called by user to
c          get a formatted list of all the variables used in the
c          integration.  Not required in usual case.
c  *IVADE  Needed only for delay differential equations.  This is called
c          by the user from the derivative subprogram.
c  *IVAG   Required only if the user has G-Stops, i.e. places to call
c          his output subroutine where certain functions have zeroes.
c  *IVAHC  Compute coefficients that depend on the step size history.
c  *IVAIN  Used to interpolate to arbitrary points.
c  *IVAOP  Used to process user option requests.
c  *IVAPR  Used to update the differences and to predict the solution
c          at the end of the current step.
c
c External Routines
c  *1MACH  Not used in the Fortran 95 version.  ("*" is "D" for double
c          and "R" for single precision.) This returns constants that
c          depend on the floating point arithmetic.  Input arguments of
c          1 to 4 give respectively:  underflow limit, overflow limit,
c          smallest relative difference between two floating point
c          numbers, and the largest relative difference between two
c          floating point numbers.
c DERIVS (formal) Name of subroutine to be called for computing

c  OPTCHK  Used in checking storage allocation.
c  *MESS   Used to output error messages and diaganostic messages.
c          (Just MESS if no floating point is output.)
c  *ZERO   Called only if *IVAG is used.  Iterates to find zeros of
c          arbitrary (continuous) functions.
c
c Common blocks -- As a left over from the distant past, some variables
c   are in common so that they would be saved.
c  *IVAEV  Holds variables that depend on the environment.
c  *IVAMC  The main common block for the package.
c  *IVASC  The secondary common block for the package.  This contains
c          variables that are required for doing interpolation and is
c          separate to simplify saving the variables that are required
c          when the solution is being dumped (saved).
c
c Common variables and local variables
c ALPHA  (*IVAMC) Array with I-th entry = (current step size) / XI(I).
c   Used in computing integration coefficients.
c B      (*IVAHC) Array used to get started on computing integration
c   coefficients.  B(K) = 1. / (K*(K+1))
c BAKMIN (*IVADE) The largest delay at the initial point.
c BETA   (*IVAMC) Array with I-th entry = product (K=1,I-1) of
c   (current (XI(K)) / XI(K) from previous step),  BETA(1)=1.  Used in
c    updating the difference tables.
c C      (*IVAIN) Array used to hold integration/interpolation coeffs.
c C0     Parameter = 0. (in *IVAA,DE,CR,A,G,HC,IN,OP,PR)
c C1     Parameter = 1. (in *IVA,A,CR,DA,HC,IN,OP)
c C10    Parameter = 10. (in *IVAA,CR,OP)
c C1000  Parameter = 1000. (in *IVACR)
c C16    Parameter = 16. (in *IVAA,OP)
c C1M3   Parameter = .001 (in *IVAA)
c C1M5   Parameter = .00001 (in *IVAA)
c C1P125 Parameter = 1.125 (in *IVAA,HC,OP)
c C1P3   Parameter = 1.3 (in *IVAA)
c C1P4   Parameter = 1.4 (in *IVACR)
c C2     Parameter = 2. (in *IVAA,DE,BU,CR,IN,OP)
c C20    Parameter = 20. (in *IVACR)
c C2P5M3 Parameter = .0025 (in *IVAA)
c C4     Parameter = 4. (in *IVACR,OP)
c C40    Parameter = 40. (in *IVACR)
c C4096  Parameter = 4096. (in *IVAA)
c C6     Parameter = 6. (in *IVAA)
c C8M3   Parameter = .008 (in *IVAA)
c CM2    Parameter = -2. (in *IVACR)
c CM8    Parameter = -8. (in *IVACR)
c CMP5   Parameter = -.5 (in *IVACR)
c CMP75  Parameter = -.75 (in *IVAOP)
c CP0625 Parameter = .0625 (in *IVAA)
c CP1    Parameter = .1 (in *IVAA,CR,DA,HC)
c CP125  Parameter = .125 (in *IVACR)
c CP25   Parameter = .25 (in *IVAA,CR,DE,OP)
c CP3    Parameter = .3 (in *IVAA,OP)
c CP4    Parameter = .4 (in *IVAA)
c CP5    Parameter = .5 (in *IVAA,CR,DA,DE,HC,OP)
c CP5625 Parameter = .5625 (in *IVAHC)
c CP625  Parameter = .625 (in *IVAOP)
c CP75   Parameter = .75 (in *IVACR,OP)
c CP8    Parameter = .8 (in *IVACR)
c CP875  Parameter = .875 (in *IVAA, OP)
c CP9    Parameter = .9 (in *IVAOP)
c CP9375 Parameter = .9375 (in *IVACR)
c CQ3125 Parameter = .03125 (in *IVACR)
c CRBQI  Parameter = .421875 (in *IVAHC)  Initial val for computing RBQ.
c CSUM   (*IVAIN) Array used to contain partial sums of the integration
c   coefficients.  This is used to corrrect for a difference table that
c   has not yet been updated.
c D      (*IVAMC) Array to be used later to store coefficients for
c   integrating stiff equations.
c   derivatives.  Not used if option 13 is set.
c DISADJ (*IVAA) Value of stepsize when discontinuity is indicated.
c DNOISE (*IVAMC) Used in determining if noise is limiting the
c   precision.  It is usually |highest difference used in correcting|
c   of the equation with the largest error estimate.
c DS     (*IVAMC) Array to be used later to store coefficients for
c   estimating errors when integrating stiff equations.
c DVC2   (*IVADB) Array used for output of variables HC to TOUT in
c   common block *IVAMC.
c E      (*IVACR) (Estimated error) / (Requested accuracy)
c EAVE   (*IVAMC) This is a weighted average of past values of EIMAX.
c   It is adjusted to account for expected changes due to step changes.
c EEPS10 (*IVAEV) = 10. * (machine epsilon).
c EEPS16 (*IVAEV) = 16. * (machine epsilon).
c EEPS2  (*IVAEV) =  2. * (machine epsilon).
c EEPT75 (*IVAEV) = (machine epsilon) ** (.75)
c EI     (*IVACR) Estimate for what E would be if step size increased.
c EIBND  (*IVACR) Array containing limits on the estimated error with
c   the stepsize increased.  This array tends to make the code a little
c   more conservative on step size increases at low order.
c EIMAX  (*IVAMC) Estimate of (error estimate / error requested) if the
c   step size should be increased.
c EIMIN  (*IVAMC) An error estimate is small enough to allow a step
c   increase if the estimate of ((error with the step size increased) /
c   (error requested)) is less than EIMIN.
c EIMINO (*IVAA) Set to C8M3 and never changed.  When step size is being
c   reduced if EIMIN .le. EIMINO then the reduction factor is set to
c   CP875.  This variable could be a parameter.
c EMAX   (*IVAMC) Largest value computed for (error estimate) / (error
c   requested).
c EOVEP2 (*IVAEV) = EEPS2 * (largest floating point number).
c EPS    (*IVACR) Current absolute error tolerance.  Also used for
c   temporary storage when computing the desired value of EPS.
c ERCOEF (*IVACR) (Error coefficient from formula) / EPS
c EREP   (*IVAMC) If EMAX > EREP, a step is repeated.  Ordinarily
c   this has the value .3.  This is set < 0 if the error tolerance is
c   specified improperly, and is set to a large value if the user
c   requests complete control over the step size.  EREP is also set
c   < 0 after a user specified discontinuity.
c EROV10 (*IVAEV) = 10. / (largest floating point number).
c ETA    (*IVAIN) Array used in computing integration/interp. coeffs.
c EVC    (*IVADB) Array used for output of variables EEPS2 to EROV10 in
c   common block *IVAEV.
c EXR    (*IVAA) Set to CP1 and never changed.  If it is estimated the
c   the (error estimate) / (error requested) on the next step will be
c   .ge. EXR then the step size is reduced.  Could be a parameter.
c F      (formal) Array used to store derivative values, the difference
c   tables, error tolerance requests, and values used by some other
c   options. (in *IVA,A,BU,CR,DA,DB,G,IN,PR)
c FDAT  (*IVAMC) Used to store data for error messages.  (Local array in
c   *IVAIN.)
c FOPT  (formal) in *IVAOP.  Passed as place to save floating point data
c   for options.  This package passes F in for FOPT when calling *IVAOP.
c G      (*IVAMC) Integration coefficients used for predicting solution.
c   G(I, J) gives the I-th coefficient for integrating a J-th order
c   differential equation.  G(1, 1) is equal to the step size.
c GAMMA  (*IVAIN) Array used in computing integration/interp. coeffs.
c GG     (*IVAHC) Array of length = max. differential equation order
c   allowed by code - 1.  GG(K) = (HH**(K+1)) / K!
c GNEW   (formal) in *IVAG.  Current value for vector function g, whose
c   zeroes are to be found.
c GOINT  (*IVACR) Used for assigned go to used in computing integration
c   coefficients.
c GOLD   (*IVAG) Previous value for element of G whose zero search is
c   active.
c GS     (*IVAMC) Integration coefficients used in estimating errors.
c GT     (formal) in *IVAG.  Previous value of GNEW.
c HC     (*IVAMC) Ratio of (new step size) / (old step size)
c HDEC   (*IVAMC) Default value to use for HC when reducing the step
c   size.  (Values closer to 1 may be used some of the time.)
c HH     Equivalenced to G(1,1) = current step size in *IVAA,CR,DA,G,HC.
c HI     (*IVAIN) Step length from the base value of the independent
c   variable for the interpolation.
c HINC   (*IVAMC) Default value to use for HC when increasing the step
c   size.  (Values closer to 1 may be used some of the time.)
c HINCC  (*IVAMC) Actual value used for default value of HC when
c   increasing the step size.  Set to HINC after start is considered
c   complete.  During the start HINCC is set to 1.125.
c HMAX   (*IVAMC) Largest value allowed for abs(step size).  Default
c   value is a very large number.
c HMAXP9 (*IVAMC) .9 * HMAX.
c HMIN   (*IVAMC) Smallest value allowed for abs(step size).  Default
c   value is 0.
c HNEW   (*IVADE) Value of step size when iterating at initial point
c   for delay differential equations.
c I      Used for temporary storage. (*IVAA,BU,CR,DA,DE,G,IN,OP,PR)
c IA     (*IVAOP) absolute value of first integer stored for an option.
c ICF    (*IVAMC) Final index for current loop in *IVACR.  Required by
c   option 18.
c ICI    (*IVAIN) Temporary index, = 0 for interpolation, 1 or 0 for
c   differentiation, and d-1, d-2, ... 0 for integration, where d is the
c   order of the differential equation.  Index of first location
c   in C() used is ICI + an offset.
c ICS    (*IVAMC) Starting index for current loop in *IVACR.
c ID     (formal) Array use to contain integer data from common.  Values
c   are returned in locations 1 to 5 as follows.
c   1    KEMAX = Index of equation with largest error estimate
c   2    KSTEP = Current step number
c   3    NUMDT = Number of differences used for each equation
c   4            Reserved for future use
c   5            Reserved for future use
c IDAT   (*IVAMC) Used to store integer for error messages.  (Also used
c   in *IVAA for temporary storage of KORD(2).  (Local array in *IVAIN.)
c IDE    (*IVADE - formal) Array used to contain past information so
c   that delays can stretch back indefinitely.  If the first location is
c   0, then any interpolations requested must be in the range of the
c   current difference tables.  At present, only the value 0 is allowed
c   in IDE(1).  This array is intended for the support of saving long
c   past histories.  IDE(2) must contain the declared dimension of WDE.
c IDEF   (*IVADE -- formal) Flag giving indicaion of what is going on.
c   = 0  User should compute derivatives and return to the main
c        integrator.
c   = 1  Code is computing additional values in order to get past data
c        necessary for starting.  User should compute derivatives and
c        call *IVADE.
c   < 0  Indicates an error condition.  If *IVADE is called without
c        changing the value of IDEF, the integration is stopped and an
c        error message printed.  Possible error flags are:
c    -1  Difference tables do not span back far enough to compute the
c        past values of Y needed.
c    -2  There is not enough space in WDE to get the required starting
c        values.
c IDIMF  (formal) Declared dimension of F().
c IDIMK  (formal) Declared dimension of KORD().
c IDIMT  (formal) Declared dimension of TSPECS().
c IDIMY  (formal) Declared dimension of Y().
c IDT    (*IVAIN) Used as a base index into the difference table.
c IFLAG  (formal in *IVAG) Used for communication with user.
c   = 1  Continue as if *IVAG was not called.
c   = 2  Check KORD(1) as one would do at start of OUTPUT if no G-Stops
c        were present. (Exit if in DERIVS.)
c   = 3  Return to the integrator.
c   = 4  Compute G and return to *IVAG.
c   = 5  A G-Stop has been found, and NSTOP gives its index.  (If NSTOP
c        < 0, the stop was an extrapolating stop.)
c   = 6  Same as 5, but requested accuracy was not met.
c   = 7  Same as 5, but there is a probable error in computing G.
c   = 8  Fatal error of some type.  (An error message has been printed.)
c IG     (*IVAG)  IG = KORD(2) on the initial entry (0 for extrapolating
c   G-Stops, and 1 for interpolating).
c IGFLG  (*IVAMC) Used primarily in *ivag, but also used in *iva to keep
c   track of the state of GSTOP calculations.
c   = -2 Extrapolatory G's initialized, but not the interpolatory.
c   = -1 Interpolatory G's initialized, but not the extrapolatory.
c   =  0 Set when integration is started or restarted, or option setting
c        GSTOP is set.
c   =  1 Iterating to find a GSTOP.
c   =  2 User told that a GSTOP was found.
c   =  3 Checking G's at point where a GSTOP was located.
c   =  4 Checking G's at a T output point.
c   =  5 Usual case, no sign change detected.
c IGSTOP (*IVAMC) IGSTOP(k) is set in *ivag to the index of the last G
c   with a 0, where k is one for an interpolatory G-Stop, and k is two
c   for an extrapolatory G-Stop.
c IGTYPE (*IVAMC) Array with two elements as for IGSTOP, but this saves
c   a flag giving the nature of convergence to the stop.
c   = 0  All known G-stops completely processed.
c   = 4  Need to compute next value while iterating.
c   = 5  Got good convergence.
c   = 6  Got convergence, but not to desired accuracy.
c   = 7  Problem in getting convergence.
c   = 8  A fatal error of some type.
c IHI    (*IVA) Last location used by the current option.
c ILGREP (*IVAMC) Used when correction to keep track of equations that
c   are to use a certain error tolerance.
c ILGROR (*IVACR) Index of last equation in the current group of
c   equations grouped for selecting integration order.
c ILOW   (*IVA) First location used by the current option.
c INCOM  (*IVADE) Array equivalenced to LDT in the common block *IVASC.
c   Used to simplify saving information in the common block.
c INCOP  (*IVAOP) Array containing data giving the amount of space in
c   IOPT used for each of the options.
c INGS   Current index for G-stop being examined in DIVAG.
c INICAS (*IVADE) Used to track the initialization for a delay equation.
c   = 1  Very beginning.
c   = 2  Getting derivative at the very beginning.
c   = 3  Getting derivatives at points prior to the initial point.
c   = 4  Getting derivative at initial point after iteration is started.
c INTCHK (*IVA) Array passed to OPTCHK containing information on storage
c   allocation.  See comments in OPTCHK for details.
c INTEG  (*IVAIN) Number of integrations being done. (<0 for
c   differentiations and =0 for interpolation.)  Also used as counter
c   when computing integration coefficients.
c        (*IVAPR) Number of integrations being done.
c INTEGS (*IVAPR) = -1 for equations that are not stiff, 0 for those
c   that are stiff.
c INTEGZ (*IVAIN) min(INTEG, 0)
c INTERP (*IVAIN) added to the usual integration order to get the order
c   to be used when interpolating: 3-KQMAXI, if HI=0; 1, if
c   |HI| > |XI(1)| and HI * XI(1) < 0; 0, otherwise -- the usual case.
c IOP10  (*IVAMC) Number of times diagnostic output is to be given when
c   leaving *ivacr (the corrector).
c IOP11  (*IVAMC) Gives current step number of the method.  Tells how
c   many of certain coefficients must be computed. (Has nothing to do
c   with options.) = min(max integ order + 1, KDIM).  Also set when
c   starting to flag that certain memory locations must be set to 0.
c IOP12  (*IVAMC) Points to location in F() where user supplied values
c   of HINC, HDEC, HMIN, and HMAX.  (0 if option 12 not used.)
c IOP13  (*IVAMC) If not zero, reverse communication will be used for
c   getting the values of derivatives.  Associated with option 13.
c IOP14  (*IVAMC) If not zero, reverse communication will be used in
c   place of calls to the output routine.  Associated with option 14.
c IOP15  (*IVAMC) If not zero, a return will be made to the user after
c   the initialization.  Associated with option 15.  This might be used
c   to overlay *iva, some of the user's code, and perhaps *ivaop.
c IOP16  (*IVAMC) Points to location in KORD() where information for
c   specifying the error tolerance is specified.  See option 16.
c IOP17  (*IVAMC) Used in initialization for option 17, afterwards this
c   cell is used by KEXIT which is equivalenced to IOP17.
c IOP18  (*IVAMC) Points to location in KORD() where information for
c   specifying a grouping of equations for derivative evaluation is
c   stored.  See option 18.
c IOP19  (*IVAMC) Points to location in KORD() where information for
c   specifying a grouping of equations for integration order control
c   is stored.  See option 19.
c IOP20  (*IVAMC) Used for option 20, gives first location in F where
c   estimated errors are to be stored.  Expected to be useful in a
c   program for solving boundary value problems using multiple shooting.
c IOP21  (*IVAMC) Was used for stiff equations option (never completely
c   coded).  The optional code still uses this (don't activate it!).
c   Now used to flag the location if F where the user has stored the
c    tolerance to use in finding G-Stops.
c IOP21S (*IVAMC) Was used for stiff equations see above.
c IOP22  (*IVAMC) Set aside for possible option for stiff equations.
c IOP3   (*IVAMC) Value set by option 3.
c   =  0 Interpolate to final point. (The default)
c   =  1 Integrate to final point.
c   = -1 Extrapolate to final point.
c IOP4   (*IVAMC) Value set by option 4.  The output routine is called
c   with KORD(1) = 4, every IOP4 steps.  (Default value for IOP4 is a
c   very large number.
c IOP5   (*IVAMC) Value provided by option 5, used to specify extra
c   output points.
c IOP6   (*IVAMC) Value provided by option 6.  If nonzero, the output
c   routine is called at the end of every step.  If > 0, there are
c   IOP6 interpolating G-Stops.
c IOP7   (*IVAMC) Value provided by option 7.  If > 0, there are K7
c   extrapolating G-Stops.
c IOP8   (*IVAMC) Value provided by option 8.  If nonzero, the output
c   routine is called with KORD(1)=8 whenever the step size is changed.
c IOP9   (*IVAMC) Value provided by option 9.  Used to specify that the
c   user wishes to save the solution.
c IOPIVA (*IVA) Used to save length of IOPT vector for error messages.
c IOPST  (*IVASC) Intended for possible use in stiff equations.
c IOPT   (formal *IVA and IVAOP) Used to specify options.
c IOPTC  (*IVAOP) In *IVAOP equivalenced so that IOPTC(3) is equivalent
c   to IOP3.
c IOPTS  (*IVAOP) Array containing the current default values to be
c   stored into IOPTC.
c IORD   (*IVACR) Index of first equation in the current group of
c   equations grouped for selecting integration order.
c IOUTKO (*IVADC) Used in *IVADI to point to KORD to keep track of
c   equation grouping for diagnostic output.
c ISVCOM (*IVADE) Used to save info. in the common block *IVASC.
c ITERS  (*IVADE) Counts iterations in starting delay differential
c   equations.  Max. value for this is arbitrarily 100.
c ITOLEP (*IVAMC) Used for temporary storage, and for the index of a
c   tolerance relative to the start of tolerances.
c IVC1   (*IVADB) Array used for output of variables IOPST to NUMDT in
c   common block *IVASC.
c IVC2   (*IVADB) Array used for output of variables ICF to NY in
c   common block *IVAMC.
c IWB    (*IVADE) Current base index for saving F values in WDE when
c   starting delay differential equations.
c IY     (*IVAMC) Used for the current index to the Y() array.  (Local
c   variable in *IVAIN used in computing IYI.)  Equivalenced to
c   IZFLAG in *IVAG.
c IYI    (*IVAIN) Y(IYI) is currently being computed.
c IYN    (*IVAIN) Y(IYN) is base Y() corresponding to Y(IYI).
c IYNI   (*IVAIN) Used as base index for computing IYN as IY is for INI.
c IYO    (*IVADE) Points to first base value of Y for current
c   interpolation when getting values for a delay differential equation.
c IZFLAG (*IVAG)  Equivalenced to IY.  Set to 0 initially, and later
c   set to the value returned by *ZERO.
c    = 0  Value set on entry at start of search.
c    = 1  Compute next g again.
c    = 2  Normal terminiation.
c    = 3  Normal termination -- error criterion not satisfied.
c    = 4  Apparent discontinuity -- no zero found.
c    = 5  Couldn't find a sign change.
c    = 6  *ZERO was called with a bad value in IZFLAG.
c J      For temporary storage. (In *IVA,A,BU,CR,DA,DB,DE,HC,IN,OP,PR)
c J1     (*IVAA & DA) Used for temporary storage.
c J2     (*IVAA) Used for temporary storage.
c JL     (*IVA) Used for checking storage.
c JLGREP (*IVACR) Contents of first location of KORD (called LGROUP in
c   *IVACR) for the current error tolerance rule.
c JLGROR (*IVACR) Contents of first location of KORD for the current
c   integration order control.
c JLIM   (*IVA) Used for checking second item in KORD list for options
c   16 and 19.
c K      For temporary storage.  (In *IVA,A,BU,CR,DA,DB,DE,HC,IN,OP,PR)
c KDIM   Parameter giving the largest number of differences supported.
c        Used in all the routines.
c KEMAX  (*IVAMC) Index associated with equation giving the largest
c   value for (estimated error) / (requested error).
c KEXIT  (*IVAMC) Equivalenced to IOP17 which is not used after
c   initialization.  Defines actions when KORD2I = -7.  (Referenced in
c   (*IVAA,DA,G).)
c   =  1  Take the step over with reduced H.
c   =  2  Take the step over.
c   =  3  Do the end of step call to OUTPUT.
c   =  4  Reset TMARK, then do same as for KEXIT = 2.
c   =  5  Reset TMARK, then do same as for KEXIT = 3.
c   =  6  Give the fatal error diagnostic.
c KFERR  (*IVA)  Temporary storage in checking for option 16.
c KGO    (*IVA)  Used to tell from whence a check is being done or an
c   error message is being written.
c   = 1 Checking an equation group for variational equations.
c   = 2 Checking an equation group for diagnostic print.
c   = 3 Checking an equation group for integration order control.
c   = 4 Checking an equation group for error control.
c   = 5 Checking an equation group for specifying ODE orders.
c   = 6 Found a problem with output type for printing.
c   = 7 Found a problem with an output group for printing.
c   = 8 Found a problem with input NEQ.
c   = 9 Order specified for the ODE's in the system is out of range.
c   =10 Option 16 was not used (an error).
c   =11 Error tolerance of 0 specified without proper flags.
c KIS    (*IVAMC) Used to check if it is time to dump the solution.
c   The check involves incrementing KIS at the end of the step, and
c   dumping the solution if KIS is 0.
c   = -1  Set in *ivacr when it is time to dump solution
c   =  0  When starting
c   =  2  After being dumped.
c   This is set to 1000 just after a user specified discontinuity, and
c   counted up from that point.
c KMARK  (*IVAMC) Identifies the type of output associated with the next
c   output point specified by TSPECS.
c KONV   (*IVADE) Counts iterations.  Test for convergence if KONV > 1.
c KORD   (formal in *IVA,A,BU,CR,DA,DB,DE,G,IN,PR) KORD(1) is used to
c   return flags for the user to test, and KORD(2) tells what routine
c   the flag is associated with.  See KORD1I and KORD2I below and the
c   write up for the program.  KORD(3) is used for communicating extra
c   information to the user in some cases.  KORD(4) to KORD(NTE+3) are
c   used for integration order for the equations, and the rest of KORD()
c   is available for user options.
c KORD1I (*IVAMC) Helps in defining the state of the integrator.
c   Frequently has the same value as KORD(1).  Meaning depends on the
c   value of KORD(2), or the value about to be assigned to KORD(2).
c   <  0  Happens when preparing to give output with extrapolation.
c   =  0  Happens when checking F at points for noise test.
c   =  1  (KORD(2)=-1)  End of integration has been reached.
c   =  1  (KORD(2)= 0)  Computing first predicted derivative.
c   =  1  (KORD(2)= 1)  Output for initial point.
c   =  2  (KORD(2)=-1)  Giving diagnostic for noise limiting precision.
c   =  2  (KORD(2)= 0)  Computing corrected derivative.
c   =  2  (KORD(2)= 1)  Output for TSPECS(3).
c   =  3  (KORD(2)=-1)  Diagnostic for step size reduction too fast.
c   =  3  (KORD(2)= 0)  Computing variational derivative.
c   =  3  (KORD(2)= 1)  Output for TSPECS(4).
c   =  4  (KORD(2)=-1)  Error, discontinuity.
c   =  4  (KORD(2)= 1)  Output for certain number of steps.
c   =  5  (KORD(2)= 0)  Get initial derivatives for stiff equations.
c   =  5  (KORD(2)= 1)  Extra output from TSPECS.
c   =  6  (KORD(2)= 1)  End of step output.
c   =  7  (KORD(2)= 0)  Evaluate G before extrapolated output point.
c   =  7  (KORD(2)= 1)  Evaluate G before extrapolated output point.
c                       (Also used when checking for other G's after
c                        finding one.)
c   =  8  (KORD(2)= 1)  Tell user step size has changed.
c   =  9  (KORD(2)= 1)  Request for user to save solution.
c   = 11  (KORD(2)=-1)  Error, step size too small at end of start.
c   = 12  (KORD(2)=-1)  Error, step size is too small.
c   = 13  (KORD(2)=-1)  Error, output points specified badly.
c   = 21  (KORD(2)=-1)  H too small to give reasonable change when added
c                       to T.
c   = 22  (KORD(2)=-1)  Error, bad tolerance.
c   = 23  (KORD(2)=-1)  Set after message for a fatal error.
c   = 24  Set on error message in *iva, along with KORD2I = -4.
c   Also used as an index into MLOC in *IVAA when an error is being
c   processsed, see MLOC below.
c KORD2I (*IVAMC) Helps in defining the state of the integrator.
c   Frequently has the same value as KORD(2).
c   = -3  Set in *ivag, to get a derivative evaluation.
c   = -2  Set in *ivag, to get another entry to OUTPUT.
c   = -1  Return to calling program, done, interrupt, or got an error.
c   =  1  Calling OUTPUT or returning to user for OUTPUT type action.
c   =  0  Calling DERIVS or returning to user for DERIVS type action.
c   = -4  Error message in *iva and in *ivaop, along with KORD1I = 24.
c   = -5  Starting
c   = -6  Starting, getting the initial derivative value or derivatives
c         for the noise test.
c   = -7  Done some extrapolation, KEXIT defines the action to take.
c         Set in *ivag to activate KEXIT action in *iva.
c   = -8  Set when user has requested adjustment of the difference
c         tables for a discontinutiy.
c KORDI  (*IVASC) Order of differential equation being integrated.  If
c   all orders are the same, this set once at the beginning.
c KOUNT   (*IVADE) Count of number of points back from the initial point
c   when solving a delay differential equation.
c KOUNTM  (*IVADE) Largest value currrently allowed for KOUNT.
c KOUNTX  (*IVADE) Largest value allowed for KOUNTM.
c KOUTKO  Used in DIVACR to track where output is wanted.
c KPRED  (*IVAMC) Value assigned to KORD1I when getting a predicted
c   derivative.  (1 used now, 5 planned for use with stiff equations.)
c KQD    (*IVACR) = max(2, integration order)
c KQDCON (*IVAMC) Number of coefficients computed with constant step
c   size for stiff equations.
c KQICON (*IVAMC) Number of coefficients computed with constant step
c   size for nonstiff equations.
c KQL    (*IVACR) Integration order at start of (*IVACR)
c KQLORD (*IVACR) Saved value of KQL when equations are grouped for
c   controlling the integration order.
c KQMAXD (*IVASC) Maximum integration order used for stiff equations.
c KQMAXI (*IVASC) Maximum integration order used for nonstiff equations.
c KQMAXS (*IVAMC) Maximum integration order for equations that have
c   some limit on the error that can be committed.
c KQMXDS (*IVAMC) Used to save KQMAXD in case step is repeated and the
c   solution must be dumped.
c KQMXI  (*IVAIN) Maximum integration order used for integration or
c   interpolation, = KQMAXI+INTERP-1.
c KQMXS  (*IVAIN) Maximum step number, = max(KQMXI, KQMAXD).
c KQMXIL (*IVAMC) Value of KQMAXI the last time integration coefficients
c   were computed.
c KQMXIP (*IVAMC) = KQMAXI + MAXINT, for computing integration coeffs.
c KQMXIS (*IVAMC) Used to save KQMAXI in case step is repeated and the
c   solution must be dumped.
c KQN    (*IVACR) Value of integration order at end of *IVACR.
c KQQ    Used for the integration order for current equation.  (Values
c   < 0 are intended for stiff equations.)  (In *IVA,BU,DA,IN,PR)
c KSC    (*IVAMC) Number of steps that have been taken with a constant
c   step size.
c KSOUT  (*IVAMC) When KSTEP reaches this value, the output routine is
c   called with KORD(1) = 4.  The default value is a very large number.
c KSSTRT (*IVAMC) Set when ending one derivative per step to KSTEP + 2.
c   Checked later in *IVAHC to decide whether to set the step changing
c   factors to their nominal values.
c KSTEP  (*IVAMC) Number of steps taken since the start of integration.
c L      Used for temporary storage.  In *IVAIN, L is the initial value
c   of LDT, except L=1 if LDT=-1, and MAXINT .ge. 0.  (Used in *IVAA,BU
c   CR,DA,DB,IN,PR.)
c LAHAG  (*IVADB) Used to get proper offset into an diagnostic message.
c LAIAG  (*IVADB) Used to get proper offset into an diagnostic message.
c LDIS   (*IVAA) Count of steps since user flagged a discontinuity.
c LDT    (*IVASC) Used to keep track of state of difference table.
c   = -5  Used only on first step to indicate that an extra iteration
c         is desired to get a firm estimate on the error.
c   = -4  Set on initialization before there is any difference table.
c   = -3  Set just after predicting, interpolation is not allowed when
c         this value is set.
c   = -2  Set when difference table is to be updated to the end of the
c         current step, but no interpolation is to be done.  (For
c         dumping the solution.)
c   =  0  Calculations for current step are complete, it is not o.k. to
c         update the difference table.
c   =  1  Difference table has been updated to the end of the current
c         step, e.g. by doing an interpolation.
c   =  2  Set when doing a special interpolation during computation of
c         derivatives.  (For delay equations.)
c LEX    (*IVAMC) Indicates how to get values at next output point:
c   = -1  Extrapolate
c   =  0  Interpolate (The usual case.)
c   =  1  Integrate to the output point, integration is not continued.
c LGO    (*IVAIN) Used as an an assigned go to.  Result is to add in
c   extra correction term when LDT has been set to 2.
c LGROUP (formal) This is a part of KORD passed into *IVACR.  The first
c   location is the start of the information on the grouping of
c   equations for error control.
c LINC   (*IVAMC) Used to indicate state of step size selection.
c   = -10 After computed derivatives at base time, after computing other
c         extra derivatives for the noise test.
c   = -9  After computed second extra derivative for noise test.
c   = -8  After computed first extra derivative for noise test.
c   = -7  Dumping the solution and then doing a user initiated restart,
c         or getting ready to compute extra derivatives for the noise
c         test.
c   = -6  Dumping the solution before a restart.
c   = -5  Set on the first step, and also set when dumping the solution
c         after a discontinuity.
c   = -4  Repeat step with no change in the step size.
c   = -3  Set when the error tolerance is set improperly.
c   = -2  User has complete control of selecting the step size.
c   = -1  Step is being repeated.
c   =  0  Step size is not to be increased on this step.
c   = k>0 Step size can be increased by HINCC**k.
c LINCD  (*IVAMC) Value of smallest k for which HINCC**k .ge. 2.
c   (=-2 if user is specifying all step size changes.)
c LINCQ  (*IVAMC) Value of smallest k for which HINCC**k .ge. 4.
c LIOPT  (*IVAOP) Value of the last index in IOPT on the last call.
c   Used so *IVA can print IOPT in error messages.
c LL     (*IVACR) Temporary variable used when equations are grouped
c   for integration order control.
c LNOTM1 (*IVAIN) Logical variable = L .ne. -1.  If LNOTM1 is true,
c   storage in Y() is different in some way lost to antiquity.  Such
c   a case can only arise in the case of stiff equations.
c LOCF1  (*IVADB) Gives packed data needed for output of tables by the
c   message processor MESS.  See comments there under METABL for defs.
c LOCF2  (*IVADB) As for LOCF1 above.
c LOCM   (*IVAA) Parameter = 32*256, used to unpack integers stored
c   in MLOC for use in error message processing.
c LPRINT (formal, *IVADB) Defines how much printing is to be done in
c   *IVADB.  Let |LPRINT| = 10*N1 + N2     (N1,N2 digits)
c    N1=1   Do not print any variables external to the integrator.
c    N1=2   Print  tspecs, current y, past y, current f, all pertinent
c           contents of KORD, and TOL.
c    N1=3   Above + difference tables up to highest difference used.
c    N1=4   Same as N1=1 + all in storage allocated for differences.
c    N2=1   Do not print any variables internal to the integrator.
c    N2=2   Print all scalar variables in interpolation common block.
c    N2=3   Above + all scalar variables in main integ. common block.
c    N2=4   Same as N1=3 + all used in arrays XI,BETA,ALPHA, first
c           column of G, GS,RBQ,SIGMA
c    N2=5   Same as N1=4 + all used in arrays G,D,DS,V
c LSC    (*IVAMC) Indicates if starting or if noise may be present.
c   =k<0 -k steps have been taken for which noise appears to be limiting
c        the precision.
c   = 0  Usual case
c   = 1  Doing 1 derivative per step after initial part of start.
c   = 2  Used as flag that it is time to set LSC=0.
c   = 3  Third step, hold the order constant.
c   = 4  Second step, increase orders from 2 to 3.
c   = 5  First step, third time through the first step (if required).
c   = 6  First step, second time through.
c   = 7  First step, first time through.
c   = 8  Set on initialization.
c LTXT?? Names of this form are used in setting up data statements for
c   error messages.  These names are generated automatically by PMESS,
c   the program that makes up these messages.
c LX     (*IVAA) Used for temporary storage in computing TMARKA().
c        ( formal *IVADE)  An integer array containing extra
c   information, as follows.
c  LX(1) Points to a location in Y beyond those already in use.  Values
c        of Y requested are computed at TSPECS(1) - Y(LX(1)) and stored
c        starting at Y(LX(1)+1).  If this index is 0, no more extra Y
c        values are to be computed.
c  LX(2) Index of the first equation for which the Y's above are to be
c        computed.  Y(LX(1)+1) will correspond to this first equation
c        index.
c  LX(3) Index of the last equation for which the Y's above are to be
c        computed.  Thus the Y's stored starting at Y(LX(1)+1) will
c        require no more space than half the space ordinarily required
c        for the array Y(), and may require significantly less.
c  LX(4) Maximum number of times to integrate F to get Y.  This should
c        be > 0, and less than or equal to the order of the highest
c        order differential equation.  (= 0 is allowed, but probably
c        not what you want.  It would give a value only for F.)  Space
c        must be set aside for all integrals of F, even if not all are
c        requested.  For a first order system, all Y's are just the
c        first integrals of the corresponding F's.  For higher order
c        equations, the first Y associated with a given F is the d-th
c        integral of the corresponding F, where d is the order of the
c        equation, and the last Y corresponding to the F is the first
c        integral of that F.
c  LX(5) As for LX(4), but gives the index for the fewest number of
c        times to integrate F.  Ordinarily this should be > 0.  If 0 is
c        requested, an estimate for the value of F at the delay point is
c        computed.  This should not be 0 more than once, for equations
c        covering the same index, since later such requests would write
c        over the earlier results.
c  LX(5i+k) , k = 1, 2, ... 5.  Treated as for the cases above.  If i
c        different cases of delayed Y's are to be computed, then
c        LX(5i+1) must be 0.
c LX2    (*IVADE) Value of LX(5i+2), when working on the i-th delay.
c MACT   Used in the programs which call the error message program.
c   This array difines the actions to be taken by that program.  (In
c   (*IVA,A,DA,DE,G,IN,OP)
c MACT0  (*IVADB) Used to call the message program, see MACT.
c MACT?  As for MACT, in (*IVA,CR,DB)
c MACTFV (*IVADB) As for MACT0.
c MAXDIF (*IVASC) Maximum differentiations required for stiff equations.
c MAXINT (*IVASC) Maximum integrations required.  (= max. order of
c   differential equations if equations are not stiff.)
c MAXKQ  (*IVA, BU)e
c MAXKQD (*IVAMC) Largest integration order allowed for stiff equations.
c MAXKQI (*IVAMC) Largest integ. order allowed for nonstiff equations.
c ME???? Parameters defining constants used for interaction with the
c   error message program MESS.  See comments there for definitions.
c   (In *IVA,A,DA,DE,G,IN,OP)
c METHOD (*IVAMC) Defines kind of methods being used.
c   = -1  Only stiff equations are being integrated.
c   =  0  Only nonstiff equations are being integrated.
c   =  1  Both kinds of methods are required.
c MLOC   (*IVA,A,DE) Contains locations in MTEXT for error messages.  In
c   *IVAA this data is packed using MLOC??, see below.
c MLOC?? (*IVAA) Parameters constructed to aid in making up packed data
c   for processing error messages.  Low two digits give the value of
c   KORD1I to use for the error index and later processing, the next two
c   give the error severity level, and the rest point to text used for
c   the message.
c MODF2  (*IVADB) Used in constructing the same kind of packed data as
c   described for LOCF1 above.
c MULTJ  Local to DIVAOP for calls not using F.
c MTEXT  (*IVA,A,CR,IN,OP) Text for error messages.
c MTXT?? (*IVA,A,CR,DA,DB,DE,G,IN,OP) Equivalenced into MTEXT.
c N      Used for temporary storage.  (In *IVAHC,IN,PR)
c NDTF   (*IVASC) Location in F() where difference table starts.
c NE     (*IVAMC) Number of equations in the first group.  (=NTE if
c   option 18 is not used.)
c NEDDIG (*IVADB) Parameter = -MEDDIG.
c NEPTOL (*IVAMC) Used for temporary storage and to save the value of
c   ITOLEP for error messages.
c NEQ    (formal) Total number of equations being integrated.
c NG     (*IVAMC) Used in *ivag for the number of g's in the current
c   context.
c NGSTOP (*IVAG) Dimension 2 array equivalenced to IOP6, and IOP7.  To
c   get the number of interpolating and extrapolating G-Stops.
c NGTOT  (*IVAMC) NGTOT(1) gives the number of interpolating G-Stops,
c   and NGTOT(2) gives the number of extrapolating G-Stops.
c NKDKO  (*IVASC) If this is nonzero (option 17), it gives the location
c   in KORD() where a vector defining the order of each equation is
c   specified.
c NLX    (*IVADE) Temporary index used to keep track of interpolations
c   being done to get Y() values for a delay differential equation.
c NOISEQ (*IVAMC) max(2, order of equation for which (error estimate)/
c   (error requested) is a maximum).
c NOUTKO (*IVAMC) If nonzero, gives the index in KORD where information
c   on what equations are to be included in the diagnostic output is
c   given.   See option 10.
c NSTOP  (formal) In *IVAG.  Index of the G-stop, see IFLAG.
c NTE    (*IVASC) Total number of equations being integrated = NEQ.
c NTEXT  (formao *IVADB) Character variable containing heading text.
c NTOLF  (*IVAMC) First location in F() where tolerance specifying
c   accuracy desired is stored.
c NUMDT  (*IVASC) Maximum allowed number of differences available for
c   doing an integration.
c NXTCHK (*IVA) Equivalenced to INTCHK(1), which gives the next
c   available location in INTCHK for storing data on storage allocation.
c NY     (*IVAMC) Total order of the system.
c NYNY   (*IVASC) Location in Y() where the base value for Y() is saved.
c NYNYO  (*IVADE) Equivalenced to the saved value from common of NYNY.
c OUTPUT (formal) Name of subroutine to be called for the output of
c   data or for computing G-Stops.  Not used if option 14 is set.
c OVD10  (*IVAEV) (largest floating point number) / 10.
c OVTM75 (*IVAEV) (largest floating point number) ** (-.75)
c RBQ    (*IVAMC) Array containing data for the preliminary noise test.
c RD     (formal *IVACO) Array use to contain floating point data from
c   common.  Values are returned in locations 1 to 3 as follows.
c   1    EMAX =  Max. ratio of estimated error to requested error
c   2            Reserved for future use
c   3            Reserved for future use
c REF    (*IVACR) Array of length 3 used for translating error tolerance
c   type into the factor used for exponential averaging for that type.
c RND    (*IVACR) Usually the current estimated error.  Used in deciding
c   if noise is limiting precision.
c RNOISE (*IVACR) Value used in comparison with RBQ() for preliminary
c   noise test.
c ROBND  (*IVAMC) Used to influence the selection of integration order.
c   The larger ROBND, the harder it is to increase the order and the
c   easier it is to decrease it.
c RVC2   (*IVADB) Array used for output of variables DNOISE to SNOISE in
c   common block *IVAMC.  These are variables that don't require a great
c   deal of precision.
c S      (*IVACR) Estimate of (step size) * eigenvalue of Jacobian.
c SIGMA  (*IVAMC) The k-th entry of this array contains a factor that
c   gives the amount the k-th difference is expected to increase if the
c   step size in increased.  These numbers get bigger it there is a past
c   history of increasing the step size.
c SIGMAS (*IVAA) Saved value of SIGMA(k) from the last step, where k =
c   integration order for equation with index KEMAX.
c SNOISE (*IVAMC) Value used in comparison with RBQ() on equation with
c   largest value for (error estimate) / (error requested).
c T      (formal) in *IVAIN. T(1) contains the point to be interpolated
c   to, and T(2) is used in a check that |HI| .le. |T(2)|.  When used by
c   other routines in this package, TSPECS is passed in for T.
c TB      (*IVADE) Base time for current interpolation.
c TC      (*IVADE) Original value of TN when getting past Y's for a
c   delay differential equation.
c TEMP   Used for temporary storage, in *IVAHC,PR
c TEMPA  (*IVACR) Array equivalenced to (TPS1,TPS2,TPS3,TPS4).
c TEMPAO (*IVACR) Array used to accumulate values in TEMPA.
c TG     (*IVAMC) TG(1) gives the last value of TSPECS(1) for which an
c   interpolatory G-Stop has been computed.  TG(2) is defined similarly
c   for extrapolatory G-Stops.
c TGSTOP (*IVAMC) TGSTOP(1) gives the value of TSPECS(1) where the last
c   0 for an interpolatory G-Stop was found.  TGSTOP(2) is defined
c   similarly for extrapolatory G-Stops.
c TMARK  (*IVAMC) Location of the next output point.
c TMARKA (*IVAA)  Array of length 2 equivalenced to TMARK (and TMARKX).
c TMARKX (*IVAMC) Location of the next output point to be found using
c   integration or extrapolation.  This variable must follow immediately
c   after TMARK in the common block.
c TN     (*IVASC) The value of TSPECS(1) at the conclusion of the last
c   step.
c TNEQ   (*IVADB) Array of dimension 1 equivalenced to TN so that an
c   array can be passed to *MESS.
c TOL    (formal) This is a part of F passed into *IVACR.  The first
c   location is the start of the information on the tolerances for error
c   control.
c TOLD   (*IVAG) Value of TSPECS(1) on one side of a zero.
c TOLG   (*IVAMC) Tolerance to pass to dzero when locating G-Stops.
c TOUT   (*IVAMC) Location of next output point defined by value of
c   TSPECS(3).  Such output is given with KORD(1) = 2.
c TP     (*IVA,A,DA,DE,HC) Used for temporary storage.
c TP1    (*IVAA,DA,HC,IN,PR) Used for temporary storage.
c TP2    (*IVAA,DA,HC,PR) Used for temporary storage.
c TP3    (*IVAA) Used for temporary storage.
c TPD    (*IVABU) Used for temporary storage.
c TPP    (*IVACR) Used for temporary storage.  Usually same as TPS3.
c TPS1   (*IVAA,CR) Used for temporary storage.  (In *IVACR is the
c   difference of order KQQ-2)
c TPS2   (*IVAA,CR) Used for temporary storage.  (In *IVACR is the
c   difference of order KQQ-1)
c TPS3   (*IVACR) Contains the difference of order KQQ.  This is the
c   last difference used in the corrector.
c TPS4   (*IVACR) Contains the difference of order KQQ+1.
c TPS5   (*IVACR) Temporary storage.
c TPS6   (*IVACR) Temporary storage.
c TPS7   (*IVACR) Temporary storage.
c TSAVE  (*IVAG) Value of TSPECS(1) before starting the search for a 0.
c TSPECS (formal *IVA,A,DB,DE,G)
c   TSPECS(1) is the current value of the independent variable.
c   TSPECS(2) is the current value of the step size.
c   TSPECS(3) is the increment to use between output points that give
c             output with KORD(1) = 2.
c   TSPECS(4) is the "final" output point.
c V      (*IVAMC) Array used in computing integration coefficients.
c XI     (*IVASC) XI(K) = TSPECS(1) - value of TSPECS(1) K steps
c   previous.
c W      (*IVAHC) Array used in computing integration coefficients.
c WDE    (formal, *IVADE)  Array used for working storage.  This storage
c   is used to save derivative values when iterating to get started.  To
c   be safe one should allow as much space as is allowed for differences
c   in F.  In most cases the start will not require this much space
c   however.  This array is also intended for the support of saving long
c   past histories.
c Y      (formal, *IVA,A,CR,DA,DB,DE,G,IN,PR) Array containing the
c   independent variable and all derivatives up to order one less than
c   the order of the differential equation.  Also use to save these
c   values at the beginning of the current step, the base values.
c YN     (formal, in *IVAPR)  Base values of y, these follow the
c   current values of the dependent variable, y, in Y().
c
c
c++S Default KDIM = 16
c++  Default KDIM = 20
c++  Default MAXORD = 2, MAXSTF = 1
c++  Default INTEGO, VAREQ, OUTPUT, DUMP, GSTOP, EXTRAP
c++  Default STIFF=.F., ARGM=.F., ERRSTO=.F.
c
      integer NEQ, IDIMT, IDIMY, IDIMF, IDIMK
      integer KORD(*), IOPT(*)
c--D Next line special: P=>D, X=>Q
      double precision TSPECS(*), Y(*)
      double precision F(*)
      external DIVAF, DIVAO
c
c *********************** Internal Variables ***************************
c
c Comments for variables used in this package can be found in the file
c   IVACOM.
c
c *********************** Type Declarations ****************************
c
      integer KDIM, MAXORD, MAXSTF
c++ Substitute for KDIM, MAXORD, MAXSTF below
      parameter (KDIM = 20, MAXORD = 2, MAXSTF = 1)
c--D Next line special: P=>D, X=>Q
      double precision TN
      double precision XI(KDIM)
c
c--D Next line special: P=>D, X=>Q
      double precision TG(2), TGSTOP(2), TMARK, TMARKX, TOUT, TOLG
      double precision ALPHA(KDIM), BETA(KDIM+1)
      double precision  D(MAXSTF+MAXORD,MAXORD), G(KDIM,MAXORD)
      double precision V(KDIM+MAXORD)
      double precision HC, HDEC, HINC, HINCC, HMAX, HMAXP9, HMIN
      double precision FDAT(11)
c
      double precision DS(MAXSTF+MAXORD, MAXORD), GS(KDIM)
      double precision SIGMA(KDIM), RBQ(KDIM), DNOISE
      double precision EAVE, EIMAX, EIMIN, EMAX, EREP, ROBND, SNOISE
c
c.    SPECIFICATION OF ENVIRONMENTAL CONSTANTS.
      double precision EEPS10, EEPS16, EROV10, EEPS2
      double precision EEPT75, EOVEP2, OVTM75, OVD10
      common / DIVAEV / EEPS2, EEPT75, EOVEP2, OVTM75, OVD10, EEPS10,
     1   EEPS16, EROV10
      save / DIVAEV /
      integer IOPST, KORDI, KQMAXD, KQMAXI, LDT, MAXDIF, MAXINT, NKDKO,
     1   NTE, NYNY, NDTF, NUMDT
      common / DIVASC / TN, XI, IOPST, KORDI, KQMAXD, KQMAXI, LDT,
     1   MAXDIF, MAXINT, NKDKO, NTE, NYNY, NDTF, NUMDT
c
      integer ICF,ICS,IGFLG,IGTYPE(2),IGSTOP(2),ILGREP,INGS,IOP3,IOP4,
     1   IOP5,IOP6,IOP7,IOP8,IOP9,IOP10,IOP11,IOP12,IOP13,IOP14,IOP15,
     2   IOP16,IOP17,IOP18,IOP19,IOP20,IOP21,IOP22,IOP21S,ITOLEP,IY,
     3   KEMAX,KIS,KMARK,KORD1I,KORD2I,KPRED,KQDCON,KQICON,KQMAXS,
     4   KQMXDS,KQMXIL,KQMXIP,KQMXIS,KSC,KSOUT,KSSTRT,KSTEP,LEX,LINC,
     5   LINCD,LINCQ,LSC,MAXKQD,MAXKQI,METHOD,NE,NEPTOL,NG,NGTOT,
     6   NOISEQ,NOUTKO,NTOLF,NY,IDAT(6)
      common /DIVAMC/ TG,TGSTOP,TMARK,TMARKX,TOUT,TOLG,HC,HDEC,HINC,
     1   HINCC,HMAX,HMAXP9,HMIN,ALPHA,BETA,D,G,V,DS,GS,SIGMA,RBQ,DNOISE,
     2   EAVE,EIMAX,EIMIN,EMAX,EREP,ROBND,SNOISE,FDAT,ICF,ICS,IGFLG,
     3   IGTYPE,IGSTOP,ILGREP,INGS,IOP3,IOP4,IOP5,IOP6,IOP7,IOP8,IOP9,
     4   IOP10,IOP11,IOP12,IOP13,IOP14,IOP15,IOP16,IOP17,IOP18,IOP19,
     5   IOP20,IOP21,IOP22,IOP21S,ITOLEP,IY,KEMAX,KIS,KMARK,KORD1I,
     6   KORD2I,KPRED,KQDCON,KQICON,KQMAXS,KQMXDS,KQMXIL,KQMXIP,KQMXIS,
     7   KSC,KSOUT,KSSTRT,KSTEP,LEX,LINC,LINCD,LINCQ,LSC,MAXKQD,MAXKQI,
     8   METHOD,NE,NEPTOL,NG,NGTOT,NOISEQ,NOUTKO,NTOLF,NY,IDAT
      save / DIVAMC / , / DIVASC /
c
      integer KGO, INTCHK(0:30), NXTCHK
      integer IHI, JL, J, ILOW, K, KQQ, JLIM, KFERR
c
      equivalence (INTCHK(1), NXTCHK)
      double precision CM1
      parameter (CM1 = (-1.D0))
      integer IOPIVA(2)
      save IOPIVA
c
c                      Declarations for error message processing.
c
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
c
      integer MACT(16), MLOC(12), MACT1(4)
c
c ********* Error message text ***************
c[Last 2 letters of Param. name]  [Text generating message.]
cAA DIVA$B
cAB The interval [1, 10**6], bounds the allowed values for NTE=$I.$E
cAC For option $I, the interval [$I, $I], bounds the allowed $C
c   values for the integration order which is set to $I.$E
cAD Option 16 must be used for error control.$E
cAE F($I) = $F, but it must be -1.0 when skipping the error check.$E
cAF For option $I, the interval [$I, $I] bounds the allowed $C
c   values for KORD($I)=$I, which is used to specify an $B
cAG output type for printing.$E
cAH output group for printing.$E
cAI equation group for variational equations.$E
cAJ order for a differential equation.$E
cAK equation group for diagnostic print.$E
cAL equation group for integration order control.$E
cAM equation group for error control.$E
cAN Option 5 argument must be .le. 0 or .gt. 4.$E
c   $
cAO KORD values for this option starting at KORD($M) are:$E
      integer LTXTAA,LTXTAB,LTXTAC,LTXTAD,LTXTAE,LTXTAF,LTXTAG,LTXTAH,
     * LTXTAI,LTXTAJ,LTXTAK,LTXTAL,LTXTAM,LTXTAN,LTXTAO
      parameter (LTXTAA=  1,LTXTAB=  7,LTXTAC= 71,LTXTAD=183,LTXTAE=226,
     * LTXTAF=290,LTXTAG=400,LTXTAH=427,LTXTAI=455,LTXTAJ=498,
     * LTXTAK=535,LTXTAL=573,LTXTAM=620,LTXTAN=655,LTXTAO=  1)
      character MTXTAA(3) * (233)
      character MTXTAB(1) * (55)
      data MTXTAA/'DIVA$BThe interval [1, 10**6], bounds the allowed val
     *ues for NTE=$I.$EFor option $I, the interval [$I, $I], bounds the$
     * allowed values for the integration order which is set to $I.$EOpt
     *ion 16 must be used for error control.$EF($I) = ','$F, but it must
     * be -1.0 when skipping the error check.$EFor option $I, the interv
     *al [$I, $I] bounds the allowed values for KORD($I)=$I, which is us
     *ed to specify an $Boutput type for printing.$Eoutput group for pri
     *nting.$Eequation gro','up for variational equations.$Eorder for a$
     * differential equation.$Eequation group for diagnostic print.$Eequ
     *ation group for integration order control.$Eequation group for err
     *or control.$EOption 5 argument must be .le. 0 or .gt. 4.$E'/
      data MTXTAB/'KORD values for this option starting at KORD($M) are:
     *$E'/

c End of automatically generated error message code.
c
c        for KGO =     1      2      3      4      5      6      7
      data MLOC / LTXTAI,LTXTAK,LTXTAL,LTXTAM,LTXTAJ,LTXTAG,LTXTAH,
     2            LTXTAB,LTXTAC,LTXTAD,LTXTAE, LTXTAN /
c           KGO        8      9     10     11      12
c
c                      1  2  3 4       5  6       7       8       9 10
      data MACT / MEEMES,38,24,0, MENTXT, 0, METEXT, MECONT, MEMDA1,0,
     1  METEXT, MEIDAT,0, MEIVEC,0, MECONT /
c           11      12 13     14 15     16
      data MACT1 / METEXT, MEIVEC, 0, MERET /
      data IOPIVA(1) / 1111 /
      data TEXT1 / 'IOPT()= $B' /
c
c ************** START OF EXECUTABLE CODE ******************
c
c     **** TEST IF CONTINUING AN INTEGRATION
      if (KORD(1) .ne. 0) go to 330
c     **** INITIALIZE VARIOUS SCALARS
      KSTEP = 0
      KQMXIS = 0
      KORD2I = -5
      KORD(2) = -1
      NTE = NEQ
      NE = NTE
      TOLG = 0.D0
c     **** SET UP OPTIONS
      if (IOPT(1) .ne. 0) call DIVAOP(IOPT, F)
      call DIVAOP(IOPIVA, F)
      if (IOPT(1) .eq. 0) IOPIVA(2) = 1
c
      if ((NE .le. 0) .or. (NE .gt. 1000000)) then
         IDAT(1) = NE
         KGO = 8
         go to 650
      end if
c                         Set up diagnostic print on storage allocation.
      INTCHK(0) = 245
      if (IOP10 .ne. 0) INTCHK(0) = 247
c
c     **** CHECK TSPECS STORAGE ALLOCATION
      INTCHK(2) = IDIMT
      INTCHK(3) = 4
      NXTCHK = 4
      if (IOP5 .ne. 0) then
         INTCHK(4) = 5
         INTCHK(5) = 5
         if (IOP5 .gt. 0) then
            INTCHK(6) = IOP5 - 4
            if (IOP5 .lt. 5) then
               KGO = 12
               go to 600
            end if
         else
            IHI = -IOP5
            JL = 4
            do 15 IHI = IHI, IDIMK-3, 3
               J = abs(KORD(IHI))
               if (J .eq. 0) go to 20
               if (abs(KORD(IHI + 2)) .gt. 1) then
                  IDAT(2) = -1
                  IDAT(3) = 1
                  KGO = 6
                  go to 600
               end if
               if ((J .le. JL) .or. (J .gt. KORD(IHI+1))) then
                  KGO = 7
                  IDAT(2) = JL + 1
                  IDAT(3) = KORD(IHI+1)
                  go to 610
               end if
               JL = KORD(IHI+1)
   15       continue
            if (KORD(IHI) .ne. 0) IHI = IHI + 3
   20       INTCHK(6) = JL - 4
         end if
         NXTCHK = 7
      end if
   25 call OPTCHK(INTCHK, IOPT, 'DIVA / TSPECS$E')
      if (NXTCHK .lt. 0) KORD2I = -4
c
c     **** CHECK KORD STORAGE ALLOCATION
      INTCHK(2) = IDIMK
      INTCHK(3) = NE + 3
      NXTCHK = 4
      if (IOP5 .lt. 0) then
         INTCHK(4) = 5
         INTCHK(5) = -IOP5
         INTCHK(6) = IHI + IOP5
         NXTCHK = 7
      end if
c
c++  Code for VAREQ is active
      if (IOP18 .ne. 0) then
         NE = abs(KORD(IOP18))
         INTCHK(NXTCHK) = 18
         ILOW = IOP18
         KGO = 1
c.       **** CHECK OPTION FOR VALID INPUT
         go to 430
      end if
c++  End
   30 continue
      if (NKDKO .ne. 0) then
c                        **** STORAGE ALLOCATED FOR ODE ORDERS
         INTCHK(NXTCHK) = 17
         INTCHK(NXTCHK+1) = NKDKO
         INTCHK(NXTCHK+2) = NTE
         NXTCHK = NXTCHK + 3
      end if
c++  Code for STIFF is inactive
c      IF (IOPST .ne. 0) then
c         INTCHK(NXTCHK) = 17
c         INTCHK(NXTCHK+1) = IOPST
c         INTCHK(NXTCHK+2) = NTE
c         NXTCHK = NXTCHK + 3
c      end if
c++  End
c
c **** SET INITIAL INTEGRATION ORDERS, TEST ODE ORDERS ****
c
      MAXINT = 0
      MAXDIF = 0
      NY = 0
      do 80 K = 1, NTE
         if (NKDKO .ne. 0) KORDI = KORD(NKDKO + K - 1)
         NY = NY + abs(KORDI)
c++  Code for STIFF is inactive
c      IF (IOPST .EQ. 0) GO TO 60
cc.    **** CHECK FOR POTENTIAL STIFF EQUATION
c      JS = abs(KORD(IOPST+K-1)) - 1
c      IF ( JS ) 52,60,54
cc.    **** EQUATION IS NOT ACTIVE
c   52 KQQ = 0
c      GO TO 56
cc.    **** EQUATION USES IMPLICIT METHOD
c   54 KQQ = -1
c      IF (JS .GT. abs(KORDI)) then
c        Set up an error message.
c      end if
c      MAXINT = max(MAXINT, abs(KORDI) - JS)
c   56 IF (KORDI .GE. 0) GO TO 70
c      KORDI = -1 - KORDI
c      JS = JS - 1
c      MAXDIF = max(MAXDIF, JS, 1)
c      GO TO 70
c++  End
c     **** EQUATION IS TO USE AN EXPLICIT METHOD
   60    KQQ = 1
         MAXINT = max(MAXINT, KORDI)
   70    if ((KORDI .gt. MAXORD) .or. (KORDI .le. 0)) then
c                    Set up error message.  KORDI is out of range.
            IDAT(1) = 17
            IDAT(2) = 1
            IDAT(3) = MAXORD
            if (NKDKO .ne. 0) then
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
c     **** SET FLAGS WHICH DEPEND ON METHOD USED
c++  Code for STIFF is inactive
c      METHOD = 1
c      IF (MAXINT .GT. 0) IF (MAXDIF) 85,90,85
c      METHOD = -1
c   85 CONTINUE
c      KPRED = 5
c      GO TO 100
c++  End
   90 METHOD = 0
      KPRED = 1
  100 continue
c
c ******* CHECK KORD FOR DIAGNOSTIC OUTPUT CONTROL *********
c
c++  Code for OUTPUT is active
      if (IOP10 .gt. 0) then
         if (NOUTKO .ne. 0) then
            INTCHK(NXTCHK) = 10
            ILOW = NOUTKO
c.    **** Check option for valid input
            KGO = 2
            go to 430
         end if
      end if
c++  End
  110 continue
c
c ********** CHECK KORD FOR INTEGRATION ORDER CONTROL ******
c
c++  Code for INTEGO is active
      if (IOP19 .ne. 0) then
c.           **** Check option for valid input
         INTCHK(NXTCHK) = 19
         ILOW = IOP19
         JLIM = -30
         KGO = 3
         go to 430
      end if
c++  End
  120 continue
c
c ********** CHECK SET UP FOR ERROR TOLERANCES *************
c
      INTCHK(NXTCHK) = 16
      ILOW = IOP16
      JLIM = -5
      KGO = 4
      if (IOP16 .ne. 0) go to 430
c.                      **** IN CURRENT CODE, IOP16=0 IS AN ERROR
      KGO = 10
      go to 650
  150 continue
c     **** CHECK KORD STORAGE ALLOCATION
      call OPTCHK(INTCHK, IOPT, 'DIVA / KORD$E')
      if (NXTCHK .lt. 0) KORD2I = -4
c
c     ******** DONE CHECKING KORD STORAGE ALLOCATION *******
c
c     **** CHECK  Y  STORAGE ALLOCATION
      INTCHK(2) = IDIMY
      INTCHK(3) = NY + NY
      NXTCHK = 4
      NYNY = NY + 1
      call OPTCHK(INTCHK, IOPT, 'DIVA / Y$E')
      if (NXTCHK .lt. 0) KORD2I = -4
c
c     **** CHECK  F  STORAGE ALLOCATION
      INTCHK(2) = IDIMF
      INTCHK(3) = NTE
      NXTCHK = 4
      if (IOP16 .ne. 0) then
c                                Error tolerance info.
         INTCHK(4) = 16
         INTCHK(5) = NTOLF
         INTCHK(6) = IHI - IOP16 + 1
         NXTCHK = 7
      end if
      if (IOP12 .gt. 0) then
         INTCHK(NXTCHK) = 12
         INTCHK(NXTCHK+1) = IOP12
         INTCHK(NXTCHK+2) = 4
         NXTCHK = NXTCHK + 3
      end if
      if (IOP21 .gt. 0) then
         INTCHK(NXTCHK) = 21
         INTCHK(NXTCHK+1) = IOP21
         INTCHK(NXTCHK+2) = 1
         NXTCHK = NXTCHK + 3
      end if
c
c++  Code for ERRSTO is inactive
c      IF (IOP20 .ne. 0) then
cc.                                Space for saving error estimates
c         INTCHK(NXTCHK) = 20
c         INTCHK(NXTCHK) = IOP20
c         INTCHK(NXTCHK) = NTE
c         NXTCHK = NXTCHK + 3
c      end if
c++  Code for STIFF is inactive
c      if (IOP21 .gt. 0) then
cc.                               Info. for stiff equations
c         INTCHK(NXTCHK) = 21
c         INTCHK(NXTCHK+1) = IOP21
c         INTCHK(NXTCHK+2) = IOP21S
c         NXTCHK = NXTCHK + 3
c      end if
c      MAXKQD = min(MAXKQI, 6)
c++  End
c                          Set aside space for the difference tables.
      INTCHK(NXTCHK) = 0
      INTCHK(NXTCHK+1) = -KDIM * NTE
      INTCHK(NXTCHK+2) = 0
      NXTCHK = NXTCHK + 3
      INTCHK(NXTCHK) = -5 * NTE
      call OPTCHK(INTCHK, IOPT, 'DIVA / F$E')
      if (NXTCHK .lt. 0) then
         KORD2I = -4
      else if (KORD2I .ne. -4) then
         do 290 K = NXTCHK+1, INTCHK(NXTCHK)
            if (INTCHK(INTCHK(K)) .eq. 0) then
               NDTF = INTCHK(INTCHK(K)+1)
               NUMDT = min(KDIM, (INTCHK(INTCHK(K)+2)-NDTF+1) / NTE)
               MAXKQI = NUMDT - 1
            else
c         Take a quick return if needed space was not specified by user.
               KORD2I = -4
            end if
  290    continue
      end if
      if (IOP9 + abs(IOP10) + IOP11 .ne. 0) then
c Insure user doesn't get in trouble with F not iniitalized.
         do 300 K = NDTF, NDTF + NTE*NUMDT - 1
            F(K) = 0.D0
  300    continue
      end if
  320 continue
      if ((KORD2I .eq. -4) .or. (IOP10 .ne. 0)) then
         MACT1(3) = IOPIVA(2)
         call MESS(MACT1, TEXT1, IOPT)
         KORD1I = 24
         KORD(1) = 24
      end if
      TMARK = TSPECS(1)
      TMARKX = TSPECS(4) + TSPECS(2)
c
c     **** DONE WITH INITIALIZATION AND CHECKING INPUTS
      if (IOP13 + IOP14 + IOP15 .ne. 0) return
  330 call DIVAA(TSPECS, Y, F, KORD, DIVAF, DIVAO)
      return
c
c ************ LOOP TO CHECK OPTION SPECIFICATIONS *********
c
  430 JL = 0
      do 560 IHI = ILOW, IDIMK
         J = KORD(IHI)
         go to (460, 480, 490, 490), KGO
c     **** CHECK ON VARIATIONAL EQUATIONS
  460    continue
c++  Code for VAREQ is active
         if (J - NTE) 470, 565, 620
  470    if (J .eq. 0) go to 560
         if (J .le. JL) go to 620
c++  End
c     **** Check on diagnostic output option
  480    continue
c++  Code for OUTPUT is active
c.    **** CHECK IF DONE
         if (J .ge. NTE) go to 565
         if (J .le. JL) go to 620
         go to 550
c++  End
  490    continue
c     **** Check integration order control (KGO=3) and
c     **** error tolerance equation grouping (KGO=4).
         if (J - NTE) 500, 565, 620
  500    if (J) 510, 530, 540
  510    if ((JL .le. 0) .and. (IHI .ne. ILOW)) go to 620
         if (J .lt. JLIM) then
c                         Output an error message.
            IDAT(2) = JLIM
            IDAT(3) = 0
            go to 630
         end if
  520    JL = -JL
         go to 560
  530    if (KGO .eq. 3) go to 520
         KFERR = NTOLF + IHI - ILOW
         if (F(KFERR) .eq. CM1) go to 510
c                         Set up error message, TOL must be -1.
            IDAT(1) = KFERR
            KGO = 11
            go to 650
  540    if (abs(JL) .ge. abs(J)) go to 620
  550    JL = J
  560    continue
  565 NXTCHK = NXTCHK + 3
      INTCHK(NXTCHK-2) = ILOW
      INTCHK(NXTCHK-1) = IHI - ILOW + 1
      go to (30, 110, 120, 150), KGO
c
c     **** AN ERROR HAS BEEN MADE
c                  Error in setting up TSPECS for extra output
  600 IHI = IHI + 2
  610 ILOW = -IOP5
      go to 630
c                  Error in KORD indices
  620 IDAT(2) = abs(JL) + 1
      IDAT(3) = NTE
c                  Set up for print of message about KORD
  630 IDAT(1) = INTCHK(NXTCHK)
  640 IDAT(4) = IHI
      IDAT(5) = KORD(IHI)
c
c ***************** Process Errors *************************************
c
  650 KORD2I = -4
      MACT(4) = LTXTAF
      if (KGO .ge. 8) MACT(4) = -1
      MACT(6) = MLOC(KGO)
c--D Next line special: P=>S, X=>D
      CALL DMESS(MACT, MTXTAA, IDAT, FDAT)
      if (KGO .lt. 8) then
         MACT(10) = ILOW
         MACT(13) = ILOW
         MACT(15) = -min(IHI+2, IDIMK)
         CALL MESS(MACT(9), MTXTAB, KORD)
         if (KGO .le. 4) go to 565
      end if
c              5   6   7    8    9   10   11  12
      go to (100, 25, 25, 320, 100, 150, 660, 25), KGO - 4
  660 KGO = 4
      go to 565
      end
c   End of DIVA

      subroutine DIVAA(TSPECS, Y, F, KORD, DIVAF, DIVAO)
c>> 1989-02-24 DIVAA  Krogh   Big error with BETA(2)=1+epsilon -- looped
c>> 1988-07-27 DIVAA  Krogh   Fixed to allow restart on a restart.
c>> 1988-03-07 DIVAA  Krogh   Initial code.
c
c  MAIN SUBROUTINE FOR VARIABLE ORDER INTEGRATION OF ORDINARY
c  DIFFERENTIAL EQUATIONS
c
      integer KORD(*)
c--D Next line special: P=>D, X=>Q
      double precision TSPECS(*), Y(*)
      double precision F(*)
      external DIVAF, DIVAO
c
      integer KDIM, MAXORD, MAXSTF
c++ Substitute for KDIM, MAXORD, MAXSTF below
      parameter (KDIM = 20, MAXORD = 2, MAXSTF = 1)
c--D Next line special: P=>D, X=>Q
      double precision TN
      double precision XI(KDIM)
c
c--D Next line special: P=>D, X=>Q
      double precision TG(2), TGSTOP(2), TMARK, TMARKX, TOUT, TOLG
      double precision ALPHA(KDIM), BETA(KDIM+1)
      double precision  D(MAXSTF+MAXORD,MAXORD), G(KDIM,MAXORD)
      double precision V(KDIM+MAXORD)
      double precision HC, HDEC, HINC, HINCC, HMAX, HMAXP9, HMIN
      double precision FDAT(11)
c
      double precision DS(MAXSTF+MAXORD, MAXORD), GS(KDIM)
      double precision SIGMA(KDIM), RBQ(KDIM), DNOISE
      double precision EAVE, EIMAX, EIMIN, EMAX, EREP, ROBND, SNOISE
c
c.    SPECIFICATION OF ENVIRONMENTAL CONSTANTS.
      double precision EEPS10, EEPS16, EROV10, EEPS2
      double precision EEPT75, EOVEP2, OVTM75, OVD10
      common / DIVAEV / EEPS2, EEPT75, EOVEP2, OVTM75, OVD10, EEPS10,
     1   EEPS16, EROV10
      save / DIVAEV /
      integer IOPST, KORDI, KQMAXD, KQMAXI, LDT, MAXDIF, MAXINT, NKDKO,
     1   NTE, NYNY, NDTF, NUMDT
      common / DIVASC / TN, XI, IOPST, KORDI, KQMAXD, KQMAXI, LDT,
     1   MAXDIF, MAXINT, NKDKO, NTE, NYNY, NDTF, NUMDT
c
      integer ICF,ICS,IGFLG,IGTYPE(2),IGSTOP(2),ILGREP,INGS,IOP3,IOP4,
     1   IOP5,IOP6,IOP7,IOP8,IOP9,IOP10,IOP11,IOP12,IOP13,IOP14,IOP15,
     2   IOP16,IOP17,IOP18,IOP19,IOP20,IOP21,IOP22,IOP21S,ITOLEP,IY,
     3   KEMAX,KIS,KMARK,KORD1I,KORD2I,KPRED,KQDCON,KQICON,KQMAXS,
     4   KQMXDS,KQMXIL,KQMXIP,KQMXIS,KSC,KSOUT,KSSTRT,KSTEP,LEX,LINC,
     5   LINCD,LINCQ,LSC,MAXKQD,MAXKQI,METHOD,NE,NEPTOL,NG,NGTOT,
     6   NOISEQ,NOUTKO,NTOLF,NY,IDAT(6)
      common /DIVAMC/ TG,TGSTOP,TMARK,TMARKX,TOUT,TOLG,HC,HDEC,HINC,
     1   HINCC,HMAX,HMAXP9,HMIN,ALPHA,BETA,D,G,V,DS,GS,SIGMA,RBQ,DNOISE,
     2   EAVE,EIMAX,EIMIN,EMAX,EREP,ROBND,SNOISE,FDAT,ICF,ICS,IGFLG,
     3   IGTYPE,IGSTOP,ILGREP,INGS,IOP3,IOP4,IOP5,IOP6,IOP7,IOP8,IOP9,
     4   IOP10,IOP11,IOP12,IOP13,IOP14,IOP15,IOP16,IOP17,IOP18,IOP19,
     5   IOP20,IOP21,IOP22,IOP21S,ITOLEP,IY,KEMAX,KIS,KMARK,KORD1I,
     6   KORD2I,KPRED,KQDCON,KQICON,KQMAXS,KQMXDS,KQMXIL,KQMXIP,KQMXIS,
     7   KSC,KSOUT,KSSTRT,KSTEP,LEX,LINC,LINCD,LINCQ,LSC,MAXKQD,MAXKQI,
     8   METHOD,NE,NEPTOL,NG,NGTOT,NOISEQ,NOUTKO,NTOLF,NY,IDAT
      save / DIVAMC / , / DIVASC /
c
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
c
      integer LDIS, KEXIT, I, J, K, J1, J2, L, LX
      double precision TP, TP1, TP2, TP3, HH, DISADJ
      double precision SIGMAS, TPS1, TPS2, EIMINO, EXR
c--D Next line special: P=>D, X=>Q
      double precision  TMARKA(2), XP, XP1
      equivalence (G(1, 1), HH), (TMARKA(1), TMARK)
      equivalence (KEXIT, IOP17)
      save EIMINO, EXR, SIGMAS, TP, TP1, TP2, TP3, TPS1, TPS2, DISADJ,
     1   LDIS
c
c                      Declarations for error message processing.
c
      integer MACT(17), MLOC(8), MENTXT, MERET, MEEMES, METEXT
      parameter (MENTXT =23)
      parameter (MERET  =51)
      parameter (MEEMES =52)
      parameter (METEXT =53)
c
c ********* Error message text ***************
c[Last 2 letters of Param. name]  [Text generating message.]
cAA DIVAA$B
cAB At: TN=$F, KSTEP=$I, with H=$F$E
cAC A previously reported error was fatal.$E
cAD Print points not properly ordered: TSPEC($I)=$F$E
cAE An error tolerance of 0 requires setting special flags.  $B
cAF Step size reduced too fast, doing a restart.  $B
cAG H is so small that TN + H = TN.  $B
cAH Error tolerance too small.  $B
cAI Step size at end of start < HMIN=$F, $B
cAJ Error estimates require a stepsize < HMIN=$F, $B
cAK (Estimated Error) / (Requested Error) for equation $I is $F.  $B
cAL Tolerance $I is F($I) = $F.$B
cAM Tolerance $I is F($I) * F($I) = $F * $F = $F.$B
cAN   Replacing F($I) with $F.$E
      integer LTXTAA,LTXTAB,LTXTAC,LTXTAD,LTXTAE,LTXTAF,LTXTAG,LTXTAH,
     * LTXTAI,LTXTAJ,LTXTAK,LTXTAL,LTXTAM,LTXTAN
      parameter (LTXTAA=  1,LTXTAB=  8,LTXTAC= 40,LTXTAD= 80,LTXTAE=129,
     * LTXTAF=189,LTXTAG=237,LTXTAH=272,LTXTAI=302,LTXTAJ=342,
     * LTXTAK=390,LTXTAL=455,LTXTAM=484,LTXTAN=531)
      character MTXTAA(3) * (186)
c
      integer LOCM, MLOCAC, MLOCAD, MLOCAE, MLOCAF, MLOCAG, MLOCAH,
     1   MLOCAI, MLOCAJ

      parameter (LOCM = 32 * 256)
c                     KORD1I   Severity   Loc. message
      parameter (MLOCAC = 23 + 32 * (99 + 256 * LTXTAC))
      parameter (MLOCAD = 13 + 32 * (38 + 256 * LTXTAD))
      parameter (MLOCAE = 22 + 32 * (38 + 256 * LTXTAE))
      parameter (MLOCAF =  3 + 32 * (14 + 256 * LTXTAF))
      parameter (MLOCAG = 21 + 32 * (38 + 256 * LTXTAG))
      parameter (MLOCAH =  2 + 32 * (25 + 256 * LTXTAH))
      parameter (MLOCAI = 11 + 32 * (38 + 256 * LTXTAI))
      parameter (MLOCAJ = 12 + 32 * (38 + 256 * LTXTAJ))
c
      data MTXTAA/'DIVAA$BAt: TN=$F, KSTEP=$I, with H=$F$EA previously r
     *eported error was fatal.$EPrint points not properly ordered: TSPEC
     *($I)=$F$EAn error tolerance of 0 requires setting special flags. $
     * ','$BStep size reduced too fast, doing a restart.  $BH is so smal
     *l that TN + H = TN.  $BError tolerance too small.  $BStep size at$
     * end of start < HMIN=$F, $BError estimates require a steps','ize <
     * HMIN=$F, $B(Estimated Error) / (Requested Error) for equation $I$
     * is $F.  $BTolerance $I is F($I) = $F.$BTolerance $I is F($I) * F(
     *$I) = $F * $F = $F.$B  Replacing F($I) with $F.$E'/
      data MLOC / MLOCAC, MLOCAD, MLOCAE, MLOCAF, MLOCAG, MLOCAH,
     1   MLOCAI, MLOCAJ /
c
c                      1 2 3 4       5 6       7       8 9      10
      data MACT / MEEMES,0,0,0, MENTXT,0, METEXT, MENTXT,0, METEXT,
     1  MENTXT, 0, METEXT, MENTXT, LTXTAN, METEXT, MERET /
c           11 12      13      14      15      16     17
c
      data EXR, EIMINO / CP1, C8M3 /
      data LDIS / 0 /
c
c ************** START OF EXECUTABLE CODE ******************************
c
  660 if (KORD2I) 670, 1380, 1840
  670 if (KORD2I .eq. -1) go to 2140
      if (KORD2I .eq. -5) go to 720
      if (KORD2I+8) 1840, 710, 1840
c     **** SPECIAL OUTPUT CASE (EXTRAPOLATION OR GSTOP)
  680 go to (1220, 1190, 830, 690, 690, 2190), KEXIT
c     **** RESET TMARK BEFORE GOING WHERE DIRECTED BY KEXIT
  690 KEXIT = KEXIT - 2
      go to 1890
c     **** USER HAS REQUESTED A RESTART
  700 KORD2I = -5
      IGFLG = 0
      if ((KORD(2) .ge. 2) .and. (LSC .lt. 3)) KORD2I = -8
      if ((KORD1I .le. 3) .or. (KORD1I .eq. 5)) go to 1890
      if (KORD2I .eq. -5) go to 720
c                Set up for a discontinuity
  710 XP = TSPECS(1)
      if (KORD(2) .eq. 3) then
c                Adjust for discontinuity in Y
         J = 0
         do 712 I = 1, NTE
            if (NKDKO .ne. 0) KORDI = KORD(NKDKO + I -1)
            K = 1
            J = J + KORDI
            XP1 = Y(J)
  711       Y(NYNY + J - K) = Y(NYNY + J - K) + XP1
            if (K .lt. KORDI) then
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
      if (XP1 .lt. CP25) then
         K = 1
         if (XP1 .lt. CMP75) K = 2
         TSPECS(1) = TN - XI(K)
         TSPECS(2) = 2.D0*(TN - TSPECS(1))
         call DIVAIN(TSPECS(1), Y, F, KORD)
         do 713 J = 1, NY
            Y(NYNY + J - 1) = Y(J)
  713    continue
c          Move difference tables back one step
         TN = TSPECS(1)
  714    call DIVABU(F, KORD)
         if (K .eq. 2) then
            KSC = max(KSC - 1, 1)
            do 715 K = max(1, KSC), IOP11-1
               BETA(K+1) = BETA(K) * (XI(K) / (XI(K+1) - XI(1)))
  715       continue
            K = 0
            go to 714
         end if
      end if
c      Take step to discontinuity.
      HH = (XP - TN)
      EREP = -abs(EREP)
      KIS = 1000
      LDIS = 1
      LSC = 0
      LINC = 0
      HINCC = C1P125
      KSSTRT = KSTEP + 2
      go to 1115
c ********
c INITIALIZE FOR STARTING AN INTEGRATION
c ********
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


c GO COMPUTE INITIAL DERIVATIVES
  730 KORD2I = -6
c++  Code for VAREQ is active
      ICS = 0
      ICF = NE
c++  End
  740 KORD1I = KPRED
      go to 1360
c   RETURN AFTER COMPUTING INITIAL (OR NOISE TEST) DERIVATIVES
  750 continue
c++  Code for VAREQ is active
      if (IOP18 .eq. 0) go to 790
c.    **** SPECIAL LOGIC TO COMPUTE VARIATIONAL DERIVATIVES
  760 if (KORD1I .eq. 3) go to 770
      KORD1I = 3
      KORD(3) = 0
  770 if (ICF .eq. NTE) go to 790
      if (KORD2I .eq. -6) then
         if (ICS .eq. 1) go to 790
      end if
      ICS = ICF + 1
  780 KORD(3) = KORD(3) + 1
      ICF = IOP18 + KORD(3)
      ICF = abs(KORD(ICF))
      if (ICF) 1360, 780, 1360
c++  End
  790 ICS = 1
      ICF = NE
c++  Code for VAREQ is active
      if (KORD2I .eq. 0) if (EREP) 2220, 2220, 1430
c++  End
      if (LINC + 5) 1490, 810, 1490
c END OF SPECIAL CODE FOR INITIALIZATION AT THE START
c ********
c UPDATE VARIABLES TO PREPARE FOR NEXT STEP
c ********
  800 LDT = 0
      EIMIN = C2P5M3+EIMIN*(C6*EAVE+EIMAX)/((EIMIN+C1)*(C6*EIMAX+EAVE))
  810 do 820 J = 1, NY
  820    Y(NYNY + J - 1) = Y(J)
      TN = TSPECS(1)
c ********
c TEST FOR VARIOUS TYPES OF OUTPUT
c ********
c     TEST FOR END OF STEP OUTPUT (OR IF DUMP OUTPUT TO BE TESTED FOR)
      if (IOP6 .eq. 0) go to 840
c   SET UP FOR END OF STEP OUTPUT
      KORD1I = 6
      go to 1810
c     SET UP AFTER OTHER COMPUTATIONS HAVE BEEN MADE (TSPECS(1).NE.TN)
  830 KORD1I = 6
      go to 860
c     TEST FOR AN OUTPUT POINT
  840 if (HH * (TMARK - TN)) 1760, 1760, 850
c     TEST FOR TOO MANY STEPS OUTPUT
  850 continue
      if (KSOUT .gt. KSTEP) go to 890
      KORD1I = 4
  860 if (TSPECS(1) .eq. TN) go to 1810
c     GO INTERPOLATE VALUES AT END OF LAST STEP
  870 TSPECS(1) = TN
      go to 1780
c     CONTINUE AFTER TOO MANY STEPS OUTPUT
  880 continue
      KSOUT = KSTEP + IOP4
  890 continue
c++  Code for DUMP is active
      if (IOP9 .eq. 0) go to 920
c++  Code for DUMP & STIFF is inactive
c      KQMXDS=KQMAXD
c++  Code for DUMP is active
      KQMXIS = KQMAXI
      KIS = KIS + 1
      if (KIS .ne. 0) go to 920
c.   TIME TO DUMP THE SOLUTION
  900 KORD1I = 9
c.    SET TO UPDATE THE DIFFERENCE TABLE
      if (LDT .eq. 1) go to 1810
c Note that doing this update can lead to very small differences in the
c results because of round off differences.
      LDT = -2
      go to 1780
c++  End
c RETURN AFTER DUMPING THE SOLUTION
  910 continue
c++  Code for DUMP is active
      KIS = 2
c.    TEST IF SOLUTION DUMP DUE TO RESTART, END, OR
c.    DROP IN INTEG. ORDER
      if (LINC .lt. 0) if (LINC + 6) 1860, 1750, 1180
c++  End
c END OF TESTING FOR VARIOUS TYPES OF OUTPUT
c   TEST IF STEPSIZE MAY BE INCREASED OR IF TESTS SHOULD BE MADE
c   FOR DECREASING THE STEPSIZE (OR IF STARTING OR USER SELECTING H)
  920 KSTEP = KSTEP + 1
      if (LINC) 930, 980, 940
  930 continue
c     **** ONLY POSSIBLE VALUES AT THIS POINT ARE
c          LINC = -2 OR -5
      if (LINC + 5) 1120, 1110, 1120
c ********
c ERROR ESTIMATES INDICATE STEPSIZE CAN BE INCREASED
c ********
  940 HC = HINCC
      if (LINC .gt. 1) HC = HC ** LINC
      HH = HC * HH
c     TEST IF NEW STEPSIZE IS TOO BIG
      if (abs(HH) .gt. HMAX) if (HMAX) 970, 970, 960
  950 EAVE = EIMAX
      ROBND = CP3 + HINCC
      go to 1110
c     NEW STEPSIZE IS TOO BIG
  960 if (abs(XI(1)) .ge. HMAXP9) go to 970
      HH = sign(HMAX, HH)
      go to 950
c     RESTORE THE OLD STEPSIZE
  970 HH = XI(1)
      LINC = 0
      go to 1150
c END OF CODE FOR CASE WHEN ERROR ESTIMATES INDICATE STEPSIZE INCREASE
c ********
c TEST IF ESTIMATED ERRORS INDICATE STEPSIZE SHOULD BE DECREASED
c ********
  980 ROBND = C1P3
      if (EIMAX .le. EAVE) go to 990
      EAVE = EAVE + CP4 * (EIMAX - EAVE)
      if ((EIMAX * EMAX) - C1M3) 1000, 1010, 1010
  990 EAVE = EIMAX
      if ((EIMAX * EMAX) .ge. EIMIN) go to 1010
 1000 ROBND = CP3 + (SIGMA(KQMAXS)/SIGMA(KQMAXS-1))
      go to 1180
c     TEST IF STEPSIZE SHOULD BE REDUCED
 1010 if (EMAX * EIMAX .lt. EXR * EAVE) go to 1180
c ********
c ERROR ESTIMATES INDICATE STEPSIZE SHOULD BE REDUCED
c ********
      HC = HDEC
      if (EIMIN .le. EIMINO) go to 1030
      EIMIN = EIMINO
      HC = CP875
 1030 HH = HC * XI(1)
      if (LSC - 1) 1040, 1080, 1090
 1040 if (abs(HH) .ge. HMIN) go to 1090
      if (abs(CP875*XI(1)) .le. HMIN) if (LINC) 1050, 970, 970
      HH = sign(HMIN, HH)
      go to 1090
c     STEPSIZE IS TOO SMALL TO BE REDUCED
c     SET UP ERROR INDICATORS AND PREPARE FOR RETURN TO USER
 1050 KORD1I = 8
      go to 2240
c     PROCEED WITH CURRENT STEPSIZE DESPITE ERROR BEING TOO BIG
 1070 HH = XI(1)
      TSPECS(2) = HH
      EMAX = C0
      LINC = 0
      if (KORD1I - 2) 1420, 1430, 1430
c     SET LSC TO END STARTING PHASE
 1080 LSC = 2
c     CHECK IF REPEATING A STEP
 1090 if (LINC .ne. -1) go to 1110
c   WHEN REPEATING A STEP, BACK UP THE DIFFERENCES AND STEPSIZE INFO.
 1100 call DIVABU(F, KORD)
c   TEST IF NOISE TEST (LINC = -7) OR IF H IS NOT
c     BEING CHANGED (LINC = -4)
      if (LINC + 4) 1780, 1180, 1110
c ********
c STEPSIZE IS TO BE CHANGED
c ********
 1110 continue
c MODIFY STEPSIZE TO REDUCE ROUNDOFF ERROR IN ACCUMULATING INDEP. VAR.
      TP = C2 * abs(TN) + C4096 * abs(HH)
      TP = (TP + abs(HH)) - TP
      if (TP .ne. C0) HH = sign(TP, HH)
c     TEST IF NEW STEPSIZE SELECTED ACTUALLY GIVES A CHANGE
      if (HH .eq. TSPECS(2)) go to 1140
 1115 TSPECS(2) = HH
      if (IOP8 .eq. 0) go to 1140
c     SETUP TO TELL USER ABOUT STEPSIZE CHANGE (OR TO CHANGE STEPSIZE)
 1120 KORD1I = 8
      go to 860
c     RETURN AFTER TELLING USER ABOUT STEPSIZE CHANGE
 1130 HH = TSPECS(2)
 1140 if (HH .ne. XI(1)) KQICON = -1
 1150 HC = min(EOVEP2, abs(HH)) / EEPS2
c ********
c PREPARE FOR BEGINNING A NEW STEP
c ********
      if (LINC .gt. 0) then
         LINC = min(LINC, LINCQ) + LINCQ
         go to 1190
      end if
 1180 LINC = LINCD
 1190 if (HC .gt. abs(TN)) go to 1200
c     **** GIVE SINGULARITY DIAGNOSTIC
      KORD1I = 5
      go to 2240
 1200 TSPECS(1) = TN + HH
      if (LEX .eq. 0) go to 1250
      if (HH * (TSPECS(1) - TMARKX) .lt. C0) go to 1250
      TSPECS(1) = TMARKX
      HH = TMARKX - TN
      LINC = 64
      if (LEX .gt. 0) go to 1240
      if ((LSC .lt. 4) .and. (HH / XI(1) .lt. CP3)) go to 1230
 1220 HH = CP875 * HH
      go to 1110
c     **** GIVE OUTPUT AT CURRENT TMARK (WITH EXTRAPOLATION)
 1230 KORD1I = -KMARK
      go to 1770
c     **** INTEGRATE TO TMARKX
 1240 KQICON = -1
c   TEST IF SUBROUTINE FOR COMPUTING INTEGRATION COEFF. SHOULD BE CALLED
 1250 continue
c++  Code for STIFF is inactive
c      IF ((KQMAXI .LT.KQICON) .OR. (KQMAXD.LT.KQDCON)) GO TO 1320
c++  Code for ~STIFF is active
      if (KQMAXI .lt. KQICON) go to 1320
c++  End
c   GO COMPUTE COEFFICIENTS REQUIRED FOR THE INTEGRATION
c     TEST IF STARTING
      if (LSC .lt. 7) go to 1310
 1260 KQMAXI = 2
c++  Code for STIFF is inactive
c      IF (METHOD) 1262,1270,1264
c 1262 KQMAXI=0
c 1264 KQMAXD=max(MAXDIF,2)
c      CALL DIVAHC
cc.  SET UP TO GO DO INITIALIZATION FOR CASE OF STIFF EQUATIONS
c      KORD1I=5
c      GO TO 1350
c++  End
c   INITIALIZE FOR EQUATIONS WHICH ARE NOT STIFF
 1270 KQMAXD = 0
      call DIVAHC
      J = NDTF
      do 1300 I = 1, NTE
c++  Code for STIFF is inactive
c         if (KORD(I + 3) .le. 0) go to 1290
c++  End
         KORD(I + 3) = 1
c     INITIALIZE THE DIFFERENCE TABLE
         if (LDT .eq. -4) F(J) = F(I)
         F(J + 1) = C0
         F(J + 2) = C0
 1290    continue
         J = J + NUMDT
 1300    continue
      if (LSC .eq. 5) go to 1340
      LSC = 7
      LDT = 1
      go to 1330
c   INTEGRATION IS NOT BEING STARTED
 1310 K = KORD(KEMAX + 3)
      SIGMAS = SIGMA(K)
      call DIVAHC
c     **** ADJUST EAVE
      TPS1 = BETA(K)
      if (TPS1 .gt. C1) TPS1 = CP5 * TPS1 + CP5
      EAVE = EAVE * TPS1 * (SIGMA(K) / SIGMAS)
c     TEST BELOW USED TO GET SAME RESULTS WITH/WITHOUT EXTRA EQUATIONS
      if (K .gt. KQICON) LSC = max(LSC, -3)
c END OF SPECIAL LOGIC FOR CASE WHEN INTEG. COEFF. ROUTINE IS CALLED
 1320 continue
c ********
c PREDICT Y
c ********
 1330 continue
c++  Code for ~ARGM is active
      call DIVAPR(Y, Y(NYNY), F, KORD)
c++  Code for ARGM is inactive
c      CALL DIVAPE
c++  End
c     GO GET PREDICTED DERIVATIVES
 1340 KORD1I = KPRED
c ********
c CALL DIVAF  (OR RETURN)
c ********
 1350 KORD2I = 0
 1360 KORD(1) = KORD1I
      KORD(2) = 0
      if (IOP13 .ne. 0) return
      call DIVAF(TSPECS(1), Y, F, KORD(1))
c     TEST FOR SPECIAL USER RETURN
 1380 if (KORD(1) .lt. 0) go to 2130
c     TEST FOR SPECIAL CASE
      if (KORD2I .ne. 0) go to 660
c ********
c TRANSFER CONTROL TO PROPER PLACE AFTER COMPUTING DERIVATIVES
c ********
      if (KORD1I - 2) 1400, 800, 1390
 1390 continue
c++  Code for VAREQ is active
      if (ICS - ICF) 1410, 1410, 760
c++  End
c ********
c PREPARE FOR CORRECTING, AND CORRECT Y
c ********
 1400 ITOLEP = 0
      ILGREP = 0
      IY = 1
      EIMAX = C1M5
      EMAX = C0
      KQMAXI = 2
      KQMAXS = 2
c++  Code for STIFF is inactive
c      IF (METHOD) 1404,1410,1406
c 1404 KQMAXI=0
c 1406 KQMAXD=2
c++  End
 1410 continue
c++  Code for ~ARGM is active
      call DIVACR(Y, F, KORD, F(NTOLF), KORD(IOP16))
c++  Code for ARGM is inactive
c      CALL DIVACE
c++  End
c     TEST IF ESTIMATED ERROR IS TOO BIG (OR IF DIAGNOSTIC CALLED FOR)
      if (EMAX .gt. EREP) if (EREP) 2210, 2210, 1670
 1420 continue
c++  Code for VAREQ is active
      if (IOP18 .ne. 0) go to 760
c++  End
 1430 KORD1I = 2
c     TEST IF NOISE APPEARS TO LIMIT PRECISION
      if (EMAX .lt. C0) go to 1470
c++  Code for ~STIFF is active
      if (LSC) 1450, 1360, 1610
c++  Code for STIFF is inactive
c      IF (LSC) 1450,1460,1610
c++  End
c     SET LSC=0 IF NOISE NO LONGER APPEARS TO LIMIT PRECISION
c     OR IF THE END OF THE STARTING PHASE HAS BEEN REACHED
 1450 LSC = 0
 1460 if (METHOD) 800, 1350, 1350
c ********
c NOISE APPEARS TO BE LIMITING PRECISION
c ********
 1470 continue
      if (LSC .le. 0) LSC = max(LSC - 1, -KQMAXS)
      if (LSC .eq. -1) go to 1460
      if (abs(EMAX) .lt. EXR) go to 1590
      LINC = -7
      TPS2 = (C1 + BETA(NOISEQ - 1)) ** NOISEQ
      if (SNOISE .lt. EEPS10 * TPS2) go to 1550
      TP = sign(EEPT75 * abs(TN) + OVTM75, HH)
      if (abs(TP) .gt. abs(HH)) go to 1550
      TSPECS(1) = TN + TP
      KORD1I = 0
c     **** GO TO BACK UP THE DIFFERENCES AND GET F(TSPECS(1))
      go to 1100
c     **** SOLUTION HAS BEEN INTERPOLATED AND F COMPUTED
 1490 continue
      KORD1I = 0
      LINC = LINC - 1
      if (LINC + 9) 1510, 1520, 1500
 1500 TSPECS(1) = TN + (TP + TP)
      TP1 = F(KEMAX)
      TP2 = F(NDTF + NUMDT * KEMAX - NUMDT)
      go to 1780
c     **** COMPUTE 2-ND DIFFERENCE AT CLOSE SPACED T VALUES
 1510 TP2 = TP3
 1520 TP3 = F(KEMAX)
      TPS1 = abs((TP3 - TP1) - (TP1 - TP2))
      if ((C16 * TPS1 * TPS2) .ge. DNOISE) if (LINC + 9) 1550, 870, 1550
 1530 continue
      TPS2 = CP25 * SNOISE / RBQ(NOISEQ)
      do 1540 K = 2, NUMDT
         TPS1 = TPS1 + TPS1
         RBQ(K) = max(TPS1, TPS2 * RBQ(K))
 1540    continue
      LINC = 0

cFTK Next two lines added 2009-10-15
      if (abs(EMAX) .lt. EREP) go to 1460
cFTK  LINC = -1  And then on 2015-03-14 commented out this line

      HH = CP875 * HH
      go to 1040
c     **** SET UP TO GIVE NOISE DIAGNOSTIC
 1550 KORD1I = 6
      go to 2240
c     **** AFTER GIVING NOISE DIAGNOSTIC
 1560 KORD1I = 2
      if (KORD(2) .ge. 0) then
        TPS1 = EEPS10
cFTK Next line added 2009-10-15
        if (TPS1 .lt. .49D0 * RBQ(2)) go to 1530
      end if
c     **** SET NEW VALUE FOR OFFENDING TOL
      F(NTOLF + ITOLEP - 1) = FDAT(7)
      if (LINC + 7) 1180, 1570, 1180
 1570 LINC = 0
 1580 if (LSC) 1460, 1460, 1610
c     **** CHANGE HINCC AND ADJUST SIGMA( )
 1590 if (LSC .ne. -4) go to 1580
      if (HINCC .eq. C1P125) go to 1580
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
c   END OF CODE FOR CASE WHEN NOISE APPEARS TO LIMIT PRECISION
c ********
c SPECIAL LOGIC FOR STARTING THE INTEGRATION
c ********
 1610 if (LSC .eq. 1) go to 800
      LSC = LSC - 1
      if (LSC - 2) 1620, 1640, 1650
 1620 if (EIMAX .le. (CP0625*EAVE*(SIGMA(KQMAXS)/SIGMAS)*(BETA(KQMAXS+
     1   1))**2)) go to 800
 1630 KSSTRT = KSTEP + 2
c   TEST IF STEPSIZE IS TOO SMALL BEFORE ENDING STARTING PHASE
      if (abs(HH) .ge. HMIN) go to 1450
c     GIVE DIAGNOSTIC FOR STEPSIZE TOO SMALL AT END OF START
      KORD1I = 7
      go to 2240
c   SET LSC TO DO ONE DERIVATIVE EVAL. PER STEP
 1640 LSC = 1
      go to 800
c     TEST IF FIRST TIME THROUGH THE FIRST STEP
 1650 if (LSC .eq. 6) go to 1340
c     END STARTING PHASE IF CONVERGENCE OF CORRECTOR ITERATES TOO SLOW
      if (LDT .eq. -5) go to 1660
      LSC = min(LSC, 4)
      go to 800
 1660 LDT = 0
      if (LSC - 4) 1260, 1630, 1260
c END OF SPECIAL LOGIC FOR STARTING THE INTEGRATION
c ********
c ESTIMATED ERROR IS TOO BIG
c ********
 1670 if (BETA(2) - C1) 1690, 1730, 1680
 1680 HC = C1 / BETA(2)
      if (BETA(2) .ge. C1P125) go to 1740
 1690 if (BETA(2) .gt. CP1) go to 1730
c   REQUIRED STEPSIZE REDUCTION IS TOO RAPID -- GIVE A DIAGNOSTIC
      KORD1I = 4
      go to 2240
c
c     TEST KORD(2) AFTER ABOVE DIAGNOSTIC OR A DISCONTINUITY DIAGNOSTIC
 1700 continue
      if (KORD(2) .eq. 0) go to 1730
c  TEST IF SOLUTION MUST BE DUMPED BEFORE A RESTART
 1710 LINC = -1
c++  Code for DUMP is active
      if (IOP9 .eq. 0) go to 1750
      if (KIS .eq. 2) go to 1750
      LINC = -6
c.    GO DUMP SOLUTION BEFORE REPEATING THE STEP
 1720 KQMAXI = KQMXIS
c++  Code for DUMP & STIFF is inactive
c      KQMAXD=KQMXDS
c++  Code for DUMP is active
      call DIVABU(F, KORD)
      go to 900
c++  End
c   SET UP TO REPEAT THE STEP
 1730 HC = CP5
 1740 LINC = -1
      if (LSC .le. 3) go to 1030
c   RESTART THE INTEGRATION IF ERROR IS TOO BIG ON FIRST OR SECOND STEP
c LOOP TO SELECT A NEW INITIAL STEPSIZE
 1750 LSC = 7
 1755 HH = HH * CP5
      EMAX = EMAX * CP25
      if (EMAX .ge. CP3) go to 1755
      go to 1090
c   END OF SELECTING A NEW INITIAL STEPSIZE
c END OF LOGIC FOR CASE WHEN ESTIMATED ERROR IS TOO BIG
c ********
c INTEGRATION HAS REACHED AN OUTPUT POINT
c ********
 1760 if (KMARK .eq. 0) go to 1920
      KORD1I = min(KMARK, 5)
      KORD(3) = KMARK
      if (TSPECS(1) .eq. TMARK) go to 1790
 1770 TSPECS(1) = TMARK
 1780 call DIVAIN(TSPECS(1), Y, F, KORD)
 1790 continue
      if (KORD1I) 1800, 730, 1810
c   OUTPUT POINT IS OBTAINED BY EXTRAPOLATION
 1800 continue
c++  Code for EXTRAP is active
      KORD1I = -KORD1I
      KORD2I = -7
      KEXIT = 4
c.  TEST IF GSTOP-S ARE PRESENT
c++  Code for EXTRAP &  GSTOP is active
      if (NGTOT .eq. 0) go to 1820
      IGFLG = 4
      KEXIT = 2
      KORD1I = 7
      if (IOP7) 740, 1820, 740
c++  End
c ********
c CALL DIVAO  (OR RETURN)
c ********
 1810 KORD2I = 1
 1820 KORD(1) = KORD1I
      KORD(2) = 1
      if (IOP14 .ne. 0) return
      call DIVAO(TSPECS(1), Y, F, KORD(1))
c     TEST FOR SPECIAL USER RETURN OR FOR A RESTART
c++  Code for ~DUMP is inactive
c 1840 IF (KORD(1)) 2130,700,1880
c++  Code for DUMP is active
 1840 if (KORD(1) .gt. 0) go to 1880
 1850 if (IOP9 .eq. 0) go to 1870
c.    **** GO DUMP THE SOLUTION
      LINC = -7
      ITOLEP = KORD(1)
      IDAT(1) = KORD(2)
      NEPTOL = KORD1I
      if (LSC .ne. 8) go to 900
 1860 LINC = min(0, LINCD)
      KORD1I = NEPTOL
      KORD(1) = ITOLEP
      KORD(2) = IDAT(1)
 1870 if (KORD(1)) 2130, 700, 2100
c++  End
 1880 if (KORD2I .lt. 0) go to (2140, 1810, 1350, 2110, 720, 750, 680,
     1   710), -KORD2I
      if (KORD2I .eq. 0) go to 1380
c ********
c TRANSFER CONTROL TO PROPER PLACE AFTER OUTPUT
c ********
 1890 if (KORD1I - 5) 1910, 1930, 1900
 1900 if (KORD1I - 8) 840, 1130, 910
 1910 if (KORD1I - 3) 1920, 1930, 880
c   GET NEW TOUT
 1920 TOUT = TSPECS(1) + TSPECS(3)
c GET NEW TMARK (NEXT INDEP. VAR. OUTPUT POINT)
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
      if (J2 .ge. J1) go to 1990
 1980 J1 = 4
      J2 = 4
      L = IOP3
c
c     **** LOOP TO SET NEW TMARK (AND TMARKX)
 1990 do 2060 J = J1, J2
c        **** TEST IF EXTRAPOLATION NOT POSSIBLE
         if (L .eq. 0) go to 2010
         LX = 2
         if (LEX) 2020, 2030, 2020
 2000    LEX = L
 2010    LX = 1
 2020    if (HH * (TSPECS(J) - TMARKA(LX))) 2030, 2060, 2060
 2030    if (J .eq. 4) go to 2050
         if (HH * (TSPECS(J) - XP)) 2060, 2040, 2050
 2040    if ((K .ge. J) .or. (K .eq. 3)) go to 2060
 2050    TMARKA(LX) = TSPECS(J)
         if (LX .eq. 2) go to 2000
         KMARK = J
 2060    continue
      if (IOP5 .lt. 0) go to 1950
      if (J1 .ne. 4) go to 1980
      if (KMARK .eq. 4) KMARK = 3
c     **** TEST IF NEW TMARK IS ACCEPTABLE
      if (HH * (XP - TMARK)) 2070, 2080, 2090
 2070 if (KORD2I - 1) 670, 840, 670
 2080 if (K .ne. KMARK) go to 2070
c++  Code for DUMP is active
      if (KORD1I .eq. 3) go to 1850
c++  Code for ~DUMP is inactive
c      IF (KORD1I .EQ. 3) GO TO 2100
c++  End
 2090 if (KORD1I .eq. 13) go to 2190
c SETUP TO INDICATE ERROR IN SPECIFICATION OF OUTPUT POINTS
      KORD1I = 2
      IDAT(2) = KMARK
      if (KMARK .le. 3) IDAT(2) = KMARK + 1
      FDAT(3) = TSPECS(IDAT(2))
      go to 2240
c     SET KORD1I=1 TO INDICATE THAT END OF INTEGRATION HAS BEEN REACHED
 2100 KORD1I = 1
c ********
c RETURN TO USER
c ********
 2110 KORD2I = -1
      KORD(1) = KORD1I
 2130 KORD(2) = -1
      return
c ********
c TRANSFER CONTROL TO PROPER PLACE AFTER RETURN TO USER
c ********
 2140 if (KORD1I - 2) 2150, 1560, 2160
 2150 KORD2I = 1
      go to 1930
 2160 if (KORD1I - 4) 1700, 2200, 2170
 2170 if (KORD1I - 13) 2180, 1930, 2190
 2180 if (abs(HH) .ge. HMIN) if (KORD1I - 11) 1030, 1450, 1030
      if (KORD(2) .eq. 0) if (KORD1I - 11) 1070, 800, 1070
c   ERROR MESSAGES HAVE BEEN IGNORED -- COMPUTATION CAN NOT CONTINUE
 2190 KORD1I = 1
      go to 2240
c
c        AFTER A DISCONTINUITY RETURN
 2200 LINC = -4
      if (KORD(2)) 1710, 1730, 1100
c ********
c PROBLEM ENCOUNTERED WHEN CORRECTING
c ********
 2210 if (LDIS .eq. 0) go to 2230
c           Extra checks when had a user specified discontinuity.
c++  Code for VAREQ is active
      if (IOP18 .ne. 0) go to 760
c++  End
 2220 KORD1I = 2
      LDIS = LDIS + 1
      TP = DISADJ / HH
      if (KIS .ge. 1000) then
         if (LDIS .eq. 2) then
            if (KQMAXS .le. 3) then
               LDIS = 0
               EREP = abs(EREP)
               TSPECS(2) = HH*min(min(TP, TP**2),
     1            (CP25 * EXR/EMAX)**.333333333D0)
               go to 720
            end if
            LINC = -5
            if (IOP9 .eq. 0) KIS = 1001
            go to 800
         end if
         if (IOP9 .eq. 0) KIS = KIS + 1
         if (KQMAXS .le. LDIS + 2) KIS = LDIS + 1
         LINC = min(LINC, LDIS-2)
      end if
      if (LDIS .gt. 2*KQMAXS) then
         EREP = abs(EREP)
         LDIS = 0
         if (EMAX .gt. EREP) go to 1670
         go to 1430
      end if
      if (TP .ge. HINCC**(LINC+2)) then
         if ((LDIS .ne. 3) .and. (TP .gt. dble(KQMAXS))) LSC = 1
         EIMIN = CP5
         EAVE = EAVE * TP**8
      end if
      if (LSC .eq. 2) go to 1630
      if (EMAX .gt. EXR) go to 1730
      go to 1430
c
 2230 EREP = abs(EREP)
c++  Code for DUMP is active
      if (LINC .lt. -3) go to 1720
c++  End
c     BAD TOL
      KORD1I = 3
c ********
c ERROR PROCESSING
c ********
 2240 FDAT(1) = TN
      FDAT(2) = HH
      IDAT(1) = KSTEP
      ITOLEP = max(NEPTOL, -NEPTOL - 1)
      J = 3
      if (KORD1I .ge. 7) then
         J = 4
         FDAT(3) = HMIN
      end if
      if (KORD1I .le. 3) then
         if (KORD1I .lt. 3) then
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
         if (KORD1I .eq. 6) then
            K = 17
            IDAT(5) = IDAT(4)
            FDAT(7) = 32.D0 * abs(EMAX) * FDAT(J+1)
            FDAT(J+2) = FDAT(7)
         end if
         MACT(12) = LTXTAL
         if (NEPTOL .lt. 0) then
            MACT(12) = LTXTAM
            IDAT(6) = IDAT(4)
            IDAT(5) = IDAT(4) + 1
            FDAT(J+2) = F(IDAT(5))
            FDAT(J+3) = FDAT(J+1) * FDAT(J+2)
         end if
      end if
c Set the location for the first part of the message that varies, set
c the error severity, and the index number, print the error and
c return or stop.
      L = MLOC(KORD1I)
      MACT(6) = L / LOCM
      MACT(2) = (L - MACT(6) * LOCM) / 32
      KORD1I = mod(L, 32)
      MACT(3) = KORD1I
      MACT(K) = MERET
c--D Next line special: P=>S, X=>D
      call DMESS(MACT, MTXTAA, IDAT, FDAT)
      MACT(K) = MENTXT
      go to 2110
c
      end
c   End of DIVAA

      subroutine DIVABU(F, KORD)
c>> 1987-12-07 DIVABU Krogh   Initial code.
c
c THIS SUBROUTINE RESTORES THE DIFFERENCE TABLE TO ITS STATE
c AT THE BEGINNING OF THE CURRENT STEP.  IF THE INTEGRATION ORDER
c WAS INCREASED, IT IS REDUCED. THE COMMON ARRAY XI IS ALSO
c RESTORED TO ITS STATE AT THE BEGINNING OF THE STEP. IF THE
c STEPSIZE IS NOT BEING CHANGED, THE ARRAY V USED TO COMPUTE
c INTEGRATION COEFFICIENTS IS RESTORED.
c
      integer KORD(*)
      double precision F(*)
c
      integer KDIM, MAXORD, MAXSTF
c++ Substitute for KDIM, MAXORD, MAXSTF below
      parameter (KDIM = 20, MAXORD = 2, MAXSTF = 1)
c--D Next line special: P=>D, X=>Q
      double precision TN
      double precision XI(KDIM)
c
c--D Next line special: P=>D, X=>Q
      double precision TG(2), TGSTOP(2), TMARK, TMARKX, TOUT, TOLG
      double precision ALPHA(KDIM), BETA(KDIM+1)
      double precision  D(MAXSTF+MAXORD,MAXORD), G(KDIM,MAXORD)
      double precision V(KDIM+MAXORD)
      double precision HC, HDEC, HINC, HINCC, HMAX, HMAXP9, HMIN
      double precision FDAT(11)
c
      double precision DS(MAXSTF+MAXORD, MAXORD), GS(KDIM)
      double precision SIGMA(KDIM), RBQ(KDIM), DNOISE
      double precision EAVE, EIMAX, EIMIN, EMAX, EREP, ROBND, SNOISE
c
      integer IOPST, KORDI, KQMAXD, KQMAXI, LDT, MAXDIF, MAXINT, NKDKO,
     1   NTE, NYNY, NDTF, NUMDT
      common / DIVASC / TN, XI, IOPST, KORDI, KQMAXD, KQMAXI, LDT,
     1   MAXDIF, MAXINT, NKDKO, NTE, NYNY, NDTF, NUMDT
c
      integer ICF,ICS,IGFLG,IGTYPE(2),IGSTOP(2),ILGREP,INGS,IOP3,IOP4,
     1   IOP5,IOP6,IOP7,IOP8,IOP9,IOP10,IOP11,IOP12,IOP13,IOP14,IOP15,
     2   IOP16,IOP17,IOP18,IOP19,IOP20,IOP21,IOP22,IOP21S,ITOLEP,IY,
     3   KEMAX,KIS,KMARK,KORD1I,KORD2I,KPRED,KQDCON,KQICON,KQMAXS,
     4   KQMXDS,KQMXIL,KQMXIP,KQMXIS,KSC,KSOUT,KSSTRT,KSTEP,LEX,LINC,
     5   LINCD,LINCQ,LSC,MAXKQD,MAXKQI,METHOD,NE,NEPTOL,NG,NGTOT,
     6   NOISEQ,NOUTKO,NTOLF,NY,IDAT(6)
      common /DIVAMC/ TG,TGSTOP,TMARK,TMARKX,TOUT,TOLG,HC,HDEC,HINC,
     1   HINCC,HMAX,HMAXP9,HMIN,ALPHA,BETA,D,G,V,DS,GS,SIGMA,RBQ,DNOISE,
     2   EAVE,EIMAX,EIMIN,EMAX,EREP,ROBND,SNOISE,FDAT,ICF,ICS,IGFLG,
     3   IGTYPE,IGSTOP,ILGREP,INGS,IOP3,IOP4,IOP5,IOP6,IOP7,IOP8,IOP9,
     4   IOP10,IOP11,IOP12,IOP13,IOP14,IOP15,IOP16,IOP17,IOP18,IOP19,
     5   IOP20,IOP21,IOP22,IOP21S,ITOLEP,IY,KEMAX,KIS,KMARK,KORD1I,
     6   KORD2I,KPRED,KQDCON,KQICON,KQMAXS,KQMXDS,KQMXIL,KQMXIP,KQMXIS,
     7   KSC,KSOUT,KSSTRT,KSTEP,LEX,LINC,LINCD,LINCQ,LSC,MAXKQD,MAXKQI,
     8   METHOD,NE,NEPTOL,NG,NGTOT,NOISEQ,NOUTKO,NTOLF,NY,IDAT
      save / DIVAMC / , / DIVASC /
c
      integer I, L, KQQ, J, K
      double precision TPD, C0, C2
      parameter (C0 = 0.D0)
      parameter (C2 = 2.D0)
c ********* START OF EXECUTABLE CODE **********
c
c ********
c BEGIN LOOP TO BACK UP DIFFERENCE TABLES
c ********
      L = NDTF - 1
      do 2410 I = 1, NTE
         KQQ = KORD(I + 3)
c++  Code for STIFF is inactive
c         IF (KQQ) 2302,2400,2310
cc.           EQUATION IS STIFF
c 2302    IF (LINC.GE.0) GO TO 2310
c         IF (F(L+1+I)) 2306,2308,2304
cc.     ORDER WAS INCREASED, AND THUS MUST BE DECREASED (KQQ.LT.0)
c 2304    KQQ=KQQ+1
c         KORD(I+3) = KQQ
c         GO TO 2308
cc.     ORDER WAS DECREASED
c 2306    KQQ=KQQ-1
c 2308    KQQ=max(2,-KQQ)
c         GO TO 2350
c++  End
c     EQUATION IS NOT STIFF
 2310    if (KQQ .gt. 2) then
            if (F(L + KQQ) .eq. C0) then
c                 ORDER WAS INCREASED, AND THUS MUST BE DECREASED
               KQQ = KQQ - 1
               KORD(I + 3) = KQQ
            end if
         end if
         J = min(KQQ, KSC)
         KQMAXI = max(KQMAXI, KQQ)
         if (KQQ .ne. 1) F(L + KQQ + 1) = 0.D0
c           BACK UP FOR BACKWARD DIFFERENCES
         do 2360 K = 1, J
            F(L + K) = F(L + K) - F(L + K + 1)
 2360    continue
         if (KQQ .gt. KSC) then
c           BACK UP FOR MODIFIED DIVIDED DIFFERENCES
            do 2390 K = J+1, KQQ
               F(L + K) = (F(L+K) - F(L+K+1)) / BETA(K)
 2390       continue
         end if
 2400    F(L + KQQ + 1) = F(L + KQQ + 1) / BETA(KQQ + 1)
         L = L + NUMDT
 2410 continue
c END OF LOOP TO BACK UP DIFFERENCE TABLES
c ********
c BACK UP XI TO BEGINNING OF THE STEP
c ********
      I = KSC + 1
      if (I - IOP11 - 1) 2420, 2440, 2450
 2420 TPD = XI(1)
c                Check below needed when starting?
      if (TPD .eq. XI(2)) go to 2450
      do 2430 K = I, IOP11
 2430    XI(K - 1) = XI(K) - TPD
 2440 XI(IOP11) = C2 * XI(IOP11 - 1)
      if (IOP11 .ne. 2) XI(IOP11) = XI(IOP11) - XI(IOP11 - 2)
 2450 KQICON = -1
      ICF = NE
      ICS = 1
      LDT = 1
      return
      end
c   End of DIVABU

      subroutine DIVACO(ID, RD)
c>> 1987-12-07 DIVACO Krogh   Initial code.
c
c THIS SUBROUTINE RETURNS THE FOLLOWING DATA FROM COMMON
c ID(1) = KEMAX  =  INDEX OF EQUATION WITH LARGEST ERROR ESTIMATE
c ID(2) = KSTEP  =  CURRENT STEP NUMBER
c ID(3) = NUMDT  =  NUMBER OF DIFFERENCES USED FOR EACH EQUATION
c ID(4) =           RESERVED FOR FUTURE USE
c ID(5) =           RESERVED FOR FUTURE USE
c RD(1) = EMAX   =  MAX. RATIO OF ESTIMATED ERROR TO REQUESTED ERROR
c RD(2) =           RESERVED FOR FUTURE USE
c RD(3) =           RESERVED FOR FUTURE USE
c
      integer ID(5)
      double precision RD(3)
c
      integer KDIM, MAXORD, MAXSTF
c++ Substitute for KDIM, MAXORD, MAXSTF below
      parameter (KDIM = 20, MAXORD = 2, MAXSTF = 1)
c--D Next line special: P=>D, X=>Q
      double precision TN
      double precision XI(KDIM)
c
c--D Next line special: P=>D, X=>Q
      double precision TG(2), TGSTOP(2), TMARK, TMARKX, TOUT, TOLG
      double precision ALPHA(KDIM), BETA(KDIM+1)
      double precision  D(MAXSTF+MAXORD,MAXORD), G(KDIM,MAXORD)
      double precision V(KDIM+MAXORD)
      double precision HC, HDEC, HINC, HINCC, HMAX, HMAXP9, HMIN
      double precision FDAT(11)
c
      double precision DS(MAXSTF+MAXORD, MAXORD), GS(KDIM)
      double precision SIGMA(KDIM), RBQ(KDIM), DNOISE
      double precision EAVE, EIMAX, EIMIN, EMAX, EREP, ROBND, SNOISE
c
      integer IOPST, KORDI, KQMAXD, KQMAXI, LDT, MAXDIF, MAXINT, NKDKO,
     1   NTE, NYNY, NDTF, NUMDT
      common / DIVASC / TN, XI, IOPST, KORDI, KQMAXD, KQMAXI, LDT,
     1   MAXDIF, MAXINT, NKDKO, NTE, NYNY, NDTF, NUMDT
c
      integer ICF,ICS,IGFLG,IGTYPE(2),IGSTOP(2),ILGREP,INGS,IOP3,IOP4,
     1   IOP5,IOP6,IOP7,IOP8,IOP9,IOP10,IOP11,IOP12,IOP13,IOP14,IOP15,
     2   IOP16,IOP17,IOP18,IOP19,IOP20,IOP21,IOP22,IOP21S,ITOLEP,IY,
     3   KEMAX,KIS,KMARK,KORD1I,KORD2I,KPRED,KQDCON,KQICON,KQMAXS,
     4   KQMXDS,KQMXIL,KQMXIP,KQMXIS,KSC,KSOUT,KSSTRT,KSTEP,LEX,LINC,
     5   LINCD,LINCQ,LSC,MAXKQD,MAXKQI,METHOD,NE,NEPTOL,NG,NGTOT,
     6   NOISEQ,NOUTKO,NTOLF,NY,IDAT(6)
      common /DIVAMC/ TG,TGSTOP,TMARK,TMARKX,TOUT,TOLG,HC,HDEC,HINC,
     1   HINCC,HMAX,HMAXP9,HMIN,ALPHA,BETA,D,G,V,DS,GS,SIGMA,RBQ,DNOISE,
     2   EAVE,EIMAX,EIMIN,EMAX,EREP,ROBND,SNOISE,FDAT,ICF,ICS,IGFLG,
     3   IGTYPE,IGSTOP,ILGREP,INGS,IOP3,IOP4,IOP5,IOP6,IOP7,IOP8,IOP9,
     4   IOP10,IOP11,IOP12,IOP13,IOP14,IOP15,IOP16,IOP17,IOP18,IOP19,
     5   IOP20,IOP21,IOP22,IOP21S,ITOLEP,IY,KEMAX,KIS,KMARK,KORD1I,
     6   KORD2I,KPRED,KQDCON,KQICON,KQMAXS,KQMXDS,KQMXIL,KQMXIP,KQMXIS,
     7   KSC,KSOUT,KSSTRT,KSTEP,LEX,LINC,LINCD,LINCQ,LSC,MAXKQD,MAXKQI,
     8   METHOD,NE,NEPTOL,NG,NGTOT,NOISEQ,NOUTKO,NTOLF,NY,IDAT
      save / DIVAMC / , / DIVASC /
c
      ID(1) = KEMAX
      ID(2) = KSTEP
      ID(3) = NUMDT
      RD(1) = EMAX
      return
      end
c   End of DIVACO

      subroutine DIVACR(Y, F, KORD, TOL, LGROUP)
c>> 1988-08-25 DIVACR Krogh   Fix bug in relative error test.
c>> 1988-01-15 DIVACR Krogh   Initial code.
c
c THIS SUBROUTINE
c   1. CORRECTS Y FOR EQUATIONS WHICH ARE NOT STIFF
c   2. ESTIMATES ERRORS
c   3. SELECTS INTEGRATION ORDERS
c   4. TESTS IF NOISE LIMITS THE PRECISION
c
c     Y = VECTOR OF PREDICTED VALUES ON ENTRY, AND OF CORRECTED
c         VALUES WHEN THE RETURN IS MADE.
c LGROUP= VECTOR INDICATING HOW ERROR TOLERANCES ARE TO BE GROUPED
c         (AND POSSIBLY HOW INTEGRATION ORDERS ARE TO BE GROUPED).
c   TOL = VECTOR CONTAINING ERROR TOLERANCES (AND POSSIBLY RELATIVE
c         ERROR FACTORS).
c     F = VECTOR GIVING PREDICTED DERIVATIVE VALUES AND DIFF. TABLES.
c    KD = VECTOR GIVING ORDERS OF THE DIFFERENTIAL EQUATIONS
c         (IF EQUATIONS HAVE DIFFERENT ORDERS).
c    KQ = VECTOR OF INTEGRATION ORDERS.
c
      integer LGROUP(*), KORD(*)
c--D Next line special: P=>D, X=>Q
      double precision Y(*)
      double precision TOL(*), F(*)
c
      integer KDIM, MAXORD, MAXSTF
c++ Substitute for KDIM, MAXORD, MAXSTF below
      parameter (KDIM = 20, MAXORD = 2, MAXSTF = 1)
c--D Next line special: P=>D, X=>Q
      double precision TN
      double precision XI(KDIM)
c
c--D Next line special: P=>D, X=>Q
      double precision TG(2), TGSTOP(2), TMARK, TMARKX, TOUT, TOLG
      double precision ALPHA(KDIM), BETA(KDIM+1)
      double precision  D(MAXSTF+MAXORD,MAXORD), G(KDIM,MAXORD)
      double precision V(KDIM+MAXORD)
      double precision HC, HDEC, HINC, HINCC, HMAX, HMAXP9, HMIN
      double precision FDAT(11)
c
      double precision DS(MAXSTF+MAXORD, MAXORD), GS(KDIM)
      double precision SIGMA(KDIM), RBQ(KDIM), DNOISE
      double precision EAVE, EIMAX, EIMIN, EMAX, EREP, ROBND, SNOISE
c
      integer IOPST, KORDI, KQMAXD, KQMAXI, LDT, MAXDIF, MAXINT, NKDKO,
     1   NTE, NYNY, NDTF, NUMDT
      common / DIVASC / TN, XI, IOPST, KORDI, KQMAXD, KQMAXI, LDT,
     1   MAXDIF, MAXINT, NKDKO, NTE, NYNY, NDTF, NUMDT
c
      integer ICF,ICS,IGFLG,IGTYPE(2),IGSTOP(2),ILGREP,INGS,IOP3,IOP4,
     1   IOP5,IOP6,IOP7,IOP8,IOP9,IOP10,IOP11,IOP12,IOP13,IOP14,IOP15,
     2   IOP16,IOP17,IOP18,IOP19,IOP20,IOP21,IOP22,IOP21S,ITOLEP,IY,
     3   KEMAX,KIS,KMARK,KORD1I,KORD2I,KPRED,KQDCON,KQICON,KQMAXS,
     4   KQMXDS,KQMXIL,KQMXIP,KQMXIS,KSC,KSOUT,KSSTRT,KSTEP,LEX,LINC,
     5   LINCD,LINCQ,LSC,MAXKQD,MAXKQI,METHOD,NE,NEPTOL,NG,NGTOT,
     6   NOISEQ,NOUTKO,NTOLF,NY,IDAT(6)
      common /DIVAMC/ TG,TGSTOP,TMARK,TMARKX,TOUT,TOLG,HC,HDEC,HINC,
     1   HINCC,HMAX,HMAXP9,HMIN,ALPHA,BETA,D,G,V,DS,GS,SIGMA,RBQ,DNOISE,
     2   EAVE,EIMAX,EIMIN,EMAX,EREP,ROBND,SNOISE,FDAT,ICF,ICS,IGFLG,
     3   IGTYPE,IGSTOP,ILGREP,INGS,IOP3,IOP4,IOP5,IOP6,IOP7,IOP8,IOP9,
     4   IOP10,IOP11,IOP12,IOP13,IOP14,IOP15,IOP16,IOP17,IOP18,IOP19,
     5   IOP20,IOP21,IOP22,IOP21S,ITOLEP,IY,KEMAX,KIS,KMARK,KORD1I,
     6   KORD2I,KPRED,KQDCON,KQICON,KQMAXS,KQMXDS,KQMXIL,KQMXIP,KQMXIS,
     7   KSC,KSOUT,KSSTRT,KSTEP,LEX,LINC,LINCD,LINCQ,LSC,MAXKQD,MAXKQI,
     8   METHOD,NE,NEPTOL,NG,NGTOT,NOISEQ,NOUTKO,NTOLF,NY,IDAT
      save / DIVAMC / , / DIVASC /
c.    SPECIFICATION OF ENVIRONMENTAL CONSTANTS.
      double precision EEPS10, EEPS16, EROV10, EEPS2, EEPT75, EOVEP2
      double precision OVTM75, OVD10
      common / DIVAEV / EEPS2, EEPT75, EOVEP2, OVTM75, OVD10, EEPS10,
     1   EEPS16, EROV10
      save / DIVAEV /
c
      integer L, I, KQL, KQN, KQD, JLGREP, J, K, ILGROR, ITOLOR, JLGROR,
     1   IORD, KOUTKO, KQLORD, LL, LKQMAX
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
c++  Code for INTEGO is active
      double precision TEMPA(4), TEMPAO(4)
c++  End
      save KOUTKO, LKQMAX
      equivalence (TPS1,TEMPA(1)), (TPS2,TEMPA(2)), (TPS3,TEMPA(3)),
     1   (TPS4, TEMPA(4))
      equivalence (G(1, 1), HH)
      integer MACT1(2), MACT2(12)
c             Parameters for Interface to MESS and DMESS
      integer MERET, METEXT, METABL
      parameter (MERET  =51)
      parameter (METEXT =53)
      parameter (METABL =55)
c ********* Error message text ***************
c[Last 2 letters of Param. name]  [Text generating message.]
cAA KSTEP=$(I6) T=$(E15.8) H=$(E12.5) LSC=$(I3) $C
c   EIMIN=$(E8.2) EAVE=$G KSC=$(I2) SIGMA($J)=$G $C
c   RQ=$(E11.5)$G$E
c   $
cAB I$HKQ$HLI$HE$HEI$HEPS$HF$H$H$H$H
c   HIGH ORDER PREDICTED DIFFERENCES$HRNOISE$HSTIFF$HBETA$E
      integer LTXTAA,LTXTAB
      parameter (LTXTAA=  1,LTXTAB=  1)
      character MTXTAA(1) * (104)
      character MTXTAB(1) * (88)
      data MTXTAA/'KSTEP=$(I6) T=$(E15.8) H=$(E12.5) LSC=$(I3) EIMIN=$(E
     *8.2) EAVE=$G KSC=$(I2) SIGMA($J)=$G RQ=$(E11.5)$G$E'/
      data MTXTAB/'I$HKQ$HLI$HE$HEI$HEPS$HF$H$H$H$HHIGH ORDER PREDICTED$
     * DIFFERENCES$HRNOISE$HSTIFF$HBETA$E'/
c
      data MACT1 / METEXT, MERET /
c (rr=repeat, t=3/5 for I/E format)  wwddtrr  wwddtrr  wwddtrr
      data MACT2 / METABL, 1, 0, 14, 0400201, 0300202, 0801503,
     1    1507501, 1103504, 0901501, 1002501, 1205501 /
c         wwddtrr  wwddtrr  wwddtrr  wwddtrr  wwddtrr
c          End of stuff for interface to message processor
c
      data REF(1), REF(2), REF(3), REF(4) / C1, CP9375, CP75, CP5 /
c++ Save data by elements if ~.C.
c++ Of next 20 lines, only the first KDIM-1 are active
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
C     data EIBND(20) / C1 /
c
c++  Code for ARGM is inactive
c      RETURN
c      ENTRY DIVACE
c++  End
c ********
c START OF CODE
c ********
      L = NDTF - 1
      if (ICS .ne. 1) L = L + (ICS - 1) * NUMDT
      do 3340 I = ICS, ICF
         if (NKDKO .ne. 0) KORDI = KORD(NKDKO + I - 1)
         IY = IY + abs(KORDI)
         KQL = KORD(I + 3)
         KQN = abs(KQL)
         KQD = max(2, KQN)
c ********
c OBTAIN ERROR TOLERANCE SPECIFIED BY THE USER
c ********
         if (I .le. ILGREP) if (KQL) 2600, 3310, 2610
         ITOLEP = abs(ITOLEP) + 1
         EPS = TOL(ITOLEP)
         ILGREP = LGROUP(ITOLEP)
c   TEST IF SIMPLE ABSOLUTE ERROR TEST IS BEING USED
         if (ILGREP .gt. 0) go to 2580
         JLGREP = ILGREP
c     GET OLD RELATIVE ERROR FACTOR
         TPS6 = TOL(ITOLEP + 1)
         ILGREP = LGROUP(ITOLEP + 1)
         ITOLEP = -ITOLEP - 1
c
         if (JLGREP + 1) 2540, 2570, 2510
c   NO CHECK ON THE ERROR ESTIMATE IS TO BE MADE
 2510    if (EPS + C1) 2520, 2590, 2520
c   ERROR TOLERANCE IS SPECIFIED IMPROPERLY
 2520    KEMAX = I
         NEPTOL = ITOLEP
         LINC = -3
         EREP = -abs(EREP)
         return
c   COMPUTE NEW RELATIVE ERROR FACTOR
 2540    continue
         TPS1 = C0
         do 2550 J = I, ILGREP
            TPS1 = TPS1 + abs(F(J))
 2550       continue
         TPS1 = abs(HH) * TPS1 / dble(ILGREP - I + 1)
         if (LSC .le. 2) go to 2560
c     ON FIRST 3 STEPS INCREASE TPS6 WHEN COMPUTING REL. ERROR FACTOR
         TPS6 = max(C4 * TPS1, TPS6)
c     ON 1-ST TIME THROUGH THE FIRST STEP, REL. ERR. FAC. IS NOT STORED
         if (LSC .eq. 7) go to 2570
 2560    continue
         TPS6 = max(TPS1, TPS6)
c   STORE NEW RELATIVE ERROR FACTOR
         TOL(-ITOLEP) = TPS6 * REF(-JLGREP - 1)
c   COMPUTE ABSOLUTE ERROR TOLERANCE
 2570    EPS = EPS * TPS6
 2580    if (EPS .le. C0) go to 2520
 2590    if (KQL) 2600, 3330, 2610
c END OF OBTAINING ERROR TOLERANCE
c ********
c OBTAIN INFORMATION USED FOR ERROR ESTIMATION, ORDER SELECTION, ETC.
c ********
c EQUATION IS STIFF
 2600    continue
c++  Code for STIFF is inactive
c      JS=abs(KORD(NJSKO+I-1))-1
c      JSI=JS
c      TPP=C0
c      TPS4=F(L+KQD+2)
c      TPS3=F(L+KQD+1)
c      TPS2=F(L+KQD)
c      TPS1=F(L+KQD-1)
c      IF (KQD.EQ.2) TPS1=Y(IY-1)
c      E=ABS(TPS3)+ABS(TPS4)
c      EI=E+ABS(TPS2)
c      RND=EI
c      IF (KORDI.GE.0) GO TO 2604
cc.    EQUATION IS IMPLICIT
c      JSI=JSI-1
c      IF (JSI.NE.0) GO TO 2604
c      IF (KORDI.EQ.-1) GO TO 2602
c      ERCOEF=GS(KQN+1)
c      GO TO 2606
c 2602 ERCOEF=.5D0*DS(KQD,1)
c      JSI=1
c      GO TO 2606
cc.    END OF SPECIAL CODE FOR IMPLICIT EQUATIONS
c 2604 ERCOEF = DS(KQD,JSI)
c 2606 ERCOEF = ABS(ERCOEF) / EPS
c      IF (LSC.LE.2)  GO TO 2710
c      IF (LSC-5) 2650,2710,2710
cc.  END OF CODE FOR STIFF EQUATIONS
c++  End
c
c EQUATION IS NOT STIFF
 2610    TPP = F(I) - F(L + 1)
         TPS3 = TPP
         TPS4 = TPP - F(L + KQD + 1)
         TPS2 = TPP + F(L + KQD)
         TPS1 = TPP + F(L + KQD - 1)
         E = abs(TPS3) + abs(TPS4)
         RND = E
         EI = E + abs(TPS2)
         ERCOEF = abs(GS(KQN + 1)) / EPS
         if (KQL .ge. 4) go to 2710
c   TEST IF STARTING OR IF INTEGRATION ORDER IS ONE
         if (LSC .le. 2) if (KQL - 2) 2660, 2710, 2710
c ********
c LOGIC ASSOCIATED WITH STARTING THE INTEGRATION
c ********
         TPS4 = C0
         if (LSC - 4) 2650, 2640, 2620
c FIRST STEP
 2620    E = E * CQ3125
         TPS3 = C0
         F(L + 4) = C0
         S = C0
c   TEST IF FIRST TIME THROUGH THE FIRST STEP
         if (LSC .eq. 7) go to 2690
c   COMPUTE S=ESTIMATE OF H * EIGENVALUE OF JACOBIAN = 2*(F(A)-F(B))/
c   (F(B)-F(C)) WHERE F(A)=CURRENT F(I), AND F(B) AND F(C) PRECEDING
c   VALUES OR ESTIMATES OF F(I)
         TPP = F(I) - F(L + 5)
         TPS4 = TPP
         E = C2 * abs(TPS4)
         if (S .ne. C0) S = (TPS4 + TPS4) / S
         if (S + CP125) 2630, 2700, 2700
c     SET LDT=-5  TO INDICATE POSSIBLE PROBLEMS DUE TO INSTABILITY
 2630    LDT = -5
         go to 2690
c   ADJUST CORRECTION MADE ON SECOND STEP
 2640    TPP = CP8 * TPP
c   ADJUST ESTIMATED ERRORS ON SECOND AND THIRD STEPS
 2650    E = abs(TPS3)
         RND = C4 * E
         go to 2710
c END OF SPECIAL LOGIC FOR STARTING
c ********
c INTEGRATION ORDER =1 IS TREATED AS A SPECIAL CASE
c ********
 2660    TPP = TPP + F(L + 2)
         if (BETA(2) .ge. C1P4) EI = EI * C1000
c   ESTIMATE NEW VALUE FOR S
         S = F(L + 4)
         if (S .eq. C0) go to 2680
         S = max(CM8, C2 * BETA(2) * (TPS1 - TPS2 - F(L + 5)) / S)
         if (S .ge. CMP5) go to 2670
c   MODIFY TPP (TO GET BETTER STABILITY CHARACTERISTICS)
         TPP = TPP * max(CP25, (CM2 - C2 * S) / (S * S))
 2670    TPS4 = TPS4 * abs(S)
 2680    E = CP25 * (E + abs(TPS4))
         EI = EI + abs(TPS4 * S)
c     STORE INFORMATION REQUIRED TO ESTIMATE S ON NEXT STEP
 2690    F(L + 4) = TPP
 2700    F(L + 5) = F(I)
c END OF SPECIAL CODE FOR INTEGRATION ORDER =1
c ********
c CODE FOR NOISE TEST AND GETTING ERROR ESTIMATE
c ********
 2710    E = E * ERCOEF
         RNOISE = C0
         if (EPS .lt. C0) go to 2810
         TPS5 = abs(F(L + 2)) + abs(F(I))
         if (TPS5 .eq. C0) go to 2760
 2720    RNOISE = RND / TPS5
         if (RNOISE .gt. RBQ(KQD)) if (RNOISE - C1) 2760, 2750, 2750
c   NOISE IS APPARENTLY SLOWING CONVERGENCE OF THE DIFFERENCES
c     REDUCE EI
         EI = RND
         TPS5 = abs(EEPS2 * Y(IY - 1)) / EPS
         if (TPS5 .lt. abs(E)) if (LSC) 2730, 2730, 2760
         E = TPS5
         RNOISE = C0
 2730    E = -abs(E)
         if (EIMIN .gt. CP1) EI = (C10 * EIMIN) * EI
c     COMPUTE REDUCTION TO BE MADE IN EI
         if (RNOISE .gt. (C20 * RBQ(KQD))) go to 2760
         K = -6 - LSC
 2740    if (K .le. 0) go to 2760
c     REDUCE EI WHEN NOISE APPARENTLY LIMITS PRECISION
         K = K - 1
         EI = CP5 * EI
         if (EI .gt. EIMIN) go to 2740
         go to 2760
 2750    TPS4 = 1.1D0 * RND
         TPS3 = RND
 2760    continue
c   TEST FOR STIFFNESS GOES HERE WHEN IMPLEMENTED
c *       INGREDIENTS OF TEST MAY INCLUDE --
c *       RNOISE, WHETHER (ABS(TPS4).GT.ABS(TPS3)),
c *       WHETHER EMAX IS INCREASING, RESULT OF TEST ON
c *       PREVIOUS STEPS, ETC.
c
c ********
c COMPUTE ERROR ESTIMATES AND INFORMATION FOR SELECTING THE STEPSIZE
c ********
         if (E .ge. abs(EMAX)) go to 2770
         if (-E .le. abs(EMAX)) go to 2780
         SNOISE = RNOISE
         DNOISE = RND
         NOISEQ = KQD
c   STORE PARAMETERS ASSOCIATED WITH LARGEST VALUE OF E
 2770    EMAX = E
         KEMAX = I
         NEPTOL = ITOLEP
c   DETERMINE HOW MUCH STEPSIZE CAN BE INCREASED
 2780    EI = EI * ERCOEF * SIGMA(KQD)
         EIMAX = max(EIMAX, EI)
         if (LINC .le. 0) go to 2810
         K = 0
 2790    if (EI .ge. min(EIMIN, EIBND(KQN))) go to 2800
         K = K + 1
         if (K .eq. LINC) go to 2810
         EI = EI * SIGMA(KQD)
         go to 2790
 2800    LINC = K
c END OF COMPUTING ERROR ESTIMATES
 2810    continue
c++  Code for ERRSTO is inactive
c      IF (IOP20 .EQ. 0) GO TO 780
cc.********
cc.STORE ERROR ESTIMATE (OPTIONAL)
cc.********
c      F(IOP20+I-1)=TPS3*GS(KQN+1)
cc.END OF STORING ERROR ESTIMATE
c++  Code for INTEGO | ERRSTO is active
         if (IOP19 .eq. 0) go to 3090
c.********
c.EQUATIONS ARE GROUPED TO USE SAME INTEGRATION METHOD (OPTIONAL)
c.********
c++  Code for INTEGO is active
         if (I .gt. 1) if (I - ILGROR) 2900, 2900, 2830
         ITOLOR = IOP19
 2830    JLGROR = KORD(ITOLOR)
         ITOLOR = ITOLOR + 1
         if (JLGROR .gt. 0) go to 2870
         ILGROR = KORD(ITOLOR)
         ITOLOR = ITOLOR + 1
         if (JLGROR + 1) 2840, 2850, 2890
 2840    if (JLGROR .lt. -2) if (KQD + JLGROR) 2850, 2880, 2880
c.INITIALIZE FOR ACCUMULATING VARIABLES USED IN ORDER SELECTION
 2850    IORD = I
         KQLORD = KQL
         do 2860 K = 1, 4
 2860       TEMPAO(K) = abs(TEMPA(K))
         go to 2930
c.ORDERS IN CURRENT GROUP CAN BE DIFFERENT
 2870    ILGROR = JLGROR
         go to 3090
c.ORDER IS NOT GOING TO BE CHANGED
 2880    JLGROR = 0
 2890    if (KQL) 3240, 3270, 3270
c.TAKE ACTION FOR EQUATION WHICH IS NOT THE FIRST IN THE GROUP
 2900    if (JLGROR) 2910, 2890, 3090
c.ACCUMULATE VARIABLES USED IN ORDER SELECTION
 2910    do 2920 K = 1, 4
 2920       TEMPAO(K) = TEMPAO(K) + abs(TEMPA(K))
c.    TEST IF THIS IS LAST EQUATION IN THE GROUP
 2930    if (I .ne. ILGROR) if (KQL) 3310, 3290, 3290
c.SET UP TO GO SELECT INTEGRATION ORDER
         KQL = 0
         do 2940 K = 1, 4
 2940       TEMPA(K) = TEMPAO(K)
         go to 3090
c.INTEGRATION ORDER HAS BEEN SELECTED
c++  Code for INTEGO | STIFF is active
 2950    continue
c++  Code for INTEGO is active
         KQL = KQLORD
         if (KQN - abs(KQL)) 2960, 2980, 3020
c.  TEST IF ORDER CAN BE DECREASED
 2960    if (JLGROR .ge. -2) if (KQL) 3010, 3040, 3040
c.    INTEGRATION ORDER WAS SELECTED OUTSIDE PERMITTED RANGE
 2970    KQN = abs(KQL)
c.    INTEGRATION ORDER IS NOT GOING TO BE CHANGED
 2980    if ((KQL .ne. 1) .or. (LSC .gt. 0)) if (KQL) 3030, 3040, 3040
c.    SET  4-TH ENTRY IN DIFFERENCE TABLES SO THAT STANDARD ADAMS
c.    METHOD IS USED WHEN KQL=1
 2990    do 3000 K = IORD, I
 3000       F(NDTF + K*NUMDT - NUMDT + 3) = C0
         go to 3270
c.  ORDER FOR STIFF EQUATION WAS REDUCED
 3010    continue
c++  Code for INTEGO & STIFF is inactive
c      IF (KQN.LT.JSI) GO TO 990
c      TPP=-C1
c      GO TO 1090
cc.  TEST IF ORDER CAN BE INCREASED
c++  Code for INTEGO is active
 3020    if (JLGROR .eq. -2) go to 2970
c++  Code for INTEGO & STIFF is inactive
c      IF (KQL.GE.0) GO TO 1140
c      IF ((JSI.NE.0).AND.(KQN.GT.(MAXKQD+JSI))) GO TO 990
c      TPP=C1
cc.  STORE RESULTS FOR STIFF EQUATIONS
c++  Code for INTEGO is active
 3030    continue
c++  Code for INTEGO & STIFF is inactive
c      DO 3035 K=IORD,I
c      KORD(K+3) = -KQN
c 3035 F(NDTF+K*NUMDT-NUMDT)=TPP
c      GO TO 3245
cc.  STORE RESULTS FOR EQUATIONS WHICH ARE NOT STIFF
c++  Code for INTEGO is active
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
c++  End
c.********
c.SELECT INTEGRATION ORDER
c.********
 3090    if (LSC .le. 0) go to 3120
c. SPECIAL ORDER SELECTION WHEN STARTING
         if (LSC - 3) 3110, 3210, 3100
 3100    if (LSC .eq. 5) if (S + .125D0) 3160, 3130, 3130
         if (LSC - 6) 3130, 3130, 3210
 3110    if (C40 * min(abs(TPS4), abs(TPS3)) .gt. abs(TPS2)) then
            if (EPS .ne. -C1) LSC = 2
         end if
         if (abs(TPS4) .lt. abs(TPS3)) if (C4 * abs(TPS4) - abs(TPS2))
     1      3130, 3130, 3210
c.  CHECK IF ORDER CAN BE INCREASED OR SHOULD BE DECREASED
 3120    TPS5 = ROBND * abs(TPS4)
         TPS6 = ROBND * (TPS5 + abs(TPS3))
         TPS7 = abs(TPS1) + abs(TPS2)
         if (TPS5 .ge. abs(TPS3)) go to 3140
         if (TPS6 .ge. TPS7) go to 3210
 3130    if (KQN .ge. MAXKQI) go to 3210
c.    INCREASE THE INTEGRATION ORDER
         KQN = KQN + 1
c++  Code for INTEGO | STIFF is active
         if (KQL) 3230, 2950, 3250
c++  Code for ~(INTEGO | STIFF) is inactive
c      GO TO 3250
c++  End
c.  CHECK IF ORDER SHOULD BE DECREASED
 3140    if (TPS6 .lt. TPS7) go to 3210
         if (TPS5 .lt. abs(TPS3 - TPS4)) go to 3210
         if ((TPS3.eq.TPS4) .and. (LSC.le.0)) go to 3210
         if (KQN - 2) 3210, 3160, 3180
 3160    KQN = 1
c++  Code for INTEGO | STIFF is active
         if (KQL) 3220, 2950, 3170
c++  End
c.    WHEN ORDER IS REDUCED TO 1 WITH ADAMS METHOD SET F(L+4)=0
 3170    F(L + 4) = C0
         go to 3260
c.    DECREASE THE INTEGRATION ORDER
 3180    KQN = KQN - 1
c++  Code for INTEGO | STIFF is active
         if (KQL) 3220, 2950, 3200
c++  End
 3200    F(L+KQD) = F(L+KQD) + TPP
         go to 3260
c   NO CHANGE IN INTEGRATION ORDER IS BEING MADE
 3210    continue
c++  Code for INTEGO is active
         if (KQL) 3240, 2950, 3270
c++  Code for ~INTEGO is inactive
c         TPS1 = EEPS10
c      GO TO 1530
c++  End
c END OF SELECTING INTEGRATION ORDER
c ********
c COMPUTE MAXIMUM INTEGRATION ORDERS AND SET NEW ONES (IF ANY)
c ********
c EQUATION IS STIFF
c     ORDER WAS DECREASED
c++  Code for INTEGO | STIFF is active
 3220    continue
c++  Code for STIFF is inactive
c      IF (KQN.LT.JSI) GO TO 3236
c      F(L+1)=-C1
c      GO TO 3233
cc.    ORDER WAS INCREASED
c++  Code for INTEGO |  STIFF  is active
 3230    continue
c++  Code for STIFF is inactive
c      IF ((JSI.NE.0).AND.(KQN.GT.(MAXKQD+JSI))) GO TO 3236
c      F(L+1)=C1
c 3233 KORD(I+3) = -KQN
c      GO TO 3245
c      ORDER WAS SET TO AN UNACCEPTABLE VALUE
c 3236 KQN=abs(KQL)
c      ORDER IS NOT BEING CHANGED
c++  Code for STIFF |  INTEGO is active
 3240    continue
c++  Code for STIFF is inactive
c      F(L+1)=C0
c 3245 IF (JSI.NE.0) KQMAXD=max(KQN,KQMAXD)
c      IF (JS.LT.abs(KORDI)) KQMAXI=max(KQN,KQMAXI)
c      GO TO 3290
c++  End
c EQUATION IS NOT STIFF
c     ORDER INCREASED
 3250    F(L + KQN + 1) = -F(L + KQD + 1)
         if (LSC .gt. 0) F(L + KQN + 1) = F(L + 1) - F(I)
c     ORDER CHANGED
 3260    KORD(I + 3) = KQN
 3270    KQMAXI = max(KQN, KQMAXI)
         if (EPS .gt. C0) KQMAXS = max(KQN, KQMAXS)
         F(L + KQD + 1) = C0
 3290    continue
         if (KQN .gt. KIS) go to 3310
c.********
c.DETERMINE IF TIME TO STORE SOLUTION (OPTIONAL)
c.********
         if (KIS .ge. 1000) then
            TP2 = max(1.5D0, dble(KQN) * C2 ** (1001 - KIS)) * abs(TPS4)
 3295       if (TP2 .gt. abs(F(L+KQN))) then
               if (KQN .le. KQL) then
                  KQN = KQN - 1
                  if (KQN .gt. 1) go to 3295
                  KQN = 1
               end if
            end if
            KORD(I+3) = KQN
            if (I .eq. 1) LKQMAX = 0
            LKQMAX = max(KQN, LKQMAX)
            KQMAXI = LKQMAX
            if (KIS .eq. 1000) then
               if (I .eq. KEMAX) EMAX = dble(8 + KQN**2) * abs(EMAX)
               go to 3325
            end if
c++  Code for DUMP is active
         else if ((E .ne. C0) .and. (EPS .gt. C0)) then
            if (IOP9 .gt. 0) if((abs(E)*dble(KIS-KQN+2)**(KQN+1))-1.D-2)
     1         3310, 3310, 3300
 3300       KIS = -1
c++  End
         end if
 3310    continue
c ********
c CORRECT
c ********
         do 3320 K = 1, KORDI
c++  Code for ~{p,x} is active
            Y(IY - K) = Y(IY - K) + G(KQL + 1, K) * TPP
c++  Code for {p,x} is inactive
Cc--D Next line special: P=>D, X=>Q
C            Y(IY - K) = Y(IY - K) + dble(G(KQL + 1, K)) * dble(TPP)
c++  END
 3320    continue
c END OF CORRECTING
 3325 continue
c++  Code for OUTPUT is active
      if (IOP10 .gt. 0) then
         if (I .eq. 1) then
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
c--D Next line special: P=>S, X=>D
            call DMESS(MACT1, MTXTAA, IDAT, FDAT)
            KOUTKO = NOUTKO
         end if
         if (KOUTKO .ne. 0) then
            if (KORD(KOUTKO) .gt. 0) then
               if (I .lt. KORD(KOUTKO)) go to 3328
               KOUTKO = KOUTKO + 1
            else
               if (I .ge. abs(KORD(KOUTKO))) KOUTKO = KOUTKO + 1
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
         if (KQL .eq. 1) FDAT(10) = S
         FDAT(11) = BETA(KQD)
c--D Next line special: P=>S, X=>D
         call DMESS(MACT2, MTXTAB, IDAT, FDAT)
 3328    if (I .eq. NTE) IOP10 = IOP10 - 1
      end if
c++  End
 3330    L = L + NUMDT
 3340    continue
      return
      end
c   End of DIVACR

      subroutine DIVAHC
c>> 1988-05-20 DIVAHC Krogh   Initial code.
c
c SUBROUTINE TO COMPUTE COEFFICIENTS REQUIRED FOR INTEGRATING
c ORDINARY DIFFERENTIAL EQUATIONS
c
      integer KDIM, MAXORD, MAXSTF
c++ Substitute for KDIM, MAXORD, MAXSTF below
      parameter (KDIM = 20, MAXORD = 2, MAXSTF = 1)
c--D Next line special: P=>D, X=>Q
      double precision TN
      double precision XI(KDIM)
c
c--D Next line special: P=>D, X=>Q
      double precision TG(2), TGSTOP(2), TMARK, TMARKX, TOUT, TOLG
      double precision ALPHA(KDIM), BETA(KDIM+1)
      double precision  D(MAXSTF+MAXORD,MAXORD), G(KDIM,MAXORD)
      double precision V(KDIM+MAXORD)
      double precision HC, HDEC, HINC, HINCC, HMAX, HMAXP9, HMIN
      double precision FDAT(11)
c
      double precision DS(MAXSTF+MAXORD, MAXORD), GS(KDIM)
      double precision SIGMA(KDIM), RBQ(KDIM), DNOISE
      double precision EAVE, EIMAX, EIMIN, EMAX, EREP, ROBND, SNOISE
c
      integer IOPST, KORDI, KQMAXD, KQMAXI, LDT, MAXDIF, MAXINT, NKDKO,
     1   NTE, NYNY, NDTF, NUMDT
      common / DIVASC / TN, XI, IOPST, KORDI, KQMAXD, KQMAXI, LDT,
     1   MAXDIF, MAXINT, NKDKO, NTE, NYNY, NDTF, NUMDT
c
      integer ICF,ICS,IGFLG,IGTYPE(2),IGSTOP(2),ILGREP,INGS,IOP3,IOP4,
     1   IOP5,IOP6,IOP7,IOP8,IOP9,IOP10,IOP11,IOP12,IOP13,IOP14,IOP15,
     2   IOP16,IOP17,IOP18,IOP19,IOP20,IOP21,IOP22,IOP21S,ITOLEP,IY,
     3   KEMAX,KIS,KMARK,KORD1I,KORD2I,KPRED,KQDCON,KQICON,KQMAXS,
     4   KQMXDS,KQMXIL,KQMXIP,KQMXIS,KSC,KSOUT,KSSTRT,KSTEP,LEX,LINC,
     5   LINCD,LINCQ,LSC,MAXKQD,MAXKQI,METHOD,NE,NEPTOL,NG,NGTOT,
     6   NOISEQ,NOUTKO,NTOLF,NY,IDAT(6)
      common /DIVAMC/ TG,TGSTOP,TMARK,TMARKX,TOUT,TOLG,HC,HDEC,HINC,
     1   HINCC,HMAX,HMAXP9,HMIN,ALPHA,BETA,D,G,V,DS,GS,SIGMA,RBQ,DNOISE,
     2   EAVE,EIMAX,EIMIN,EMAX,EREP,ROBND,SNOISE,FDAT,ICF,ICS,IGFLG,
     3   IGTYPE,IGSTOP,ILGREP,INGS,IOP3,IOP4,IOP5,IOP6,IOP7,IOP8,IOP9,
     4   IOP10,IOP11,IOP12,IOP13,IOP14,IOP15,IOP16,IOP17,IOP18,IOP19,
     5   IOP20,IOP21,IOP22,IOP21S,ITOLEP,IY,KEMAX,KIS,KMARK,KORD1I,
     6   KORD2I,KPRED,KQDCON,KQICON,KQMAXS,KQMXDS,KQMXIL,KQMXIP,KQMXIS,
     7   KSC,KSOUT,KSSTRT,KSTEP,LEX,LINC,LINCD,LINCQ,LSC,MAXKQD,MAXKQI,
     8   METHOD,NE,NEPTOL,NG,NGTOT,NOISEQ,NOUTKO,NTOLF,NY,IDAT
      save / DIVAMC / , / DIVASC /
c.    SPECIFICATION OF ENVIRONMENTAL CONSTANTS.
      double precision EEPS10, EEPS16, EROV10, EEPS2, EEPT75, EOVEP2
      double precision OVTM75, OVD10
      common / DIVAEV / EEPS2, EEPT75, EOVEP2, OVTM75, OVD10, EEPS10,
     1   EEPS16, EROV10
      save / DIVAEV /
c                 K - 1 + 1 / K  is equivalent to max(1, K-1)
      double precision GG(MAXORD - 1 + 1/MAXORD), B(KDIM+MAXORD),
     1   W(KDIM+MAXORD)
      integer  K, N, J
      double precision C0, CP1, CRBQI, CP5, CP5625, C1, C1P125
      parameter (C0 = 0.D0)
      parameter (CP1 = .1D0)
      parameter (CRBQI = .421875D0)
      parameter (CP5 = .5D0)
      parameter (CP5625 = .5625D0)
      parameter (C1 = 1.D0)
      parameter (C1P125 = 1.125D0)
c++  Code for STIFF is inactive
c      INTEGER          GODIF
c++  End
      double precision TP1, TP2, HH, TEMP, TP
      equivalence (G(1, 1), HH)
c
      save GG, W
c
c  B(K)= 1/(K*(K+1))
c++ Save data by elements if ~.C.
c++ Of next 23 lines, only the first KDIM+MAXORD are active
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
C     data B(23) / 1.811594202898550724637681159420289855072D-3 /
c
c ********
c START OF CODE
c ********
c     SET STEP NUMBER OF METHOD
c++  Code for STIFF is inactive
c      IOP11 = MIN(max(KQMAXI,KQMAXD) + 1), KDIM)
c++  Code for ~STIFF is active
      IOP11 = MIN(KQMAXI + 1, KDIM)
c++  End
c     TEST IF STEPSIZE WAS CHANGED
      if (KQICON .ge. 0) go to 3510
c ********
c STEPSIZE JUST CHANGED
c ********
c     SET CONSTANTS DEPENDING ON NEW STEPSIZE
      KQMXIL = KQMAXI
      TP1 = HH
      GG(1) = TP1 * TP1
      G(1, 2) = GG(1) * CP5
      if (MAXINT .le. 2) go to 3450
      do 3440 K = 3, MAXINT
         GG(K - 1) = G(1, K - 1) * TP1
         G(1, K) = GG(K - 1) / dble(K)
 3440    continue
c     SET CONSTANTS INDICATING STEP CHANGE
 3450 KQICON = 0
c++  Code for STIFF is inactive
c      KQDCON=0
c++  End
      KQMXIP = 1
      KSC = 1
      if (LSC .lt. 7) go to 3490
c     SPECIAL SET-UP OF CONSTANTS ON THE VERY FIRST STEP
      HINCC = C1P125
      LINCD = 6
      LINCQ = 12
      if (HINC .gt. C0) go to 3460
      LINCD = -2
      LINC = -2
      ROBND = C1
 3460 SIGMA(1) = 1.0D0
      BETA(1) = C1
      do 3470 N = 1, IOP11
c++  Code for STIFF is inactive
c      D(1,N)=C0
c++  End
         XI(N) = TP1
         ALPHA(N) = C1
         BETA(N + 1) = C1
         SIGMA(N + 1) = dble(N + 1) * SIGMA(N) * HINCC
 3470    continue
      TEMP = EEPS16
      RBQ(1) = C1
      RBQ(2) = CP1
      TP = CRBQI
c     **** IN THE LOOP BELOW RBQ(K) IS COMPUTED TO BE
c          APPROXIMATELY (3/4 ** ((K-1) ** 2 - 1) / 10
c          .5625 = (3/4) ** 2    TP = (3/4) ** (2*K -3)
      do 3480 K = 3, KDIM
         TEMP = TEMP + TEMP
         RBQ(K) = max(TEMP, RBQ(K - 1) * TP)
 3480    TP = TP * CP5625
      go to 3560
c     SET-UP AFTER THE FIRST STEP
 3490 TP2 = XI(1)
      XI(1) = TP1
      BETA(2) = TP1 / TP2
      K = 2
      if (HINCC .eq. HINC) go to 3540
      if ((LSC .ne. 0) .or. ((KSTEP-KSSTRT-KQMAXS) .lt. 10)) go to 3540
      HINCC = C1
      LINCD = 0
 3500 LINCD = LINCD + 1
      HINCC = HINCC * HINC
      if (HINCC .lt. 2.D0) go to 3500
      LINC = (LINC * (LINCD + LINCD)) / LINCQ
      LINCQ = LINCD + LINCD
      HINCC = HINC
      go to 3540
c END OF LOGIC FOR CASE WHEN STEPSIZE JUST CHANGED
c     TEST IF MAXIMUM INTEGRATION ORDER DID NOT INCREASE
 3510 if (KQMAXI .gt. KQMXIL) then
c ********
c INTEGRATION ORDER WAS INCREASED -- GET NEW V'S
c ********
         KQMXIL = KQMAXI
         KQMXIP = KQMXIL + MAXINT
         K = KQMXIP
         V(K) = B(K)
         if (KQICON .eq. 1) go to 3530
c     if (KQICON .eq. K) KQICON = KQICON - 1 --- Removed 1999-08-19
         do 3520 N = 2, KQICON
            K = K - 1
 3520       V(K) = V(K) - ALPHA(N) * V(K + 1)
c END OF GETTING NEW V'S
      else
         IOP11 = max(IOP11, KQMXIL+1)
      end if
 3530 if (IOP11 .le. KSC) go to 3560
c ********
c COMPUTE PARAMETERS WHICH ARE STILL CHANGING AS A RESULT OF
c A CHANGE IN THE STEPSIZE
c ********
      TP2 = XI(KSC)
c     UPDATE CONSTANT STEP COUNTER
      KSC = KSC + 1
      K = KSC
      BETA(K) = C1
 3540 continue
      TEMP = HINCC
c
c   LOOP TO COMPUTE NEW VALUES OF PARAMETERS
      do 3550 N = K, IOP11
         TP1 = TP2 + HH
         TP2 = XI(N)
         XI(N) = TP1
         ALPHA(N) = HH / TP1
         BETA(N + 1) = BETA(N) * (TP1 / TP2)
         TEMP = max(TEMP, dble(N) * (ALPHA(N) * HINCC))
         SIGMA(N) = SIGMA(N - 1) * TEMP
 3550    continue
      if (IOP11 .ne. KDIM) XI(IOP11 + 1) = TP2 + HH
c END OF CODE FOR COMPUTING PARAMETERS WHICH ARE STILL CHANGING
c
 3560 if (KQICON .ge. KQMXIP) go to 3690
c ********
c COMPUTE INTEGRATION COEFFICIENTS WHICH ARE STILL CHANGING
c ********
      KQMXIL = max(KQMAXI, KQMXIL)
      KQMXIP = KQMXIL + MAXINT
      J = KQMXIP - KQICON
      N = KQICON + 1
      KQICON = N
      if (N .ne. 1) go to 3580
c INITIALIZE V AND W
      do 3570 K = 1, J
         V(K) = B(K)
 3570    W(K) = V(K)
      go to 3600
c UPDATE V AND INITIALIZE W
 3580 if (N .eq. KDIM) go to 3690
      do 3590 K = 1, J
         V(K) = V(K) - ALPHA(N) * V(K + 1)
 3590    W(K) = V(K)
c SET TRANSFER FOR LOOP BELOW DEPENDING ON VALUE OF MAXINT
 3600 continue
      go to 3660
c
 3640 J = J - 1
c INNER LOOP FOR COMPUTING INTEGRATION COEFFICIENTS
      do 3650 K = 1, J
 3650    W(K) = W(K) - ALPHA(N) * W(K + 1)
c     STORE INTEGRATION COEFFICIENTS
 3660 G(N + 1, 1) = HH * W(1)
      GS(N + 1) = G(N + 1, 1) - G(N, 1)
c++  Code for MAXORD >= 2 is active
      if (MAXINT .ge. 2) then
         G(N + 1, 2) = GG(1) * W(2)
c++  Code for MAXORD >= 3 is inactive
c        if (MAXINT .gt. 2) then
c           DO 3665 K=3,MAXINT
c3665          G(N+1,K)=GG(K-1)*W(K)
c        end if
c++  Code for MAXORD >= 2 is active
      end if
c++  End
      N = N + 1
      if (N .le. KQMXIL) go to 3640
c END OF COMPUTING INTEGRATION COEFFICIENTS
c
 3690 continue
c++  Code for STIFF is inactive
c      IF (KQDCON.GT.KQMAXD) GO TO 4662
cc.********
cc.COMPUTE DIFFERENTIATION COEFFICIENTS WHICH ARE STILL CHANGING
cc.********
cc.SET TRANSFER FOR LOOP BELOW, DEPENDING ON VALUE OF MAXDIF
c++  Code for STIFF & MAXORD >= 2 is inactive
c      IF (MAXDIF-2) 3693,3692,3691
c 3691 ASSIGN 3696 TO GODIF
c      GO TO 3694
c 3692 ASSIGN 3698 TO GODIF
c      GO TO 3694
c 3693 ASSIGN 3699 TO GODIF
c++  Code for STIFF is inactive
c 3694 KQDCON=KQDCON+1
cc.LOOP FOR COMPUTING DIFFERENTIATION COEFFICIENTS
c      DO 3699 N=KQDCON,KQMAXD
c      DS(N+1,2)=C1/XI(N)
c      D(N+1,1)=DS(N+1,2)+D(N,1)
c      DS(N+1,1)=DS(N+1,2)/D(N+1,1)
c++  Code for STIFF & MAXORD >= 2 is inactive
c      GO TO GODIF, (3696,3698,3699)
c 3696 CONTINUE
c++  Code for STIFF & MAXORD >= 3 is inactive
c      DO 3697 K=3,MAXDIF
c      DS(N+1,K)=D(N,K-2) * (K-1)/XI(N)
c 3697 D(N+1,K-1)=DS(N+1,K) + D(N,K-1)
c++  Code for STIFF is inactive
c 3698 CONTINUE
c++  Code for STIFF & MAXORD >= 2 is inactive
c      D(N+1,MAXDIF)=D(N,MAXDIF) + D(N,MAXDIF-1) * (MAXDIF)/XI(N)
c++  Code for STIFF is inactive
c 3699 CONTINUE
c++  End
c
c END OF COMPUTING DIFFERENTIATION COEFFICIENTS
      return
      end
c   End of DIVAHC

      subroutine DIVAIN(T, Y, F, KORD)
c>> 1988-01-14 DIVAIN Krogh   Initial code.
c
c  SUBROUTINE TO DO INTERPOLATION FOR VARIABLE ORDER INTEG. ROUTINE
c
      integer KORD(*)
c--D Next line special: P=>D, X=>Q
      double precision T(*), Y(*)
      double precision F(*)
      integer KDIM, MAXORD
c++ Substitute for KDIM, MAXORD below
      parameter (KDIM = 20, MAXORD = 2)
c--D Next line special: P=>D, X=>Q
      double precision TN
      double precision XI(KDIM)
c
      integer IOPST, KORDI, KQMAXD, KQMAXI, LDT, MAXDIF, MAXINT, NKDKO,
     1   NTE, NYNY, NDTF, NUMDT
      common / DIVASC / TN, XI, IOPST, KORDI, KQMAXD, KQMAXI, LDT,
     1   MAXDIF, MAXINT, NKDKO, NTE, NYNY, NDTF, NUMDT
      save / DIVASC /
      integer I, ICI, IDT, INTERP, INTEG, INTEGZ, IY, IYI, IYN, IYNI, J,
     1    K, KQMXI, KQMXS, KQQ, L, N
      double precision C0, C1, C2
      parameter (C0 = 0.D0)
      parameter (C1 = 1.D0)
      parameter (C2 = 2.D0)
      double precision C(KDIM+MAXORD-1), ETA(KDIM)
      double precision GAMMA(KDIM)
      double precision TP1, HI
      double precision CSUM(KDIM+MAXORD-1)
c--D Next line special: P=>D, X=>Q
      double precision XP1
      logical LNOTM1
c
c              Stuff for processing error messages
      integer IDAT(1)
      double precision FDAT(6)
      integer MENTXT, MERET, MEEMES, METEXT
      parameter (MENTXT =23)
      parameter (MERET  =51)
      parameter (MEEMES =52)
      parameter (METEXT =53)
      integer MACT(8)
c ********* Error message text ***************
c[Last 2 letters of Param. name]  [Text generating message.]
cAA DIVAIN$B
cAB Interpolating at T(1)=$F with $B
cAC TN=$F, T(2)=$F and H=$F.  T(1) must be in [$F, $F].$E
cAD internal variable LDT = $I.  Interpolation not allowed now.$E
      integer LTXTAA,LTXTAB,LTXTAC,LTXTAD
      parameter (LTXTAA=  1,LTXTAB=  9,LTXTAC= 41,LTXTAD= 94)
      character MTXTAA(1) * (154)
      data MTXTAA/'DIVAIN$BInterpolating at T(1)=$F with $BTN=$F, T(2)=$
     *F and H=$F.  T(1) must be in [$F, $F].$Einternal variable LDT = $I
     *.  Interpolation not allowed now.$E'/
c
c                      1 2 3 4       5 6       7      8
      data MACT / MEEMES,0,0,0, MENTXT,0, METEXT, MERET /
c
c++  Code for ARGM is inactive
c      ENTRY DIVAIE
c++  End
c ********
c START OF CODE -- CHECK ON STATE OF DIFFERENCE TABLE
c ********
      L = LDT
      if (L) 3710, 3730, 3780
 3710 if (L + 2) 4170, 3730, 3720
 3720 if (MAXINT .ge. 0) L = 1
      go to 3840
c ********
c UPDATE DIFFERENCE TABLE TO START OF NEXT STEP
c ********
 3730 K = NDTF
      do 3770 I = 1, NTE
         KQQ = KORD(I + 3)
         if (KQQ .le. 0) go to 3760
c EQUATION IS NOT STIFF
         TP1 = F(I) - F(K)
c LOOP TO DO UPDATING
         N = K + max(abs(KQQ), 2)
         do 3750 J = K, N
 3750       F(J) = F(J) + TP1
 3760    continue
         K = K + NUMDT
 3770    continue
      LDT = 1
      if (L .ne. 0) return
c END OF UPDATING DIFFERENCE TABLE
c ********
c INITIALIZE FOR COMPUTATION OF COEFFICIENTS
c ********
 3780 INTERP = 0
      HI = T(1) - TN
      GAMMA(1) = HI / XI(1)
      if (GAMMA(1)) 3790, 3800, 3810
 3790 if (GAMMA(1) .ge. -C1) go to 3820
      INTERP = 1
      if (abs(HI) - abs(T(2))) 3820, 3820, 4180
 3800 INTERP = 2 - KQMAXI
      go to 3820
 3810 if (GAMMA(1) .gt. C2) if (LDT - 2) 4180, 3820, 3820
 3820 KQMXI = KQMAXI + INTERP - 1
c++  Code for STIFF is inactive
c      KQMXS=max(KQMXI,KQMAXD)
c++  Code for ~STIFF is active
      KQMXS = KQMXI
c++  End
      do 3830 N = 2, KQMXS
 3830    GAMMA(N) = (HI + XI(N-1)) / XI(N)
 3840 LNOTM1 = L .ne. -1
      INTEG = MAXINT
      if (INTEG .le. 0) if (INTEG + MAXDIF) 4160, 3950, 3950
c ********
c COMPUTE INTEGRATION COEFFICIENTS
c ********
c     INITIAL SET-UP
c         COMPUTE INITIAL C VALUES
      do 3850 N = 1, INTEG
         C(N) = HI / dble(N)
 3850 continue
      I = INTEG + 1
      INTEG = INTEG + KQMXI
      do 3860 N = I, INTEG
         C(N) = C(N - 1) * (dble(N - MAXINT) / dble(N))
 3860    continue
c         COMPUTE ETA'S
      do 3870 N = 1, KQMXI
 3870    ETA(N) = HI / XI(N)
c         COMPUTE C(K)'S TO CORRESPOND TO G(K-MAXINT+1,MAXINT),
c         K=MAXINT, MAXINT+1,..., MAXINT+KQMXI-1
      I = INTEG
 3880 J = INTEG
      INTEG = J - 1
      if (INTEG .le. MAXINT) go to 3900
      do 3890 N = J, I
 3890    C(N) = ETA(N - INTEG) * C(N) + C(N - 1)
      go to 3880
 3900 do 3910 N = J, I
 3910    C(N) = ETA(N - INTEG) * C(N)
c         END OF COMPUTING  G(---,MAXINT)
      INTEGZ = 0
      go to 3940
c         COMPUTE C(K)-S TO CORRESPOND TO G(K-INTEG+1,INTEG),
c         K=INTEG+1,INTEG+2,..., INTEG+KQMXI
 3920 do 3930 N = 1, KQMXI
 3930    C(INTEG+N) = GAMMA(N)*C(INTEG+N-1) - ETA(N)*C(INTEG+N)
 3940 ICI = INTEG - 1
      go to 4020
c END OF COMPUTING INTEGRATION COEFFICIENTS
c ********
c COMPUTE COEFFICIENTS FOR INTERPOLATION
c ********
 3950 C(1) = C1
      ICI = 0
      do 3960 N = 1, KQMXS
 3960    C(N + 1) = GAMMA(N) * C(N)
      if (INTEG + 1) 3970, 3990, 4010
c END OF COMPUTING INTERPOLATION COEFFICIENTS
c
c     SET-UP TO COMPUTE DIFFERENTIATION COEFFICIENTS REQUIRED
c     IN ORDER TO GET COEFFICIENTS ACTUALLY USED
 3970 INTEG = 0
      ICI = 1
 3980 INTEG = INTEG - 1
      if (INTEG .eq. MAXINT) ICI = 0
c ********
c COMPUTE DIFFERENTIATION COEFFICIENTS
c ********
 3990 INTERP = max(INTERP, 0)
      TP1 = dble(-INTEG)
      C(1) = TP1 * C(1) / XI(-INTEG)
      J = KQMAXD + INTEG
      do 4000 N = 1, J
 4000    C(N + 1) = (TP1*C(N)) / XI(N - INTEG) + GAMMA(N - INTEG) * C(N)
c     C(N) NOW CORRESPONDS TO THE DIFFERENTIAL COEFFICIENT
c          D(N-INTEG,-INTEG)
 4010 INTEGZ = INTEG
      if (ICI .ne. 0) go to 3980
c END OF COMPUTING DIFFERENTIATION COEFFICIENTS
c ********
c BEGINNING OF LOOP TO DO
c         INTEGRATION       (INTEG.GT.0)
c         INTERPOLATION     (INTEG.EQ.0)
c         DIFFERENTIATION   (INTEG.LT.0)
c TO THE POINT INDICATED BY T.
c ********
c     SET UP INITIAL INDICES
 4020 if (NYNY .lt. 0) then
         IY = -NYNY
         IYNI = NYNY + ICI + 1
         if (LDT .eq. 2) then
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
         if (NKDKO .ne. 0) KORDI = KORD(NKDKO + I - 1)
         IY = IY + abs(KORDI)
         KQQ = KORD(I + 3)
c         GET INDEX OF HIGHEST ORDER DIFFERENCE TO BE USED
         K = max(abs(KQQ) + INTERP, 2)
         IYI = -INTEG
         if (KQQ) 4030, 4130, 4040
c EQUATION IS STIFF
 4030    continue
c++  Code for STIFF is inactive
c      JS=abs(KORD(NJSKO+I-1))-1
c      IYI=IYI-JS
c      IF(LNOTM1) IF (IYI) 4034,4032,4130
c      IF (KORDI.LT.0) IYI=IYI+1
c      IYI=IYI+MAXINT-abs(KORDI)
c      IF (IYI) 4034,4130,4130
cc.      IF EQUATION IS IMPLICIT DO NOT COMPUTE AN F
c 4032 IF (KORDI.LT.0) GO TO 4130
cc.      TEST IF INTEG TOO BIG FOR THIS EQUATION
c 4034 IF (abs(KORDI).LT.-IYI) GO TO 4130
c      IYI=IYI+IY
c      IYN=IYI+IYNI
cc. COMPUTE INNER PRODUCT FOR STIFF EQUATIONS
c      IF (INTEGZ.EQ.0) GO TO ???
cc.    DIFFERENTIATING
c      TP1 = C0
c      DO 4036 J = K+INTEGZ, 1, -1
c         TP1 = TP1 + C(J) * F(IDT+J-1)
c 4036 CONTINUE
cc.    TEST WHETHER TO STORE RESULT IN Y OR F
c      IF (IYI-IY) 4080, 4090, 4080
cc.    INTEGRATING OR INTERPOLATING
c      TP1 = C0
c      DO 4037 J = ICI + K, ICI + 2, -1
c         TP1 = TP1 + C(J) * F(IDT+J-ICI-1)
c 4037 CONTINUE
c      IF (INTEG.EQ.0) GO TO 4120
c      TP1=TP1 + C(ICI+1)*Y(IYN+1)
c++  End
         go to 4100
c END OF SPECIAL CODE FOR STIFF EQUATIONS
c
c EQUATION IS NOT STIFF
 4040    if (LNOTM1) if (IYI) 4050, 4060, 4130
         IYI = IYI + MAXINT - KORDI
         if (IYI .ge. 0) go to 4130
c       TEST IF INTEG TOO BIG FOR THIS EQUATION
 4050    if (KORDI .lt. -IYI) go to 4130
 4060    IYI = IYI + IY
         IYN = IYI + IYNI
c  COMPUTE INNER PRODUCT FOR EQUATION WHICH IS NOT STIFF
         XP1 = C0
         if (LDT .eq. 2) then
            if (KQQ .ne. KQMAXI) XP1 = CSUM(K+INTEGZ+ICI) *
     1         F(IDT+INTEGZ+NUMDT-1)
         end if
         do 4070 J = K + INTEGZ + ICI, ICI + 1, -1
            XP1 = XP1 + C(J) * F(IDT - ICI - 1 + J)
 4070       continue
         if (INTEG) 4080, 4090, 4100
c STORE FINAL RESULT IN Y WHEN DIFFERENTIATING
 4080    continue
         Y(IYI) = XP1
         go to 4130
c STORE INTERPOLATED VALUE IN F (OR STIFF DIFFERENTIATION)
 4090    F(I) = XP1
         go to 4130
c PICK UP EXTRA STUFF TO ADD TO INNER PRODUCT WHEN INTEGRATING
 4100    K = ICI
         if (K .eq. 0) go to 4120
 4110    continue
         XP1 = C(K) * (XP1 + Y(IYN))
         IYN = IYN - 1
         K = K - 1
         if (K .ne. 0) go to 4110
c STORE FINAL RESULT IN Y WHEN INTEGRATING (OR STIFF INTERPOLATION)
 4120    Y(IYI) = XP1 + Y(IYN)
 4130    continue
         IDT = IDT + NUMDT
 4140    continue
c
      INTEG = INTEG - 1
      if (INTEG .ge. -MAXDIF) if (INTEG) 3990, 3950, 3920
 4160 return
c ********
c ERROR PROCESSING
c ********
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
      if (XI(1) .lt. 0) then
         FDAT(5) = FDAT(6)
         FDAT(6) = TN - T(2)
      end if
 4190 FDAT(1) = T(1)
c--D Next line special: P=>S, X=>D
      call DMESS(MACT, MTXTAA, IDAT, FDAT)
      if (MACT(2) .lt. 50) go to 3820
      return
      end
c   End of DIVAIN

      subroutine DIVAOP(IOPT, FOPT)
c>> 1987-12-07 DIVAOP Krogh   Initial code.
c
c  SUBROUTINE TO SET UP OPTIONS FOR DIFFERENTIAL EQUATION  PACKAGE -IVA
      double precision FOPT(*)
      integer IOPT(*)
c
      integer KDIM, MAXORD, MAXSTF
c++ Substitute for KDIM, MAXORD, MAXSTF below
      parameter (KDIM = 20, MAXORD = 2, MAXSTF = 1)
c--D Next line special: P=>D, X=>Q
      double precision TN
      double precision XI(KDIM)
c
c--D Next line special: P=>D, X=>Q
      double precision TG(2), TGSTOP(2), TMARK, TMARKX, TOUT, TOLG
      double precision ALPHA(KDIM), BETA(KDIM+1)
      double precision  D(MAXSTF+MAXORD,MAXORD), G(KDIM,MAXORD)
      double precision V(KDIM+MAXORD)
      double precision HC, HDEC, HINC, HINCC, HMAX, HMAXP9, HMIN
      double precision FDAT(11)
c
      double precision DS(MAXSTF+MAXORD, MAXORD), GS(KDIM)
      double precision SIGMA(KDIM), RBQ(KDIM), DNOISE
      double precision EAVE, EIMAX, EIMIN, EMAX, EREP, ROBND, SNOISE
c
      integer IOPST, KORDI, KQMAXD, KQMAXI, LDT, MAXDIF, MAXINT, NKDKO,
     1   NTE, NYNY, NDTF, NUMDT
      common / DIVASC / TN, XI, IOPST, KORDI, KQMAXD, KQMAXI, LDT,
     1   MAXDIF, MAXINT, NKDKO, NTE, NYNY, NDTF, NUMDT
c
      integer ICF,ICS,IGFLG,IGTYPE(2),IGSTOP(2),ILGREP,INGS,IOP3,IOP4,
     1   IOP5,IOP6,IOP7,IOP8,IOP9,IOP10,IOP11,IOP12,IOP13,IOP14,IOP15,
     2   IOP16,IOP17,IOP18,IOP19,IOP20,IOP21,IOP22,IOP21S,ITOLEP,IY,
     3   KEMAX,KIS,KMARK,KORD1I,KORD2I,KPRED,KQDCON,KQICON,KQMAXS,
     4   KQMXDS,KQMXIL,KQMXIP,KQMXIS,KSC,KSOUT,KSSTRT,KSTEP,LEX,LINC,
     5   LINCD,LINCQ,LSC,MAXKQD,MAXKQI,METHOD,NE,NEPTOL,NG,NGTOT,
     6   NOISEQ,NOUTKO,NTOLF,NY,IDAT(6)
      common /DIVAMC/ TG,TGSTOP,TMARK,TMARKX,TOUT,TOLG,HC,HDEC,HINC,
     1   HINCC,HMAX,HMAXP9,HMIN,ALPHA,BETA,D,G,V,DS,GS,SIGMA,RBQ,DNOISE,
     2   EAVE,EIMAX,EIMIN,EMAX,EREP,ROBND,SNOISE,FDAT,ICF,ICS,IGFLG,
     3   IGTYPE,IGSTOP,ILGREP,INGS,IOP3,IOP4,IOP5,IOP6,IOP7,IOP8,IOP9,
     4   IOP10,IOP11,IOP12,IOP13,IOP14,IOP15,IOP16,IOP17,IOP18,IOP19,
     5   IOP20,IOP21,IOP22,IOP21S,ITOLEP,IY,KEMAX,KIS,KMARK,KORD1I,
     6   KORD2I,KPRED,KQDCON,KQICON,KQMAXS,KQMXDS,KQMXIL,KQMXIP,KQMXIS,
     7   KSC,KSOUT,KSSTRT,KSTEP,LEX,LINC,LINCD,LINCQ,LSC,MAXKQD,MAXKQI,
     8   METHOD,NE,NEPTOL,NG,NGTOT,NOISEQ,NOUTKO,NTOLF,NY,IDAT
      save / DIVAMC / , / DIVASC /
c.    SPECIFICATION OF ENVIRONMENTAL CONSTANTS.
      double precision EEPS10, EEPS16, EROV10, EEPS2, EEPT75, EOVEP2
      double precision OVTM75, OVD10
      common / DIVAEV / EEPS2, EEPT75, EOVEP2, OVTM75, OVD10, EEPS10,
     1   EEPS16, EROV10
      save / DIVAEV /
c
      integer IOPTS(23), INCOP(22), IOPTC(23), I, IA, J, K, LIOPT, MULTJ
      double precision CMP75, C0, CP25, CP3, CP5, CP625, CP75, CP875,
     1   CP9, C1, C1P125, C2, C4, C10, C16
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
      external D1MACH
      double precision D1MACH
      equivalence (IOPTC(3), IOP3)
      save IOPTS, LIOPT
c
c                      Declarations for error message processing.
c
      integer MECONT, MERET, MEEMES, MEIVEC
      parameter (MECONT =50)
      parameter (MERET  =51)
      parameter (MEEMES =52)
      parameter (MEIVEC =57)
      integer MACT(7), MACT1(5)
c
c ********* Error message text ***************
c[Last 2 letters of Param. name]  [Text generating message.]
cAA DIVAOP$B
cAB Error in IOPT() specifications: IOPT =$E
cAC HMIN = $F is > HMAX = $F.$E
      integer LTXTAA,LTXTAB,LTXTAC
      parameter (LTXTAA= 1,LTXTAB= 9,LTXTAC=49)
      character MTXTAA(1) * (76)
      data MTXTAA/'DIVAOP$BError in IOPT() specifications: IOPT =$EHMIN$
     * = $F is > HMAX = $F.$E'/
c **** End of text generated by pmess
c                      1   2  3   4       5  6      7
      data MACT / MEEMES, 88, 24, 0, MEIVEC, 0, MERET /
      data MACT1 / MEEMES, 28, 24, LTXTAC, MERET /
c
c                      IOP4       IOP17
      data IOPTS / 3*0, 500000, 12*0, 1, 6*0 /
c
c                  1  2    3  8  9 10 11 12  13  16   17 21 22
      data INCOP / 1, 3, 5*2, 1, 2, 3, 1, 2, 3*1, 3, 4*2, 2, 2 /
c
c ********* START OF EXECUTABLE CODE ***********************
c
      MULTJ = 1
      K = 1
 4200 I = IOPT(K)
      IA = abs(I)
c 1 and 6 lines below need 21 changed if more options are added.
      if (IA .le. 21) if (I) 4220, 4520, 4280
      if (IA .ne. 1111) go to 4490
      if (I .lt. 0) then
        MULTJ = -1
        K = K + 1
        go to 4200
      end if
      IOPT(2) = LIOPT
c
c     ****  INITIALIZE FOR STARTING A NEW INTEGRATION
      do 4210 J = 3, 23
 4210    IOPTC(J) = IOPTS(J)
      KSOUT = IOPTS(4)
      KMARK = 1 - IOPTS(1)
      KORDI = IOPTS(17)
      NKDKO = max(-KORDI, 0)
      IOPST = IOPTS(22)
      go to 4260
c
c     **** SET A NOMINAL VALUE
 4220 IOPTS(IA) = 0
      if (IA .eq. 12) go to 4420
      if (IA - 2) 4400, 4240, 4230
 4230 if (IA .eq. 4) IOPTS(4) = 500000
      if (IA .eq. 21) TOLG = 0.D0
      go to 4390
c
c     **** SET ALL OPTIONS TO THEIR NOMINAL VALUES
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
      if (IOPTS(12) .eq. 0) go to 4420
 4270 return
c
c     **** SET SPECIFIED OPTION
 4280 J = IOPT(K + 1)
      if (INCOP(IA) - 2) 4290, 4330, 4300
c     **** OPTION INVOLVES NO EXTRA PARAMETERS
 4290 IOPTS(IA) = 1
      if (IA - 2) 4400, 4400, 4390
c     **** TAKE CARE OF SECOND EXTRA PARAMETER
 4300 if (IA .ne. 10) go to 4310
      NOUTKO = IOPT(K + 2)
      if (NOUTKO) 4500, 4350, 4350
 4310 if (IA .ne. 16) go to 4320
      NTOLF = IOPT(K + 2)
      if (NTOLF) 4500, 4500, 4350
 4320 if (J .eq. 3) then
        if (KMARK .ne. 3) then
           if (XI(1)*(FOPT(IOPT(K+2)) - TMARK) .ge. C0) go to 4400
        end if
      end if
      TMARK = FOPT(IOPT(K+2))
      KMARK = J
      go to 4400
c     **** TAKE CARE OF FIRST EXTRA PARAMETER
 4330 continue
      if (IA .eq. 12) go to 4410
      if (IA .eq. 4) KSOUT = J
      if (IA .eq. 21) TOLG = FOPT(J)
 4350 IOPTS(IA) = J * MULTJ
      if (abs(IA - 7) .gt. 1) go to 4360
c     **** SET SPECIAL PARAMETERS FOR GSTOP-S
      IGFLG = 0
      NGTOT = IOPTS(7) + max(IOPTS(6), 0)
c     **** TEST FOR ERROR
      if (J .gt. 500) go to 4500
 4360 if (J .gt. 0) go to 4390
      if ((IA .eq. 5) .or. (IA .eq. 17)) go to 4390
      if (J + 1) 4500, 4370, 4380
 4370 if (IA .eq. 7) go to 4500
 4380 if ((IA .eq. 4) .or. (IA .eq. 11) .or. (IA .ge. 16)) go to 4500
c     **** STORE SAVED VALUE IN COMMON
 4390 IOPTC(IA) = IOPTS(IA)
c
c     **** INCREMENT K TO GET NEXT OPTION
 4400 K = K + INCOP(IA)
      go to 4200
c
c ******* SET UP INFORMATION FOR CHANGING STEPSIZE *********
c
c     **** TEST IF VALUES ARE ALREADY SET
 4410 if (IOPTS(12) .ne. 0) go to 4430
c     **** SET NOMINAL VALUES FOR VARIABLES ONLY SET ONCE
 4420 EREP = CP3
c     **** SET NOMINAL VALUES FOR STEPSIZE CONTROL AND ENV. CONSTANTS
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
      if (I .ne. 12) go to 4470
 4430 IOPTS(12) = J * MULTJ
      if (J) 4450, 4470, 4460
c     **** SET UP TO GIVE USER COMPLETE STEPSIZE CONTROL
 4450 EREP = C1 / EROV10
      HINC = -C2
      IOP8 = 1
c## Recent code 12/16/94
      LINCD = -2
      LINC = -2
      ROBND = C1
c## Recent code 9/6/2001
      TOLG = 0.D0
c## End of recent code
      go to 4480
c     **** SET USER VALUES FOR STEPSIZE CONTROL
 4460 if (FOPT(J) .ne. C0) HINC = max(C1P125, min(FOPT(J),C4))
      if (FOPT(J + 1) .ne. C0) HDEC = min(CP875, max(FOPT(J + 1), CP25))
      if (FOPT(J + 2) .ne. C0) HMIN = FOPT(J + 2)
      if (FOPT(J + 3) .ne. C0) HMAX = FOPT(J + 3)
      if ((HMIN .gt. HMAX) .or. (HMAX .le. 0.D0)) then
         call DMESS(MACT1, MTXTAA, IOPT, FOPT(J+2))
         KORD2I = -4
      end if
 4470 KQICON = -1
 4480 HMAXP9 = HMAX * CP9
      if (I - 1111) 4400, 4270, 4400
c
c ***************** ERROR  IN  IOPT ************************
c
 4490 IA = 1
 4500 MACT(6) = K + INCOP(IA) - 1
      call MESS(MACT, MTXTAA, IOPT)
      KORD1I = 24
      KORD2I = -4
c Usual return with no error is here.
 4520 LIOPT = K
      return
      end
c   End of DIVAOP

      subroutine DIVAPR(Y, YN, F, KORD)
c>> 1988-01-13 DIVAPR Krogh   Initial code.
c
c THIS SUBROUTINE
c   1. UPDATES THE DIFFERENCE TABLE FROM THE PREVIOUS STEP (IF NOT
c      DONE ALREADY).
c   2. PREDICTS WHAT THE VALUES OF THE DEPENDENT VARIABLES, Y, AND
c      THE DIFFERENCE TABLE, DT, WILL BE AT THE END OF THE CURRENT STEP.
c
c   Y = VECTOR OF PREDICTED VALUES COMPUTED BY THIS SUBROUTINE.
c   YN= VECTOR OF VALUES OF Y COMPUTED ON THE LAST STEP.
c   F = VECTOR OF DERIVATIVE VALUES.
c   DT= ARRAY CONTAINING DIFFERENCE TABLES.
c   KD= VECTOR GIVING ORDERS OF THE DIFFERENTIAL EQUATIONS (IF
c       EQUATIONS HAVE DIFFERENT ORDERS).
c   KQ= VECTOR OF INTEGRATION ORDERS.
c
      integer KORD(*)
c--D Next line special: P=>D, X=>Q
      double precision Y(*), YN(*)
      double precision F(*)
c
      integer KDIM, MAXORD, MAXSTF
c++ Substitute for KDIM, MAXORD, MAXSTF below
      parameter (KDIM = 20, MAXORD = 2, MAXSTF = 1)
c--D Next line special: P=>D, X=>Q
      double precision TN
      double precision XI(KDIM)
c
c--D Next line special: P=>D, X=>Q
      double precision TG(2), TGSTOP(2), TMARK, TMARKX, TOUT, TOLG
      double precision ALPHA(KDIM), BETA(KDIM+1)
      double precision  D(MAXSTF+MAXORD,MAXORD), G(KDIM,MAXORD)
      double precision V(KDIM+MAXORD)
      double precision HC, HDEC, HINC, HINCC, HMAX, HMAXP9, HMIN
      double precision FDAT(11)
c
      double precision DS(MAXSTF+MAXORD, MAXORD), GS(KDIM)
      double precision SIGMA(KDIM), RBQ(KDIM), DNOISE
      double precision EAVE, EIMAX, EIMIN, EMAX, EREP, ROBND, SNOISE
c
c.    SPECIFICATION OF ENVIRONMENTAL CONSTANTS.
      double precision EEPS10, EEPS16, EROV10, EEPS2
      double precision EEPT75, EOVEP2, OVTM75, OVD10
      common / DIVAEV / EEPS2, EEPT75, EOVEP2, OVTM75, OVD10, EEPS10,
     1   EEPS16, EROV10
      save / DIVAEV /
      integer IOPST, KORDI, KQMAXD, KQMAXI, LDT, MAXDIF, MAXINT, NKDKO,
     1   NTE, NYNY, NDTF, NUMDT
      common / DIVASC / TN, XI, IOPST, KORDI, KQMAXD, KQMAXI, LDT,
     1   MAXDIF, MAXINT, NKDKO, NTE, NYNY, NDTF, NUMDT
c
      integer ICF,ICS,IGFLG,IGTYPE(2),IGSTOP(2),ILGREP,INGS,IOP3,IOP4,
     1   IOP5,IOP6,IOP7,IOP8,IOP9,IOP10,IOP11,IOP12,IOP13,IOP14,IOP15,
     2   IOP16,IOP17,IOP18,IOP19,IOP20,IOP21,IOP22,IOP21S,ITOLEP,IY,
     3   KEMAX,KIS,KMARK,KORD1I,KORD2I,KPRED,KQDCON,KQICON,KQMAXS,
     4   KQMXDS,KQMXIL,KQMXIP,KQMXIS,KSC,KSOUT,KSSTRT,KSTEP,LEX,LINC,
     5   LINCD,LINCQ,LSC,MAXKQD,MAXKQI,METHOD,NE,NEPTOL,NG,NGTOT,
     6   NOISEQ,NOUTKO,NTOLF,NY,IDAT(6)
      common /DIVAMC/ TG,TGSTOP,TMARK,TMARKX,TOUT,TOLG,HC,HDEC,HINC,
     1   HINCC,HMAX,HMAXP9,HMIN,ALPHA,BETA,D,G,V,DS,GS,SIGMA,RBQ,DNOISE,
     2   EAVE,EIMAX,EIMIN,EMAX,EREP,ROBND,SNOISE,FDAT,ICF,ICS,IGFLG,
     3   IGTYPE,IGSTOP,ILGREP,INGS,IOP3,IOP4,IOP5,IOP6,IOP7,IOP8,IOP9,
     4   IOP10,IOP11,IOP12,IOP13,IOP14,IOP15,IOP16,IOP17,IOP18,IOP19,
     5   IOP20,IOP21,IOP22,IOP21S,ITOLEP,IY,KEMAX,KIS,KMARK,KORD1I,
     6   KORD2I,KPRED,KQDCON,KQICON,KQMAXS,KQMXDS,KQMXIL,KQMXIP,KQMXIS,
     7   KSC,KSOUT,KSSTRT,KSTEP,LEX,LINC,LINCD,LINCQ,LSC,MAXKQD,MAXKQI,
     8   METHOD,NE,NEPTOL,NG,NGTOT,NOISEQ,NOUTKO,NTOLF,NY,IDAT
      save / DIVAMC / , / DIVASC /
      double precision C0
      parameter (C0 = 0.D0)
c
      integer I, INTEG, INTEGS, J, K, KQQ, L, N
      double precision TEMP(KDIM)
      double precision TP1
c--D Next line special: P=>D, X=>Q
      double precision XP
      data INTEGS / -1 /
c
c++  Code for ARGM is inactive
c      RETURN
c      ENTRY DIVAPE
c++  End
c ********
c START OF CODE
c ********
      IY = 0
      L = NDTF - 1
      do 4680 I = 1, NTE
         INTEG = KORDI
         if (NKDKO .ne. 0) INTEG = KORD(NKDKO + I - 1)
         KQQ = KORD(I + 3)
         K = max(abs(KQQ), 2)
         if (KQQ) 4530, 4520, 4540
 4520    IY = IY + abs(INTEG)
         go to 4670
c ********
c EQUATION IS STIFF, OR IMPLICIT
c ********
 4530    continue
c++  Code for STIFF is inactive
c      KQQ=-KQQ
c      N=KQQ-1
c      JS=abs(KORD(NJSKO+I-1))-1
c      IMPLIC=INTEG
c      INTEG=abs(IMPLIC)-JS
cc.    SET INTEGS FOR STIFF EQUATIONS
c      INTEGS=0
c      IF (K-KSC) 160,160,140
cc.END OF SET-UP FOR STIFF EQUATIONS
c++  End
c ********
c EQUATION IS NOT STIFF
c ********
 4540    N = KQQ
         if (LDT .ne. 0) if (K - KSC) 4570, 4570, 4550
c     DIFFERENCE TABLE HAS NOT BEEN UPDATED
         TP1 = F(I) - F(L + 1)
         if (K - KSC) 4610, 4610, 4590
c END OF SET-UP FOR EQUATIONS WHICH ARE NOT STIFF
c ********
c GET PREDICTED DIFFERENCES FROM UPDATED DIFFERENCE TABLE
c ********
 4550    F(L + K + 1) = F(L + K + 1) * BETA(K + 1)
         TEMP(K) = F(L + K) * BETA(K)
         F(L + K) = TEMP(K)
c LOOP FOR MODIFIED DIVIDED DIFFERENCES
 4560    K = K - 1
         if (K .le. KSC) go to 4580
         TEMP(K) = F(L + K) * BETA(K)
         F(L + K) = TEMP(K) + F(L + K + 1)
         go to 4560
c CODE FOR BACKWARD DIFFERENCES
 4570    F(L + K + 1) = F(L + K + 1)
         TEMP(K) = F(L + K)
         K = K - 1
c
 4580    TEMP(K) = F(L + K)
         F(L + K) = TEMP(K) + F(L + K + 1)
         K = K - 1
         if (K .ne. 0) go to 4580
         go to 4630
c ********
c UPDATE DIFFERENCE TABLE AND GET PREDICTED DIFFERENCES
c ********
c CODE FOR MODIFIED DIVIDED DIFFERENCES
 4590    F(L + K + 1) = (F(L+K+1) + TP1) * BETA(K + 1)
         TEMP(K) = (F(L + K) + TP1) * BETA(K)
         F(L + K) = TEMP(K)
 4600    K = K - 1
         if (K .le. KSC) go to 4620
         TEMP(K) = (F(L + K) + TP1) * BETA(K)
         F(L + K) = TEMP(K) + F(L + K + 1)
         go to 4600
c CODE FOR BACKWARD DIFFERENCES
 4610    F(L + K + 1) = (F(L+K+1) + TP1)
         TEMP(K) = F(L + K) + TP1
         F(L + K) = TEMP(K)
         K = K - 1
c
 4620    TEMP(K) = F(L + K) + TP1
         F(L + K) = TEMP(K) + F(L + K + 1)
         K = K - 1
         if (K .ne. 0) go to 4620
c ********
c COMPUTE Y-S OBTAINED USING INTEGRATION
c ********
c     TEST IF NEXT Y TO BE OBTAINED BY INTERPOLATION
 4630    continue
c++  Code for STIFF is inactive
c      IF (INTEG.EQ.0) GO TO 4662
c++  End
         IY = IY + 1
c     FORM INNER PRODUCT
         XP = C0
         do 4650 J = INTEGS + N + 1, INTEGS + 2, -1
c++  Code for ~{p,x} is active
            XP = XP + G(J, INTEG) * TEMP(J)
c++  Code for {p,x} is inactive
Cc--D Next line special: P=>D, X=>Q
C            XP = XP + dble(G(J, INTEG)) * dble(TEMP(J))
c++  END
 4650    continue
         K = INTEG + INTEGS
         do 4660 J = K, 1, -1
c++  Code for ~{p,x} is active
            XP = XP + G(1, J) * YN(IY + J)
c++  Code for {p,x} is inactive
Cc--D Next line special: P=>D, X=>Q
C            XP = XP + dble(G(1, J)) * dble(YN(IY + J))
c++  END
 4660    continue
         Y(IY) = YN(IY) + XP
         INTEG = INTEG - 1
         if (K) 4670, 4670, 4630
c END OF COMPUTING Y-S OBTAINED BY INTEGRATION
c ********
c COMPUTE Y-S OBTAINED USING INTERPOLATION AND DIFFERENTIATION
c ********
c++  Code for STIFF is inactive
cc.    RESTORE INTEGS FOR EQUATIONS WHICH ARE NOT STIFF
c 4662 INTEGS=-1
c      IY=IY+1
cc.    COMPUTE Y USING INTERPOLATION
c      Y(IY)=YN(IY) + F(L+2)
c      IF (KQQ.EQ.1) Y(IY)=YN(IY)
c 4663 INTEG=INTEG+1
c      IF (INTEG.EQ.JS) IF (IMPLIC) 4680,4680,4664
cc.    COMPUTE INTEG-TH DERIVATIVE
c      XP = C0
c 4664 DO 4666 J = KQQ+1, INTEG+1, -1
c         XP = XP + D(J, INTEG) * TEMP(J)
c 4666 CONTINUE
c      IF (INTEG.EQ.JS) GO TO 4667
c      IY=IY+1
c      Y(IY)=XP
c      GO TO 4663
cc.STORE PREDICTED VALUE FOR F
c 4667 CONTINUE
c      F(L+NUMDT)=XP
c++  End
 4670    L = L + NUMDT
 4680    continue
      LDT = -3
      return
      end
      subroutine DIVADB(LPRINT, TSPECS, Y, F, KORD, TEXT)
c Copyright (c) 1996 California Institute of Technology, Pasadena, CA.
c ALL RIGHTS RESERVED.
c Based on Government Sponsored Research NAS7-03001.
C>> 2009-11-04 DIVADB Krogh Included TOLG, initilized the unitialized.
C>> 2000-12-01 DIVADB Krogh Removed unused parameter METXTF.
C>> 1996-07-02 DIVADB Krogh Transpose flag for matrix output in C.
C>> 1996-03-25 DIVADB Krogh Introduced TEXT1-TEXT4 to comply with F77.
C>> 1996-01-19 DIVADB Krogh Changed NTEXT to TEXT to agree with doc.
C>> 1995-04-26 DIVADB Krogh Fixed print of V & G's for high order eqs.
C>> 1994-11-11 DIVADB Krogh Declared all vars.
c>> 1994-10-20 DIVADB Krogh Changes to use M77CON
c>> 1994-09-12 DIVADB Krogh Added CHGTYP code.
c>> 1994-03-07 DIVADB Krogh Allow larger order in single precision.
c>> 1993-05-03 DIVADB Krogh Additions for Conversion to C.
c>> 1993-04-14 DIVADB Krogh Changes for new MESS usage.
c>> 1992-04-08 DIVADB Krogh Unused labels 10 and 60 removed.
c>> 1992-03-10 DIVADB Krogh Fixed value for KDIM in single p. version.
c>> 1992-02-17 DIVADB Krogh Made tabs depend on # digits output.
c>> 1991-11-04 DIVADB Krogh Switched to use MESS, DMESS
c>> 1990-03-08 DIVADB Krogh Unused stiff vars. set to 0.
c>> 1989-07-21 DIVADB Krogh Code for integrating discontinuities
c>> 1988-06-07 DIVADB Krogh Dim. of IVC2 and DVC2 upped by 1 (old bug)
c>> 1987-12-07 DIVADB Krogh Initial code.
c
c--D replaces "?": ?IVADB, ?IVAEV, ?IVAMC, ?IVASC, ?MESS
c
c SUBROUTINE TO GIVE DEBUG PRINT FOR VARIABLE ORDER INTEGRATOR
c
c  LET ABS(LPRINT)=  10*N1 + N2     (N1,N2 DIGITS)
c    N1=1   DO NOT PRINT ANY VARIABLES EXTERNAL TO THE INTEGRATOR
c    N1=2   PRINT  TSPECS, CURRENT Y, PAST Y, CURRENT F,
c           ALL PERTINENT CONTENTS OF KORD, AND TOL.
c    N1=3   ABOVE + DIFFERENCE TABLES UP TO HIGHEST DIFFERENCE USED
c    N1=4   SAME AS N1=1 + ALL IN STORAGE ALLOCATED FOR DIFFERENCES
c
c    N2=1   DO NOT PRINT ANY VARIABLES INTERNAL TO THE INTEGRATOR
c    N2=2   PRINT ALL SCALAR VARIABLES IN INTERPOLATION COMMON BLOCK
c    N2=3   ABOVE + ALL SCALAR VARIABLES IN MAIN INTEG. COMMON BLOCK
c    N2=4   SAME AS N1=3 + ALL USED IN ARRAYS XI,BETA,ALPHA, FIRST
c           COLUMN OF G, GS,RBQ,SIGMA
c    N2=5   SAME AS N1=4 + ALL USED IN ARRAYS G,D,DS,V
c
      integer LPRINT, KORD(*)
      character TEXT*(*)
      character TEXT1(1)*11, TEXT2(1)*4, TEXT3(1)*5, TEXT4(1)*4
      integer IVC1(12), IVC2(65), J, K, L, N1, N2
      double precision  DVC2(7), RVC2(8), EVC(8)
c--D Next line special: P=>D, X=>Q
      double precision TSPECS(*), Y(*), TNEQ(1), DVC1(7)
      double precision F(*)
c
c++S Default KDIM = 16
c++  Default KDIM = 20
c++  Default MAXORD = 2, MAXSTF = 1
c++  Default STIFF=.F.
      integer KDIM, MAXORD, MAXSTF
c++ Substitute for KDIM, MAXORD, MAXSTF below
      parameter (KDIM = 20, MAXORD = 2, MAXSTF = 1)
c--D Next line special: P=>D, X=>Q
      double precision TN
      double precision XI(KDIM)
c
c--D Next line special: P=>D, X=>Q
      double precision TG(2), TGSTOP(2), TMARK, TMARKX, TOUT, TOLG
      double precision ALPHA(KDIM), BETA(KDIM+1)
      double precision  D(MAXSTF+MAXORD,MAXORD), G(KDIM,MAXORD)
      double precision V(KDIM+MAXORD)
      double precision HC, HDEC, HINC, HINCC, HMAX, HMAXP9, HMIN
      double precision FDAT(11)
c
      double precision DS(MAXSTF+MAXORD, MAXORD), GS(KDIM)
      double precision SIGMA(KDIM), RBQ(KDIM), DNOISE
      double precision EAVE, EIMAX, EIMIN, EMAX, EREP, ROBND, SNOISE
c
c.    SPECIFICATION OF ENVIRONMENTAL CONSTANTS.
      double precision EEPS10, EEPS16, EROV10, EEPS2
      double precision EEPT75, EOVEP2, OVTM75, OVD10
      common / DIVAEV / EEPS2, EEPT75, EOVEP2, OVTM75, OVD10, EEPS10,
     1   EEPS16, EROV10
      save / DIVAEV /
      integer IOPST, KORDI, KQMAXD, KQMAXI, LDT, MAXDIF, MAXINT, NKDKO,
     1   NTE, NYNY, NDTF, NUMDT
      common / DIVASC / TN, XI, IOPST, KORDI, KQMAXD, KQMAXI, LDT,
     1   MAXDIF, MAXINT, NKDKO, NTE, NYNY, NDTF, NUMDT
c
      integer ICF,ICS,IGFLG,IGTYPE(2),IGSTOP(2),ILGREP,INGS,IOP3,IOP4,
     1   IOP5,IOP6,IOP7,IOP8,IOP9,IOP10,IOP11,IOP12,IOP13,IOP14,IOP15,
     2   IOP16,IOP17,IOP18,IOP19,IOP20,IOP21,IOP22,IOP21S,ITOLEP,IY,
     3   KEMAX,KIS,KMARK,KORD1I,KORD2I,KPRED,KQDCON,KQICON,KQMAXS,
     4   KQMXDS,KQMXIL,KQMXIP,KQMXIS,KSC,KSOUT,KSSTRT,KSTEP,LEX,LINC,
     5   LINCD,LINCQ,LSC,MAXKQD,MAXKQI,METHOD,NE,NEPTOL,NG,NGTOT,
     6   NOISEQ,NOUTKO,NTOLF,NY,IDAT(6)
      common /DIVAMC/ TG,TGSTOP,TMARK,TMARKX,TOUT,TOLG,HC,HDEC,HINC,
     1   HINCC,HMAX,HMAXP9,HMIN,ALPHA,BETA,D,G,V,DS,GS,SIGMA,RBQ,DNOISE,
     2   EAVE,EIMAX,EIMIN,EMAX,EREP,ROBND,SNOISE,FDAT,ICF,ICS,IGFLG,
     3   IGTYPE,IGSTOP,ILGREP,INGS,IOP3,IOP4,IOP5,IOP6,IOP7,IOP8,IOP9,
     4   IOP10,IOP11,IOP12,IOP13,IOP14,IOP15,IOP16,IOP17,IOP18,IOP19,
     5   IOP20,IOP21,IOP22,IOP21S,ITOLEP,IY,KEMAX,KIS,KMARK,KORD1I,
     6   KORD2I,KPRED,KQDCON,KQICON,KQMAXS,KQMXDS,KQMXIL,KQMXIP,KQMXIS,
     7   KSC,KSOUT,KSSTRT,KSTEP,LEX,LINC,LINCD,LINCQ,LSC,MAXKQD,MAXKQI,
     8   METHOD,NE,NEPTOL,NG,NGTOT,NOISEQ,NOUTKO,NTOLF,NY,IDAT
      save / DIVAMC / , / DIVASC /
      equivalence (IVC1(1), IOPST), (IVC2(1), ICF)
      equivalence (TNEQ, TN)
      equivalence (RVC2(1), DNOISE), (DVC1(1), TG), (DVC2(1), HC),
     1   (EVC(1), EEPS2)
c
c                      Declarations for error message processing.
      integer MEDDIG, NEDDIG, METDIG, METABS, MERET, METEXT,
     1   METABL, MEIVEC, MEFVEC, MEFMAT
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
      integer MACT0(3), MACT1(2), MACT2(7), MACT3(7), MACT4(8),
     1   MACT5(11), MACT6(3), MACT7(14), MACTFV(4)
      integer KPI, KPE
c                      wddtrr        wwddtrr
      parameter (KPI = 400301, KPE = 1305501)
c
c ********* Error message text ***************
c[Last 2 letters of Param. name]  [Text generating message.]
cAA $NKORD:    $B
cAB Int. Ord.: $B
c   $
cAC D.E. Ord.: $B
c   $
cAD Meth.Type: $B
c   $
cAE Tolerance Groups: $B
cAF Tolerances: $B
c   $
cAG $NDifferences$B
cAH Eq. $#
cAI Ord. $#
c   $
cAJ $NTN=$F$E
c   $
cAK $NIOPST=$I$TKORDI=$I$TKQMAXD=$I$TKQMAXI=$I$TLDT=$I$T$C
c   MAXDIF=$I$TMAXINT=$I$TNKDKO=$I$TNTE=$I$TNYNY=$I$TNDTF=$I$C
c   $TNUMDT=$I$E
c   $
cAL $NICF=$I$TICS=$I$TIGFLG=$I$TIGTYPE(1)=$I$TIGTYPE(2)=$I$T$C
c   IGSTOP(1)=$I$TIGSTOP(2)=$I$TILGREP=$I$TINGS=$I$TIOP3=$I$T$C
c   IOP4=$I$TIOP5=$I$TIOP6=$I$TIOP7=$I$TIOP8=$I$TIOP9=$I$T$C
c   IOP10=$I$TIOP11=$I$TIOP12=$I$TIOP13=$I$TIOP14=$I$TIOP15=$I$T$C
c   IOP16=$I$TIOP17=$I$TIOP18=$I$TIOP19=$I$TIOP20=$I$TIOP21=$I$T$C
c   IOP22=$I$TIOP21S=$I$TITOLEP=$I$TIY=$I$TKEMAX=$I$TKIS=$I$T$C
c   KMARK=$I$TKORD1I=$I$TKORD2I=$I$TKPRED=$I$TKQDCON=$I$T$C
c   KQICON=$I$TKQMAXS=$I$TKQMXDS=$I$TKQMXIL=$I$TKQMXIP=$I$T$C
c   KQMXIS=$I$TKSC=$I$TKSOUT=$I$TKSSTRT=$I$TKSTEP=$I$TLEX=$I$T$C
c   LINC=$I$TLINCD=$I$TLINCQ=$I$TLSC=$I$TMAXKQD=$I$TMAXKQI=$I$T$C
c   METHOD=$I$TNE=$I$TNEPTOL=$I$TNG=$I$TNGTOT=$I$TNOISEQ=$I$T$C
c   NOUTKO=$I$TNTOLF=$I$TNY=$I$E
cAM $NDNOISE=$F$TEAVE=$F$TEIMAX=$F$TEIMIN=$F$TEMAX=$F$T$C
c   EREP=$F$TROBND=$F$E
c   $
cAN $NTG(1)=$F$TTG(2)=$F$TTGSTOP(1)=$F$TTGSTOP(2)=$F$C
c   $TTMARK=$F$TTMARKX=$F$TTOUT=$F$E
c   $
cAO HC=$F$THDEC=$F$THINC=$F$THINCC=$F$THMAX=$F$T$C
c   HMAXP9=$F$THMIN=$F$T$N$E
c   $
cAP K$HXI(K)$HBETA(K)$HALPHA(K)$HG(K,1)$HRBQ(K)$HSIGMA(K)$H
c   GS(K)$HV(K)$HG(K,2..MAXINT)$E
c   $
cAQ $NEEPS2=$F$TEEPT75=$F$TEOVEP2=$F$TOVTM75=$F$TOVD10=$F$T$C
c   EEPS10=$F$TEEPS16=$F$TEROV10=$F$E
      integer LTXTAA,LTXTAB,LTXTAC,LTXTAD,LTXTAE,LTXTAF,LTXTAG,LTXTAH,
     * LTXTAI,LTXTAJ,LTXTAK,LTXTAL,LTXTAM,LTXTAN,LTXTAO,LTXTAP,LTXTAQ
      parameter (LTXTAA=  1,LTXTAB= 14,LTXTAC=  1,LTXTAD=  1,LTXTAE=  1,
     * LTXTAF= 21,LTXTAG=  1,LTXTAH= 16,LTXTAI= 22,LTXTAJ=  1,
     * LTXTAK=  1,LTXTAL=  1,LTXTAM=655,LTXTAN=  1,LTXTAO=  1,
     * LTXTAP=  1,LTXTAQ=  1)
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
      data MTXTAG/'$NIOPST=$I$TKORDI=$I$TKQMAXD=$I$TKQMAXI=$I$TLDT=$I$TM
     *AXDIF=$I$TMAXINT=$I$TNKDKO=$I$TNTE=$I$TNYNY=$I$TNDTF=$I$TNUMDT=$I$
     *E'/
      data MTXTAH/'$NICF=$I$TICS=$I$TIGFLG=$I$TIGTYPE(1)=$I$TIGTYPE(2)=$
     *I$TIGSTOP(1)=$I$TIGSTOP(2)=$I$TILGREP=$I$TINGS=$I$TIOP3=$I$TIOP4=$
     *I$TIOP5=$I$TIOP6=$I$TIOP7=$I$TIOP8=$I$TIOP9=$I$TIOP10=$I$TIOP11=$I
     *$TIOP12=$I$TIOP13=$I$TIOP14=$I$TIOP15=$I$TIOP16=$I$TIOP17','=$I$TI
     *OP18=$I$TIOP19=$I$TIOP20=$I$TIOP21=$I$TIOP22=$I$TIOP21S=$I$TITOLEP
     *=$I$TIY=$I$TKEMAX=$I$TKIS=$I$TKMARK=$I$TKORD1I=$I$TKORD2I=$I$TKPRE
     *D=$I$TKQDCON=$I$TKQICON=$I$TKQMAXS=$I$TKQMXDS=$I$TKQMXIL=$I$TKQMXI
     *P=$I$TKQMXIS=$I$TKSC=$I$TKSOUT=$I$TKSS','TRT=$I$TKSTEP=$I$TLEX=$I$
     *TLINC=$I$TLINCD=$I$TLINCQ=$I$TLSC=$I$TMAXKQD=$I$TMAXKQI=$I$TMETHOD
     *=$I$TNE=$I$TNEPTOL=$I$TNG=$I$TNGTOT=$I$TNOISEQ=$I$TNOUTKO=$I$TNTOL
     *F=$I$TNY=$I$E$NDNOISE=$F$TEAVE=$F$TEIMAX=$F$TEIMIN=$F$TEMAX=$F$TER
     *EP=$F$TROBND=$F$E  '/
      data MTXTAI/'$NTG(1)=$F$TTG(2)=$F$TTGSTOP(1)=$F$TTGSTOP(2)=$F$TTMA
     *RK=$F$TTMARKX=$F$TTOUT=$F$E'/
      data MTXTAJ/'HC=$F$THDEC=$F$THINC=$F$THINCC=$F$THMAX=$F$THMAXP9=$F
     *$THMIN=$F$T$N$E'/
      data MTXTAK/'K$HXI(K)$HBETA(K)$HALPHA(K)$HG(K,1)$HRBQ(K)$HSIGMA(K)
     *$HGS(K)$HV(K)$HG(K,2..MAXINT)$E'/
      data MTXTAL/'$NEEPS2=$F$TEEPT75=$F$TEOVEP2=$F$TOVTM75=$F$TOVD10=$F
     *$TEEPS10=$F$TEEPS16=$F$TEROV10=$F$E'/
c
      data MACT0 / METABS, 10, MERET /
      data MACT1 / METEXT, MERET /
c                       1       2  3       4       5  6      7
      data MACT2 / METEXT, MEIVEC, 3, METEXT, MEIVEC, 0, MERET /
c                       1       2  3       4       5  6      7
      data MACT3 / METEXT, MEIVEC, 0, METEXT, MEFVEC, 0, MERET /
c                       1       2  3  4  5       6       7
      data MACT4 / METEXT, MEFMAT, 0, 0, 0, LTXTAI, LTXTAH, MERET /
c                       1   2       3       4   5       6  7       8
      data MACT5 / METABS, 12, METEXT, METABS, 18, METDIG, 5, METEXT,
     1             METABS, 0, MERET /
c                       9 10     11
      data MACT6 / NEDDIG, 0, MERET /
c                         2 3 4   5 6 7 8 9  10  11  12 13 14
      data MACT7 / METABL,0,0,0,KPI,0,0,0,0,KPE,KPE,KPE, 0, 0 /
c                        1       2  3      4
      data MACTFV / METEXT, MEFVEC, 3, MERET /
c
      data TEXT1 / '$NTSPECS:$B' /
      data TEXT2 / 'Y:$B' /
      data TEXT3 / 'YN:$B' /
      data TEXT4 / 'F:$B' /
c
c ********
c START OF CODE -- PRINT TEXT AND SET INDEX FOR F
c ********
c    Getting variables that are not yet assigned some values.
c++  Code for ~STIFF is active
      KQDCON = 0
      KQMXDS = 0
      MAXKQD = 0
c++  End
      GS(1) = 1.D0
      if (IOP6 .eq. 0) then
        IGTYPE(1) = 0
        IGSTOP(1) = 0
        TG(1) = 0.D0
        TGSTOP(1) = 0.D0
      end if
      if (IOP7 .eq. 0) then
        IGTYPE(2) = 0
        IGSTOP(2) = 0
        TG(2) = 0.D0
        TGSTOP(2) = 0.D0
      end if
      if (IOP6 + IOP7 .eq. 0) then
        INGS = 0
        NG = 0
      end if
      if (IOP10 .eq. 0) NOUTKO = 0
      J = 0
      call MESSFT(MACT0, TEXT)
c
      N1 = LPRINT / 10
      N2 = LPRINT - 10 * N1
      if (N1 .le. 1) go to 80
c ********
c PRINT ALL EXTERNAL VARIABLES EXCEPT FOR THE DIFFERENCES
c ********
      MACTFV(3) = max(IOP5, 4)
c--D Next line special: P=>D, X=>Q
      call DMESS(MACTFV, TEXT1, KORD, TSPECS)
      MACTFV(3) = NY
c--D Next line special: P=>D, X=>Q
      call DMESS(MACTFV, TEXT2, KORD, Y)
c--D Next line special: P=>D, X=>Q
      call DMESS(MACTFV, TEXT3, KORD, Y(NYNY))
      MACTFV(3) = NTE
c--D Next line special: P=>S, X=>D
      call DMESS(MACTFV, TEXT4, KORD, F)
      MACT2(6) = NTE
      call MESS(MACT2, MTXTAA, KORD)
      if (NKDKO .gt. 0) call MESS(MACT2(4), MTXTAB, KORD(NKDKO))
      if (IOPST .gt. 0) call MESS(MACT2(4), MTXTAC, KORD(IOPST))
c WRITE TOL
      K = IOP16
   70 if (KORD(K) .lt. 0) K = K + 1
      K = K + 1
      if (KORD(K - 1) .lt. NTE) go to 70
      MACT3(3) = K - IOP16
      MACT3(6) = MACT3(3)
c--D Next line special: P=>S, X=>D
      call DMESS(MACT3, MTXTAD, KORD(IOP16), F(NTOLF))
      if (N1 .eq. 2) go to 80
c ********
c WRITE THE DIFFERENCE TABLES
c ********
      K = NUMDT
      if (N1 .eq. 3) K = KQMAXS
      MACT4(3) = NUMDT
      MACT4(4) = -K
      MACT4(5) = NTE
c--D Next line special: P=>S, X=>D
      call DMESS(MACT4, MTXTAE, KORD, F(NDTF))
c
   80 if (N2 .le. 1) return
c ********
c WRITE SCALARS IN COMMON
c ********
c--D Next line special: P=>D, X=>Q
      call DMESS(MACT1, MTXTAF, KORD, TNEQ)
c
c ===== COMMON 1  -- INTEGER
c
      call MESS(MACT1, MTXTAG, IVC1)
      if (N2 .eq. 2) return
      call MESS(MACT6, MTXTAA, IDAT)
      MACT5(10) = MACT6(2) + 14
c
c ===== COMMON 2  -- INTEGER AND FLOATING POINT
c
c--D Next line special: P=>S, X=>D
      call DMESS(MACT5, MTXTAH, IVC2, RVC2)
c--D Next line special: P=>D, X=>Q
      call DMESS(MACT1, MTXTAI, IVC2, DVC1)
c--D Next line special: P=>S, X=>D
      call DMESS(MACT1, MTXTAJ, IVC2, DVC2)
      if (N2 .eq. 3) return
c         wddtrr              wddtrr
      J = 101000 * MACT6(2) + 800501
      MACT7(2) = 1
      MACT7(3) = KQMAXS
      MACT7(4) = 8
      do 90 K = 6, 9
         MACT7(K) = J
   90 continue
      if (N2 .gt. 0) then
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
         if (N2 .ge. 4) then
            FDAT(8) = V(K)
            do 95 J = 2, L
               FDAT(7+J) = G(K, J)
   95       continue
         end if
c--D Next line special: P=>S, X=>D
         call DMESS(MACT7, MTXTAK, IDAT, FDAT)
  100 continue
c++  Code for STIFF is inactive
c     if (MAXDIF .le. 0) return
c        Need to define MACT8 and set values
cc--D Next line special: P=>S, X=>D
c     call DMESS(MACT8, 'D$B', IDAT, D)
cc--D Next line special: P=>S, X=>D
c     call DMESS(MACT8, 'DS$B', IDAT, DS)
c++  End
c
c--D Next line special: P=>S, X=>D
      call DMESS(MACT1, MTXTAL, IDAT, EVC)
      return
      end
      subroutine DIVAG(TSPECS, Y, F, KORD, IFLAG, NSTOP, GNEW, GT)
c Copyright (c) 1996 California Institute of Technology, Pasadena, CA.
c ALL RIGHTS RESERVED.
c Based on Government Sponsored Research NAS7-03001.
c>> 2001-09-07 DIVAG  Krogh  Changes to allow user tol on G-Stops.
c>> 1995-06-20 DIVAG  Krogh  Fixed problem introduced with last change.
c>> 1995-05-09 DIVAG  Krogh  Fixed G-Stop/discontinuity code interaction
C>> 1994-11-11 DIVAG  Krogh  Declared all vars.
c>> 1994-10-20 DIVAG  Krogh  Changes to use M77CON
c>> 1994-09-12 DIVAG  Krogh  Added CHGTYP code.
c>> 1994-08-17 DIVAG  Krogh  Modified internal comment.
c>> 1994-03-07 DIVAG  Krogh  Allow larger order in single precision.
c>> 1993-04-27 DIVAG  Krogh  Additions for Conversion to C.
c>> 1992-10-12 DIVAG  Krogh  Fixed G-Stop/discontinuity code interaction
c>> 1992-09-17 DIVAG  Krogh  Slight change in check for sign change.
c>> 1992-04-08 DIVAG  Krogh  Unused labels 140,150,230, and 250 removed.
c>> 1992-03-10 DIVAG  Krogh  Fixed value for KDIM in single p. version.
c>> 1990-01-29 DIVAG  Krogh  Added arg to call to DERMN.
c>> 1988-03-04 DIVAG  Krogh  Initial code.
c
c--D replaces "?": ?IVAG,?IVABU,?IVAEV,?IVAIN,?IVAMC,?IVASC,?MESS,?ZERO
c
c     -XXXG(TSPECS, Y, F, KORD, IFLAG, NSTOP, GNEW, GT)
c
c  SUBROUTINE TO LOCATE OUTPUT POINTS AT ZEROS OF ARBITRARY
c  FUNCTIONS  **** GSTOPS **** FOR DIFFERENTIAL EQUATION
c  INTEGRATOR -ODE (OR -IVA).
c
      integer KORD(*), IFLAG, NSTOP
c--D Next line special: P=>D, X=>Q
      double precision TSPECS(*), Y(*), TOLD, TSAVE
      double precision F(*), GNEW(*), GT(*), GOLD
      save GOLD, TOLD, TSAVE
c
c++SP Default KDIM = 16
c++  Default KDIM = 20
c++  Default MAXORD = 2, MAXSTF = 1
      integer KDIM, MAXORD, MAXSTF
c++ Substitute for KDIM, MAXORD, MAXSTF below
      parameter (KDIM = 20, MAXORD = 2, MAXSTF = 1)
c--D Next line special: P=>D, X=>Q
      double precision TN
      double precision XI(KDIM)
c
c--D Next line special: P=>D, X=>Q
      double precision TG(2), TGSTOP(2), TMARK, TMARKX, TOUT, TOLG
      double precision ALPHA(KDIM), BETA(KDIM+1)
      double precision  D(MAXSTF+MAXORD,MAXORD), G(KDIM,MAXORD)
      double precision V(KDIM+MAXORD)
      double precision HC, HDEC, HINC, HINCC, HMAX, HMAXP9, HMIN
      double precision FDAT(11)
c
      double precision DS(MAXSTF+MAXORD, MAXORD), GS(KDIM)
      double precision SIGMA(KDIM), RBQ(KDIM), DNOISE
      double precision EAVE, EIMAX, EIMIN, EMAX, EREP, ROBND, SNOISE
c
c.    SPECIFICATION OF ENVIRONMENTAL CONSTANTS.
      double precision EEPS10, EEPS16, EROV10, EEPS2
      double precision EEPT75, EOVEP2, OVTM75, OVD10
      common / DIVAEV / EEPS2, EEPT75, EOVEP2, OVTM75, OVD10, EEPS10,
     1   EEPS16, EROV10
      save / DIVAEV /
      integer IOPST, KORDI, KQMAXD, KQMAXI, LDT, MAXDIF, MAXINT, NKDKO,
     1   NTE, NYNY, NDTF, NUMDT
      common / DIVASC / TN, XI, IOPST, KORDI, KQMAXD, KQMAXI, LDT,
     1   MAXDIF, MAXINT, NKDKO, NTE, NYNY, NDTF, NUMDT
c
      integer ICF,ICS,IGFLG,IGTYPE(2),IGSTOP(2),ILGREP,INGS,IOP3,IOP4,
     1   IOP5,IOP6,IOP7,IOP8,IOP9,IOP10,IOP11,IOP12,IOP13,IOP14,IOP15,
     2   IOP16,IOP17,IOP18,IOP19,IOP20,IOP21,IOP22,IOP21S,ITOLEP,IY,
     3   KEMAX,KIS,KMARK,KORD1I,KORD2I,KPRED,KQDCON,KQICON,KQMAXS,
     4   KQMXDS,KQMXIL,KQMXIP,KQMXIS,KSC,KSOUT,KSSTRT,KSTEP,LEX,LINC,
     5   LINCD,LINCQ,LSC,MAXKQD,MAXKQI,METHOD,NE,NEPTOL,NG,NGTOT,
     6   NOISEQ,NOUTKO,NTOLF,NY,IDAT(6)
      common /DIVAMC/ TG,TGSTOP,TMARK,TMARKX,TOUT,TOLG,HC,HDEC,HINC,
     1   HINCC,HMAX,HMAXP9,HMIN,ALPHA,BETA,D,G,V,DS,GS,SIGMA,RBQ,DNOISE,
     2   EAVE,EIMAX,EIMIN,EMAX,EREP,ROBND,SNOISE,FDAT,ICF,ICS,IGFLG,
     3   IGTYPE,IGSTOP,ILGREP,INGS,IOP3,IOP4,IOP5,IOP6,IOP7,IOP8,IOP9,
     4   IOP10,IOP11,IOP12,IOP13,IOP14,IOP15,IOP16,IOP17,IOP18,IOP19,
     5   IOP20,IOP21,IOP22,IOP21S,ITOLEP,IY,KEMAX,KIS,KMARK,KORD1I,
     6   KORD2I,KPRED,KQDCON,KQICON,KQMAXS,KQMXDS,KQMXIL,KQMXIP,KQMXIS,
     7   KSC,KSOUT,KSSTRT,KSTEP,LEX,LINC,LINCD,LINCQ,LSC,MAXKQD,MAXKQI,
     8   METHOD,NE,NEPTOL,NG,NGTOT,NOISEQ,NOUTKO,NTOLF,NY,IDAT
      save / DIVAMC / , / DIVASC /
c
      integer I, IG, IGFLGS, IZFLAG, KEXIT, NGSTOP(2)
      double precision HH
      equivalence (G(1,1), HH), (NGSTOP(1), IOP6), (KEXIT, IOP17),
     1   (IZFLAG, IY), (IGFLGS, ITOLEP)
c
c                      Declarations for error message processing.
c
      integer MEMDA1, MERET, MEEMES
      parameter (MEMDA1 =27)
      parameter (MERET  =51)
      parameter (MEEMES =52)
      integer MACT(7)
c ********* Error message text ***************
c[Last 2 letters of Param. name]  [Text generating message.]
cAA DIVAG$B
cAB Call with bad values of KORD.  KORD(1)=$I, KORD(2)=$I, when $C
c   TSPECS(1)=$F and KSTEP=$M.$E
      integer LTXTAA,LTXTAB
      parameter (LTXTAA= 1,LTXTAB= 8)
      character MTXTAA(1) * (95)
      data MTXTAA/'DIVAG$BCall with bad values of KORD.  KORD(1)=$I, KOR
     *D(2)=$I, when TSPECS(1)=$F and KSTEP=$M.$E'/
      data MACT / MEMDA1, 0, MEEMES, 68, 24, 0, MERET /
c
c ****************** START OF EXECUTABLE CODE **************
c
      IFLAG = 1
      IG = KORD(2)
      if ((IG .ne. 0) .and. (IG .ne. 1)) go to 500
      if (IGFLG - 3) 10, 80, 70
   10 if (IGFLG - 1) 20, 210, 60
c
c ******************** INITIAL POINT ***********************
c
   20 if (IGFLG .eq. -3) then
         IGFLG = 5
         return
      end if
      IGTYPE(IG + 1) = 0
      if ((NGSTOP(IG + 1) .le. 0) .or. (IG + IGFLG .eq. -1)) go to 30
      IGFLG = IG - 2
      go to 40
   30 IGFLG = 5
   40 NG = NGSTOP(2 - IG)
      do 50 I = 1, NG
   50    GT(I) = GNEW(I)
      go to 480
c
c     **** USER HAS BEEN TOLD THAT A GSTOP WAS FOUND
c          TEST IF CALLED FROM OUTPUT WHEN CALL SHOULD
c          BE FROM DERIVS
   60 if ((IGTYPE(1) .eq. 0) .and. (IG .ne. 0)) go to 420
c     **** PROTECT AGAINST NOISEY G NEAR THE ZERO
      if (GNEW(INGS) .eq. 0.D0) GNEW(INGS) = GT(INGS)
c
c ****** TEST FOR CHANGE IN THE SIGN OF A G ****************
c
   70 NG = NGSTOP(2 - IG)
      INGS = 0
   80 INGS = INGS + 1
      if (INGS .gt. NG) if (IGFLG - 4) 400, 380, 480
      if (GNEW(INGS)) 90, 100, 110
   90 if (GT(INGS)) 120, 120, 130
  100 if (GT(INGS)) 130, 80, 130
  110 if (GT(INGS)) 130, 120, 120
  120 GT(INGS) = GNEW(INGS)
      go to 80
c
c ********* A SIGN CHANGE HAS BEEN FOUND *******************
c
  130 NSTOP = INGS
      if (IG .eq. 0) NSTOP = -INGS
      if (IGFLG .ne. 5) go to 200
c     **** USUAL CASE -- TEST IF OUTPUT POINT PRECEDES THE
c          SIGN CHANGE, OR IF PREDICTING, CORRECTING, OR
c          NOT THROUGH THE FIRST STEP.
c     **** TEST IF AN INTERPOLATED G WAS WHAT CHANGED SIGN
      if (IG .ne. 0) go to 180
c     **** BACK UP DIFFERENCES AND OTHER STUFF TO BEGINNING
c          OF THE STEP
      call DIVABU(F, KORD)
c     **** TEST IF CORRECTING
      if (KORD1I .eq. 2) go to 170
c     **** TEST IF THROUGH THE FIRST STEP
      if (LSC .lt. 4) go to 180
c     **** IF FIRST DERIVATIVE EVALUATION OF THE FIRST
c          STEP, FIND THE GSTOP, AND USE IT TO GET A NEW
c          INITIAL STEPSIZE
      if (LSC .eq. 7) go to 200
c     **** SET NEW STEPSIZE AFTER SIGN CHANGE WHEN STARTING
  160 HH = TSPECS(1) - TN
c     **** SET KEXIT TO TRY NEW STEPSIZE
  170 KEXIT = 1
      TSPECS(1) = TN
      go to 460
c     **** SET KEXIT FOR USUAL CASE
  180 KEXIT = IG + 2
c     **** TEST IF SIGN CHANGE IN G PRECEDES NEXT OUTPUT PT.
      if (HH * (TSPECS(1) - TMARK)) 200, 200, 190
c     **** SET UP TO EVALUATE G AT OUTPUT POINT
  190 IGFLG = 4
      TSPECS(1) = TMARK
      NSTOP = 0
      go to 240
c
c ***************** FIND THE ZERO OF G *********************
c
c     **** INITIALIZE ZERO FINDER
  200 TOLD = TG(2 - IG)
      GOLD = GT(INGS)
      TSAVE = TSPECS(1)
      IGFLGS = IGFLG
      IGFLG = 1
      IZFLAG = 0
      go to 220
c     **** TEST IF ZERO ALREADY FOUND
  210 if (IZFLAG - 1) 350, 220, 310
  220 continue
      call DZERO(TSPECS(1), GNEW(INGS), TOLD, GOLD, IZFLAG, TOLG)
c     **** TEST FOR CONVERGENCE
      if (IZFLAG .ne. 1) go to 260
c     **** INTERPOLATE NEW Y, AND GO COMPUTE G AGAIN
  240 call DIVAIN(TSPECS(1), Y, F, KORD)
      IFLAG = 4
      KORD2I = IG - 3
      return
c     **** CONVERGENCE -- CHOOSE TOLD TO GIVE A CHANGE
c          IN SIGN
  260 if (GNEW(INGS) .eq. 0.D0) go to 290
      if (TSPECS(1) - TOLD) 270, 300, 280
  270 if (HH) 290, 300, 300
  280 if (HH) 300, 300, 290
  290 TOLD = TSPECS(1)
c     **** CHECK IF SIGN CHANGE DUE TO NOISE
  300 TSPECS(1) = TOLD + XI(1)
      go to 240
  310 TSPECS(1) = TOLD
      if (GNEW(INGS)) 320, 340, 330
  320 if (GT(INGS)) 340, 360, 360
  330 if (GT(INGS)) 360, 360, 340
c     **** ZERO WAS EVIDENTLY DUE TO NOISE
  340 TSPECS(1) = TSAVE
      IZFLAG = 0
      go to 370
  350 IGFLG = IGFLGS
c     SET KORD2I TO INITIAL VALUE TO AVOID LOOP
      KORD2I = IG
      go to 80
c     **** SAVE INFORMATION ABOUT THE STOP
  360 IGFLG = 3
      TGSTOP(2 - IG) = TSPECS(1)
      IGTYPE(2 - IG) = IZFLAG + 3
      IGSTOP(2 - IG) = NSTOP
  370 NSTOP = 0
      go to 240
c
c ************** AFTER SEARCH FOR A SIGN CHANGE ************
c
c     **** NO SIGN CHANGE AT A T OUTPUT POINT
c     TEST IF CALLED FROM OUTPUT
  380 if (IG .ne. 0) go to 390
c     SET UP FOR CALL TO OUTPUT
      KORD1I = 7
      KORD2I = -2
      IFLAG = 3
      go to 480
c     **** ADJUST KEXIT AND SET UP TO GIVE OUTPUT
  390 KEXIT = KEXIT + 2
      KORD1I = min(KMARK, 5)
      KORD(3) = KMARK
      KORD(1) = KORD1I
      IGFLG = 5
      IFLAG = 2
      go to 470
c     **** TEST IF USER HAS BEEN TOLD OF GSTOP
  400 if (IGFLG .eq. 2) go to 450
c     **** A GSTOP HAS BEEN FOUND
c     TEST IF STARTING
      if (LSC .eq. 7) go to 160
      IFLAG = IGTYPE(2 - IG)
      NSTOP = IGSTOP(2 - IG)
      INGS = abs(NSTOP)
      if (INGS .eq. 0) go to 410
      GT(INGS) = -GT(INGS)
      if (IG .eq. 0) go to 430
      IGFLG = 2
c     If interpolated GSTOP was found set to check again in case of
c     multiple stops at exactly the same point.
      if (IGTYPE(1) .ne. 0) go to 440
c     **** TELL USER OF AN EXTRAPOLATED GSTOP
  410 IFLAG = IGTYPE(2)
      NSTOP = IGSTOP(2)
      INGS = abs(NSTOP)
c     **** SET SO DERIVS IS CALLED WITH KORD(1) = KPRED
  420 KORD1I = KPRED
      KORD2I = -3
      return
c     **** AN EXTRAPOLATED GSTOP WAS FOUND, SET UP TO CHECK
c          INTERPOLATED STOPS (IF ANY)
  430 NG = NGSTOP(1)
      INGS = 0
      IFLAG = 3
      NSTOP = IGSTOP(2)
c     **** SET SO OUTPUT IS CALLED WITH KORD(1) = 7
  440 KORD1I = 7
      KORD2I = -2
      go to 490
c     **** CHECK THAT AN EXTRAPOLATED G-STOP IS NOT MISSED
  450 if ((IG .eq. 0) .or. (IGTYPE(2) .eq. 0)) go to 460
c     SET TO CHECK FOR INTERPOLATED G-S.
      TG(1) = TSPECS(1)
      IGTYPE(1) = 0
      TSPECS(1) = TGSTOP(2)
      INGS = 0
      IGFLG = 3
      go to 240
c     **** SET SO INTEGRATOR GOES TO PLACE DIRECTED BY KEXIT
  460 NSTOP = 0
      IGFLG = 5
      IFLAG = 3
  470 KORD2I = -7
c     **** STORE INFO. ON LAST G COMPUTED
  480 IGTYPE(2 - IG) = 0
  490 TG(2 - IG) = TSPECS(1)
      return
c
c ********************** ERROR PROCESSING ******************
c
  500 MACT(2) = KSTEP
c--D Next line special: P=>S, X=>D
      call DMESS(MACT, MTXTAA, KORD, TSPECS)
      IFLAG = 8
      KEXIT = 6
      go to 470
      end
      subroutine DMESS (MACT, TEXT, IDAT, FDAT)
c Copyright (c) 1996 California Institute of Technology, Pasadena, CA.
c ALL RIGHTS RESERVED.
c Based on Government Sponsored Research NAS7-03001.
c>> 2009-09-27 DMESS Krogh  Same as below, in another place.
c>> 2009-07-23 DMESS Krogh  Changed ,1x to :1x in write to FMTF.
c>> 2008-06-13 DMESS Krogh  Changed -0's to 0.
c>> 2007-09-08 DMESS Krogh  Fixed definitions of MEVLAS.
c>> 2006-10-08 DMESS Krogh  Another try, see 2005-05-26
c>> 2006-10-08 DMESS Krogh  Fixed minor problem in matrix/vector output.
c>> 2006-10-01 DMESS Krogh  Print NaN's and infity (at least for g77).
c>> 2006-07-01 DMESS Krogh  messxc => dmessxc (and not static) (for C)
c>> 2006-04-07 DMESS Krogh  Major rewrite of code for F.Pt. format.
c>> 2006-04-04 DMESS Krogh  Fixes in C code for vectors & matrices.
c>> 2006-04-02 DMESS Krogh  Added code for output of sparse vector.
c>> 2005-07-10 DMESS Krogh  Small adjustment for last correction.
c>> 2005-05-26 DMESS Krogh  Fixed "*****" output in boundary case.
c>> 2002-05-16 DMESS Krogh  Added way for user to get error count.
c>> 2002-03-27 DMESS Krogh  Fixed crash when number is -INF.
c>> 2001-06-08 DMESS Krogh  Eliminated Hollerith in formats.
c>> 2001-05-25 DMESS Krogh  Added a couple of commas in formats.
c>> 1997-06-17 DMESS Krogh  In C code made messxc, static.
c>> 1996-07-12 DMESS Krogh  Changes to use .C. and C%%.
c>> 1996-03-30 DMESS Krogh  Added external statement.
c>> 1994-10-20 DMESS Krogh  Changes to use M77CON
c>> 1994-09-21 DMESS Krogh  Added CHGTYP code.
c>> 1994-09-08 DMESS Krogh  Added new matrix/vector capabilities.
c>> 1994-08-17 DMESS Krogh  Removed duplicate save statement.
c>> 1994-04-19 DMESS Krogh  Removed blank line from DMESS.
c>> 1993-05-14 DMESS Krogh  Changed TEXT to array of character strings.
c>> 1993-04-14 DMESS Krogh  Fixes for conversion to C. (C%% comments.)
c>> 1992-07-12 DMESS Krogh  Fixed so negative KDFDEF works.
c>> 1992-05-27 DMESS Krogh  Initialized LDFDEF in a data statement.
c>> 1992-05-14 DMESS Krogh  Put common blocks in save statement.
c>> 1992-04-28 DMESS Krogh  Corrected minor error in floating pt. format
c>> 1992-02-28 DMESS Krogh  Initial Code.
c
c--D replaces "?": ?MESS,?MESSXC
c
c Processes Messages -- Actions are controlled by MACT().  See
c comment is subroutine MESS.  This program is for the extra
c argument of type real.
c
c BUF    In common CMESSC, see MESS.
c DOLS   In common for intitialization, not used here.  See MESS.
c EUNIT  In common for intitialization, not used here.  See MESS.
c FDAT   Formal argument -- gives floating point data to print.
c FBIG   Largest magnitude of floating point number to output.
c FSMA   Smalles magnitude floating point number to be printed.
c FMTF   In common CMESSC, format for printing floating point number.
c FMTG   In common CMESSC, user format to use in place of FMTF.
c FMTSP  Format for printing sparse vectors.
c FOUT   Floating point number to be output.
c FSMA   Smallest postitive floating point number.
c ICOL   In common CMESSI, see MESS.
c ID     Number of decimal digits for floating point format statement.
c IDAT   Integer data -- passed to MESS.
c IVAR   In common CMESSI, see MESS.
c IWF    In common CMESSI, see MESS.
c IWG    In common CMESSI, see MESS.
c J      Temporary index.
c K      Temporary index.
c KSMA   Number of leading 0's in "F" format.  If < 0, -KSMA gives the
c    number of extra digits to the left of the decimal point.
c    KSMA depends on abs(smallest number printed).
c KBIG   Number of extra digits before the decimal point required by the
c    largest number to be printed if "F" format is used.
c KDF    In common CMESSI, see MESS.
c KDFDEF In common CMESSI, see MESS.
c KDIAG  In common CMESSI, not used here, see MESS.
c KEXE   Extra space required for E format.
c KF     In common CMESSI, see MESS.
c KLINE  In common CMESSI, see MESS.
c KSCRN  In common CMESSI, see MESS.
c KRES1  In common CMESSI, see MESS.
c KSPEC  In common CMESSI, see MESS.
c LASTER In common CMESSI, not used here, see MESS.
c LASTI  In common CMESSI, see MESS.
c LBUF   In common CMESSI, see MESS.
c LDFDEF Value of KDFDEF for this routine.  (Saved)
c LENBUF In common CMESSI, see MESS.
c LENLIN In common CMESSI, not used here, see MESS.
c LENTRY In common CMESSI, see MESS.
c LHEAD  In common CMESSI, not used here, see MESS.
c LINERR In common CMESSI, not used here, see MESS.
c LINMSG In common CMESSI, not used here, see MESS.
c LOCBEG In common CMESSI, see MESS.
c LPRINT In common CMESSI, not used here, see MESS.
c LSTOP  In common CMESSI, not used here, see MESS.
c LSTRT  In common CMESSI, see MESS.
c MACT   Formal argument, see MESS.
c MDAT   In common CMESSI, not used here, see MESS.
c MEMDA5 In common CMESSI, see MESS.
c MESS   Program called for most of the message processing.
c MPT    In common CMESSI, see MESS.
c MUNIT  In common CMESSI, not used here, see MESS.
c NCOL   In common CMESSI, see MESS.
c NDIM   In common CMESSI, see MESS.
c NEG    1 if any number is < 0, else it is 0.
c NFDAT  In common CMESSI, see MESS.
c NIDAT  In common CMESSI, not used here, see MESS.
c NMDAT  In common CMESSI, not used here, see MESS.
c NTEXT  In common CMESSI, not used here, see MESS.
c OUNIT  In common CMESSI, not used here, see MESS.
c D1MACH External func. giving floating pt. info. about the environment.
c SUNIT  In common CMESSI, not used here, see MESS.
c TEXT   Formal argument, passed to MESS, see there.
c XARGOK In common CMESSI, see MESS.
c
C%% void dmessxc(long int);
c
      external         D1MACH
      integer          MACT(*), IDAT(*)
      double precision FDAT(*)
      character        TEXT(*)*(*)
      character        FMTSP*29
      integer          ICOL, ID, J, K, KSMA, KBIG, KEXE, LDFDEF, NEG
      double precision FBIG, FOUT, FSMA, D1MACH
      save LDFDEF, FMTSP
      save /CMESSI/, /CMESSC/
c++ CODE for .C. is inactive
c      integer  kciwid, kccwid, kcrwid, lbeg, lend, lfprec, lgprec
c      common /MESSCC/ kciwid,kccwid,kcrwid,lbeg,lend,lfprec,lgprec
c++ END
c
c ************************** Data from common block ********************
c
c For comments on these variables, see the listing for MESS.
c
      integer   LENBUF, MEVBAS, MEVLAS
      parameter (LENBUF=250)
      parameter (MEVBAS=10)
      parameter (MEVLAS=33)
      logical          XARG, GOTFMT, XARGOK
      integer          ERRCNT, EUNIT, ICHAR0, IRC, IVAR(MEVBAS:MEVLAS),
     1   IMAG, INC, ITEXT, IWF, IWG, KDF, KDFDEF, KDI, KDIAG, KDJ,
     2   KLINE, KSCRN, KSHIFT, KSPEC, KT, MAXERR, LASTI, LBUF, LENLIN,
     3   LENOUT, LENTRY, LENTXT, LHEAD, LINERR, LINMSG, LOCBEG, LPRINT,
     4   LSTOP, LSTRT, LTEXT, MAXWID(2), MDAT(5), MPT, MUNIT, NCOL,
     5   NDIM, NFDAT, NIDAT, NMDAT, NROW, NTEXT, OUNIT, SUNIT, TABSPA
c
      character BUF*(LENBUF), DOLS*72, FMTC*7, FMTF*20, FMTG*15,
     1  FMTI*7, FMTIM(2)*7, FMTJ*7, FMTR*7, FMTT*15
      common /CMESSI/ SUNIT, LHEAD, KDFDEF, LINMSG, LINERR, MUNIT,
     1   EUNIT, KSCRN, KDIAG, MAXERR, LSTOP, LPRINT, KDF, NTEXT, NIDAT,
     2   NFDAT, NMDAT, MDAT, TABSPA, ERRCNT, ICHAR0, IMAG, INC, IRC,
     3   ITEXT, IWF, IWG, KDI, KDJ, KLINE, KSHIFT, KSPEC, KT, LASTI,
     4   LBUF, LENLIN, LENOUT, LENTRY, LENTXT, LOCBEG, LSTRT, LTEXT,
     5   MAXWID, MPT, NROW, NCOL, NDIM, OUNIT, GOTFMT, XARG, XARGOK
      common /CMESSC / BUF, DOLS, FMTF, FMTG, FMTI, FMTJ, FMTT, FMTIM
      equivalence (IVAR(MEVBAS), SUNIT)
      equivalence (FMTIM(1), FMTR), (FMTIM(2), FMTC)
c
      data LDFDEF / 0 /
c

c ************************* Start of Executable Code *******************
c
      XARGOK = .true.
      if (LDFDEF .eq. 0) then
         LDFDEF = 1 - int(log10(d1mach(3)))
      end if
      KDFDEF = LDFDEF
      KDF = KDFDEF
   10 call MESS (MACT, TEXT, IDAT)
c             4    5    6    7    8    9   10   11
      go to (20, 100, 100, 200, 300, 400, 100, 500), LENTRY-3
      XARGOK = .false.
      LDFDEF = KDFDEF
      return
c                                      Print from FDAT
   20 J = LBUF + 1
      FOUT = FDAT(NFDAT)
      NFDAT = NFDAT + 1
      if (KSPEC .ge. 8) then
         LBUF = LBUF + IWG
C%% messcc.lend = cmessi.lbuf;
C%% cmessc.buf[messcc.lend] = ' ';
C%% if ((j > 1) && (cmessc.buf[j-2] >= '0') &&
C%%    (cmessc.buf[j-2] <= '9')) cmessc.buf[j++ - 1] = ' ';
C%% sprintf(&cmessc.buf[j-1], cmessc.fmtg, cmessi.iwg,
C%%    messcc.lgprec, fout);
C%% if (cmessc.buf[messcc.lend] != 0) {messcc.lbeg=j; dmessxc(kexe);}
         write (BUF(J:LBUF), FMTG) FOUT
         go to 10
      end if
      if (FOUT .le. 0.D0) then
        if (FOUT .eq. 0.D0) then
          FDAT(NFDAT-1) = 0.D0
          FOUT = 0.D0
        else
          NEG = 1
        end if
      else if (FOUT .gt. 0.D0) then
        NEG = 0
      else
c               Must be a Nan
        NEG = 0
        FBIG = 1.0
        FSMA = 1.0
        IWF = 2
        go to 40
      end if

      FBIG = abs(FOUT)
      FSMA = FBIG
      IWF = 2
c                                      Get the format.
   40 continue
      if (KDF .eq. 0) KDF = KDFDEF
      KEXE = 0
      if (FBIG .ne. 0.D0) then
        if (FSMA .eq. 0.D0) FSMA = 1.D0
        FBIG = FBIG * (1.D0 + .5D0 * .1D0**abs(KDF))
        IWF = IWF + NEG
        if (KDF .lt. 0) then
          KSMA = 0
        else
          KSMA = -log10(FSMA)
          if (FSMA .lt. 1.D0) KSMA = KSMA + 1
        end if
        KBIG = log10(FBIG)
        if (FBIG .lt. 1.D0) then
          KBIG = KBIG - 1
          if (FBIG .gt. 1.D0 - .1D0**abs(KDF-1)) KBIG = 0
        end if
c         This is to get ininities out (at least with g77)
        if ((KBIG .lt. -1000) .or. (KBIG .gt. 1000)) KBIG = 8
        if ((KSMA .lt. -1000) .or. (KSMA .gt. 1000)) KSMA = 8
        if (max(KBIG, 0) + max(KSMA,0) .ge. 4) then
c               Want to use an "E" format
          KEXE = 3 + max(0, int(log10(dble(max(KBIG,abs(KSMA))+1.D-5))))
          if (KDF .lt. 0) then
            ID = -KDF
          else
            ID = KDF - 1
          end if
          IWF = IWF + ID + KEXE
c++ CODE for ~.C. is active
          if (LENTRY .eq. 10) IWF = IWF - 1
          write (FMTF, '(''(1P,99(E'',I2,''.'',I2,''E'',I1,'':1X))'')')
     1      IWF, ID, KEXE - 2
c++ CODE for .C. is inactive
cc WATCOM C and others (?) man need an extra 1 here??
c%%    strcpy(cmessc.fmtf, "%*.*E ");
c      lfprec = id
c++ END
          go to 60
        end if
      else
        KSMA = 1
        KBIG = -1
      end if
c               Want to use an "F" format
      if (KDF .lt. 0) then
        ID = -KDF
      else
        ID = KDF + KSMA - 1
      end if
c++ CODE for ~.C. is active
      IWF = IWF + ID + max(KBIG, -1)
      write (FMTF, '(''(0P,99(F'',I2,''.'',I2,'':1X))'')') IWF,ID
c++ CODE for .C. is inactive
c      IWF = IWF + ID + max(KBIG, 0)
c%%    strcpy(cmessc.fmtf, "%*.*f ");
c      lfprec = id
c++ END
   60 if (LENTRY .ne. 4) then
        IWF = IWF + 1
        if (LENTRY .ne. 10) go to 10
c               Take care of setup for sparse vector
        IMAG = 0
        do 70 J = LOCBEG, LASTI
          IMAG = max(abs(IMAG), IDAT(J))
 70     continue
        call MESSFI
c  Format forms:     12345678901234567890   123456789012345678  1234567
c                    (1P,99(Edd.ddEd:1X))   (0P,99(Fxx.xx:1X))  (99Idd)
c++ CODE for ~.C. is active
        if (FMTF(8:8) .eq. 'F') then
          FMTSP=
     1      '(99(' // FMTI(4:6) // ''') '',0P,' // FMTF(8:18)
        else
          FMTSP=
     1     '(99(' // FMTI(4:6) // ''')'' ,1P,' // FMTF(8:20)
        end if
c++ CODE for .C. is inactive
Cc Using cmessc.fmtf in place of fmtsp
C%%      if (cmessc.fmtf[4] == 'f') {
C%%      strcpy(cmessc.fmtf, "%*ld) %*.*f ");
C%%      kexe = 0;
C%%      }
c%%      else strcpy(cmessc.fmtf, "%*ld) %*.*E ");
c%%      cmessi.iwf++;
c++ END
        IWF = IWF + KDI + 1
        go to 10
      end if
c
      LBUF = LBUF + IWF
C%% messcc.lend = cmessi.lbuf;
C%% cmessc.buf[messcc.lend] = ' ';
C%% if ((j > 1) && (cmessc.buf[j-2] >= '0') &&
C%%    (cmessc.buf[j-2] <= '9')) cmessc.buf[j++ - 1] = ' ';
C%% sprintf(&cmessc.buf[j-1], cmessc.fmtf, cmessi.iwf,
C%%    messcc.lfprec, fout);
C%% if (cmessc.buf[messcc.lend] != 0) {messcc.lbeg=j; dmessxc(kexe);}
      write (BUF(J:LBUF),FMTF) FOUT
      go to 10
c                                     Get format for a vector or matrix
  100 ICOL = 1
      if (FDAT(LOCBEG) .lt. 0.D0) then
        NEG = 1
      else if (FDAT(LOCBEG) .ge. 0.D0) then
        NEG = 0
      else
c               Must be a Nan
        NEG = 0
        FBIG = 1.0
        FSMA = 1.0
        go to 110
      end if

      FBIG = abs(FDAT(LOCBEG))
      FSMA = FBIG
  110 do 120 J = LOCBEG, LASTI, INC
        if (FDAT(J) .le. 0.D0) then
          if (FDAT(J) .eq. 0.D0) then
            FDAT(J) = 0.D0
          else
            NEG = 1
          end if
        end if
        FBIG = max(FBIG, abs(FDAT(J)))
        if (FSMA .eq. 0.D0) then
          FSMA = abs(FDAT(J))
        else if (FDAT(J) .ne. 0.D0) then
          FSMA = min(FSMA, abs(FDAT(J)))
        end if
  120 continue
      if (NCOL .ne. 0) then
         ICOL = ICOL + 1
         LOCBEG = LOCBEG + NDIM
         LASTI = LASTI + NDIM
         if (ICOL .le. NCOL) go to 110
      end if
      IWF = 2
      go to 40
c                                    Floating point vector output
  200 continue
C%% messcc.lend = cmessi.lstrt-1;
C%% neg = 0;
C%% for (j=cmessi.mpt; j<cmessi.mpt+cmessi.kline; j++){
C%%   messcc.lbeg = messcc.lend;
C%%   messcc.lend = messcc.lend + cmessi.iwf;
C%%   if (kexe) {
C%%     if (kexe == 5)
C%%        neg = ((fabs(fdat[cmessi.inc*j-1]) < 1.e100)) ? -1: 0;
C%%     else if (kexe ==3) neg = 1;
C%%   }
C%%   sprintf(&cmessc.buf[messcc.lbeg], cmessc.fmtf,
C%%    cmessi.iwf+neg, messcc.lfprec, fdat[cmessi.inc*j-1]);
C%%   if ((kexe == 3) || ((kexe == 5) && neg)) dmessxc(kexe); }
      write(BUF(LSTRT:LBUF),FMTF)(FDAT(K),K=MPT,MPT+INC*(KLINE-1),INC)


c      print '(/A/)', BUF(1:LBUF)

      MPT = MPT + KLINE * INC
      go to 10
c                                    Floating point matrix output
  300 continue
C%% messcc.lend = cmessi.lstrt-1;
C%% neg = 0;
C%% for (j=cmessi.mpt; j<=cmessi.lasti; j+=cmessi.ndim){
C%%    messcc.lbeg = messcc.lend;
C%%    messcc.lend = messcc.lend + cmessi.iwf;
C%%   if (kexe) {
C%%     if (kexe == 5) neg = ((fabs(fdat[j-1]) < 1.e100)) ? -1: 0;
C%%     else if (kexe ==3) neg = 1;
C%%   }
CC%%    if ((messcc.lbeg > 1) && (cmessc.buf[messcc.lbeg-1] >= '0') &&
CC%%       (cmessc.buf[messcc.lbeg-1] <= '9'))
CC%%       cmessc.buf[messcc.lbeg++] = ' ';
C%%    sprintf(&cmessc.buf[messcc.lbeg],
C%%       cmessc.fmtf, cmessi.iwf+neg, messcc.lfprec, fdat[j-1]);
C%%    if ((kexe == 3) || ((kexe == 5) && neg)) dmessxc(kexe); }
      write (BUF(LSTRT:LBUF), FMTF) (FDAT(K), K = MPT, LASTI, NDIM)
      go to 10
c                                    Table output
  400 continue
C%% messcc.lend = cmessi.lstrt-1;
C%% neg=0;
C%% for (j=cmessi.mpt; j<cmessi.mpt+cmessi.kline; j++){
C%%    messcc.lbeg = messcc.lend;
C%%    messcc.lend = messcc.lend + cmessi.iwf;
CC%%    if ((messcc.lbeg > 1) && (cmessc.buf[messcc.lbeg-1] >= '0') &&
CC%%       (cmessc.buf[messcc.lbeg-1] <= '9'))
CC%%       cmessc.buf[messcc.lbeg++] = ' ';
C%%    sprintf(&cmessc.buf[messcc.lbeg], cmessc.fmtf,
C%%       cmessi.iwf+neg, messcc.lfprec, fdat[j-1]); }
      write (BUF(LSTRT:LBUF), FMTF) (FDAT(K), K = MPT, MPT+KLINE-1)
      go to 10


c                                   Sparse vector output
 500  continue
C%%  messcc.lend = -1;
C%%  neg = 0;
C%%  for (j=cmessi.mpt; j<cmessi.mpt+cmessi.kline; j++) {
C%%    messcc.lbeg = messcc.lend + 1;
C%%     messcc.lend = messcc.lend + cmessi.iwf;
C%%   if (kexe) {
C%%     if (kexe == 5) neg = ((fabs(fdat[j-1]) < 1.e100)) ? -1: 0;
C%%     else if (kexe == 3) neg = 1;
C%%   }
C%%   sprintf(&cmessc.buf[messcc.lbeg], cmessc.fmtf, cmessi.kdi,
C%%     idat[j-1],cmessi.iwf-cmessi.kdi-2+neg,messcc.lfprec,fdat[j-1]);
C%%     if ((kexe == 3) || ((kexe == 5) && neg)) dmessxc(kexe); }
      write (BUF(1:LBUF), FMTSP) (IDAT(K),FDAT(K),K=MPT,MPT+KLINE-1)
      MPT = MPT + KLINE
      go to 10

      end

c%%  void dmessxc(long int kexe)
c%%{
c%%  /* Adjust for lack of control on digits in exponent */
c%%  char c;
c%% if (cmessc.fmtf[4] == 'f') return;
c%% if (kexe == 4) return;
c%% if (kexe == 3) { // Should only be one digit in the exponent
c%%   cmessc.buf[messcc.lend-1] = cmessc.buf[messcc.lend];
c%%   cmessc.buf[messcc.lend] = ' ';
c%% }
c%% else { // Should be at least 3 digits in the exponent.
c%%   c =cmessc.buf[messcc.lend-4];
c%%   if ((c < '0') || (c > '9')) {
c%%     cmessc.buf[messcc.lend-1] = cmessc.buf[messcc.lend-2];
c%%     cmessc.buf[messcc.lend-2] = cmessc.buf[messcc.lend-3];
c%%     cmessc.buf[messcc.lend-3] = '0';
c%%     cmessc.buf[messcc.lend] = ' ';
c%%   }
c%% }
c%% return;
c%%} /* end of function */
      subroutine DZERO(X1, F1, X2, F2, MODE, TOL)
c Copyright (c) 1996 California Institute of Technology, Pasadena, CA.
c ALL RIGHTS RESERVED.
c Based on Government Sponsored Research NAS7-03001.
c>> 2010-04-14 DZERO  Krogh  No discontinuity message if |F1-F2| small.
c>> 2010-04-12 DZERO  Krogh  Fixed KNKP to get discontinuity diagnostic.
c>> 2010-02-20 DZERO  Krogh  $G => $F for print out of iterations
c>> 2008-03-01 DZERO  Krogh  Minor change in diagnostic print.
c>> 2000-12-01 DZERO  Krogh  Removed unused variable C1P01.
c>> 1998-11-01 DZERO  Krogh  Set so errors stop less easily.
c>> 1998-11-01 DZERO  Krogh  For "mangle", INDIC replaced with MACT(3).
c>> 1996-03-30 DZERO  Krogh  Added external statement.
c>> 1995-11-09 DZERO  Krogh  Fixed so char. data at col. 72 is not ' '.
C>> 1994-11-11 DZERO  Krogh  Declared all vars.
c>> 1994-10-20 DZERO  Krogh  Changes to use M77CON
c>> 1994-09-08 DZERO  Krogh  Added CHGTYP code.
c>> 1993-04-27 DZERO  Krogh  Additions for Conversion to C.
c>> 1993-04-13 DZERO  Krogh  Minor change for new MESS.
c>> 1992-04-08 DZERO  Krogh  Unused label 400 removed.
c>> 1992-01-09 DZERO  Krogh  Moved calc. of XXMXO up (for error msg.)
c>> 1991-11-26 DZERO  Krogh  Converted to new error processor.
c>> 1988-08-14 DZERO  Krogh  Labels runumbered.
c>> 1988-03-07 DZERO  Krogh  Initial code.
c
c--D replaces "?": ?ZERO, ?MESS
c
c SUBROUTINE TO FIND A BOUNDED ZERO
c
c analysis and coding by Fred T.Krogh at the Jet Propulsion
c Laboratory, Pasadena, Calif.  April 25, 1972.
c Modified for portability, April 1984 by Krogh.
c Algorithmic changes, vars. added to save stmt., Sept. 1987 by Krogh
c
c Parameters in the calling sequence are defined as follows:
c
c  X1  = independent variable
c  F1  = dependent variable --  initially   F1=F(X1).
c        When MODE=1 (or 5) the user is to compute F(X1) given X1
c  X2  = second value of independent variable
c  F2  = F(X2) on the initial entry.  When MODE = 2-4, F2=F(X2) and
c        F1*F2 .le. 0.
c  MODE  is a parameter used for communication between this
c        subroutine and the user. (The user should set MODE
c        only to initialize it to 0 before the first call)
c      =1  compute F(X1) and call $ZERO
c      =2  F(X1) is approximately 0, the iteration is finished
c          and the error criterion is satisfied.
c      =3  same as MODE=2, except the error criterion can
c          not be satisfied.
c      =4  apparently the function has a discontinuity
c          between X1 and X2 -- No zero can be found
c      =5  F1*F2 was greater than zero on the first call, and an attempt
c          to bound the zero on both sides have failed.
c      =6  fatal error -- $ZERO was called after mode was set .ge.2.
c          If $ZERO is called again, the program will be stopped.
c          (Unless MODE is set to 0)
c      <0  If MODE is set <0 and $ZERO is called, no action is taken
c          except that print is turned on for -MODE calls to $ZERO.
c          This print gives all values of X and F used in the iteration.
c  TOL    is the error tolerance
c     TOL.GT.0  Iterate until values of X1 and X2 are known
c              for which abs(X1-X2) .le. tol and F1*F2 .le. 0.
c     TOL.LT.0  Iterate until a value of X1 is found for which
c              abs(F1) .le. abs(TOL).
c     TOL  = 0  Iterate until the zero is determined as
c              precisely as possible.  MODE = 3 is impossible
c              in this case.
c
c Parameters in the calling sequence have the following types
c
      integer MODE
      double precision X1, X2, F1, F2, TOL
c
c Usage is as follows (of course, variations are possible.)
c         Somehow one has available X1, F1, X2, and F2 such
c         that F1 = F(X1), F2 = F(X2) and F1*F2 .le. 0.
c         In addition, one should assign a value to TOL.
c     MODE = 0
c***  In the statement below, $ is replaced by an 'S' for single
c***  precision and a 'D' for double.
c XXX call $ZERO(X1,F1,X2,F2,MODE,TOL)
c     go to  (N1,N2,N3,N4,N5,N6), MODE
c  N1 COMPUTE  F1=F(X1)
c     go to XXX
c
c  N4 continue
c  N5 continue
c  N6 stop
c  N3 If you want to -- print results to note that error
c                       is too big.
c  N2 zero is found, do whatever you want to with it.
c
c End of comments explaining usage.
c
c ************************* Usage of internal variables ****************
c
c C0     Parameter = 0.
c C1     Parameter = 1.
c C1P25  Parameter = 1.25
c C2     Parameter = 2.
c C4     Parameter = 4.
c CP01   Parameter = 0.01
c CP125  Parameter = 1.25
c CP25   Parameter = 0.25
c CP5    Parameter = 0.5
c CP75   Parameter = 0.75
c CP99   Parameter = 0.99
c D1MACH Gets constants associated with floating point arithmetic.
c DFDXXX = (XXMXL/FFMFL) * (est. deriv. of f w.r.t. x at x = XX).  All
c   derivatives are base on a second degree polynonial that interpolates
c   the last three points generated.
c DFDXXX = (XXMXO/FFMFL) * (est. deriv. of f w.r.t. x at x = X0).
c DIV    If artificial subdivision of the interval is used, determines
c   the amount of the sudivision.  (-XXMXOO * DIV / (1. + DIV))
c DMESS  Prints error messages.
c DXDFFF = (FFMFL/XXMXL) * (est. deriv. of x w.r.t. f at f = FF).
c DXDFFO = (FFMFO/XXMXL) * (est. deriv. of x w.r.t. f at f = FO).
c F1     (formal arg.) The last value of F computed, on return the value
c   of F(X1).
c F2     (formal arg.) The other initial value provided for F.  Set to
c   the value of F(X2) on returns.
c FDAT   Temporary storage for floating point values for messages.
c FF     Value of F1 after F is computed.
c FFDFO  FF / FO
c FFMFB  FFMFL + FLMFB = FF - value of FF 2 iterations back.
c FFMFL  FF - FL
c FL     Value of FF from the previous iteration.
c FLMFB  Value of FFMFL from the previous iteration
c FO     F(XO)
c I      Comments for LCHNG define how I is set.
c IDAT   Temporary storage for integer values for messages.
c J      This is 1 if FF .le. 0., and is 2 if FF > 0.
c KNKP   KNKP(J) (see J above) is decreased by 3 when there are signs of
c   decent convergence.  It is counted up when convergence is slow.
c KS     =-1  initially,
c        = 0  whenever F changes sign, otherwise
c        = number of times the sign of F has remained the same
c KTYP   = 1 if interpolation was used to get the last iterate, = 0 if
c   an artificial subdivision was used.
c LCHG  the J-th continuation in the data statement for LCHG below gives
c new states for the case when the state number is J-1.  State 0 is the
c initial state.  The I-th entry on a row gives the state for case on I
c as follows:  (TP is the ratio (new f) / (last f of same sign)
c    I = 1   TP < 0.01
c    I = 2   .01 <= TP < 1
c    I = 3   TP = 1
c    I = 4   1 < TP <= 4
c    I = 5   TP > 4.
c States are as follows:
c    0   initial state, or big increase, or small increase in state 0
c    1   after big decrease, perhaps followed by small decreases
c    2   after one or more small decreases from state 0
c    3   one or more small increases from state 2
c    4   one or more small decreases from state 3
c    5   decision made that noise is a problem on this side of zero.
c LINIT  - the maximum number of iterations that can be taken with out
c   getting a sign change in F.
c LMODE  The value of MODE the last time in this routine.
c LNLP
c LTXTxx Names of this form are used in setting up data statements for
c   error messages.  These names are generated automatically by PMESS,
c   the program that makes up these messages.
c MACT   This array difines the actions to be taken by the error message
c   program.  See comments in MESS for details.  MODE is set to MACT(3)
c   on exit.
c MACT1  As for MACT except used for the diagnostic print.
c MExxxx Parameters defining constants used for interaction with the
c   error message program MESS.  See comments there for definitions.
c MLOC   Contains locations in MTXTAA for error messages.
c MODE   (formal) See comments above.
c MTXTAA Text for error messages.
c MTXTAB Text for diagnostic message.
c MTXTAC Text for diagnostic message.
c NP     If > 0, gives number of iterations till diagnostic print stops.
c QFM
c QXM
c RND    Largest relative difference between succesive floating point
c   numbers.
c SMALL  .5 / (RND * largest floating point number)
c TOL    (Formal) See description above.
c TOLX   Actually tolerance required for accuracy in X.  Usually =
c   max(TOL, XRND).  It can be cut by a factor of 2 for use in setting
c   bounds on an acceptable interval.
c TP     Ordinarily the ratio (FF / prev. FF of the same sign.
c TP1    Used for temporary storage.
c X1     (Formal) Value of x where F is to be computed, and value
c   returned for the zero after convergence.
c X2     (Formal) Initially other value of x where F is given.  After
c   convergence gives the other closest x which gives an F of opposite
c   sign from that given by x1.
c XL     Value of XX from the previous iteration.
c XLMXB  Value of XXMXL from the previous iteration.
c XO     Value of x on opposite side of the zero from the current x.
c XRND   Best accuracy that one could hope for based on the finite
c   precision of floating point numbers.
c XX     Current x, the last value of X1 where F was computed.
c XXMXL  XX - XL
c XXMXO  XX - XO = length of interval in which 0 lies.
c XXMXOL Value of XXMXO from a previous iteration.
c
      external D1MACH
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
      double precision D1MACH
c
      parameter (C0 = 0.D0, C1 = 1.D0, C2 = 2.D0, C4 = 4.D0)
      parameter (C8 = 8.D0)
      parameter (CP125 = 0.125D0, CP25 = 0.25D0, CP75 = 0.75D0)
      parameter (CP5 = 0.5D0)
      parameter (C1P25 = 1.25D0)
      parameter (CP01 = 0.01D0)
      parameter (CP001 = 0.001D0)
      parameter (CP99 = 0.99D0)
      parameter (C1P031 = 1.03125D0)
c
c                      Declarations for error message processing.
c
      integer MERET, MEEMES, METEXT
      double precision FDAT(4)
      integer MACT(5), MACT1(2), MLOC(4), IDAT(2)
      save DIV, FL, FLMFB, FO, KNKP, KS, KTYP, LCHG, LMODE,
     1   LNLP, MACT, NP, RND, SMALL, XL, XLMXB, XO, XX, XXMXOL
      parameter (MERET  =51)
      parameter (MEEMES =52)
      parameter (METEXT =53)
c ********* Error message text ***************
c[Last 2 letters of Param. name]  [Text generating message.]
cAA DZERO$B
cAB Best bound for zero is [$F, $F], but tolerance is $F.$E
cAC Apparent discontinuity in function near X = $F.$E
cAD Can not find a sign change: X1=$F, X2=$F, F1=$F, F2=$F$E
cAE Called with MODE = $I.$E
c   $
cAF In DZERO -- X1=$F F1=$F KTYP=$I DIV=$G KS=$I$E
c   $
cAG             X2=$F F2=$F$E
      integer LTXTAA,LTXTAB,LTXTAC,LTXTAD,LTXTAE,LTXTAF,LTXTAG
      parameter (LTXTAA=  1,LTXTAB=  8,LTXTAC= 63,LTXTAD=112,LTXTAE=169,
     * LTXTAF=  1,LTXTAG=  1)
      character MTXTAA(1) * (193)
      character MTXTAB(1) * (46)
      character MTXTAC(1) * (25)
      data MTXTAA/'DZERO$BBest bound for zero is [$F, $F], but tolerance
     * is $F.$EApparent discontinuity in function near X = $F.$ECan not$
     * find a sign change: X1=$F, X2=$F, F1=$F, F2=$F$ECalled with MODE$
     * = $I.$E'/
      data MTXTAB/'In DZERO -- X1=$F F1=$F KTYP=$I DIV=$F KS=$I$E'/
      data MTXTAC/'            X2=$F F2=$F$E'/
c **** End of automatically generated text
c                      1  2  3  4      5
      data MACT / MEEMES, 0, 0, 0, MERET /
      data MACT1 / METEXT, MERET /
      data MLOC / LTXTAB, LTXTAC, LTXTAD, LTXTAE /
c
      data RND / C0 /
      data KS, KTYP, LMODE, DIV / 0, 0, 2, C0 /
      data LCHG /
     1   1, 2, 0, 0, 0,
     2   1, 1, 4, 5, 0,
     3   1, 2, 3, 3, 0,
     4   1, 4, 4, 3, 0,
     5   1, 4, 5, 5, 0,
     6   1, 5, 5, 5, 0 /
      data NP / 0 /

c
c INITIALIZE
c
      if (MODE .lt. 0) then
         NP = -1 - MODE
         return
      end if
      if (NP .gt. 0) then
         NP = NP - 1
         FDAT(1) = X1
         FDAT(2) = F1
         FDAT(3) = DIV
         IDAT(1) = KTYP
         IDAT(2) = KS
         call DMESS(MACT1, MTXTAB, IDAT, FDAT)
         if (MODE .ne. 0) if (LMODE - 1) 70, 80, 450
         FDAT(1) = X2
         FDAT(2) = F2
         call DMESS(MACT1, MTXTAC, IDAT, FDAT)
      else if (MODE .ne. 0) then
         if (LMODE - 1) 70, 80, 450
      end if
c
      if (RND .eq. C0) then
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
c             Take care of points on same side of zero.
   70 FF = F1
      XX = X1
      TP = FF / FL
      if (TP .lt. C0) go to 30
      LMODE = LMODE - 1
      if (LMODE .lt. LINIT) then
         MACT(3) = 5
         FDAT(1) = XX
         FDAT(2) = XL
         FDAT(3) = FF
         FDAT(4) = FL
         go to 250
      end if
      if (TP .gt. C1) then
         FF = FL
         XX = XL
         FL = F1
         XL = X1
      end if
      if (abs(FF) .ge. C8 * abs(FL-FF)) then
         TP = C8
      else
         TP = max(-CP25*dble(LMODE), FF / (FL - FF))
      end if
      FL = FF
      XO = XL
      XL = XX
      if (XX .eq. XO) XO = C1P031 * XX + sign(CP001, XX)
      XX = XX + TP * (XX - XO)
      X1 = XX
      MODE = 1
      return
c
   75 X1 = XL
      F1 = FL
      go to 250
c END OF INITIALIZATION
c
c
c ENTRY AFTER COMPUTING F FOR THE LAST ITERATE
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
      if (FF .gt. C0) J = 2
      if (TP - C1) 150, 140, 130
  130 I = 4
      if (TP .gt. C4) I = 5
      go to 160
  140 I = 3
      go to 160
  150 I = 2
      if (TP .lt. CP01) I = 1
      if (TP .lt. CP99) go to 170
  160 KNKP(J) = KNKP(J) + 1
      go to 180
  170 KNKP(J) = 0
  180 XXMXO = XX - XO
      LNLP(J) = LCHG(5*LNLP(J) + I)
      if (LNLP(J) .ge. 4) then
         if (LNLP(3 - J) .ge. 4) go to 210
      end if
c XXMXO GIVES THE LENGTH OF THE INTERVAL INSIDE WHICH
c THE ZERO IS KNOWN TO LIE.
      if (C2 * abs(XXMXO) .lt. abs(XXMXOL)) then
         KNKP(J) = max(0, KNKP(1) - 3)
      end if
      XXMXOL = XXMXO
      XRND = RND * (abs(XX) + abs(XO) + SMALL)
c
c TEST FOR CONVERGENCE
      if (TOL) 190, 200, 200
  190 continue
      if (abs(FF) .le. abs(TOL)) go to 220
  200 continue
      TOLX = max(TOL, XRND)
      if (abs(XXMXO) .gt. TOLX) go to 310
c
c CONVERGENCE -- PREPARE FOR FINAL EXIT
  210 if ((abs(XXMXO) .gt. TOL) .and. (TOL .ne. C0)) then
         MACT(3) = 3
         FDAT(3) = TOL
         if (XXMXO .gt. 0) then
            FDAT(2) = XX
            FDAT(1) = XO
         else
            FDAT(1) = XX
            FDAT(2) = XO
         end if
      end if
c SET FINAL VALUES FOR X1,F1,X2,AND F2
  220 continue
      if (abs(FF) .le. abs(FO)) go to 240
      F1 = FO
      X1 = XO
  230 FO = FF
      XO = XX
  240 X2 = XO
      F2 = FO
c TEST FOR DISCONTINUITY
      if ((KNKP(1) .gt. 5) .or. (KNKP(2) .gt. 5)) then
        if (abs(F1 - F2) .gt. RND * max(X1, 1.D0)) then
          MACT(3) = 4
          FDAT(1) = XX
        end if
      end if
  250 MODE = MACT(3)
      if (MACT(3) - 2) 420, 420, 430
c END OF CODE FOR FINAL EXIT
c
c F NOT DECREASING (OR THE FIRST ITERATE)
c PREPARE TO DIVIDE THE INTERVAL
  260 TP = C1
      if (KS) 370, 280, 270
  270 if (KTYP .eq. 0) go to 290
  280 DIV = C2
  290 continue
      DIV = max(DIV, FFDFO)
c KTYP=0 IF AND ONLY IF THE INTERVAL WAS DIVIDED (USING DIV)
c ON THE LAST ITERATION
      if (KTYP .eq. 0) DIV = DIV * (C1P25 / (C1P25 - TP))
c DIVIDE THE INTERVAL AS SPECIFIED BY DIV
  300 TP1 = -XXMXO * (DIV/(DIV+C1))
      KTYP = 0
      go to 410
c
  310 continue
      XXMXL = XX - XL
      FFMFL = FF - FL
      FFDFO = abs(FF / FO)
      TOLX = CP5 * TOLX
      if (TP .ge. C1) go to 260
c DIVIDE THE INTERVAL IF F HAS HAD THE SAME SIGN FOR
c FOUR OR MORE TIMES IN SUCCESSION
      if (KS - 4) 320, 340, 290
  320 continue
      if (FLMFB .eq. C0) go to 340
c BEGINNING OF CODE TO DETERMINE IF INVERSE QUADRATIC
c INTERPOLATION IS TO BE USED.
      FFMFB = FFMFL + FLMFB
      if (FFMFB .eq. C0) go to 330
      QFM = C1 - (FFMFL / FLMFB) * (XLMXB / XXMXL)
      QXM = C1 - (XXMXL / XLMXB) * (FLMFB / FFMFL)
      DXDFFF = C1 + (FFMFL / FFMFB) * QFM
      DXDFFO = DXDFFF + C2 * ((FO - FF) / FFMFB) * QFM
      TP1 = XXMXL + XLMXB
      DFDXXX = C1 + (XXMXL / TP1) * QXM
      DFDXXO = DFDXXX + C2 * ((XO - XX) / TP1) * QXM
      TP1 = DXDFFF * DFDXXX
      if ((TP1 .le. CP25) .or. (TP1 .ge. C4)) go to 330
      TP1 = DXDFFO * DFDXXO
      if ((TP1 .gt. CP25) .and. (TP1 .lt. C4)) go to 380
c
c DERIVATIVES DO NOT MATCH WELL ENOUGH
  330 continue
      if (KS .eq. 0) if (FFDFO - C1) 350, 370, 360
  340 continue
      if ((KTYP .eq. 0) .and. (TP .ge. CP75)) go to 290
      continue
      TP = C1 - TP
      if (TP .le. FFDFO) go to 280
      FFDFO = FFDFO / TP
      DIV = CP125
      go to 290
  350 continue
      DIV = CP5 * max(max(CP25, FFDFO), TP / (C1P25 - min(TP, C1)))
      go to 300
  360 continue
      DIV = min(C4, CP5 * FFDFO)
      go to 300
c INTERPOLATE WITH SECANT METHOD
  370 TP1 = -XXMXL
      go to 390
c
c DERIVATIVES MATCH UP PRETTY WELL.
  380 continue
c INTERPOLATE USING THE INVERSE QUADRATIC
      TP1 = XXMXL * (QFM * (FL / FFMFB) - C1)
  390 TP1 = (FF/FFMFL) * TP1
      KTYP = 1
c
c EXIT TO GET F(X)
  410 continue
      FL = FF
      FLMFB = FFMFL
      XLMXB = XXMXL
      XL = XX
c COMPUTE X1, INSURING THAT IT IS NOT TOO CLOSE TO THE
c ENDS OF THE INTERVAL
      XX = min(max(XL + TP1, min(XL, XO) + TOLX), max(XL, XO) - TOLX)
      X1 = XX
  420 LMODE = MODE
      return
c
  430 MACT(2) = 11*MACT(3)  - 19
  440 MACT(4) = MLOC(MACT(3)-2)
      call DMESS(MACT, MTXTAA, IDAT, FDAT)
      go to 420
c
c A CALL TO THE SUBROUTINE HAS BEEN MADE WITH MODE.NE.1
  450 IDAT(1) = MODE
      MACT(3) = 6
      MODE = 6
      if (LMODE .ne. 6) go to 430
      MACT(2) = 99
      go to 440
      end
c++ CODE for .C. is inactive
C%% static FILE *c_handle[2], *scratch_file;
C%% static char *c_fname[2]={"MESSF-xx", "MESSF-xx"};
C%% char *ctmp;
c++ END
      subroutine MESS(MACT, TEXT, IDAT)
c Copyright (c) 1996 California Institute of Technology, Pasadena, CA.
c ALL RIGHTS RESERVED.
c Based on Government Sponsored Research NAS7-03001.
c>> 2010-02-22 MESS  Krogh  Moved NSKIP=0 to start of code.
c>> 2009-10-30 MESS  Krogh  Defined DSCRN.
c>> 2009-02-28 MESS  Krogh  Added FMTT = ' ' for NAG compiler.
c>> 2009-02-28 MESS  Krogh  Fixed "f" format for C code.
c>> 2007-09-08 MESS  Krogh  Fixed definitions of MEVLAS.
c>> 2006-07-27 MESS  Krogh  Fixed boundary case in printing long text.
c>> 2006-03-20 MESS  Krogh  Added code for output of sparse vector.
c>> 2005-04-07 MESS  Krogh  Declared LFLGDB integer in MESSMH.
c>> 2004-12-15 MESS  Krogh  Added " - 1" at end of line on label 410.
c>> 2002-05-17 MESS  Krogh  Added way for user to get error count.
c>> 2001-12-28 MESS  Krogh  Added NSKIP for more flexible output values.
c>> 2000-12-30 MESS  Krogh  Fixed some types/casts in C code.
c>> 1997-12-12 MESS  Krogh  Prefixed 0P edit descriptor to F format.
c>> 1996-07-11 MESS  Krogh  Transpose matrix output for C.
c>> 1996-06-27 MESS  Krogh  fprintf(stdout, => printf( & memset now used
c>> 1996-06-18 MESS  Krogh  "Saved" NTEXTR.
c>> 1996-05-15 MESS  Krogh  Changes to use .C. and C%%.
c>> 1996-03-30 MESS  Krogh  Added external statement.
C>> 1996-01-24 MESS  Krogh  Fixed minor bug introduced with "$ " stuff.
C>> 1996-01-23 MESS  Krogh  Minor changes for C conversion.
C>> 1995-11-10 MESS  Krogh  Add code to change "$ " to " " in headings.
C>> 1995-08-11 MESS  Krogh  Made code default to not using UMESS.
C>> 1995-01-20 MESS  Krogh  Fixed unusual case in matrix output.
C>> 1994-12-15 MESS  Krogh  Removed block data for Cray T3D.
C>> 1994-11-11 MESS  Krogh  Declared all vars.
c>> 1994-09-14 MESS  Krogh  Fixed to get 1 more "$" in C output.
c>> 1994-09-08 MESS  Krogh  Added new matrix/vector capabilities.
c>> 1994-08-22 MESS  Krogh  Fix for conversion to C for new converter.
c>> 1994-07-05 MESS  Krogh  Fixed bug, KDI and FMTI could be inconsist.
c>> 1994-05-20 MESS  Krogh  Changes to MESSFT so line 1 can go to file.
c>> 1994-05-20 MESS  Krogh  Changes to setting output unit.
c>> 1994-05-09 MESS  Krogh  Integer vectors had overflow & space probs.
c>> 1993-05-19 MESS  Krogh  Changed TEXT to array of character strings.
c>> 1993-04-14 MESS  Krogh  Fixes for conversion to C. (C%% comments.)
c>> 1993-03-10 MESS  Krogh  Broke into smaller pieces.
c>> 1992-12-02 MESS  Krogh  Added save statement to block data subpr.
c>> 1992-07-13 MESS  Krogh  Add checks in heading set up.
c>> 1992-07-12 MESS  Krogh  Fixed so $$ prints a single $ in TEXT.
c>> 1992-07-12 MESS  Krogh  Set out of bound inputs to limit values.
c>> 1992-07-12 MESS  Krogh  Fixed so output works to alternate files.
c>> 1992-07-12 MESS  Krogh  Added integer declarations for parameters.
c>> 1992-06-24 MESS  Krogh  More blanks allowed on break of long lines.
c>> 1992-06-10 MESS  Krogh  Minor fix to vector output.
c>> 1992-05-27 MESS  Krogh  Fixed bug on line width setting.
c>> 1992-05-14 MESS  Krogh  Put common blocks in save statement.
c>> 1992-05-11 MESS  Krogh  Added label to assigned go to & a comment.
c>> 1992-04-08 MESS  Krogh  Unused labels 60, 220 and 320 removed.
c>> 1992-03-20 MESS  Krogh  Changed status on open to SCRATCH.
c>> 1992-03-10 MESS  Krogh  1 line below label 690 changed max to min.
c>> 1992-02-05 MESS  Krogh  Fixed bugs in printing matrix labels.
c>> 1992-01-29 MESS  Krogh  Added UMESS and multiple print option.
c>> 1991-12-09 MESS  Krogh  Fine tuning of vector output.
c>> 1991-10-10 MESS  Krogh  Insure no stop if stop level = 9.
c>> 1991-06-26 MESS  Krogh  Initial Code.
c Processes Messages -- Actions are controlled by MACT().
c This routine is intended for use primarily by other library routines.
c Users of library routines may want to use values of MACT from MERET-
c MESUNI, and may have an interest in using it to print messages
c from their own software.
c This routine has companion routines that are called with the same
c three arguments, plus one additional argument.  This argument is
c referred to here as FDAT since actions specified here can result
c in returns to print data from FDAT.  The name FDAT is used because
c this other routine will ordinarily print floating point data, but
c it could also print other kinds of data, e.g. logical.  At present
c only SMESS and DMESS are defined which are for single and double
c precision floating point data.
c MACT is a vector specifying sequentially the actions desired.
c Some of these actions require more than one location, in which
c case the next action follows the last datum required by the
c previous action.  Internal variables together with default
c values in parentheses which are used to keep track of locations
c are as follows:
c  NTEXT  (1)   The next text output starts at TEXT(NTEXT).
c  NIDAT  (1)   The next output from IDAT starts at IDAT(NIDAT).
c  NFDAT  (1)   The next output from FDAT starts at FDAT(NFDAT).
c  NMDAT  (1)   The next output from MDAT starts at MDAT(NMDAT), where
c               MDAT is defined by actions MEMDA1-MEMDA5 below, and
c               NMDAT is set to one at end of every text output.
c An action which uses data pointed to by one of the above will cause
c the pointer to be incremented to one past the last location used.  An
c exception is NMDAT which when it reaches 5 is not incremented and the
c value pointed to is incremented instead.
c Actions are encoded by values starting in MACT(1) as follows.
c (Arguments required are given in parentheses at the start of
c description.  These arguments follow the action index.  The next
c action follows the last argument for the preceding action.  Action
c indices have been selected so that it is easy to add new functionality
c without affecting codes using an earlier version.  Where bounds are
c indicated for an argument, if the argument is outside the bounds it is
c treated as if it had the value for the bound violated.)
c MESUNI=10  (0 .le. K10 .le. 99) Set the unit to use for a scratch
c            file.  The default unit for a scratch file is 30.  If a
c            scratch file is needed, (only needed here if a table
c            exceeds the line length), and unit 30 can not be opened as
c            a new scratch file, then units 29, 28, ..., will be tried
c            until an acceptable unit is found.  Library routines may
c            use this file but must be sure that the use does not
c            conflict with the printing of tables here, or the use by
c            any other library routines.  If K10 is 0, a scratch unit is
c            assumed not to be available, and tables with long lines
c            will be printed with each line on multiple lines.
c MEHEAD=11  (0 .le. K11 .le. 1) Defines the print that surrounds an
c            error message.  K11=0 gives nothing, and 1 gives the first
c            4 characters in TEXT repeated 18 times.  If this is not
c            used, one gets 72 $'s.  (To get a blank line use 1 with
c            TEXT = '    '.)
c MEDDIG=12  (-50 .le. K12 .le. 50) Set default digits to print for
c            floating point.  If K12 > 0 then K12 significant digits
c            will be printed, if K12 < 0, then -K12 digits will be
c            printed after the decimal point, and if K12 = 0, the
c            default will be used, which is the full machine precision.
c            Setting or getting this value will only work properly if
c            the action is taken by calling SMESS or DMESS as
c            appropriate.
c MEMLIN=13  (39 .le. K13 .le. 500) Set message line length to K13.
c            (Default is 128.)
c MEELIN=14  (39 .le. K14 .le. 500) Set error message line length to
c            K14. (Default is 79)
c MEMUNI=15  (-99 .le. K15 .le. 99) Messages go to unit K15.  If K15 = 0
c            (default), 'print' is used.  If K15 < 0, messages go to
c            both 'print' and to unit abs(K15).  If a write can not be
c            done to unit abs(K15), this unit will be opened with file
c            name MESS_Fxx.tmp, where xx is the value of abs(K15).
c MEEUNI=16  (-99 .le. K16 .le. 99) As for MEMUNI, except for Error
c            Messages.
c MESCRN=17  (0 .le. K17 .le. 100000000) Set number of lines to print to
c            standard output before pausing for "go" from user.  Default
c            is 0, which never stops.
c MEDIAG=18  (0 .le. K18 .le. 1000000000) Set the diagnostic level
c            desired.  This routine makes no use of K18.  It merely
c            serves as a place to set it and to answer inquiries on its
c            value.  It is intended to be set by users of library
c            software.  Library packages that make use of this number
c            are expected to use it as described below.  If K18 = 0 (the
c            default), no diagnostic print is being requested.  Else m =
c            mod(K18, 256) determines whether a package will do
c            diagnostic printing.  Associated with a library package is
c            a number L which must be a power of 2 < 129, and which
c            should be mentioned in the documentation for the package.
c            If the bit logical or(m,L) = L then diagnostic output for
c            the routine with the associated value of L is activated.
c            The value of L should have been selected by the following
c            somewhat vague rules.  Let base 2 log(L) = 2*i + j, where j
c            is 0 or 1.  Select i = level of the library package, where
c            the level is 0 if no other library routine that is likely
c            to be used with the package could reasonably be expected to
c            want any embedded diagnostics, and otherwise is
c            min(4, I+1), where I is the maximum level for any library
c            routine which is likely to be used with the package.
c            Select j = 0 if the user is relatively unlikely to want
c            diagnostics, and j = 1, if this is a routine for which
c            considering its level the user is relatively likely to want
c            diagnostic output.  The next 8 bits, mod(K18/256, 256), may
c            be used by the library routine to select the actual output
c            that is to be given.  These bits may be ignored,  but if
c            they are used, the lowest order bits should correspond to
c            less voluminous output that is more likely to be requested.
c            Finally, K18 / (2**16) may be used to give a count on how
c            many times to print the diagnostics that are given.  This
c            count may be interpreted by library routines in slightly
c            different ways, but when used it should serve to turn off
c            all output after a certain limit is reached.  By
c            convention, if this is 0 there is no upper bound on the
c            count.
c MEMAXE=19  (0 .le. K19 .le. 1000000000) Set the maximum error value.
c            When retrieving this value, it is the maximum value seen
c            for 10000*s + 1000*p + i, where s, p, and i are the stop
c            and print levels, and the index on the last error message
c            processed, respectively.  See MEEMES below.
c MESTOP=20  (0 .le. K20 .le. 8) Set the stop level for error messages.
c            If an error message has a stop index > min(K20, 8), the
c            program is stopped after processing the message.  The
c            default value is K20=3.
c MEPRNT=21  (0 .le. K21 .le. 8) Set the print level for error messages.
c            If an error message has a print index > K21, or the message
c            is going to stop when finished, information in an error
c            message is processed, else all the actions including
c            printing are skipped.  (MESTOP controls stopping.)  The
c            default value is MEPRNT = 3.
c An action index of -i, for i < METDIG, will return in the location
c ordinarily used for Ki the current default value for the internal
c variable set by Ki.  In the case of MESUNI, if the scratch unit has
c not been opened, it will be opened before returning the unit number.
c
c METDIG=22  (-50 .le. K22 .le. 50) As for MEDDIG, except the value here
c            is temporary, lasting until the return, or next use of this
c            action.  If 0, the internal value for K12 is used instead.
c MENTXT=23  (1 .le. K23 .le. 10000000) Set value of NTEXT to K23.
c MEIDAT=24  (1 .le. K24 .le. 1000000000) Set value of NIDAT to K24.
c MEFDAT=25  (1 .le. K25 .le. 1000000000) Set value of NFDAT to K25.
c MEMDAT=26  (1 .le. K26 .le. 5) Set value of NMDAT to K26.
c MEMDA1=27  (K27) set MDAT(1) to K27.  See description of NMDAT above.
c MEMDA2=28  (K28) set MDAT(2) to K28.
c MEMDA3=29  (K29) set MDAT(3) to K29.
c MEMDA4=30  (K30) set MDAT(4) to K30.
c MEMDA5=31  (K31) set MDAT(5) to K31.
c METABS=32  (1 .le. K32 .le. 100) set spacing for tabs to K32.
c MECONT=50  Exit, but no print of current print buffer.  The error or
c            diagnostic message is to be continued immediately.
c MERET=51   All done with diagnostic or error message, complete
c            processing and return, or for some error messages stop.
c MEEMES=52  (K52, L52, M52) Start an error message with severity level
c            K52,index for the error of L52, and message text starting
c            at TEXT(M52).  If M52 is 0, message text starts at
c            TEXT(NTEXT), and if M52 < 0, no message text is
c            printed as part of this action.  Library routines should
c            set K52 = 10*s + p, where s is the stop level desired, and
c            p the print level, and should have 10 > p .ge. s .ge. 0.
c            We offer the following guidelines as a yardstick for
c            setting the value of s.
c   = 9  User has ignored warning that program was going to be stopped.
c   = 8  Program has no way to continue.
c   = 7  User has given no indication of knowing that functionality of
c        results is reduced.  (E.g. not enough space for some result.)
c   = 6  Program could continue but with reduced functionality.
c   = 5  Results far worse than user expected to want.
c   = 4  User has given no indication of knowing that results do not
c        meet requested or expected accuracy.
c   = 3  Warning is given that program will be stopped without some
c        kind of response from the calling program.
c   = 2  Program is not delivering requested or expected accuracy.
c   = 1  Some kind of problem that user could correct with more care in
c        coding or in problem formulation.
c   = 0  Message is for information of uncritical nature.
c            Print levels might be counted down so that warnings given
c            several times are no longer given, or be counted up so
c            that a warning is only given after a certain threshold is
c            reached.  Levels should be selected with the understanding
c            that the default is to print only levels above 3.
c METEXT=53  Print TEXT, starting at TEXT(NTEXT).  Print ends
c            with the last character preceding the first '$'.  Special
c            actions are determined by the character following the '$'.
c            Except as noted, the '$' and the single character which
c            follows are not printed.  In the text below, "to continue",
c            means to continue print of TEXT with the next character
c            until the next "$".  Except for the one case noted, NTEXT
c            is set to point to the second character after the "$".
c            Note, letters must be in upper case.  Possibilities are:
c      B  Break text, but don't start a new line.
c      E  End of text and line.
c      R  Break text, don't change the value of NTEXT.  Thus next
c         text Repeats the current.
c      N  Start a New line, and continue.
c      I  Print IDAT(NIDAT), set NIDAT=NIDAT+1, and continue.
c      J  As for I above, except use the last integer format
c         defined by a "$(", see below.
c      F  Print FDAT(NFDAT), set NFDAT=NFDAT+1, and continue.
c      G  As for F above, except use the last floating format
c         defined by a "$(", see below.
c      M  Print MDAT(NMDAT), set NMDAT=NMDAT+1, and continue.
c      H  Marks terminator for column and row Headings, see table,
c         vector, and matrix output below.  This causes enough blanks to
c         be generated to keep column headings centered over their
c         columns.  After the blanks are generated, text is continued
c         until the next '$'.  This is not to be used except inside
c         column or row headings.  The last row or column should be
c         terminated with a '$E' or if appropriate, a '$#' for a row or
c         column label.
c      (  Starts the definition of a format for integer or floating
c         point output.  The format may not contain a "P" field, and
c         must require no more than 12 characters for floating point
c         (e.g. "(nnEww.ddEe)", where each of the lower case letters
c         represents a single digit), and no more than 7 characters for
c         integer output.  Following the ")" that ends the format, if
c         the next character is not a "$" then "$J" or "$G" type output
c         is done, see above.  In either case processing of TEXT then
c         continues.
c      T  Tab.
c      #  Used in matrix row or column labels this prints the current
c         row or column index, respectively, ends the text for the
c         current row or column, and resets the text pointer to where
c         it started.
c      $  a single '$' is printed, continue till the next '$'.
c      -  Start a negative number for skipping.
c     0-9 Digits for skipping.
c      C  Only used by PMESS which deletes it and the preceding '$'.
c         Used at the end of a line to indicate continued text.
c   other Don't use this -- the '$' is ignored, but new features may
c         change the action.  (E.g. $P might be added to get a prompt.)
c ME????=54  Not used.
c METABL=55  (K55, L55, M55, N55)  Note this action automatically
c            returns when done, further locations in MACT are not
c            examined.  This action prints a heading and/or data that
c            follows a heading.  If K55 is 1, then the heading text
c            starting in TEXT(NTEXT) is printed.  This text
c            should contain embedded "$H"'s to terminate columns of the
c            heading.  If there is no heading on a column, use " $H".
c            Note the leading blank.  If the heading is to continue
c            over k columns, begin the text with "$H" repeated k-1
c            times with no other embedded characters.  The very last
c            column must be terminated with "$E" rather than "$H".
c            After tabular data is printed, K55 is incremented by 1,
c            and compared with L55.  If K55 > L55, K55 is reset to 1,
c            and if the data that was to be printed had lines that were
c            too long, data saved in the scratch file is printed using
c            the headings for the columns that would not fit on the
c            first pass.  Note that only one line of tabular data can
c            be printed on one call to this subroutine.
c            M55 gives the number of columns of data associated with the
c            heading.
c            N55 is a vector containing M55 entries.  The k-th integer
c            in N55 defines the printing action for the k-th column
c            of the table.  Let such an integer have a value defined by
c            rr + 100 * (t + 10 * (dd + 100 * ww)), i.e. wwddtrr, where
c            0 .le. rr,dd,ww < 100, and 0 .le. t < 10.
c      rr    The number of items to print.
c      t     The type of output.
c            1  Print text starting at TEXT(NTEXT), rr = 01.
c            2  Print the value of K55, rr = 01.
c            3  Print integers starting at IDAT(NIDAT).
c            4  Print starting at FDAT(NFDAT), using an F format.
c            5  Print starting at FDAT(NFDAT), using an E format.
c            6  Print starting at FDAT(NFDAT), using an G format.
c      dd    Number of digits after the decimal point.
c      ww    The total number of column positions used by the column,
c            including the space used to separate this column from the
c            preceding one.  This must be big enough so that the column
c            headings will fit without overlap.
c MEIVEC=57  (K57) Print IDAT as a vector with K57 entries.  The vector
c            output starts on the current line even if the current line
c            contains text.  This is useful for labeling the vector.
c            The vector starts at IDAT(NIDAT).
c            If K57 < 0,  indices printed in the labels for the vector
c            start at at NIDAT, and entries from NIDAT to -K57 are
c            printed.
c MEIMAT=58  (K58, L58, M58, I58, J58) Print IDAT as a matrix with K58
c            declared rows, L58 actual rows, and M58 columns.  If K58<0,
c            instead of using 1 for the increment between rows, and K58
c            for the increment between columns, -K58 is used for the
c            increment between rows, and 1 is used for the increment
c            between columns.  If L58<0, the number of actual rows is
c            mod(-L58, 100000), and the starting row index is -L58 /
c            100000.  Similarly for M58<0. TEXT(I58) starts the text for
c            printing row labels.  If I58 < 0, no row labels are
c            printed.  If I58 = 0, it is as if it pointed to text
c            containing "Row $E".  Any "$" in a row or column label must
c            be followed by "H" or "E" which terminates the text for the
c            label.  In the case of $H, text for the next label follows
c            immediately, in the case of $E the current row index is
c            printed in place of the $E and the next label uses the same
c            text.  J58 is treated similarly to I58, except for column
c            labels, and with "Row $E" replaced with "Col $E".  The
c            matrix starts at IDAT(NIDAT), and NIDAT points one past the
c            end of the matrix when finished.
c MEJVEC=59  (K59) As for MEIVEC, except use format set with $(.
c MEJMAT=60  (K60, L60, M60, I60, J60) As for MEIMAT, except use the
c            format set with $(.
c MEFVEC=61  (K61) As for MEIVEC, except print FDAT as a vector with
c            K61 entries.  The vector starts at FDAT(NFDAT).
c MEFMAT=62  (K62, L62, M62, I62, J62) As for action MEIMAT, but
c            instead print FDAT, and use NFDAT in place of NIDAT.
c MEGVEC=63  (K63) As for MEFVEC, except use format set with $(.
c MEGMAT=64  (K64, L64, M64, I64, J64) As for MEIMAT, except use the
c            format set with $(.
c MEIVCI=65  (K65, L65) As for MEIVEC, except the vector entries have a
c            spacing of K65, and there are L65 entries in the vector.
c MEJVCI=66  (K66) As for MEIVCI, except use format set with $(.
c MEFVCI=67  (K67, L67) As for MEFVEC, except the vector entries have a
c            spacing of K67, and there are L67 entries in the vector.
c MEGVCI=68  (K68) As for MEFVCI, except use format set with $(.
c MEFSPV=69  (K69) Output IDAT, FDAT as a sparse vector.
c
c
c ************************** Internal Variables ************************
c
c BUF    Character string holding characters to be output.
c C      Used for temp. storage of a character.
c DOLS   A heading/trailing string of characters, default = $'s.
c ERMSG  Default part of error message.
c ERRCNT Used to keep a count of error messages.
c EUNIT  Unit number for output of error messages.
c FDAT   Formal array, containing floating point data to output.  Only
c   appears external to this subroutine.
c FIRST  Set = .true. initially, then .false. after MESS is called.
c FMTC   Format for integer output of matrix column headings.
c FMTF   Format for floating point or other output.
c FMTG   Format set by user for floating point or other output.
c FMTI   Character string holding format for integer output.
c FMTIM  Equivalenced to FMTR, FMTC.
c FMTJ   Format set by user for integer output.
c FMTR   Value of FMTI for format of row indices.
c FMTT   Format to be stored in FMTJ or FMTG.
c GETW   Set true if still need to get width for a format.
c GOTFMT Set .true. if format has been set by user.  When printing
c   tables, set true when heading has been output.
c I      Index of current action from MACT.
c ICHAR0 Value of ICHAR('0')
c ICOL   Current column index in matrix output.
c IDAT   Formal array, containing integer data to output.
c IMAG   Magnitude of integer to output, with negative sign if integer
c   is < 0.
c INC    Increment between successive elements in a vector or in the
c    column of a matrix.
c INCM   Array giving amount of space used by the options.
c INERR  0 if not processing an error message, 1 if printing an error
c   message, -1 if in an error message that is not being printed, and >1
c   if printing an error message that stops.  Set to -2 when the error
c   message is supposed to stop.
c IOUT   Integer to be output.
c IRC    = 1 for rows, = 2 for columns when determining labels for
c   matrix output.
c IROW   Row index for matrix output.  Also used in table output to
c    count lines for printing line index on read from scratch unit.
c IROW1  Starting row index for matrix output.
c ITEXT  Index of the element of TEXT use for the next text output.
c ITXTSV Saved value of NTEXT when doing matrix output.
c IVAR   Integer array, that is equivalenced to a number of integer
c   variables that can be set by the user.
c IWF    Width to be used in a floating pt. (or other) format.
c IWG    Value like IWF for user set format.
c J      Used as a temporary index.
c JJ      Used as a temporary index.
c K      Used as a temporary index.
c K      Used as a temporary index.
c K1     Used as a temporary index.
c K2     Used as a temporary index.
c KDF    Current number of digits to print for floating point.  See
c   description of MACT(*) = METDIG above.
c KDFDEF Current default for KDF, see description of MACT(*) = MEDDIG.
c KDI    Number of digits used to print last integer.
c KDIAG  Not directly referenced.  Holds place in IVAR for reference
c   from outside.  See comments above.
c KDILAB Length for printing index in vector output.
c KDJ    As for KDI, except for format set by user.
c KK     Temporary index.
c KLINE  Count of number of things to print on a line.  (In table print
c   is the number to print for one spec field.)
c KNT    In vector output gives the current index for output.
c KOLWID Length to use for printing a column heading.  Equivalenced to
c   MAXWID(2).
c KP     Index from error action input for the print action.
c KRES1  Holds place in common block for future use.
c KS     Index from error action input for the stop action.
c KSCRN  Number of lines to "print" before pausing.
c KSHIFT Amount to shift column heading before printing.
c KSPEC  Defines action after looking for character in TEXT.  (Also
c   used as a temporary index.)
c   1.  $B   Break the text here continue on same line.
c   2.  $E   Break text, print what is in BUF.
c   3.  $R   Break text, continue on same line, NTEXT set to repeat the
c            current text.
c   4.  $N   Print BUF, continue with following text.
c   5.  $I   Print IDAT(NIDAT), continue TEXT.
c   6.  $F   Print FDAT(NFDAT), continue TEXT.
c   7.  $M   Print MDAT(NMDAT), continue TEXT.
c   8.  $J   As for $I, but with user format.
c   9.  $G   As for $F, but with user format.
c  10.  $(   Set a user format.
c  11.  $T   Tab.
c  12.       Set when done with an action.
c  13.       Set when done with boiler plate text for an error message.
c   0. Other Ignore the "$", continue with TEXT.
c KT     Used for logic in output of headings.
c        = 1 Output table headings.
c        = 2 Get row/column widths for matrix output.
c        = 3 Output column headings for matrix output.
c LASKNT In vector output value index for last element to print.
c LASTI  Last index for matrix output, or for finding values that
c   determine format.
c LBUF   Position of characters in BUF, usually the last to print.
c LBUF1  Start of text to shift when shifting text in BUF to the left.
c LBUF2  End of text to shift when shifting text in BUF to the left.
c LENBUF Parameter giving the number of character in BUF.
c LENLIN Gives number of character in output lines.
c LENOUT Length of output for table or vector/matrix output.
c LENTXT Length of character array elements in TEXT.
c LENTRY Tells what to do on entry (and sometimes other places.)
c   = 1  The value on first entry.
c   = 2  A previous entry is to be continued.
c   = 3  A non printing error message is to be continued
c   = 4  Just done output from inside a METEXT action.
c   = 5  Got "maximum" value for entries in a vector.
c   = 6  Got "maximum" value for entries in a matrix.
c   = 7  Vector (either print or get format for label indices.)
c   = 8  Matrix (either print or get format for label indices.)
c   = 9  Output of data in a table.
c   =10  Get "maximum" valur for entries in a sparse vector.
c   =11  Output a sparse vector.
c LHEAD  If = 0 no print of DOLS, else DOLS printed in error messages.
c LINERR Gives LENLIN for error messages.
c LINMSG Gives LENLIN for diagnostic messages.
c LINSTR Space at start of line for label in vector and matrix output.
c   Equivalenced to MAXWID(1).
c LNERR  Parameter giving the default value for LENERR, only in XMESS.
c LNMSG  Parameter giving the default value for LINMSG, only in XMESS.
c LOCBEG Index of first item in vector and matrix output.
c LPRINT For error messages with a print level .le. LPRINT nothing is
c   printed (unless the message would result in a stop).
c LSTOP  As for LPRINT, except with the stop level, and stopping.
c LSTRT  Starting character position in BUF for storing next characters.
c LTEXT  Length of heading text in TEXT.
c M      Index for the current action.
c MACT   Formal integer array specifying the actions, see above.
c MAXERR Value save in IVAR for user to get indication of last (really
c   the maximum error seen so far = 1000 * (10*stop + print) + index.
c MAXWID Equivalenced to  LINSTR and KOLWID.
c MBNDHI Upper bounds for inputs to IVAR.
c MBNDLO Lower bounds for inputs to IVAR.
c MDAT   Array where user can store integers in IVAR for later output.
c   Also used to store column indices when tables are saved on scratch
c   unit.
c
c The following parameter names starting with ME define actions
c   which the user can request.  They have been documented in the
c   description above, except for the ones defined just below.
c MEGBAS is 1 less than the smallest action that involves something
c   other than just storing or retrieving a value.
c MEMAXI is the largest action which a user can request.
c MEVBAS is the smallest action index, used to set the starting index in
c   IVAR.
c MEVLAS is the largest index for a variable in IVAR.
c MECONT,  MEDDI,  MEELI,  MEEME, MEEUNI, MEFDAT, MEFMAT, MEFSPV,
c MEFVCI, MEFVEC, MEGBAS, MEGMAT, MEGVCI, MEGVEC, MEHEAD, MEIDAT,
c MEIMAT, MEIVCI, MEIVEC, MEJMAT, MEJVCI, MEJVEC, MEMAXE, MEMAXI,
c MEMDA1, MEMDA2, MEMDA3, MEMDA4, MEMDA5, MEMDAT, MEMLIN, MEMUNI,
c MENTXT, MEPRNT, MESCRN, MERES1, MERES2, MERES3,  MERET, MESTOP,
c MESUNI,  METAB, METDIG, METEXT
c MPT  Current pointer to data for matrix or vector output.
c MTEXT  Equivalenced to MTEXTR and MTEXTC.
c MTEXTC TEXT(MTEXTC) starts text for printing column labels.
c MTEXTR TEXT(MTEXTR) starts text for printing row labels.
c MUNIT  Output unit used for messages that aren't in an error message.
c NCOL   Number of columns for matrix output, 0 for vector output,
c   count of column left for table output.
c NDIM   Distance between columns for matrix output.
c NFDAT  Index of next item in FDAT to print.
c NIDAT  Index of next item in IDAT to print.
c NLINE  Maximum number of data items to print on a line for matrix and
c   vector output.  If the scratch file is used for table output,
c   NLINE gives the original end of the buffer.
c NMDAT  Pointer to next thing to print from MDAT.
c NROCO  Equivalenced to (NROW, NCOL).  Used in matrix output.
c NROW   Number of rows for matrix output.    When printing tables,
c   MDAT(NROW) gives place where line was split.  (=0 if not split)
c NSKIP  The amount to skip ahead on the next floating or integer
c   output.
c NTEXT  Index inside an element of TEXT for the next text output.
c NTEXTR Value of NTEXT to use if get a $R.
c NTXTSV Saved value of NTEXT when doing matrix output.
c OUNIT  Index of the current output unit.
c SC     Parameter for special character used to introduce actions.
c   Default value is '$'.  If this is changed the "$"'s in comments
c   should be changed to the new value of the character.  (Note that
c   SC = '\' is not portable.)
c SCRNAM Name of file constructed for error output or message output.
c SUNIT  Index for the scratch unit, -1 if not yet assigned.
c TEXT   Formal argument giving the character string from which all text
c   is taken.
c UMESS  Name of subroutine called that does nothing, but which may be
c   modified by the user to cause different actions to be taken.
c   The usual version of MESS has the call to UMESS commented out.
c XARG   If .true., output data is not integer, and a return is made to
c   print data from FDAT.
c XARGOK Set .true. if call is from program that will print data from
c   FDAT.
c
c++ CODE for .C. is inactive
c      integer  kciwid, kccwid, kcrwid, lbeg, lend, lfprec, lgprec
c      common /MESSCC/ kciwid,kccwid,kcrwid,lbeg,lend,lfprec,lgprec
C%%    long int kc;
c++ END
      integer LNMSG, LNERR
      parameter (LNMSG=128)
      parameter (LNERR=79)
c
c ************** Parameters Defining Actions (See Above) ***************
c
      integer   MESUNI, MEHEAD, MEDDIG, MEMLIN, MEELIN, MEMUNI, MEEUNI,
     1  MESCRN, MEDIAG, MEMAXE, MESTOP, MEPRNT, METDIG, MENTXT, MEIDAT,
     2  MEFDAT, MEMDAT, MEMDA1, MEMDA2, MEMDA3, MEMDA4, MEMDA5, METABS,
     3  MEERRS, MECONT, MERET , MEEMES, METEXT, METABL, MERES3, MEIVCI,
     4  MEIVEC, MEIMAT, MEJVCI, MEJVEC, MEJMAT, MEFVCI, MEFVEC, MEFMAT,
     5  MEGVCI, MEGVEC, MEGMAT, MEMAXI, MEGBAS, MEVBAS, MEVLAS, MEFSPV
c Parameters for changing the environment.
      parameter (MESUNI=10,MEHEAD=11,MEDDIG=12,MEMLIN=13,MEELIN=14,
     1 MEMUNI=15,MEEUNI=16,MESCRN=17,MEDIAG=18,MEMAXE=19,MESTOP=20,
     2 MEPRNT=21,METDIG=22,MENTXT=23,MEIDAT=24,MEFDAT=25,MEMDAT=26,
     3 MEMDA1=27,MEMDA2=28,MEMDA3=29,MEMDA4=30,MEMDA5=31,METABS=32,
     4 MEERRS=33)
c Parameters for actions.
      parameter (MECONT=50, MERET=51,MEEMES=52,METEXT=53,MEFSPV=54,
     1 METABL=55,MERES3=56,MEIVEC=57,MEIMAT=58,MEJVEC=59,MEJMAT=60,
     2 MEFVEC=61,MEFMAT=62,MEGVEC=63,MEGMAT=64,MEIVCI=65,MEJVCI=66,
     2 MEFVCI=67,MEGVCI=68)
c Parameter derived from those above.
      parameter (MEMAXI=68,MEGBAS=49,MEVBAS=10,MEVLAS=33)
c
c ************************** Variable Declarations *********************
c
      external MESSGS
      integer    MACT(*), IDAT(*)
      character  TEXT(*)*(*)
c
      integer    I, ICOL, INCM(MECONT:MEIMAT), INERR, IOUT, IROW, IROW1,
     1    ITEXTR, ITXTSV, J, JJ, K, K1, K2, KDILAB, KK, KNT, KOLWID, KP,
     2    KS, LASKNT, LBUF1, LBUF2, LENBUF, LINSTR, M,
     3    MBNDHI(MEVBAS:MEVLAS), MBNDLO(MEVBAS:MEVLAS), MTEXT(2),
     4    MTEXTC, MTEXTR, NLINE, NROCO(2), NSKIP, NTEXTR, NTXTSV
      integer MESSGS
      logical   GETW, FIRST
      character ERMSG*63, ERMSG1*27
      character SC, C
      parameter (SC='$')
      save  FIRST, I, ICOL, INERR, IROW, IROW1, ITXTSV, KDILAB, KNT,
     1   LASKNT, M, MTEXT, NLINE, NSKIP, NTEXTR, NTXTSV
      save /CMESSI/, /CMESSC/
      equivalence (MTEXT(1), MTEXTR), (MTEXT(2), MTEXTC)
c
c ************************** Data from common block ********************
c
      parameter (LENBUF=250)
      logical          XARG, GOTFMT, XARGOK
      integer          ERRCNT, EUNIT, ICHAR0, IRC, IVAR(MEVBAS:MEVLAS),
     1   IMAG, INC, ITEXT, IWF, IWG, KDF, KDFDEF, KDI, KDIAG, KDJ,
     2   KLINE, KSCRN, KSHIFT, KSPEC, KT, MAXERR, LASTI, LBUF, LENLIN,
     3   LENOUT, LENTRY, LENTXT, LHEAD, LINERR, LINMSG, LOCBEG, LPRINT,
     4   LSTOP, LSTRT, LTEXT, MAXWID(2), MDAT(5), MPT, MUNIT, NCOL,
     5   NDIM, NFDAT, NIDAT, NMDAT, NROW, NTEXT, OUNIT, SUNIT, TABSPA
c
      character BUF*(LENBUF), DOLS*72, FMTC*7, FMTF*20, FMTG*15,
     1  FMTI*7, FMTIM(2)*7, FMTJ*7, FMTR*7, FMTT*15
      common /CMESSI/ SUNIT, LHEAD, KDFDEF, LINMSG, LINERR, MUNIT,
     1   EUNIT, KSCRN, KDIAG, MAXERR, LSTOP, LPRINT, KDF, NTEXT, NIDAT,
     2   NFDAT, NMDAT, MDAT, TABSPA, ERRCNT, ICHAR0, IMAG, INC, IRC,
     3   ITEXT, IWF, IWG, KDI, KDJ, KLINE, KSHIFT, KSPEC, KT, LASTI,
     4   LBUF, LENLIN, LENOUT, LENTRY, LENTXT, LOCBEG, LSTRT, LTEXT,
     5   MAXWID, MPT, NROW, NCOL, NDIM, OUNIT, GOTFMT, XARG, XARGOK
      common /CMESSC / BUF, DOLS, FMTF, FMTG, FMTI, FMTJ, FMTT, FMTIM
      equivalence (IVAR(MEVBAS), SUNIT)
      equivalence (FMTIM(1), FMTR), (FMTIM(2), FMTC)
c
      equivalence (NROCO, NROW)
      equivalence (MAXWID(1), LINSTR), (MAXWID(2), KOLWID)
c ************************** End of stuff from common block ************
c
      data INERR, FIRST / 0, .true. /
      data ERMSG /
     1' reports error: Stop level = x, Print level = y, Error index = '/
      data ERMSG1 / ': Print level = y, Index = ' /
c                 50  51, 52  53  54 55 56 57 58
      data INCM /  1,  1,  4,  1,  2, 0, 0, 2, 6 /
      data MBNDLO /  0, 0, -50,  39,  39, -99, -99,         0,
     1       0,          0, 0, 0, -50,        1,          1,
     2       1, 1, -1000000000, -1000000000, -1000000000, -1000000000,
     3       -1000000000,   1,  0 /
      data MBNDHI / 99, 1,  50, 500, 500,  99,  99, 100000000,
     1  1000000000, 1000000000, 8, 8,  50, 10000000, 1000000000,
     2  1000000000, 5, 1000000000, 1000000000, 1000000000, 1000000000,
     3  1000000000, 100, 1000000000 /
c
c ************************* Start of Executable Code *******************
c
c
      NSKIP = 0
      if (FIRST) then
         FIRST = .false.
c Initialize common block
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
c++ CODE for ~.C. is active
         DOLS(1:40) = '$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$'
         DOLS(41:72) ='$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$'
         FMTI = '(99I01)'
         FMTJ = '(99I06)'
         FMTG = '(1P,99Exx.xx)  '
c++ CODE for .C. is inactive
C%%    memset(cmessc.dols,'$',72);
C      FMTI = '%*d'
C      FMTJ = '%*d\0'
C      FMTG = '%*.*E\0'
c++ END
      else
c               1  2  3   4    5    6    7    8   9   10   11
         go to (5,10,20,850,1160,1620,1130,1530,960,1210,1220), LENTRY
      end if
c                             First entry for a message
    5 LBUF = 0
c                             Usual continuation entry
   10 I = 1
      NTEXT = 1
      ITEXT = 1
      LENTXT = len(TEXT(1))
      NIDAT = 1
      NFDAT = 1
      NMDAT = 1
      go to 120
c                     Continuation entry when have a non printing error
c Skip all actions -- Inside non-printing error message.
   20 I = 1
   30 K = MACT(I)
      if (K .le. MERET) then
         if (K .eq. MERET) go to 120
         if (K .eq. MECONT) return
         if (K .le. -MEGBAS) go to 180
         I = I + 2
      else
         if (K .gt. MEIMAT) then
            if (K .gt. MEMAXI) go to 180
            K = MEIVEC + mod(K - MEIVEC, 2)
         end if
         I = I + INCM(K)
      end if
      go to 30
c
c Print BUF
   40 call MESSPR
c                             Usual place to end an action request.
  100 I = I + INCM(M)
c                             Pick up the next action request
  120 M = MACT(I)
      if (M .gt. MEGBAS) go to 140
      I = I + 2
      if (abs(M) .gt. MEVLAS) go to 180
      if (M .gt. 0) then
         IVAR(M) = MACT(I-1)
         if (IVAR(M) .lt. MBNDLO(M)) then
            IVAR(M) = MBNDLO(M)
         else if (IVAR(M) .gt. MBNDHI(M)) then
            IVAR(M) = MBNDHI(M)
         end if
c            MEHEAD, MEDDIG, MEMLIN, MEELIN, MEMUNI, MEEUNI
         go to (122,    124,    126,    126,    128,    128), M - MESUNI
         if (M .ne. MENTXT) go to 120
         ITEXT = (NTEXT-1) / LENTXT
         NTEXT = NTEXT - LENTXT*ITEXT
         ITEXT = ITEXT + 1
         go to 120
  122    if (LHEAD .ne. 0) then
         end if
         go to 120
  124    KDF = KDFDEF
         go to 120
  126    LENLIN = LINMSG
         go to 120
  128    if (IVAR(M) .ne. 0) then
C%%          k = labs(cmessi.ounit);
C%%          c_fname[m-15][6] = k / 10 + '0';
C%%          c_fname[m-15][7] = k % 10 + '0';
C%%          if (strcmp(&c_fname[16-m][6], &c_fname[m-15][6]))
C%%             c_handle[m-15] = fopen(c_fname[m-15],"w");
C%%          else
C%%             c_handle[m-15] = c_handle[16-m];
            K = abs(IVAR(M))
         end if
         OUNIT = MUNIT
         go to 120
      end if
      if (M .eq. -MESUNI) then
C%%      if (cmessi.sunit == -1L) {
C%%          scratch_file = tmpfile();
C%%          cmessi.sunit = 1L;}
         if (SUNIT .le. 0) SUNIT = MESSGS()
      end if
C
      MACT(I-1) = IVAR(-M)
      go to 120
c  ME ..    CONT  RET EMES ETXT  FSPV TABL
  140 go to (170, 200, 310, 400, 1200, 910, 180), M-MEGBAS
      if (M .le. MEGVCI) go to 1000
      go to 180
c
c Action MECONT -- Continue message on next entry
  170 LENTRY = 2
      return
c
c Some kind of error in message specification.
  180 continue
c++ CODE for ~.C. is active
      BUF(1:57) =
     1   'Actions in MESS terminated due to error in usage of MESS.'
c++ CODE for .C. is inactive
C%%   memcpy(cmessc.buf,
C%%   "Actions in MESS terminated due to error in usage of MESS.",57);
c++ END
      LBUF = 57
c
c Action MERET -- Finish a message.
  200 LENTRY = 1
      J = INERR
      INERR = 0
      if (J .ge. 2) INERR = -2
      if (J .gt. 0) go to 330
c                       Finish print before exit.
      call MESSPR
      return
c
c Action MEEMES -- Start an error message
  310 LENTRY = 3
      ERRCNT = ERRCNT + 1
c++  Code for UMESS is inactive
C      call UMESS(TEXT, MACT(I+1), IVAR)
c++  End
      IMAG = max( 0, min(999, MACT(I+2)))
      K = MACT(I+1)
      MAXERR = max(MAXERR, 1000*K + IMAG)
      KS = K / 10
      KP = K - 10 * KS
      if (KS .le. min(LSTOP, 8)) then
         if (KP .le. LPRINT) then
            INERR = -1
            go to 20
         end if
         INERR = 1
      else
         INERR = 2
      end if
      OUNIT = EUNIT
      LENLIN = LINERR
c                        Output a blank line.
      BUF(1:1) = ' '
      LBUF = 1
  330 call MESSPR
c                        Put out line of $'s
      if (LHEAD .ne. 0) then
         LBUF = min(len(DOLS), LENLIN)
c++ CODE for ~.C. is active
         BUF(1:LBUF) = DOLS(1:LBUF)
         if (INERR.lt.0) BUF(5:37)=' Fatal error -- Program stopped. '
c++ CODE for .C. is inactive
C%%      memcpy(cmessc.buf, cmessc.dols, cmessi.lbuf);
C%%      if (inerr < 0L)
C%%      memcpy(&cmessc.buf[4]," Fatal error -- Program stopped. ",34);
c++ END
         call MESSPR
      end if
      if (INERR .le. 0) then
c                                 Just finished an error message
         if (INERR .ne. 0) stop
         OUNIT = MUNIT
         LENLIN = LINMSG
         return
      end if
c                     Just starting an error message get program name
      NTEXTR = 0
      go to 410
c                     Got the program name in BUF.
  370 LBUF = min(LBUF, 40)
      if (KS .eq. 0) then
         ERMSG1(17:17) = char(KP + ICHAR0)
C%%       memcpy(&cmessc.buf[cmessi.lbuf], ermsg1, strlen(ermsg1));
         BUF(LBUF+1:LBUF+len(ERMSG1)) = ERMSG1
         LBUF = LBUF + len(ERMSG1)
      else
         ERMSG(30:30) = char(KS + ICHAR0)
         ERMSG(47:47) = char(KP + ICHAR0)
C%%       memcpy(&cmessc.buf[cmessi.lbuf], ermsg, strlen(ermsg));
         BUF(LBUF+1:LBUF+len(ERMSG)) = ERMSG
         LBUF = LBUF + len(ERMSG)
      end if
      LSTRT = LBUF + 1
      call MESSFI
      LBUF = LBUF + KDI
C%%   sprintf(&cmessc.buf[cmessi.lstrt-1L], "%*ld",
C%%           (int)messcc.kciwid, cmessi.imag);
      write (BUF(LSTRT:LBUF), FMTI) IMAG
c          Finish up the start error message action.
      if (MACT(I+3) .lt. 0) go to 40
      if (MACT(I+3) .ne. 0) then
         ITEXT = (MACT(I+3)-1) / LENTXT
         NTEXT = MACT(I+3) - LENTXT*ITEXT
         ITEXT = ITEXT + 1
      end if
      KSPEC = 13
      go to 480
c                  Take care of any left over print from error header
  390 if (LBUF .ne. 0) call MESSPR
c
c Action METEXT -- Print string from TEXT
  400 LENTRY = 4
      NTEXTR = NTEXT
      ITEXTR = ITEXT
c                  Continue with print from TEXT
c K     take at most K-1 chars., but if 0 take max number
c K1    is last loc. used from TEXT if LENTXT is BIG.
c NEXT  is first character location in TEXT(ITEXT)
c K2    is last character location in TEXT(ITEXT)
c LSTRT is first character position in BUF
c LBUF  is last used character position in BUF

 410  LSTRT = LBUF + 1
      K2 = min(LENTXT, NTEXT + (LENBUF - LSTRT))
C%%       if ((ctmp=memchr(TEXT(cmessi.itext-1L,cmessi.ntext-1), SC,
C%%          k2 - cmessi.ntext + 1)) == NULL)
C%%             k = 0;
C%%       else
C%%             k = ctmp - TEXT(cmessi.itext-1L,cmessi.ntext-1) + 1;
      K = index(TEXT(ITEXT)(NTEXT:K2), SC)
      if (K .eq. 0) then
c Want to take all that we can.
         LBUF = LSTRT + K2 - NTEXT
C%%       memcpy(&cmessc.buf[cmessi.lstrt-1L], TEXT(cmessi.itext-1L,
C%%         cmessi.ntext-1), k2 - cmessi.ntext + 1L);
         BUF(LSTRT:LBUF) = TEXT(ITEXT)(NTEXT:K2)
         if (K2 .eq. LENTXT) then
           ITEXT = ITEXT + 1
           NTEXT = 1
           if (LBUF .le. LENLIN) go to 410
         else
           NTEXT = K2 + 1
         end if
         KSPEC = 12
         if (ITEXT - ITEXTR .lt. 4000) go to 480
         KSPEC = 2
         go to 430
      end if
      LBUF = LBUF + K - 1
C%%   if (k >= 2) memcpy(&cmessc.buf[cmessi.lstrt-1],
C%%     TEXT(cmessi.itext-1L, cmessi.ntext-1), k - 1L);
      if (K .ge. 2) BUF(LSTRT:LBUF) = TEXT(ITEXT)(NTEXT:NTEXT+K-2)
c        Jump to location below if get $ after computing an NSKIP.
  415 continue
      NTEXT = NTEXT + K + 1
      if (NTEXT .gt. LENTXT) then
         ITEXT = ITEXT + 1
         if (NTEXT .eq. LENTXT + 1) then
            C = TEXT(ITEXT-1)(LENTXT:LENTXT)
            NTEXT = 1
         else
            C = TEXT(ITEXT)(1:1)
            NTEXT = 2
         end if
      else
         C = TEXT(ITEXT)(NTEXT-1:NTEXT-1)
      end if
      if (C .eq. ' ') then
c                Special code to take care of " " following "$".
         NTEXT = NTEXT - 1
         if (NTEXT .eq. 0) then
            ITEXT = ITEXT - 1
            NTEXT = LENTXT
         end if
         go to 410
      end if
      if (NTEXTR .eq. 0) then
         if (LENTRY .eq. 3) go to 370
         go to 1510
      end if
      KSPEC = index('BERNIMFJG(T', C)
  430 if (LBUF .gt. LENLIN) go to 480
c              1   2   3   4   5   6   7   8   9  10  11  12, 13
c              B   E   R   N   I   M   F   J   G   (   T done end err
      go to (455,480,450,460,700,680,900,700,900,600,690,410,390), KSPEC
c               No match  -- Check for setting NSKIP
      if (((C .ge. '0') .and. (C .le. '9')) .or. (C .eq. '-')) then
         NSKIP = 0
         K1 = 1
         if (C .ne. '-') go to 436
         K1 = -1
  433    C = TEXT(ITEXT)(NTEXT:NTEXT)
         NTEXT = NTEXT + 1
         if (NTEXT .ge. LENTXT) then
            ITEXT = ITEXT + 1
            NTEXT = 1
         end if
  436    if ((C .ge. '0') .and. (C .le. '9')) then
            NSKIP = 10 * NSKIP + K1 * (ICHAR(C) - ICHAR0)
            go to 433
         end if
         if (C .eq. '$') then
            K = 0
            go to 415
         end if
      end if
c
c Continue with the text.
  440 LBUF = LBUF + 1
      BUF(LBUF:LBUF) = C
      go to 410
c                        Reset NTEXT for $R
  450 NTEXT = NTEXTR
      ITEXT = ITEXTR
c                             Done with METEXT action.
  455 NMDAT = 1
      go to 100
c           At this point want to output all in BUF
  460 do 470 LBUF = LBUF, 1, -1
         if (BUF(LBUF:LBUF) .ne. ' ') go to 480
  470 continue
  480 LBUF2 = LBUF
      if (LBUF2 .eq. 0) then
         LBUF = 1
         BUF(1:1) = ' '
      else if (LBUF .gt. LENLIN) then
         do 485 K = LENLIN+1, LENLIN/3, -1
            if (BUF(K:K) .eq. ' ') then
               LBUF = K - 1
               go to 490
            end if
  485    continue
         LBUF = LENLIN
      end if
  490 LBUF1 = LBUF
      call MESSPR
      if (LBUF1 .ge. LBUF2) then
c                       The entire buffer has been printed.
         if (KSPEC .le. 2) go to 455
         if (KSPEC .ne. 4) go to 430
         go to 410
      end if
c                       Remove trailing blanks
      do 510 LBUF1 = LBUF1+1, LBUF2
         if (BUF(LBUF1:LBUF1) .ne. ' ') go to 520
  510 continue
c                       Shift the contents of the buffer.
  520 LBUF = LBUF2-LBUF1+1
      LSTRT = 1
  530 if (LBUF .ge. LBUF1) then
c                              Take care of overlap.
         K = 2*LBUF1 - LSTRT
C%%memcpy(&cmessc.buf[cmessi.lstrt-1],&cmessc.buf[lbuf1-1],k-lbuf1);
         BUF(LSTRT:LBUF1-1) = BUF(LBUF1:K-1)
         LSTRT = LBUF1
         LBUF1 = K
         go to 530
      end if
C%% if (cmessi.lbuf>=cmessi.lstrt) memcpy(&cmessc.buf[cmessi.lstrt-1],
C%%       &cmessc.buf[lbuf1-1L], lbuf2-lbuf1+1);
      if (LBUF .ge. LSTRT) BUF(LSTRT:LBUF) = BUF(LBUF1:LBUF2)
      go to 430
c
c Get information on user format
  600 KSPEC = 8
c              I,   i,   F,   f,   E,   e,   G,   g
      go to (604, 604, 601, 601, 602, 602, 602, 602),
     1   index('IiFfEeGg',TEXT(ITEXT)(NTEXT:NTEXT))
      go to 180
  601 continue
c++ CODE for ~.C. is active
      FMTG='(0P,99F'
c++ CODE for .C. is inactive
C%%   strcpy(cmessc.fmtg, "%*.*f\0");
C%%   messcc.lgprec = 0;
c++ END
      go to 603
  602 continue
c++ CODE for ~.C. is active
      FMTG='(1P,99'//TEXT(ITEXT)(NTEXT:NTEXT)
c++ CODE for .C. is inactive
C%%   strcpy(cmessc.fmtg, "%*.*E\0");
C      FMTG(5:5) = TEXT(ITEXT)(NTEXT:NTEXT)
C%%   messcc.lgprec = 0;
c++ END
  603 KSPEC = 9
  604 IMAG = 0
      GETW = .true.
      K = NTEXT
      FMTT = ' '
  606 continue
         NTEXT = NTEXT + 1
         if (NTEXT .gt. LENTXT) then
            ITEXT = ITEXT + 1
            NTEXT = 1
         end if
c++ CODE for ~.C. is active
         FMTT(NTEXT-K:NTEXT-K) = TEXT(ITEXT)(NTEXT:NTEXT)
c++ END
         JJ = ichar(TEXT(ITEXT)(NTEXT:NTEXT)) - ICHAR0
         if (GETW) then
            if ((JJ .ge. 0) .and. (JJ .le. 9)) then
               IMAG = 10*IMAG + JJ
            else
               if (TEXT(ITEXT)(NTEXT:NTEXT) .eq. ')')  go to 610
               if (TEXT(ITEXT)(NTEXT:NTEXT) .ne. '.')  go to 180
               GETW = .false.
            end if
         else
            if (TEXT(ITEXT)(NTEXT:NTEXT) .eq. ')') go to 610
            if ((JJ .lt. 0) .or. (JJ .gt. 9)) go to 180
c++ CODE for .C. is inactive
C%%         messcc.lgprec = 10*messcc.lgprec + jj;
c++ END
         end if
      go to 606
c
  610 NTEXT = NTEXT + 1
      if (NTEXT .gt. LENTXT) then
         ITEXT = ITEXT + 1
         NTEXT = 1
      end if
c++ CODE for ~.C. is active
      if (KSPEC .eq. 8) then
         KDJ = IMAG
         FMTJ(5:7) = FMTT
      else
         IWG = IMAG
         FMTG(8:15) = FMTT
      end if
c++ CODE for .C. is inactive
C%%   if (cmessi.kspec == 8)
C%%       cmessi.kdj = cmessi.imag;
C%%   else
C%%       cmessi.iwg = cmessi.imag;
c++ END
      if (TEXT(ITEXT)(NTEXT:NTEXT) .eq. SC) go to 410
      if (KSPEC .eq. 8) go to 700
      if (XARGOK) return
      go to 440
c
c                         Print from MDAT
  680 IOUT = MDAT(NMDAT)
      if (NMDAT .ge. 6) then
         MDAT(NMDAT) = MDAT(NMDAT) + 1
      else
         NMDAT = NMDAT + 1
      end if
      go to 720
c
c                         Process a tab
  690 LSTRT = LBUF + 1
      LBUF = min(LBUF + TABSPA - mod(LBUF, TABSPA), LENLIN+1)
C%%  for (kc=cmessi.lstrt-1; kc<cmessi.lbuf; kc++) cmessc.buf[kc]=' ';
      BUF(LSTRT:LBUF) = ' '
      go to 850
c                         Print from IDAT
  700 NIDAT = NIDAT + NSKIP
      NSKIP = 0
      IOUT = IDAT(NIDAT)
      NIDAT = NIDAT + 1
  720 LSTRT = LBUF + 1
      IMAG = IOUT
      if (KSPEC .ge. 8) then
         LBUF = LBUF + KDJ
C%%   sprintf(&cmessc.buf[cmessi.lstrt-1],"%*ld",(int)cmessi.kdj, iout);
      write (BUF(LSTRT:LBUF), FMTJ) IOUT
         go to 850
      end if
c
c                Get format for integer output.
      call MESSFI
      LBUF = LBUF + KDI
C%% sprintf(&cmessc.buf[cmessi.lstrt-1],"%*ld",(int)messcc.kciwid,iout);
      write (BUF(LSTRT:LBUF), FMTI) IOUT
c                         Entry here to check line after numeric output.
  850 if (LBUF .le. LENLIN) go to 410
      KSPEC = 12
      go to 480
c
c                          Take care of output for extra argument.
  900 if (XARGOK) return
      go to 180
c
c Action METABL -- Start a table
  910 GOTFMT = MACT(I+1) .ne. 1
      if (.not. GOTFMT) then
         IROW = 0
         KOLWID = 0
      end if
      LENTRY = 9
      if (LBUF .ne. 0) call MESSPR
  920 continue
C%%   memset(cmessc.buf,' ',LENBUF);
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
c                                 Print the data
         LSTRT = LBUF + 1
         LBUF = min(LBUF + KLINE * LENOUT, LENBUF)
         JJ = JJ / 100
         KK = mod(JJ, 10)
c              Text,   I   I',   F    E    G
         go to (948, 941, 941, 943, 945, 944), KK
         go to 180
c                             Integer output
  941    continue
c++ CODE for ~.C. is active
         KDI = LENOUT
         FMTI(5:5) = char(LENOUT / 10 + ichar0)
         FMTI(6:6) = char(mod(LENOUT, 10) + ichar0)
c++ END
         if (KK .eq. 3) then
C%%         sprintf(&cmessc.buf[cmessi.lstrt-1], "%*ld",
C%%            (int)cmessi.lenout, mact[i]);
            write (BUF(LSTRT:LBUF), FMTI) MACT(I+1)
            go to 960
         end if
c                            Regular integer output
         NIDAT = NIDAT + NSKIP
         NSKIP = 0
c++ CODE for ~.C. is active
         write (BUF(LSTRT:LBUF), FMTI) (IDAT(K), K = NIDAT,
     1      NIDAT+KLINE-1)
         NIDAT = NIDAT + KLINE
c++ CODE for .C. is inactive
C%%  kk = cmessi.nidat;
C%%  for (cmessi.nidat=kk; cmessi.nidat<kk+cmessi.kline; cmessi.nidat++)
C%%     sprintf(&cmessc.buf[cmessi.lstrt+cmessi.lenout*(cmessi.nidat
C%%       - kk) - 1], "%*ld", (int)cmessi.lenout, idat[cmessi.nidat-1]);
c++ END
         go to 960
c                           Various floating point output
  943    continue
c++ CODE for ~.C. is active
         FMTF = '(0P,99F  .  )'
c++ END
         go to 946
  944    continue
c++ CODE for ~.C. is active
         FMTF = '(1P,99G  .  )'
c++ END
         go to 946
  945    continue
c++ CODE for ~.C. is active
         FMTF = '(1P,99E  .  )'
c++ END
  946    JJ = mod(JJ/10, 100)
c++ CODE for ~.C. is active
         FMTF(8:8) = char(ICHAR0 + LENOUT / 10)
         FMTF(9:9) = char(ICHAR0 + mod(LENOUT, 10))
         FMTF(11:11) = char(ICHAR0 + JJ / 10)
         FMTF(12:12) = char(ICHAR0 + mod(JJ, 10))
c++ CODE for .C. is inactive
C%%      strcpy(cmessc.fmtf, "%*.*E\0");
C        IWF = LENOUT
C        lfprec = JJ
c++ END
         if (.not. XARGOK) go to 180
         MPT = NFDAT
         NFDAT = NFDAT + KLINE
         return
c                           Text output
  948    K1 = NTEXT + LBUF - LSTRT
C%%    memcpy(&cmessc.buf[cmessi.lstrt-1], TEXT(cmessi.itext-1,
C%%       cmessi.ntext -1), k1 - cmessi.ntext);
         BUF(LSTRT:LBUF) = TEXT(ITEXT)(NTEXT:K1-1)
         NTEXT = K1
      else
c                                 Print the heading
         KT = 1
         call MESSMH(TEXT)
         if (KT .lt. 0) go to 180
      end if
  960 if ((LBUF .le. MDAT(NROW)) .and. (NCOL .gt. 0)) go to 940
      if (NROW .eq. 1) then
         JJ = LBUF
         LBUF = MDAT(1)
         call MESSPR
         LBUF = JJ
      else
         if (IROW .eq. 0) then
            if (NROW .eq. 2) then
c++ CODE for ~.C. is active
               if (SUNIT .le. 0) SUNIT = MESSGS()
               rewind(SUNIT)
c++ CODE for .C. is inactive
C%%        if (cmessi.sunit == -1) {
C%%           scratch_file = tmpfile();
C%%           cmessi.sunit = 1;}
C%%        rewind(scratch_file);
c++ END
            end if
         end if
C%%       fwrite(&cmessc.buf[4], cmessi.mdat[cmessi.nrow-1]-4, 1,
C%%          scratch_file);
         write(SUNIT) BUF(5:MDAT(NROW))
      end if
      if (LBUF .gt. MDAT(NROW)) then
C%%  memcpy(&cmessc.buf[4], &cmessc.buf[cmessi.mdat[cmessi.nrow-1]],
C%%     cmessi.lbuf - cmessi.mdat[cmessi.nrow-1]);
         BUF(5:LBUF - MDAT(NROW) + 4) = BUF(MDAT(NROW)+1:LBUF)
         LBUF = LBUF - MDAT(NROW) + 4
         NROW = NROW + 1
         if (.not. GOTFMT) then
            if (NROW .gt. 5) go to 180
            MDAT(NROW) = LBUF
         end if
         if (NCOL .eq. 0) go to 960
         go to 940
      end if
      LBUF = 0
      if (.not. GOTFMT) then
         GOTFMT = .true.
         IROW = IROW - 1
         go to 920
      end if
      MACT(I+1) = MACT(I+1) + 1
      if (MACT(I+1) .le. MACT(I+2)) go to 999
      MACT(I+1) = 1
      if (NROW .eq. 1) go to 999
C%%    fputc(EOF, scratch_file);
      endfile SUNIT
      KK = 1
  994 KK = KK + 1
      if (KK .gt. NROW) go to 999
C%%   rewind(scratch_file);
      rewind(SUNIT)
      IROW = -1
      K = KK
  995 LBUF = 5
      IROW = IROW + 1
      if (IROW .ne. 0) then
C%%      sprintf(cmessc.buf, "%4ld",  irow%10000);
         write(BUF(1:4), '(I4)') mod(IROW, 10000)
      else
C%%    memset(cmessc.buf,' ',4);
         BUF(1:4) = ' '
      end if
      do 996 J = 2, K
         if (J .eq. K) LBUF = MDAT(KK)
C%%       if (fread(&cmessc.buf[4], cmessi.lbuf-4, 1,
C%%         scratch_file) == 0) goto L_994;
         read(SUNIT, END = 994) BUF(5:LBUF)
  996 continue
      K = NROW
      call MESSPR
      go to 995
  999 LENTRY = 1
      return
c
c                          Get started with vector or matrix output
 1000 INC = 1
      LOCBEG = NIDAT
      if (M .gt. MEGMAT) then
c Have a user set increment between entries of a vector.
        M = MEIVEC + 2 * (M - MEIVCI)
        I = I + 1
        INC = MACT(I)
      end if
      XARG = M .gt. MEJMAT
      if (XARG) then
         M = M - 4
         LOCBEG = NFDAT
         if (.not. XARGOK) go to 40
      end if
      GOTFMT = M .gt. MEIMAT
      if (GOTFMT) M = M - 2
      LOCBEG = LOCBEG + NSKIP
      NSKIP = 0
      MPT = LOCBEG
      if (M .eq. MEIMAT) go to 1300
c                           Take care of setup for vector output
      KNT = 0
      LASKNT = MACT(I+1)
      if (LASKNT .le. 0) then
         LASKNT = -LASKNT
         KNT = LOCBEG - 1
         if (LASKNT .le. KNT) go to 40
      end if
      IMAG = LASKNT
      LASTI = LOCBEG + INC * (LASKNT - 1 - KNT)
      NCOL = 0
c                          Get format for label output.
      call MESSFI
C%%   messcc.kcrwid = messcc.kciwid;
      FMTR = FMTI
      KDILAB = KDI+1
      LINSTR = 2*KDILAB+2
      if (XARG) then
         if (.not. GOTFMT) go to 1150
         IWF = IWG
         FMTF = FMTG
c++ CODE for .C. is inactive
C%%      cmessi.iwf = cmessi.iwg;
C%%      messcc.lfprec = messcc.lgprec;
c++ END
         go to 1160
      end if
      call MESSFD(IDAT)
c                                          After integer format
      LENOUT = KDI
      NIDAT = LASTI + 1
c                                          Common code continues here
 1080 NLINE = (LENLIN - LINSTR + 1) / LENOUT
      if (LBUF .eq. 0) go to 1090
      K = max(LINSTR, LBUF+1)
      if (((K-LINSTR)/LENOUT + (LENLIN-K+1)/LENOUT) .lt. NLINE) K = K +
     1   LENOUT - mod(K-LINSTR, LENOUT)
      KLINE = (LENLIN - K + 1) / LENOUT
      if (KLINE .lt. min(LASKNT-KNT, NLINE/2)) go to 1085
      LINSTR = K - LENOUT * ((K - LINSTR) / LENOUT)
      if (KLINE .ge. LASKNT-KNT)  then
         KLINE = LASKNT - KNT
         K = LBUF + 1
      end if
      KNT = KNT + KLINE
C%%    for (kc=cmessi.lbuf; kc < k; kc++) cmessc.buf[kc] = ' ';
      BUF(LBUF+1:K) = ' '
      LBUF = K
      go to 1110
 1085 call MESSPR
 1090 continue
c++ CODE for ~.C. is active
      BUF = ' '
      write (BUF(1:KDILAB), FMTR) KNT+1
c++ CODE for .C. is inactive
C%%   memset(cmessc.buf,' ',LENBUF);
C%%   sprintf(cmessc.buf, "%*ld", (int)messcc.kcrwid, knt+1);
c++ END
      BUF(KDILAB:KDILAB) = '-'
      KLINE = min(NLINE, LASKNT - KNT)
      KNT = KNT + KLINE
C%%    sprintf(&cmessc.buf[kdilab], "%*ld", (int)messcc.kcrwid, knt);
      write (BUF(KDILAB+1:2*KDILAB), FMTR) KNT
C%%    cmessc.buf[kdilab*2L-1] = ':';
C%%    for (kc=kdilab*2L; kc < *linstr-1; kc++) cmessc.buf[kc] = ' ';
      BUF(2*KDILAB:LINSTR-1) = ':'
      LBUF = LINSTR
 1110 LSTRT = LBUF
      LBUF = LBUF + LENOUT * KLINE - 1
      if (XARG) return
c                                    Integer output
c++ CODE for ~.C. is active
      write (BUF(LSTRT:LBUF), FMTI) (IDAT(K), K = MPT,
     1    MPT+INC*(KLINE-1), INC)
c++ CODE for .C. is inactive
C%%   for (k=cmessi.mpt; k<=cmessi.mpt+cmessi.kline-1; k++)
C%%  sprintf(&cmessc.buf[cmessi.lstrt+messcc.kciwid*(k-cmessi.mpt)-1],
C%%      "%*ld", (int)messcc.kciwid, idat[cmessi.inc*k-1]);
c++ END
      MPT = MPT + KLINE * INC
c
c                                     Entry here after vector output.
 1130 if (MPT .le. LASTI) go to 1085
      go to 40
c                                          Get other format
 1150 LENTRY = 5
      return
c                                          After other format
 1160 LENOUT = IWF
      LENTRY = 7
      NFDAT = LASTI + 1
      go to 1080


c                         Sparse vector output.
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

c                Entry after getting format for sparse data output.
 1210 LENOUT = IWF
      LENTRY = 11
      NLINE = LENLIN / IWF

 1220 call MESSPR
      KLINE = min(LASKNT - MPT + 1, NLINE)
      if (KLINE .le. 0) go to 40
      LBUF = LENOUT * KLINE
      return

c
c                           Take care of setup for matrix output
 1300 continue
      NDIM = MACT(I+1)
      if (NDIM .le. 0) then
         if (NDIM .eq. 0) go to 40
         INC = -NDIM
         NDIM = 1
      end if
      ICOL = 1
      IROW1 = 1
      NROW = MACT(I+2)
      if (NROW .le. 0) then
         if (NROW .eq. 0) go to 40
         IROW1 = -NROW / 100000
         NROW = -NROW - 99999 * IROW1 - 1
      end if
      NCOL = MACT(I+3)
      if (NCOL .le. 0) then
         if (NCOL .eq. 0) go to 40
         ICOL = -NCOL / 100000
         NCOL = -NCOL - 99999 * IROW1 - 1
      end if
      NTXTSV = NTEXT
      ITXTSV = ITEXT
      IRC = 1
c                        Compute widths for row and column labels
 1320 MAXWID(IRC) = 0
      MTEXT(IRC) = MACT(I+IRC+3)
      IMAG = NROCO(IRC)
      KLINE = IMAG
 1330 NTEXT = MTEXT(IRC)
      if (NTEXT .ge. 0) then
         if (NTEXT .eq. 0) then
            LTEXT = 5
         else
c                        Go get row/column widths
            KT = 2
            call MESSMH(TEXT)
            if (KT .lt. 0) then
               MTEXT(IRC) = 0
               go to 1330
            end if
         end if
         call MESSFI
         MAXWID(IRC) = max(MAXWID(IRC), LTEXT + KDI+1)
C%%      if (cmessi.irc == 1)
C%%         messcc.kcrwid = cmessi.kdi;
C%%      else
C%%         messcc.kccwid = cmessi.kdi;
         FMTIM(IRC) = FMTI
      end if
      IRC = IRC + 1
      if (IRC .eq. 2) go to 1320
c                 Widths for Row and column titles have been computed.
      KSHIFT = 1
      LASTI = LOCBEG + INC * (NROW - IROW1)
      if (XARG) then
         if (.not. GOTFMT) go to 1610
c++ CODE for ~.C. is active
         IWF = IWG
         FMTF = FMTG
c++ CODE for .C. is inactive
C%%      cmessi.iwf = cmessi.iwg;
C%%      messcc.lfprec = messcc.lgprec;
c++ END
         go to 1620
      end if
      call MESSFD(IDAT)
c
      If (KDI .ge. KOLWID) then
         LENOUT = KDI
      else
         KSHIFT = (KOLWID - KDI + 2) /2
         LENOUT = KOLWID
c++ CODE for ~.C. is active
         KDI = KOLWID
         FMTI(5:5) = char(ICHAR0 + KOLWID / 10)
         FMTI(6:6) = char(ICHAR0 + mod(KOLWID, 10))
c++ CODE for .C. is inactive
C%%  messcc.kciwid = *kolwid;
c++ END
      end if
      NIDAT = NIDAT + NDIM*NCOL
c                              Continue with commmon code
 1390 NLINE = (LENLIN - LINSTR) / LENOUT
      if (LBUF .le. LINSTR) go to 1420
 1400 call MESSPR
 1420 IROW = IROW1
      KLINE = min(NLINE, NCOL-ICOL+1)
c                       Output column labels (if any)
      if (MTEXTC .lt. 0) go to 1480
      NTEXT = MTEXTC
      IMAG = ICOL
      KT = 3
      call MESSMH(TEXT)
      if (KT .lt. 0) go to 180
c                       Return from output of column labels.
      MTEXTC = NTEXT
 1480 ICOL = ICOL + KLINE
 1490 call MESSPR
c
c                      Output row labels (if any)
      if (MTEXTR .lt. 0) go to 1520
      if (MTEXTR .eq. 0) then
C%%       memcpy(&cmessc.buf[cmessi.lbuf],"Row ", 4);
         BUF(LBUF+1:LBUF+4) = 'Row '
         LBUF = LBUF + 4
         go to 1515
      end if
      NTEXT = MTEXTR
      ITEXT = (NTEXT-1) / LENTXT
      NTEXT = NTEXT - ITEXT * LENTXT
      ITEXT = ITEXT + 1

c                     Go get text for row label
      NTEXTR = 0
      go to 410
c                     Return from getting text for row label
 1510 if (C .ne. '#') then
         MTEXTR = NTEXT + LENTXT * (ITEXT-1)
C%%    for (kc=cmessi.lbuf; kc < *linstr; kc++) cmessc.buf[kc] = ' ';
         BUF(LBUF+1:LINSTR) = ' '
         go to 1520
      end if
 1515 continue
C%%   sprintf(&cmessc.buf[cmessi.lbuf],"%*ld",(int)messcc.kcrwid,irow);
C%%    for (kc=cmessi.lbuf+messcc.kcrwid;
C%%       kc < *linstr; kc++) cmessc.buf[kc] = ' ';
      write (BUF(LBUF+1:LINSTR), FMTR) IROW
 1520 LSTRT = LINSTR + 1
      LBUF = LINSTR + LENOUT*KLINE
      LASTI = MPT + NDIM * KLINE - 1
      if (XARG) return
c                                    Integer output
C%% for (k=cmessi.mpt; k<=cmessi.lasti; k+=cmessi.ndim)
C%%  sprintf(&cmessc.buf[cmessi.lstrt + messcc.kciwid*(k-cmessi.mpt)/
C%%     cmessi.ndim - 1], "%*ld", (int)messcc.kciwid, idat[k-1]);
         write (BUF(LSTRT:LBUF), FMTI) (IDAT(K), K=MPT,LASTI,NDIM)
c
c                                     Entry here after matrix output.
 1530 MPT = MPT + INC
      IROW = IROW + 1
c
      if (IROW .le. NROW) go to 1490
      if (ICOL .gt. NCOL) then
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
c                                Need to get format for matrix print.
 1610 LENTRY = 6
      return
c                                Entry after got format for matrix print
 1620 If (IWF .ge. KOLWID) then
         LENOUT = IWF
      else
         KSHIFT = (KOLWID - IWF + 2) /2
         LENOUT = KOLWID
C%%      cmessi.iwf = *kolwid;
C%%      strcpy(cmessc.fmtf, "%*.*E\0");
         write (FMTF(7:8), '(I2)') KOLWID
      end if
      NFDAT = NFDAT + NDIM*NCOL
      LENTRY = 8
      go to 1390
      end

      subroutine MESSFD(IDAT)
c Get the format for data to be printed in vectors and arrays.
c
c ************** Variable only used here *******************************
c
c K      Temporary index.
c J      Temporary index.
c IDAT   Input array to MESS
c IMAX   Used when computing largest integer in array.
c IMIN   Used when computing smallest integer in array.
c
      integer J, K, IDAT(*), IMAX, IMIN
c
c For comments on other variables, see the listing for MESS.
      integer   LENBUF, MEVBAS, MEVLAS
      parameter (LENBUF=250)
      parameter (MEVBAS=10)
      parameter (MEVLAS=33)
      logical          XARG, GOTFMT, XARGOK
      integer          ERRCNT, EUNIT, ICHAR0, IRC, IVAR(MEVBAS:MEVLAS),
     1   IMAG, INC, ITEXT, IWF, IWG, KDF, KDFDEF, KDI, KDIAG, KDJ,
     2   KLINE, KSCRN, KSHIFT, KSPEC, KT, MAXERR, LASTI, LBUF, LENLIN,
     3   LENOUT, LENTRY, LENTXT, LHEAD, LINERR, LINMSG, LOCBEG, LPRINT,
     4   LSTOP, LSTRT, LTEXT, MAXWID(2), MDAT(5), MPT, MUNIT, NCOL,
     5   NDIM, NFDAT, NIDAT, NMDAT, NROW, NTEXT, OUNIT, SUNIT, TABSPA
c
      character BUF*(LENBUF), DOLS*72, FMTC*7, FMTF*20, FMTG*15,
     1  FMTI*7, FMTIM(2)*7, FMTJ*7, FMTR*7, FMTT*15
      common /CMESSI/ SUNIT, LHEAD, KDFDEF, LINMSG, LINERR, MUNIT,
     1   EUNIT, KSCRN, KDIAG, MAXERR, LSTOP, LPRINT, KDF, NTEXT, NIDAT,
     2   NFDAT, NMDAT, MDAT, TABSPA, ERRCNT, ICHAR0, IMAG, INC, IRC,
     3   ITEXT, IWF, IWG, KDI, KDJ, KLINE, KSHIFT, KSPEC, KT, LASTI,
     4   LBUF, LENLIN, LENOUT, LENTRY, LENTXT, LOCBEG, LSTRT, LTEXT,
     5   MAXWID, MPT, NROW, NCOL, NDIM, OUNIT, GOTFMT, XARG, XARGOK
      common /CMESSC / BUF, DOLS, FMTF, FMTG, FMTI, FMTJ, FMTT, FMTIM
      equivalence (IVAR(MEVBAS), SUNIT)
      equivalence (FMTIM(1), FMTR), (FMTIM(2), FMTC)
c
      save /CMESSI/, /CMESSC/
c
      if (GOTFMT) then
         KDI = KDJ
C%%      messcc.kciwid = cmessi.kdj;
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
      if (NCOL .ne. 0) then
         K = K + 1
         LOCBEG = LOCBEG + NDIM
         LASTI = LASTI + NDIM
         if (K .le. NCOL) go to 10
      end if
      IMAG = IMAX
      if ((IMAG/10) + IMIN .lt. 0) IMAG = IMIN
      KDI = -KDI
      call MESSFI
      return
      end


      subroutine MESSFI
c Get the format for the integer IMAG.
c
c ************** Variable only used here *******************************
c
c I, K, KD are used in determining number of characters needed to
c          represent IMAG.
c
      integer I, K, KD
c
c For comments on other variables, see the listing for MESS.
      integer   LENBUF, MEVBAS, MEVLAS
      parameter (LENBUF=250)
      parameter (MEVBAS=10)
      parameter (MEVLAS=33)
      logical          XARG, GOTFMT, XARGOK
      integer          ERRCNT, EUNIT, ICHAR0, IRC, IVAR(MEVBAS:MEVLAS),
     1   IMAG, INC, ITEXT, IWF, IWG, KDF, KDFDEF, KDI, KDIAG, KDJ,
     2   KLINE, KSCRN, KSHIFT, KSPEC, KT, MAXERR, LASTI, LBUF, LENLIN,
     3   LENOUT, LENTRY, LENTXT, LHEAD, LINERR, LINMSG, LOCBEG, LPRINT,
     4   LSTOP, LSTRT, LTEXT, MAXWID(2), MDAT(5), MPT, MUNIT, NCOL,
     5   NDIM, NFDAT, NIDAT, NMDAT, NROW, NTEXT, OUNIT, SUNIT, TABSPA
c
      character BUF*(LENBUF), DOLS*72, FMTC*7, FMTF*20, FMTG*15,
     1  FMTI*7, FMTIM(2)*7, FMTJ*7, FMTR*7, FMTT*15
      common /CMESSI/ SUNIT, LHEAD, KDFDEF, LINMSG, LINERR, MUNIT,
     1   EUNIT, KSCRN, KDIAG, MAXERR, LSTOP, LPRINT, KDF, NTEXT, NIDAT,
     2   NFDAT, NMDAT, MDAT, TABSPA, ERRCNT, ICHAR0, IMAG, INC, IRC,
     3   ITEXT, IWF, IWG, KDI, KDJ, KLINE, KSHIFT, KSPEC, KT, LASTI,
     4   LBUF, LENLIN, LENOUT, LENTRY, LENTXT, LOCBEG, LSTRT, LTEXT,
     5   MAXWID, MPT, NROW, NCOL, NDIM, OUNIT, GOTFMT, XARG, XARGOK
      common /CMESSC / BUF, DOLS, FMTF, FMTG, FMTI, FMTJ, FMTT, FMTIM
      equivalence (IVAR(MEVBAS), SUNIT)
      equivalence (FMTIM(1), FMTR), (FMTIM(2), FMTC)
c
      save /CMESSI/, /CMESSC/
c
      KD = 1
      if (KDI .lt. 0) then
c              KDI < 0 to flag need for extra space -- avoids overflows
         KDI = -KDI
         KD = 2
      end if
      K = 1
      if (IMAG .lt. 0) then
         IMAG = -IMAG
         KD = KD + 1
      end if
      I = IMAG / 10
      if (I .ne. 0) then
   10    K = 10 * K
         KD = KD + 1
         if (I .ge. K) go to 10
      end if
      if (KD .ne. KDI) then
         KDI = KD
c++ CODE for ~.C. is active
         FMTI(5:5) = char(ICHAR0 + KDI / 10)
         FMTI(6:6) = char(ICHAR0 + mod(KDI, 10))
c++ CODE for .C. is inactive
C%%      messcc.kciwid = cmessi.kdi;
c++ END
      end if
      return
      end

c++ CODE for ~.C. is active
      integer function MESSGS()
c                                 Get a scratch unit assigned.
      integer J
c
      MESSGS = 31
   10 MESSGS = MESSGS - 1
      if (MESSGS .eq. 0) stop 'Could not assign scratch unit in MESS.'
      open (MESSGS, STATUS='SCRATCH', ACCESS='SEQUENTIAL',
     1    FORM='UNFORMATTED', IOSTAT=J)
      if (J .ne. 0) go to 10
      return
      end
c++ END

      subroutine MESSMH(TEXT)
c Processing of multiple headings:
c
c J     Used as a temporary index.
c K     Used as a temporary index.
c KB    Number of initial blanks
c KK    Used as a temporary index.
c KT    Used for logic in output of headings.  Set <0 on exit if there
c       is an input error.
c     KT = 1 Output table headings.  (Set to -1 on fatal error.)
c     KT = 2 Get row/column widths for matrix output.  (Set to -2 if
c            error results in no headings.)
c     KT = 3 Output column headings.  (Set to -1 on fatal error.)
c L     Used as a temporary index.
c LFLGDB 2 in usual case, 3 if see a "$ ".
c LSTRDB Value of LSTRT when see a "$ ".
c LTXTDB Value of LTEXT when see a "$ ".
c TEXT  Original input character vector.
c
      integer J, K, KB, KK, L, LFLGDB, LSTRDB, LTXTDB
      character*(*)  TEXT(*)
      character SC, C
      parameter (SC='$')
c For comments on other variables, see the listing for MESS.
      integer   KOLWID, LINSTR, LENBUF, MEVBAS, MEVLAS
      parameter (LENBUF=250)
      parameter (MEVBAS=10)
      parameter (MEVLAS=33)
      logical          XARG, GOTFMT, XARGOK
      integer          ERRCNT, EUNIT, ICHAR0, IRC, IVAR(MEVBAS:MEVLAS),
     1   IMAG, INC, ITEXT, IWF, IWG, KDF, KDFDEF, KDI, KDIAG, KDJ,
     2   KLINE, KSCRN, KSHIFT, KSPEC, KT, MAXERR, LASTI, LBUF, LENLIN,
     3   LENOUT, LENTRY, LENTXT, LHEAD, LINERR, LINMSG, LOCBEG, LPRINT,
     4   LSTOP, LSTRT, LTEXT, MAXWID(2), MDAT(5), MPT, MUNIT, NCOL,
     5   NDIM, NFDAT, NIDAT, NMDAT, NROW, NTEXT, OUNIT, SUNIT, TABSPA
c
      character BUF*(LENBUF), DOLS*72, FMTC*7, FMTF*20, FMTG*15,
     1  FMTI*7, FMTIM(2)*7, FMTJ*7, FMTR*7, FMTT*15
      common /CMESSI/ SUNIT, LHEAD, KDFDEF, LINMSG, LINERR, MUNIT,
     1   EUNIT, KSCRN, KDIAG, MAXERR, LSTOP, LPRINT, KDF, NTEXT, NIDAT,
     2   NFDAT, NMDAT, MDAT, TABSPA, ERRCNT, ICHAR0, IMAG, INC, IRC,
     3   ITEXT, IWF, IWG, KDI, KDJ, KLINE, KSHIFT, KSPEC, KT, LASTI,
     4   LBUF, LENLIN, LENOUT, LENTRY, LENTXT, LOCBEG, LSTRT, LTEXT,
     5   MAXWID, MPT, NROW, NCOL, NDIM, OUNIT, GOTFMT, XARG, XARGOK
      common /CMESSC / BUF, DOLS, FMTF, FMTG, FMTI, FMTJ, FMTT, FMTIM
      equivalence (IVAR(MEVBAS), SUNIT)
      equivalence (FMTIM(1), FMTR), (FMTIM(2), FMTC)
c
      equivalence (MAXWID(1), LINSTR), (MAXWID(2), KOLWID)

      save /CMESSI/, /CMESSC/
c++ CODE for .C. is inactive
C%%      long int kc;
c++ END
      save LFLGDB
      data LFLGDB / 2 /
c
      if (NTEXT .ne. 0) then
         ITEXT = (NTEXT-1) / LENTXT
         NTEXT = NTEXT - ITEXT * LENTXT
         ITEXT = ITEXT + 1
      end if
      do 300 J = 1, max(1,KLINE)
         if (NTEXT .eq. 0) then
            K = KOLWID
            go to 210
         end if
         LFLGDB = 2
         LTEXT = 0
  110    continue
C%%       ctmp=memchr(TEXT(cmessi.itext-1L,cmessi.ntext-1), SC,
C%%          cmessi.lentxt - cmessi.ntext + 1);
C%%       if (ctmp == NULL)
C%%             l = 0;
C%%       else
C%%             l = ctmp - TEXT(cmessi.itext-1L,cmessi.ntext-1) + 1;
         L = index(TEXT(ITEXT)(NTEXT:LENTXT), SC)
         if (L .eq. 0) then
            LTEXT = LTEXT + LENTXT - NTEXT + 1
            if (LTEXT .lt. 80) then
               ITEXT = ITEXT + 1
               NTEXT = 1
               go to 110
            end if
            LTEXT = 0
            if (KT .eq. 3) go to 310
            go to 160
         end if
         NTEXT = NTEXT + L + 1
         LTEXT = L + LTEXT - 1
         if (NTEXT .gt. LENTXT) then
            ITEXT = ITEXT + 1
            if (NTEXT .eq. LENTXT + 1) then
               C = TEXT(ITEXT-1)(LENTXT:LENTXT)
               NTEXT = 1
            else
               C = TEXT(ITEXT)(1:1)
               NTEXT = 2
            end if
         else
            C = TEXT(ITEXT)(NTEXT-1:NTEXT-1)
         end if
         if (C .eq. 'H') go to (180, 190, 200), KT
         if (C .eq. 'E') go to (180, 310, 200), KT
         if (C .eq. '#') go to (140, 310, 200), KT
         if (C .eq. ' ') then
c  Special code to set for removing the "$" preceding a blank.
            LSTRDB = LSTRT
            LFLGDB = 3
            LTXTDB = LTEXT
            LTEXT = LTEXT + 1
            go to 110
         end if
         if (KT .ne. 1) go to 160
  140    LTEXT = LTEXT + 2
         go to 110
  160    KT = -KT
         go to 310
c
  180    KOLWID = KOLWID + LENOUT
         if (LTEXT .eq. 0) go to 300
         KB = KOLWID-LTEXT
         if (KB .lt. 0) stop
     1   'Stopped in MESS -- Column width too small in a heading.'
         if (XARG)  KB = 1 + KB/2
         LSTRT = LBUF + KB + 1
         LBUF = LBUF + KOLWID
         if (LBUF .le. LENLIN) MDAT(NROW) = LBUF
         KOLWID = 0
         go to 220
c
c                                  Set up column widths
  190    MAXWID(IRC) = max(MAXWID(IRC), LTEXT)
         go to 300
c
c                                  Output matrix column
  200    K = KOLWID
         if (C .ne. '#') K = LTEXT
  210    KB = LENOUT - KOLWID
         if (J .eq. 1) then
c                        Special setup for the first column.
            if (XARG) KB = (KB + 1) / 2
            KB = KB + KSHIFT + LINSTR - LBUF
         end if
         KB = KB + KOLWID - K
         LSTRT = LBUF + KB + 1
         LBUF = LSTRT + K - 1
c                                  Set initial blanks
  220    continue
C%%      if (kb > 0) for (kc=cmessi.lstrt-kb-1; kc<cmessi.lstrt-1; kc++)
C%%         cmessc.buf[kc] = ' ';
         if (KB .gt. 0) BUF(LSTRT-KB:LSTRT-1) = ' '
c                                  Move characters
         if (NTEXT .eq. 0) then
C%%       memcpy(&cmessc.buf[cmessi.lstrt-1],"Col ", 4);
            BUF(LSTRT:LSTRT+3) = 'Col '
            C = '#'
            LSTRT = LSTRT+4
         else
            K = NTEXT - LTEXT - LFLGDB
            if (K .le. 0) then
               KK = max(0, 3-NTEXT)
C%%       memcpy(&cmessc.buf[cmessi.lstrt-1L], TEXT(cmessi.itext-2L,
C%%         cmessi.lentxt+k-1), -k-kk+1L);
               BUF(LSTRT:LSTRT-K-KK)=TEXT(ITEXT-1)(LENTXT+K:LENTXT-KK)
               LSTRT = LSTRT-K-KK+1
               K = 1
            end if
            if (NTEXT .gt. 3) then
C%%       memcpy(&cmessc.buf[cmessi.lstrt-1L], TEXT(cmessi.itext-1L,
C%%         k-1), cmessi.ntext-k-2L);
               BUF(LSTRT:LSTRT+NTEXT-K-3) = TEXT(ITEXT)(K:NTEXT-3)
               LSTRT = LSTRT + NTEXT - K - 2
            end if
         end if
         if (LFLGDB .eq. 3) then
c  Special code to remove the "$" preceding a blank.  Only works for 1.
            do 250 L = LSTRDB + LTXTDB + max(0, KB), LSTRT
               BUF(L:L) = BUF(L+1:L+1)
  250       continue
            LFLGDB = 2
            LSTRT = LSTRT - 1
         end if
         if (C .eq. '#') then
c                                  Output column index
C%%         sprintf(&cmessc.buf[cmessi.lstrt-1], "%*ld ",
C%%           (int)(cmessi.lbuf-cmessi.lstrt), cmessi.imag+j-1);
            write (BUF(LSTRT:LBUF), FMTC) IMAG + J - 1
            if (NTEXT .ne. 0) NTEXT = K
            go to 300
         end if
c                                  Set trailing blanks
C%%      if (cmessi.lstrt <= cmessi.lbuf)
C%%           for (kc=cmessi.lstrt-1; kc < cmessi.lbuf; kc++)
C%%              cmessc.buf[kc] = ' ';
         if (LSTRT .le. LBUF) BUF(LSTRT:LBUF) = ' '
  300 continue
  310 return
      end

      subroutine MESSPR
c Prints the buffer for MESS
c
c ************** Variable only used here *******************************
c
c NSCRN  Number of lines currently on CRT from messages.
c
      integer   NSCRN, K
      character SCRNAM*12
      save      NSCRN
c
c For comments on other variables, see the listing for MESS.
      integer   LENBUF, MEVBAS, MEVLAS
      parameter (LENBUF=250)
      parameter (MEVBAS=10)
      parameter (MEVLAS=33)
      logical          XARG, GOTFMT, XARGOK
      integer          ERRCNT, EUNIT, ICHAR0, IRC, IVAR(MEVBAS:MEVLAS),
     1   IMAG, INC, ITEXT, IWF, IWG, KDF, KDFDEF, KDI, KDIAG, KDJ,
     2   KLINE, KSCRN, KSHIFT, KSPEC, KT, MAXERR, LASTI, LBUF, LENLIN,
     3   LENOUT, LENTRY, LENTXT, LHEAD, LINERR, LINMSG, LOCBEG, LPRINT,
     4   LSTOP, LSTRT, LTEXT, MAXWID(2), MDAT(5), MPT, MUNIT, NCOL,
     5   NDIM, NFDAT, NIDAT, NMDAT, NROW, NTEXT, OUNIT, SUNIT, TABSPA
c
      character BUF*(LENBUF), DOLS*72, FMTC*7, FMTF*20, FMTG*15,
     1  FMTI*7, FMTIM(2)*7, FMTJ*7, FMTR*7, FMTT*15
      common /CMESSI/ SUNIT, LHEAD, KDFDEF, LINMSG, LINERR, MUNIT,
     1   EUNIT, KSCRN, KDIAG, MAXERR, LSTOP, LPRINT, KDF, NTEXT, NIDAT,
     2   NFDAT, NMDAT, MDAT, TABSPA, ERRCNT, ICHAR0, IMAG, INC, IRC,
     3   ITEXT, IWF, IWG, KDI, KDJ, KLINE, KSHIFT, KSPEC, KT, LASTI,
     4   LBUF, LENLIN, LENOUT, LENTRY, LENTXT, LOCBEG, LSTRT, LTEXT,
     5   MAXWID, MPT, NROW, NCOL, NDIM, OUNIT, GOTFMT, XARG, XARGOK
      common /CMESSC / BUF, DOLS, FMTF, FMTG, FMTI, FMTJ, FMTT, FMTIM
      equivalence (IVAR(MEVBAS), SUNIT)
      equivalence (FMTIM(1), FMTR), (FMTIM(2), FMTC)
c
      save /CMESSI/, /CMESSC/
      data NSCRN / 0 /
c
      if (LBUF .ne. 0) then
 10     if (BUF(LBUF:LBUF) .eq. ' ') then
          if (LBUF .gt. 1) then
            LBUF = LBUF - 1
            go to 10
          end if
        end if
        if (OUNIT .le. 0) then
          if (KSCRN .gt. 0) then
            if (NSCRN .ge. KSCRN) then
C%%               printf( " Type 'Enter' to continue\n" );
              print '('' Type "Enter" to continue'')'
C%%               scanf( "%*[^\n]%*c" );
              read (*, *)
              NSCRN = 0
            end if
            NSCRN = NSCRN + 1
          end if
C%%      printf( "%.*s\n", (int)cmessi.lbuf, cmessc.buf);
          print '(1X, A)', BUF(1:LBUF)
          if (OUNIT .eq. 0) go to 20
        end if
c++ CODE for ~.C. is active
        K = abs(OUNIT)
        write (K, '(A)', ERR=30) BUF(1:LBUF)
c++ CODE for .C. is inactive
C%%      fprintf(c_handle[labs(cmessi.ounit)-1], "%.*s\n",
C%%      (int)cmessi.lbuf, cmessc.buf);
c++ END
 20     LBUF = 0
      end if
      return
c++ CODE for ~.C. is active
c              See if opening fixes the error
30    write(SCRNAM, '(A, I2.2, A)') 'MESSF_', K, '.tmp'
      open (UNIT=K, STATUS='UNKNOWN', FILE=SCRNAM)
      write (K, '(A)') BUF(1:LBUF)
      return
c++ END
      end

      subroutine MESSFT(MACT, FTEXT)
c  Prints FTEXT, which contains a Fortran character string, and then
c  call MESS to do the actions in MACT.  Actions in MACT can not do
c  anything other than actions that reference MACT.
c  This routine intended for use by library subroutines getting text in
c  the form of a Fortran character string.
c
      integer MACT(*)
      character FTEXT*(*)
c
      integer J, K, IDAT(1), MECONT, MEPRNT, MESUNI
c++ CODE for ~.C. is active
      character TEXT(1)*1
c++ CODE for .C. is inactive
C      character TEXT(1)*2
c++ END
      parameter (MESUNI=10, MEPRNT=21, MECONT=50)
c
      INTEGER LENBUF, MEVBAS, MEVLAS
      parameter (LENBUF=250)
      parameter (MEVBAS=10)
      parameter (MEVLAS=33)
      logical          XARG, GOTFMT, XARGOK
      integer          ERRCNT, EUNIT, ICHAR0, IRC, IVAR(MEVBAS:MEVLAS),
     1   IMAG, INC, ITEXT, IWF, IWG, KDF, KDFDEF, KDI, KDIAG, KDJ,
     2   KLINE, KSCRN, KSHIFT, KSPEC, KT, MAXERR, LASTI, LBUF, LENLIN,
     3   LENOUT, LENTRY, LENTXT, LHEAD, LINERR, LINMSG, LOCBEG, LPRINT,
     4   LSTOP, LSTRT, LTEXT, MAXWID(2), MDAT(5), MPT, MUNIT, NCOL,
     5   NDIM, NFDAT, NIDAT, NMDAT, NROW, NTEXT, OUNIT, SUNIT, TABSPA
c
      character BUF*(LENBUF), DOLS*72, FMTC*7, FMTF*20, FMTG*15,
     1  FMTI*7, FMTIM(2)*7, FMTJ*7, FMTR*7, FMTT*15
      common /CMESSI/ SUNIT, LHEAD, KDFDEF, LINMSG, LINERR, MUNIT,
     1   EUNIT, KSCRN, KDIAG, MAXERR, LSTOP, LPRINT, KDF, NTEXT, NIDAT,
     2   NFDAT, NMDAT, MDAT, TABSPA, ERRCNT, ICHAR0, IMAG, INC, IRC,
     3   ITEXT, IWF, IWG, KDI, KDJ, KLINE, KSHIFT, KSPEC, KT, LASTI,
     4   LBUF, LENLIN, LENOUT, LENTRY, LENTXT, LOCBEG, LSTRT, LTEXT,
     5   MAXWID, MPT, NROW, NCOL, NDIM, OUNIT, GOTFMT, XARG, XARGOK
      common /CMESSC / BUF, DOLS, FMTF, FMTG, FMTI, FMTJ, FMTT, FMTIM
      equivalence (IVAR(MEVBAS), SUNIT)
      equivalence (FMTIM(1), FMTR), (FMTIM(2), FMTC)
c
      do 10 J = 1, 100, 2
         K = abs(MACT(J))
         if ((K .gt. MEPRNT) .or. (K .lt. MESUNI)) go to 20
  10  continue
  20  K = MACT(J)
      MACT(J) = MECONT
      call MESS(MACT, TEXT, IDAT)
      MACT(J) = K
C%%      k = strlen(ftext);
      K = len(FTEXT)
      NTEXT = 1
      if (K .ne. 0) then
         if (FTEXT(1:1) .eq. '0') then
            NTEXT = 2
            K = K - 1
            if (LBUF .eq. 0) then
               BUF(1:1) = ' '
               LBUF = 1
            end if
         end if
         call MESSPR
         LBUF = K
C%%      memcpy(cmessc.buf, &ftext[cmessi.ntext-1], k);
         BUF(1:K) = FTEXT(NTEXT:NTEXT+K-1)
      end if
      ICHAR0 = ICHAR('0')
      if (MACT(J) .ne. MECONT) call mess(MACT(J), TEXT, IDAT)
      return
      end
      subroutine OPTCHK(INTCHK, IOPT, ETEXT)
c Copyright (c) 1996 California Institute of Technology, Pasadena, CA.
c ALL RIGHTS RESERVED.
c Based on Government Sponsored Research NAS7-03001.
c>> 1998-11-01 OPTCHK  Krogh  ERRSEV => MACT(2) for "mangle".
c>> 1996-05-13 OPTCHK  Krogh  Changes to use .C. and C%%.
C>> 1995-03-10 OPTCHK  Krogh  Added "abs(.) just below "do 140 ..."
C>> 1994-11-11 OPTCHK  Krogh  Declared all vars.
c>> 1993-05-17 OPTCHK  Krogh  Additions for Conversion to C.
c>> 1991-11-25 OPTCHK  Krogh  More comments, little clean up of code.
c>> 1991-10-09 OPTCHK  Krogh  More comments, little clean up of code.
c>> 1991-06-27 OPTCHK  Krogh  Initial Code.
c
c OPTCHK -- Fred T. Krogh, Jet Propulsion Lab., Pasadena, CA.
c This subroutine is intended for the use of other library routines.
c It is used to check the storage used by options for some array.
c
c INTCHK is an array that provides information on how storage has been
c   allocated.  Information is provided in groups of three words, after
c   an initial group of 4 that contain:
c    0. INTCHK(0) / 10 is the index to use for error messages, and
c       mod(INTCHK(0), 10) defines actions as follows:
c       Value  Action on Errors  Action if no errors
c         0    Prints & Stops               Returns
c         1    Prints & Returns             Returns
c         2    Prints & Stops      Prints & Returns
c         3    Prints & Returns    Prints & Returns
c        >3  Any error message will be continued, subtract 4 from the
c            value to get one of the actions above.
c    1. Contains LAST, 1 more than the index of the last location in
c       INTCHK used to store the information as described below.
c    2. The amount of storage available in the array being checked.  If
c       this is 0, it is assumed that the user would like to know what
c       size to declare for the array, and if it is -1, it is assumed
c       that a library routine is doing the check without knowing what
c       the declared size is.
c    3. The amount of storage required in this array if no options were
c       used.  (This + 1 is first loc. available for option storage.)
c   The rest should be set as follows, for k = 5, LAST, 3:
c    k-1 Option index, 0 indicates storage required by default.  (This
c        may depend on input parameters, but doesn't depend explicitly
c        on an option.)
c    k   If >0, gives the start of space used by the option.
c        If <0, gives -(amount of space required by an option), for
c               which a starting location is to be determined.
c    k+1 If preceding entry is >0, gives space required by the option.
c        Else it is assumed that the space requested is to be found,
c        and a diagnostic will be given if space can not be found.  Else
c        INTCHK(K+1) IOPT(INTCHK(k+1)) Diagnostic on successful alloc.?
c            0           ----          No
c            <0          ----          Yes
c            >0          .ne. -1       Yes
c            >0          .eq. -1       No, and IOPT(INTCHK(k+1)) is set
c                                      to starting loc. of space found.
c        When this program finds the space for an option, values of
c        INTCHK(j) for j .ge. LAST will be used.  INTCHK(k+1) is set
c        temporarily to -(space required) and INTCHK(k-1) is reduced by
c        2000000 to flag the fact that the location index must be saved
c        after INTCHK(LAST).  INTCHK(LAST) is assumed to contain 1 if
c        the largest space requested is required, and to contain
c        -(minimal space needed) if the large amount requested is not
c        essential.
c   On exit, INTCHK(1) is set = -LAST if there was some kind of error.
c   (Unless the user has called MESS this return can't happen.)
c   INTCHK(2) is set to suggest a value for the storage to be declared.
c   The remaining locations are changed as follows:
c    k-1 Negated if there was a problem with this option, or if the
c        following location was initially <0.
c    k   Starting location used or suggested for the option.
c    k+1 Last location used or suggested for the option.
c
c        In addition if values of INTCHK(j) for j .ge. LAST are used,
c        INTCHK(j), for j = LAST+1, LAST+2, ..., will be set so that
c        INTCHK(j) is equal to (one of the k's above) - 1 for an option
c        that has had a starting location assigned, and INTCHK(LAST) is
c        set to the last index value for j in the above list.
c
c IOPT  This is stored to if INTCHK(k) is < 0, see above.
c ETEXT Input text of the form 'Package_name / Argument_name$E'.  Thus
c   for argument KORD of the package DIVA, this would = 'DIVA / KORD$E'.
c
c ************************** Variable definitions **********************
c
c ERRBAD Flag to use for error severity if get an error.  17 if errors
c        are to print but not stop, 57 otherwise.
c ETEXT  (Input) Used to construct error messages, see above.
c I      Induction variable for accessing INTCHK().
c INTCHK (InOut) See above.
c IOPT   (Out) If space information is being obtained results are saved
c        here.  See above.
c ISTRT  Starting index for placing options with unknown locations.
c KEX    Points to the last place in INTCHK where indices of entries
c        that have had locations determined here are stored.
c L      Temporary index.
c LAST   First free location in INTCHK.  (=INTCHK(1))
c LNEG   Index of the first INTCHK entry that is to be positioned.
c LOPT   Location in IOPT to get location of INTCHK entry that is
c        positioned.
c LTXTAx Variables setup by PMESS defining the locations in MTXTAA where
c        various error messages start.
c LTXTEE Location in MTXTAA where data in ETEXT is stored.
c LTXTZZ Number of characters available for saving ETEXT in MTXTAA.
c LWANT  -Number of locations wanted by INTCHK entry being positioned.
c MACT   Vector used to specify error printing actions, see MESS.
c        MACT(2) flags error/diagnostic level.  = 0 none; 07 is used to
C        get diagnostics only; and ERRBAD otherwise.
c MEEMES Parameter specifying that error message is to be printed.
c MEIDAT Parameter specifying location in INTCHK to start printing
c        integer in MESS.
c MEIMAT Parameter specifying an integer matrix is to be printed.
c MENTXT Parameter specifying location in INTCHK to start printing
c        text from MTXTAA in MESS.
c MERET  Parameter specifying the end of an error message.
c MESS   Routine to print error and diagnostc messages.
c METEXT Parameter specifying that text is to be printed.
c MI     Temporary index used to save in acceptable location.
c MTXTAA Used to contain error message text and instructions, see MESS.
c MTXTAx Character variables setup by PMESS and equivalenced into ETEXT
c        used to contain parts of the error messages.
c MTXTZZ As for MTXTAx, except not setup by PMESS.  Used to hold text
c        from ETEXT.
c MV     Temporary, which contains value associated with INTCHK(MI).
c N      Temporary value.
c NERBAD Array telling what to do concerning errrors.  ERRBAD is set
c        from NERBAD(mod(INTCHK(0), 10)), and the default value for
c        MACT(2) is set from NERBAD(INTCHK(0)+4).
c
c ************************** Variable Declarations *********************
c
      integer INTCHK(0:*), IOPT(*)
      character ETEXT*(*)
      integer I, ISTRT, KEX, L, LAST, LNEG, LOPT, LWANT, MI, MV, N
c Declarations for error messages.
      integer MENTXT, MEIDAT, MECONT, MERET, MEEMES, METEXT, MEIMAT,
     1   LTXTEE, LTXEND
      parameter (MENTXT =23)
      parameter (MEIDAT =24)
      parameter (MECONT =50)
      parameter (MERET =51)
      parameter (MEEMES =52)
      parameter (METEXT =53)
      parameter (MEIMAT =58)
      integer MACT(16), ERRBAD, NERBAD(0:7)
c
c ********* Error message text ***************
c[Last 2 letters of Param. name]  [Text generating message.]
cAA OPTCHK$B
cAB "Option #" is negated if option needs attention.$N
c   "Option 0" is for space not associated with a specific option.$N
c   "First Loc." is negated if user did not set value.$N
c   Space avail. = $I; all options must have first loc. > $I$E
cAC Option #$HFirst Loc.$HLast Loc.$E
cAD From subprogram/argument: $B
cAE Space for ETEXT.$E
      integer LTXTAA,LTXTAB,LTXTAC,LTXTAD,LTXTAE
      parameter (LTXTAA=  1,LTXTAB=  9,LTXTAC=233,LTXTAD=267,LTXTAE=295)
      character MTXTAA(2) * (156)
c                          Next 4 lines not automatically generated
c%%    #define LTXTEE  137
      parameter (LTXTEE = LTXTAE - 156 - 2)
      parameter (LTXEND = 156)
C
      data MTXTAA/'OPTCHK$B"Option #" is negated if option needs attenti
     *on.$N"Option 0" is for space not associated with a specific option
     *.$N"First Loc." is negated if user di','d not set value.$NSpace av
     *ail. = $I; all options must have first loc. > $I$EOption #$HFirst$
     * Loc.$HLast Loc.$EFrom subprogram/argument: $BSpace for ETEXT.$E'/
c
c                      1 2 3      4       5  6       7      8       9
      data MACT / MEEMES,0,1,LTXTAD, MEIDAT, 2, MENTXT,LTXTAB, METEXT,
     1   MEIMAT,3,3,0,LTXTAC,-1, MERET /
c            10    13     14 15     16
      data NERBAD / 57, 17, 57, 17, 0, 0, 7, 7 /
c
c *************************** Start of Executable Code *****************
c
      MACT(3) = INTCHK(0) / 10
      MACT(16)=MERET
      I = INTCHK(0) - 10*MACT(3)
      if (I .gt. 3) then
         I = I - 4
         MACT(16) = MECONT
      end if
      ERRBAD = NERBAD(I)
      MACT(2) = NERBAD(I+4)
      LAST = INTCHK(1)
      KEX = LAST
   20 LNEG = 0
      do 100 I = 5, LAST, 3
c       Loop to sort on the low indices -- Inefficient algorithm to keep
c       code short -- LAST should never be very big.
         MI = I
         MV = INTCHK(I)
         do 50 L = I+3, LAST, 3
c                                    Find mimimum from here down.
            if (INTCHK(L) .lt. MV) then
               MI = L
               MV = INTCHK(L)
            end if
   50    continue
         if (MI .ne. I) then
c                                   Interchange to get low at top.
            do 70 L = -1, 1
               N = INTCHK(I+L)
               INTCHK(I+L) = INTCHK(MI+L)
               INTCHK(MI+L) = N
   70       continue
         end if
         if (MV .lt. 0) then
c                Save loc. of last entry that needs space to be found.
            LNEG = I
         else if (LNEG .eq. 0) then
c        Entry I and previous entries are in their correct sorted order.
            if (INTCHK(I+1) .lt. 0) then
               if (INTCHK(I-1) .lt. -1000000) then
                  INTCHK(I-1) = INTCHK(I-1) + 2000000
                  INTCHK(I+1) = -INTCHK(I+1)
c                            Save INTCHK index defining allocated space.
                  KEX = KEX + 1
                  INTCHK(KEX) = I - 1
               else
c                   Error -- Got request for a negative amount of space.
                  MACT(2) = ERRBAD
                  INTCHK(I-1) = -abs(INTCHK(I-1))
               end if
            end if
c                Save final location used by the option.
            INTCHK(I+1) = INTCHK(I) + INTCHK(I+1) - 1
            if (INTCHK(I) .le. INTCHK(I-2)) then
c                                           Error -- options overlap.
               INTCHK(I-1) = -abs(INTCHK(I-1))
               MACT(2) = ERRBAD
            end if
         end if
  100 continue
      if (LNEG .ne. 0) then
c     Find spaces that need to be allocated, starting with the smallest.
         ISTRT = LNEG
         I = LNEG
  120    LWANT = INTCHK(LNEG)
         LOPT = INTCHK(LNEG+1)
         if (I .eq. LNEG) then
c                         Make fake entry to get started.
            INTCHK(LNEG) = 1
            INTCHK(LNEG+1) = INTCHK(3)
         end if
         do 140 ISTRT = ISTRT, LAST-3, 3
            if(INTCHK(I)+abs(INTCHK(I+1))-LWANT .lt. INTCHK(ISTRT+3))
     1         go to 150
            I = ISTRT + 3
  140    continue
  150    INTCHK(LNEG) = INTCHK(I) + abs(INTCHK(I+1))
         if (LOPT .ne. 0) then
            if (LOPT .gt. 0) then
               if (IOPT(LOPT) .eq. -1) then
                  IOPT(LOPT) = INTCHK(LNEG)
                  go to 160
               end if
            end if
c                     Error -- IOPT not expecting starting loc.
            INTCHK(LNEG-1) = -abs(INTCHK(LNEG-1))
            MACT(2) = ERRBAD
         end if
  160    INTCHK(LNEG+1) = LWANT
         INTCHK(LNEG-1) = INTCHK(LNEG-1) - 2000000
         if (LNEG .lt. 8) go to 20
         I = LNEG
         LNEG = LNEG - 3
         go to 120
      end if
      if (INTCHK(LAST-1) .gt. INTCHK(2)) then
         if (INTCHK(2) .lt. 0) go to 180
         if (LAST .ne. KEX) then
            if (INTCHK(KEX) .eq. LAST - 3) then
               if (INTCHK(LAST) .le. 0) then
                  if (INTCHK(LAST-2)-INTCHK(LAST)-1 .le. INTCHK(2)) then
                     INTCHK(LAST-1) = INTCHK(2)
                     go to 180
                  end if
               end if
            end if
         end if
         INTCHK(LAST-3) = -abs(INTCHK(LAST-3))
         MACT(2) = ERRBAD
      end if
  180 if (LAST .ne. KEX) INTCHK(LAST) = KEX
      if (MACT(2) .gt. 0) then
  190    if (LAST .ne. KEX) then
            do 200 I = LAST+1, abs(KEX)
               INTCHK(INTCHK(I)+1) = -INTCHK(INTCHK(I)+1)
  200       continue
            if (KEX .lt. 0) go to 210
            KEX = -KEX
         end if
         MACT(13) = (LAST - 4) / 3
c%%       strcpy(&mtxtaa[1][LTXTEE-1], etext);
         MTXTAA(2)(LTXTEE:LTXEND)=ETEXT(1:)
         call MESS(MACT, MTXTAA, INTCHK(1))
         if (MACT(2) .gt. 10) INTCHK(1) = -LAST
         if (LAST .ne. KEX) go to 190
      end if
  210 INTCHK(2) = INTCHK(LAST-1)
      return
      end