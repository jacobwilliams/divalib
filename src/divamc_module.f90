!*************************************************************************
!>
!  The main common block for the package.

    module divamc_module

    use diva_constants

    implicit none

    public

    integer  :: icf !! Final index for current loop in DIVACR.  Required by
                    !! option 18.
    integer  :: ics !! Starting index for current loop in DIVACR.
    integer  :: idat(6) !! Used to store integer for error messages.  (Also used
                        !! in DIVAA for temporary storage of KORD(2).  (Local array in DIVAIN.)
    integer  :: igflg !! Used primarily in DIVAg, but also used in DIVA to keep
                      !! track of the state of GSTOP calculations:
                      !!
                      !!  * -2 Extrapolatory G's initialized, but not the interpolatory.
                      !!  * -1 Interpolatory G's initialized, but not the extrapolatory.
                      !!  *  0 Set when integration is started or restarted, or option setting GSTOP is set.
                      !!  *  1 Iterating to find a GSTOP.
                      !!  *  2 User told that a GSTOP was found.
                      !!  *  3 Checking G's at point where a GSTOP was located.
                      !!  *  4 Checking G's at a T output point.
                      !!  *  5 Usual case, no sign change detected.
    integer  :: igflgs
    integer  :: igstop(2) !! IGSTOP(k) is set in DIVAg to the index of the last G
                          !! with a 0, where k is one for an interpolatory G-Stop, and k is two
                          !! for an extrapolatory G-Stop.
    integer  :: igtype(2) !! Array with two elements as for IGSTOP, but this saves
                          !! a flag giving the nature of convergence to the stop:
                          !!
                          !!  * 0  All known G-stops completely processed.
                          !!  * 4  Need to compute next value while iterating.
                          !!  * 5  Got good convergence.
                          !!  * 6  Got convergence, but not to desired accuracy.
                          !!  * 7  Problem in getting convergence.
                          !!  * 8  A fatal error of some type.
    integer  :: ilgrep !! Used when correction to keep track of equations that
                       !! are to use a certain error tolerance.
    integer  :: ings
    integer  :: iop10 !! Number of times diagnostic output is to be given when
                      !! leaving DIVAcr (the corrector).
    integer  :: iop11 !! Gives current step number of the method.  Tells how
                      !! many of certain coefficients must be computed. (Has nothing to do
                      !! with options.) = min(max integ order + 1, KDIM).  Also set when
                      !! starting to flag that certain memory locations must be set to 0.
    integer  :: iop12 !! Points to location in F() where user supplied values
                      !! of HINC, HDEC, HMIN, and HMAX.  (0 if option 12 not used.)
    integer  :: iop13 !! If not zero, reverse communication will be used for
                      !! getting the values of derivatives.  Associated with option 13.
    integer  :: iop14 !! If not zero, reverse communication will be used in
                      !! place of calls to the output routine.  Associated with option 14.
    integer  :: iop15 !! If not zero, a return will be made to the user after
                      !! the initialization.  Associated with option 15.  This might be used
                      !! to overlay DIVA, some of the user's code, and perhaps DIVAop.
    integer  :: iop16 !! Points to location in KORD() where information for
                      !! specifying the error tolerance is specified.  See option 16.
    integer  :: iop17 !! Used in initialization for option 17, afterwards this
                      !! cell is used by KEXIT which is equivalenced to IOP17.
    integer  :: iop18 !! Points to location in KORD() where information for
                      !! specifying a grouping of equations for derivative evaluation is
                      !! stored.  See option 18.
    integer  :: iop19 !! Points to location in KORD() where information for
                      !! specifying a grouping of equations for integration order control
                      !! is stored.  See option 19.
    integer  :: iop20 !! Used for option 20, gives first location in F where
                      !! estimated errors are to be stored.  Expected to be useful in a
                      !! program for solving boundary value problems using multiple shooting.
    integer  :: iop21 !! Was used for stiff equations option (never completely
                      !! coded).  The optional code still uses this (don't activate it!).
                      !! Now used to flag the location if F where the user has stored the
                      !! tolerance to use in finding G-Stops.
    integer  :: iop21s !! Was used for stiff equations see above.
    integer  :: iop22 !! Set aside for possible option for stiff equations.
    integer  :: iop3 !! Value set by option 3:
                     !!
                     !!  *  0 Interpolate to final point. (The default)
                     !!  *  1 Integrate to final point.
                     !!  * -1 Extrapolate to final point.
    integer  :: iop4 !! Value set by option 4.  The output routine is called
                     !! with KORD(1) = 4, every IOP4 steps.  (Default value for IOP4 is a
                     !! very large number.
    integer  :: iop5 !! Value provided by option 5, used to specify extra
                     !! output points.
    integer  :: iop6 !! Value provided by option 6.  If nonzero, the output
                     !! routine is called at the end of every step.  If > 0, there are
                     !! IOP6 interpolating G-Stops.
    integer  :: iop7 !! Value provided by option 7.  If > 0, there are K7
                     !! extrapolating G-Stops.
    integer  :: iop8 !! Value provided by option 8.  If nonzero, the output
                     !! routine is called with KORD(1)=8 whenever the step size is changed.
    integer  :: iop9 !! Value provided by option 9.  Used to specify that the
                     !! user wishes to save the solution.
    integer  :: ioptc(23) !! In DIVAOP equivalenced so that IOPTC(3) is equivalent
                          !! to IOP3.
    integer  :: itolep !! Used for temporary storage, and for the index of a
                       !! tolerance relative to the start of tolerances.
    integer  :: ivc2(65) !! Array used for output of variables ICF to NY in
                         !! common block DIVAMC.
    integer  :: iy !! Used for the current index to the Y() array.  (Local
                   !! variable in DIVAIN used in computing IYI.)  Equivalenced to
                   !! IZFLAG in DIVAG.
    integer  :: izflag
    integer  :: kemax !! Index associated with equation giving the largest
                      !! value for (estimated error) / (requested error).
    integer  :: kexit !! Equivalenced to IOP17 which is not used after
                      !! initialization.  Defines actions when KORD2I = -7.
                      !! (Referenced in (DIVAA,DA,G).).
                      !!
                      !!  * 1  Take the step over with reduced H.
                      !!  * 2  Take the step over.
                      !!  * 3  Do the end of step call to OUTPUT.
                      !!  * 4  Reset TMARK, then do same as for KEXIT = 2.
                      !!  * 5  Reset TMARK, then do same as for KEXIT = 3.
                      !!  * 6  Give the fatal error diagnostic.
    integer  :: kis !! Used to check if it is time to dump the solution.
                    !! The check involves incrementing KIS at the end of the step, and
                    !! dumping the solution if KIS is 0.
                    !!
                    !!  * -1  Set in DIVAcr when it is time to dump solution
                    !!  *  0  When starting
                    !!  *  2  After being dumped.
                    !!
                    !! This is set to 1000 just after a user specified discontinuity, and
                    !! counted up from that point.
    integer  :: kmark !! Identifies the type of output associated with the next
                      !! output point specified by TSPECS.
    integer  :: kord1i !! Helps in defining the state of the integrator.
                       !! Frequently has the same value as KORD(1).  Meaning depends on the
                       !! value of KORD(2), or the value about to be assigned to KORD(2).
                       !!
                       !!  * <  0  Happens when preparing to give output with extrapolation.
                       !!  * =  0  Happens when checking F at points for noise test.
                       !!  * =  1  (KORD(2)=-1)  End of integration has been reached.
                       !!  * =  1  (KORD(2)= 0)  Computing first predicted derivative.
                       !!  * =  1  (KORD(2)= 1)  Output for initial point.
                       !!  * =  2  (KORD(2)=-1)  Giving diagnostic for noise limiting precision.
                       !!  * =  2  (KORD(2)= 0)  Computing corrected derivative.
                       !!  * =  2  (KORD(2)= 1)  Output for TSPECS(3).
                       !!  * =  3  (KORD(2)=-1)  Diagnostic for step size reduction too fast.
                       !!  * =  3  (KORD(2)= 0)  Computing variational derivative.
                       !!  * =  3  (KORD(2)= 1)  Output for TSPECS(4).
                       !!  * =  4  (KORD(2)=-1)  Error, discontinuity.
                       !!  * =  4  (KORD(2)= 1)  Output for certain number of steps.
                       !!  * =  5  (KORD(2)= 0)  Get initial derivatives for stiff equations.
                       !!  * =  5  (KORD(2)= 1)  Extra output from TSPECS.
                       !!  * =  6  (KORD(2)= 1)  End of step output.
                       !!  * =  7  (KORD(2)= 0)  Evaluate G before extrapolated output point.
                       !!  * =  7  (KORD(2)= 1)  Evaluate G before extrapolated output point.
                       !!    (Also used when checking for other G's after
                       !!    finding one.)
                       !!  * =  8  (KORD(2)= 1)  Tell user step size has changed.
                       !!  * =  9  (KORD(2)= 1)  Request for user to save solution.
                       !!  * = 11  (KORD(2)=-1)  Error, step size too small at end of start.
                       !!  * = 12  (KORD(2)=-1)  Error, step size is too small.
                       !!  * = 13  (KORD(2)=-1)  Error, output points specified badly.
                       !!  * = 21  (KORD(2)=-1)  H too small to give reasonable change when added
                       !!    to T.
                       !!  * = 22  (KORD(2)=-1)  Error, bad tolerance.
                       !!  * = 23  (KORD(2)=-1)  Set after message for a fatal error.
                       !!  * = 24  Set on error message in DIVA, along with KORD2I = -4.
                       !!
                       !! Also used as an index into MLOC in DIVAA when an error is being
                       !! processsed, see MLOC below.
    integer  :: kord2i !! Helps in defining the state of the integrator.
                       !! Frequently has the same value as KORD(2).
                       !!
                       !!  * -3  Set in DIVAg, to get a derivative evaluation.
                       !!  * -2  Set in DIVAg, to get another entry to OUTPUT.
                       !!  * -1  Return to calling program, done, interrupt, or got an error.
                       !!  *  1  Calling OUTPUT or returning to user for OUTPUT type action.
                       !!  *  0  Calling DERIVS or returning to user for DERIVS type action.
                       !!  * -4  Error message in DIVA and in DIVAop, along with KORD1I = 24.
                       !!  * -5  Starting
                       !!  * -6  Starting, getting the initial derivative value or derivatives
                       !!    for the noise test.
                       !!  * -7  Done some extrapolation, KEXIT defines the action to take.
                       !!     Set in DIVAg to activate KEXIT action in DIVA.
                       !!  * -8  Set when user has requested adjustment of the difference
                       !!    tables for a discontinutiy.
    integer  :: kpred !! Value assigned to KORD1I when getting a predicted
                      !! derivative.  (1 used now, 5 planned for use with stiff equations.)
    integer  :: kqdcon !! Number of coefficients computed with constant step
                       !! size for stiff equations.
    integer  :: kqicon !! Number of coefficients computed with constant step
                       !! size for nonstiff equations.
    integer  :: kqmaxs !! Maximum integration order for equations that have
                       !! some limit on the error that can be committed.
    integer  :: kqmxds !! Used to save KQMAXD in case step is repeated and the
                       !! solution must be dumped.
    integer  :: kqmxil !! Value of KQMAXI the last time integration coefficients
                       !! were computed.
    integer  :: kqmxip !! = KQMAXI + MAXINT, for computing integration coeffs.
    integer  :: kqmxis !! Used to save KQMAXI in case step is repeated and the
                       !! solution must be dumped.
    integer  :: ksc !! Number of steps that have been taken with a constant
                    !! step size.
    integer  :: ksout !! When KSTEP reaches this value, the output routine is
                      !! called with KORD(1) = 4.  The default value is a very large number.
    integer  :: ksstrt !! Set when ending one derivative per step to KSTEP + 2.
                       !! Checked later in DIVAHC to decide whether to set the step changing
                       !! factors to their nominal values.
    integer  :: kstep !! Number of steps taken since the start of integration.
    integer  :: lex !! Indicates how to get values at next output point:
                    !!
                    !!  * -1  Extrapolate
                    !!  *  0  Interpolate (The usual case.)
                    !!  *  1  Integrate to the output point, integration is not continued.
    integer  :: linc !! Used to indicate state of step size selection.
                     !!
                     !!  * -10 After computed derivatives at base time, after computing other
                     !!    extra derivatives for the noise test.
                     !!  * -9  After computed second extra derivative for noise test.
                     !!  * -8  After computed first extra derivative for noise test.
                     !!  * -7  Dumping the solution and then doing a user initiated restart,
                     !!    or getting ready to compute extra derivatives for the noise
                     !!    test.
                     !!  * -6  Dumping the solution before a restart.
                     !!  * -5  Set on the first step, and also set when dumping the solution
                     !!    after a discontinuity.
                     !!  * -4  Repeat step with no change in the step size.
                     !!  * -3  Set when the error tolerance is set improperly.
                     !!  * -2  User has complete control of selecting the step size.
                     !!  * -1  Step is being repeated.
                     !!  *  0  Step size is not to be increased on this step.
                     !!  * k>0 Step size can be increased by HINCC**k.
    integer  :: lincd !! Value of smallest k for which HINCC**k >= 2.
                      !! (=-2 if user is specifying all step size changes.)
    integer  :: lincq !! Value of smallest k for which HINCC**k >= 4.
    integer  :: lsc !! Indicates if starting or if noise may be present.
                    !!
                    !!  * k<0 -k steps have been taken for which noise appears to be limiting the precision.
                    !!  * 0   Usual case
                    !!  * 1   Doing 1 derivative per step after initial part of start.
                    !!  * 2   Used as flag that it is time to set LSC=0.
                    !!  * 3   Third step, hold the order constant.
                    !!  * 4   Second step, increase orders from 2 to 3.
                    !!  * 5   First step, third time through the first step (if required).
                    !!  * 6   First step, second time through.
                    !!  * 7   First step, first time through.
                    !!  * 8   Set on initialization.
    integer  :: maxkqd !! Largest integration order allowed for stiff equations.
    integer  :: maxkqi !! Largest integ. order allowed for nonstiff equations.
    integer  :: method !! Defines kind of methods being used.
                       !!  * -1  Only stiff equations are being integrated.
                       !!  *  0  Only nonstiff equations are being integrated.
                       !!  *  1  Both kinds of methods are required.
    integer  :: ne !! Number of equations in the first group.  (=NTE if
                   !! option 18 is not used.)
    integer  :: neptol !! Used for temporary storage and to save the value of
                       !! ITOLEP for error messages.
    integer  :: ng !! Used in DIVAg for the number of g's in the current context.
    integer  :: ngstop(2) !! Dimension 2 array equivalenced to IOP6, and IOP7.  To
                          !! get the number of interpolating and extrapolating G-Stops.
    integer  :: ngtot !! NGTOT(1) gives the number of interpolating G-Stops,
                      !! and NGTOT(2) gives the number of extrapolating G-Stops.
    integer  :: noiseq !! max(2, order of equation for which (error estimate)/
                       !! (error requested) is a maximum).
    integer  :: noutko !! If nonzero, gives the index in KORD where information
                       !! on what equations are to be included in the diagnostic output is
                       !! given.   See option 10.
    integer  :: ntolf !! First location in F() where tolerance specifying
                      !! accuracy desired is stored.
    integer  :: ny !! Total order of the system.
    real(wp) :: alpha(kdim) !! Array with I-th entry = (current step size) / XI(I).
                            !! Used in computing integration coefficients.
    real(wp) :: beta(kdim+1) !! Array with I-th entry = product (K=1,I-1) of
                             !! (current (XI(K)) / XI(K) from previous step), BETA(1)=1. Used in
                             !! updating the difference tables.
    real(wp) :: d(maxstf+maxord,maxord) !! Array to be used later to store coefficients for
                                        !! integrating stiff equations.
                                        !! derivatives.  Not used if option 13 is set.
    real(wp) :: dnoise !! Used in determining if noise is limiting the
                       !! precision.  It is usually |highest difference used in correcting|
                       !! of the equation with the largest error estimate.
    real(wp) :: ds(maxstf+maxord, maxord) !! Array to be used later to store coefficients for
                                          !! estimating errors when integrating stiff equations.
    real(wp) :: dvc1(7)
    real(wp) :: dvc2(7) !! Array used for output of variables HC to TOUT in
                        !! common block DIVAMC.
    real(wp) :: eave !! This is a weighted average of past values of EIMAX.
                     !! It is adjusted to account for expected changes due to step changes.
    real(wp) :: eimax !! Estimate of (error estimate / error requested) if the
                      !! step size should be increased.
    real(wp) :: eimin !! An error estimate is small enough to allow a step
                      !! increase if the estimate of ((error with the step size increased) /
                      !! (error requested)) is less than EIMIN.
    real(wp) :: emax !! Largest value computed for (error estimate) / (error
                     !! requested).
    real(wp) :: erep !! If EMAX > EREP, a step is repeated.  Ordinarily
                     !! this has the value .3.  This is set < 0 if the error tolerance is
                     !! specified improperly, and is set to a large value if the user
                     !! requests complete control over the step size.  EREP is also set
                     !! < 0 after a user specified discontinuity.
    real(wp) :: fdat(11) !! Used to store data for error messages.  (Local array in
                         !! DIVAIN.)
    real(wp) :: g(kdim,maxord) !! Integration coefficients used for predicting solution.
                               !! G(I, J) gives the I-th coefficient for integrating a J-th order
                               !! differential equation.  G(1, 1) is equal to the step size.
    real(wp) :: gs(kdim) !! Integration coefficients used in estimating errors.
    real(wp) :: hc !! Ratio of (new step size) / (old step size)
    real(wp) :: hdec !! Default value to use for HC when reducing the step
                     !! size.  (Values closer to 1 may be used some of the time.)
    real(wp) :: hh
    real(wp) :: hinc !! Default value to use for HC when increasing the step
                     !! size.  (Values closer to 1 may be used some of the time.)
    real(wp) :: hincc !! Actual value used for default value of HC when
                      !! increasing the step size.  Set to HINC after start is considered
                      !! complete.  During the start HINCC is set to 1.125.
    real(wp) :: hmax !! Largest value allowed for abs(step size).  Default
                     !! value is a very large number.
    real(wp) :: hmaxp9 !! .9 * HMAX.
    real(wp) :: hmin !! Smallest value allowed for abs(step size).  Default
                     !! value is 0.
    real(wp) :: rbq(kdim) !! Array containing data for the preliminary noise test.
    real(wp) :: robnd !! Used to influence the selection of integration order.
                      !! The larger ROBND, the harder it is to increase the order and the
                      !! easier it is to decrease it.
    real(wp) :: rvc2(8) !! Array used for output of variables DNOISE to SNOISE in
                        !! common block DIVAMC.  These are variables that don't require a great
                        !! deal of precision.
    real(wp) :: sigma(kdim) !! The k-th entry of this array contains a factor that
                            !! gives the amount the k-th difference is expected to increase if the
                            !! step size in increased.  These numbers get bigger it there is a past
                            !! history of increasing the step size.
    real(wp) :: snoise !! Value used in comparison with RBQ() on equation with
                       !! largest value for (error estimate) / (error requested).
    real(wp) :: tg(2) !! TG(1) gives the last value of TSPECS(1) for which an
                      !! interpolatory G-Stop has been computed.  TG(2) is defined similarly
                      !! for extrapolatory G-Stops.
    real(wp) :: tgstop(2) !! TGSTOP(1) gives the value of TSPECS(1) where the last
                          !! 0 for an interpolatory G-Stop was found.  TGSTOP(2) is defined
                          !! similarly for extrapolatory G-Stops.
    real(wp) :: tmark !! Location of the next output point.
    real(wp) :: tmarka(2) !! Array of length 2 equivalenced to TMARK (and TMARKX).
    real(wp) :: tmarkx !! Location of the next output point to be found using
                       !! integration or extrapolation.  This variable must follow immediately
                       !! after TMARK in the common block.
    real(wp) :: tolg !! Tolerance to pass to dzero when locating G-Stops.
    real(wp) :: tout !! Location of next output point defined by value of
                     !! TSPECS(3).  Such output is given with KORD(1) = 2.
    real(wp) :: v(kdim+maxord) !! Array used in computing integration coefficients.

    common /divamc/ tg,tgstop,tmark,tmarkx,tout,tolg,hc,hdec,hinc,    &
    &   hincc,hmax,hmaxp9,hmin,alpha,beta,d,g,v,ds,gs,sigma,rbq,dnoise,&
    &   eave,eimax,eimin,emax,erep,robnd,snoise,fdat,icf,ics,igflg,    &
    &   igtype,igstop,ilgrep,ings,iop3,iop4,iop5,iop6,iop7,iop8,iop9,  &
    &   iop10,iop11,iop12,iop13,iop14,iop15,iop16,iop17,iop18,iop19,   &
    &   iop20,iop21,iop22,iop21s,itolep,iy,kemax,kis,kmark,kord1i,     &
    &   kord2i,kpred,kqdcon,kqicon,kqmaxs,kqmxds,kqmxil,kqmxip,kqmxis, &
    &   ksc,ksout,ksstrt,kstep,lex,linc,lincd,lincq,lsc,maxkqd,maxkqi, &
    &   method,ne,neptol,ng,ngtot,noiseq,noutko,ntolf,ny,idat

    equivalence (tmarka(1), tmark)
    equivalence (kexit,     iop17)
    equivalence (g(1, 1),   hh)
    equivalence (ngstop(1), iop6)
    equivalence (izflag,    iy)
    equivalence (igflgs,    itolep)
    equivalence (ivc2(1),   icf)
    equivalence (ioptc(3),  iop3)

    end module divamc_module
!*************************************************************************
