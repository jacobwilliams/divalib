      program DRDIVA
!>> 2010-06-09 DRDIVA  Krogh Used parameters for all dimenssions.
!>> 2001-05-25 DRDIVA  Krogh Minor change for making .f90 version.
!>> 1996-06-14 DRDIVA  Krogh  Small change in output format
!>> 1994-11-02 DRDIVA  Krogh  Changes to use M77CON
!>> 1994-07-18 DRDIVA  Krogh   Last change.
!--D replaces "?": DR?IVA, ?IVA, ?IVAF, ?IVAO
! Sample driver for DIVA --  Set up to solve two second order equations.
!
      use diva_module
      integer INEQ, IFDIM, IKDIM, ITDIM, IYDIM
      parameter (INEQ=2,IFDIM=16*INEQ+1,IKDIM=6,ITDIM=4,IYDIM=4*INEQ)
      integer          NEQ, KORD(IKDIM), IOPT(6)
!--D Next line special: P=>D, X=>Q
      double precision TSPECS(ITDIM), Y(IYDIM), T, H, DELT, TFINAL
      integer NDIG
      double precision F(IFDIM), TOL
      equivalence (TSPECS(1), T), (TSPECS(2), H), (TSPECS(3), DELT),    &
     &  (TSPECS(4), TFINAL)
      external DIVAO, DIVAF
!++SP Default NDIG = 4
!++  Default NDIG = 10
!++ Substitute for NDIG below
      parameter (NDIG = 10)
      parameter (TOL = 10.D0 **(-NDIG))
!
      data   NEQ,    T,    H,                   DELT, TFINAL /          &
     &         2, 0.D0, 1.D0, 6.283185307179586477D0,   2.D1 /,         &
     &       Y(1), Y(2), Y(3), Y(4) /                                   &
     &       1.D0, 0.D0, 0.D0, 1.D0 /
!
! Set option for error control, local absolute error < TOL.
      data   IOPT(1), IOPT(2), IOPT(3) /                                &
     &            16,       6,   3 /,                                   &
     &                KORD(6), F(3)/                                    &
     &                      2,  TOL/
! Group the system to be treated as a single unit, set tolerance value
! Set option for second order equations
      data   IOPT(4), IOPT(5) /                                         &
     &            17,       2 /
!  Flag end of options
      data   IOPT(6) / 0/
!
! Do the integration
!
      KORD(1) = 0
  100 continue
      call DIVA(TSPECS, Y, F, KORD, NEQ, DIVAF, DIVAO,                  &
     &  ITDIM, IYDIM, IFDIM, IKDIM, IOPT)
      if (KORD(1) .NE. 1) go to 100
      stop
      END

      subroutine DIVAF(T, Y, F, KORD)
! Sample derivative subroutine for use with DIVA
! This evaluates derivatives for a simple two body problem.
!
      integer          KORD
!--D Next line special: P=>D, X=>Q
      double precision T, Y(4), TP
      double precision F(2)
!
! Evaluate the derivatives
!
      TP = Y(1)*Y(1) + Y(3)*Y(3)
      TP = 1.D0 / (TP * SQRT(TP))
      F(1) = -Y(1) * TP
      F(2) = -Y(3) * TP
      return
      END

      subroutine DIVAO(TSPECS, Y, F, KORD)
! Sample output subroutine for use with DIVA.
! This subroutine gives output for a simple two body problem.
!
      integer          KORD
!--D Next line special: P=>D, X=>Q
      double precision TSPECS(4), Y(4)
      double precision F(2)
 1000 format (12X,                                                      &
     &   'RESULTS FOR A SIMPLE 2-BODY PROBLEM'//                        &
     &   8X, 'T', 13X, 'U/V', 11X, 'UP/VP', 9X, 'UPP/VPP')
 1001 format (1P,SP,4E15.6 / 15X, 3E15.6/' ')
!
! Do the output
!
      if (KORD .EQ. 1) then
        write (*, 1000)
      end if
      write (*, 1001) TSPECS(1), Y(1), Y(2), F(1), Y(3), Y(4), F(2)
!
      return
      END
