      program DRDIVDB
!>> 2009-11-03 DRDIVDB  Krogh -- Used option 11
!>> 2001-07-16 DRDIVDB  Krogh -- Declared type for paramter TOL.
!>> 2001-05-25 DRDIVDB  Krogh -- Added comma to format.
!>> 1996-06-14 DRDIVDB  Krogh  Small change in output format
!>> 1994-09-12 DRDIVDB  Krogh  Fixed for CHGTYP.
!>> 1993-07-18 DRDIVDB  Krogh  Fixed to get same output in S.P. and D.P.
!>> 1993-05-05 DRDIVDB  Krogh  Adjusted to simplify conversion to C.
!>> 1989-05-04 Fred T. Krogh
!
!--D replaces "?": DR?IVDB, ?IVA, ?IVADB, ?IVAF, ?IVAO, ?MESS
!
! Sample driver for DIVA --  Set up to solve two second order equations.
! Tests debug output
      use diva_module

      integer INEQ, IFDIM, IYDIM
      parameter (INEQ=2, IFDIM=16*INEQ+1, IYDIM=4*INEQ)
      integer          NEQ, KORD(6), IOPT(7)
!++ Code for SHORT_LINE is inactive
!      integer          MACT(5)
!++ END
      double precision TSPECS(4), T, H, DELT, TFINAL
      double precision F(IFDIM), Y(IYDIM)
      external DIVAO, DIVAF
      equivalence (TSPECS(1), T), (TSPECS(2), H), (TSPECS(3), DELT),    &
     &   (TSPECS(4), TFINAL)
      double precision TOL
      integer NDIG
!++S Default NDIG = 4
!++  Default NDIG = 10
!++ Substitute for NDIG below
      parameter (NDIG = 10)
      parameter (TOL = 10.D0 **(-NDIG))
!
! Parameters for setting up message processor to specify the number of
! digits, and the line length.
      integer MEDDIG, MEMLIN, MERET
      parameter (MEDDIG=12, MEMLIN=13, MERET=51)
!
      data   NEQ,    T,    H,                   DELT, TFINAL /          &
     &         2, 0.D0, 1.D0, 6.283185307179586477D0,   2.D1 /,         &
     &       Y(1), Y(2), Y(3), Y(4) /                                   &
     &       1.D0, 0.D0, 0.D0, 1.D0 /
!
! Set option for error control, local absolute error < 1.D-10.
      data   IOPT(1), IOPT(2), IOPT(3) /                                &
     &            16,       6,   3 /,                                   &
     &                KORD(6), F(3)/                                    &
     &                      2, TOL/
! Group the system to be treated as a single unit, set tolerance value
! Set option for second order equations
      data   IOPT(4), IOPT(5) /                                         &
     &            17,       2 /
! Set option to initialize some space to 0, and end of options.
      data IOPT(6), IOPT(7) / 11, 0 /

!
!
! Do the integration
!
      KORD(1) = 0
  100 continue
         call DIVA(TSPECS,Y,F,KORD,NEQ,DIVAF,DIVAO,4,IYDIM,IFDIM,6,IOPT)
         if (KORD(1) /= 1) go to 100
!++ Code for SHORT_LINE is inactive
!      MACT(1) = MEDDIG
!      MACT(2) = 7
!      MACT(3) = MEMLIN
!      MACT(4) = 79
!      MACT(5) = MERET
!      call DMESS(MACT, ' ', MACT, F)
!++ END
      call DIVADB(44, TSPECS, Y, F, KORD, ' Test of DIVADB')
      stop
      END

      subroutine DIVAF(T, Y, F, KORD)
!
! Sample derivative subroutine for use with DIVA
! This evaluates derivatives for a simple two body problem.
!
      integer          KORD
      double precision T, F(2), Y(4)
      double precision TP
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
!
! Sample output subroutine for use with DIVA.
! This subroutine gives output for a simple two body problem.
!
      integer          KORD
      double precision TSPECS(4), Y(4), F(2)
 1000 format (12X,                                                      &
     &   'RESULTS FOR A SIMPLE 2-BODY PROBLEM'//                        &
     &   8X, 'T', 13X, 'U/V', 11X, 'UP/VP', 9X, 'UPP/VPP')
 1001 format (1P,SP,4E15.6 / 15X, 3E15.6/' ')
!
! Do the output
!
      if (KORD == 1) then
        write (*, 1000)
      end if
      write (*, 1001) TSPECS(1), Y(1), Y(2), F(1), Y(3), Y(4), F(2)
!
      return
      END
