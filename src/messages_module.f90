!*************************************************************************
!>
!  Processes Messages

    module messages_module

    use diva_constants

    implicit none

    public

    contains
!*************************************************************************

!*************************************************************************
!>
! Processes Messages -- Actions are controlled by MACT().  See
! comment is subroutine MESS.  This program is for the extra
! argument of type real.
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

      integer          :: MACT(*)
      integer          :: IDAT(*)
      double precision :: FDAT(*)
      character        :: TEXT(*)*(*)

      character        FMTSP*29
      integer          ICOL, ID, J, K, KSMA, KBIG, KEXE, LDFDEF, NEG
      double precision FBIG, FOUT, FSMA
      save LDFDEF, FMTSP
      save /CMESSI/, /CMESSC/

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

! ************************** Data from common block ********************
!
! For comments on these variables, see the listing for MESS.
!
      integer, parameter :: LENBUF = 250
      integer, parameter :: MEVBAS = 10
      integer, parameter :: MEVLAS = 33

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
!++ END
        IWF = IWF + KDI + 1
        go to 10
      end if
!
      LBUF = LBUF + IWF
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
      write(BUF(LSTRT:LBUF),FMTF)(FDAT(K),K=MPT,MPT+INC*(KLINE-1),INC)

!      print '(/A/)', BUF(1:LBUF)

      MPT = MPT + KLINE * INC
      go to 10
!                                    Floating point matrix output
  300 continue
      write (BUF(LSTRT:LBUF), FMTF) (FDAT(K), K = MPT, LASTI, NDIM)
      go to 10
!                                    Table output
  400 continue
      write (BUF(LSTRT:LBUF), FMTF) (FDAT(K), K = MPT, MPT+KLINE-1)
      go to 10

!                                   Sparse vector output
  500 continue
      write (BUF(1:LBUF), FMTSP) (IDAT(K),FDAT(K),K=MPT,MPT+KLINE-1)
      MPT = MPT + KLINE
      go to 10

    end subroutine dmess
!*************************************************************************


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

    subroutine MESS(MACT, TEXT, IDAT)

        integer, parameter :: LNMSG = 128
        integer, parameter :: LNERR = 79
  !
  ! ************** Parameters Defining Actions (See Above) ***************
  !

  ! Parameters for changing the environment.
      integer, parameter :: MESUNI=10
      integer, parameter :: MEHEAD=11
      integer, parameter :: MEDDIG=12
      integer, parameter :: MEMLIN=13
      integer, parameter :: MEELIN=14
      integer, parameter :: MEMUNI=15
      integer, parameter :: MEEUNI=16
      integer, parameter :: MESCRN=17
      integer, parameter :: MEDIAG=18
      integer, parameter :: MEMAXE=19
      integer, parameter :: MESTOP=20
      integer, parameter :: MEPRNT=21
      integer, parameter :: METDIG=22
      integer, parameter :: MENTXT=23
      integer, parameter :: MEIDAT=24
      integer, parameter :: MEFDAT=25
      integer, parameter :: MEMDAT=26
      integer, parameter :: MEMDA1=27
      integer, parameter :: MEMDA2=28
      integer, parameter :: MEMDA3=29
      integer, parameter :: MEMDA4=30
      integer, parameter :: MEMDA5=31
      integer, parameter :: METABS=32
      integer, parameter :: MEERRS=33

  ! Parameters for actions.
      integer, parameter :: MECONT=50
      integer, parameter ::  MERET=51
      integer, parameter :: MEEMES=52
      integer, parameter :: METEXT=53
      integer, parameter :: MEFSPV=54
      integer, parameter :: METABL=55
      integer, parameter :: MERES3=56
      integer, parameter :: MEIVEC=57
      integer, parameter :: MEIMAT=58
      integer, parameter :: MEJVEC=59
      integer, parameter :: MEJMAT=60
      integer, parameter :: MEFVEC=61
      integer, parameter :: MEFMAT=62
      integer, parameter :: MEGVEC=63
      integer, parameter :: MEGMAT=64
      integer, parameter :: MEIVCI=65
      integer, parameter :: MEJVCI=66
      integer, parameter :: MEFVCI=67
      integer, parameter :: MEGVCI=68

  ! Parameter derived from those above.
      integer, parameter :: MEMAXI=68
      integer, parameter :: MEGBAS=49
      integer, parameter :: MEVBAS=10
      integer, parameter :: MEVLAS=33
  !
  ! ************************** Variable Declarations *********************
  !
        !external MESSGS
        integer    MACT(*), IDAT(*)
        character  TEXT(*)*(*)
  !
        integer    I, ICOL, INCM(MECONT:MEIMAT), INERR, IOUT, IROW, IROW1,&
       &    ITEXTR, ITXTSV, J, JJ, K, K1, K2, KDILAB, KK, KNT, KOLWID, KP,&
       &    KS, LASKNT, LBUF1, LBUF2, LINSTR, M,                          &
       &    MBNDHI(MEVBAS:MEVLAS), MBNDLO(MEVBAS:MEVLAS), MTEXT(2),       &
       &    MTEXTC, MTEXTR, NLINE, NROCO(2), NSKIP, NTEXTR, NTXTSV
        logical   GETW, FIRST
        character ERMSG*63, ERMSG1*27
        character C

        character, parameter :: SC = '$'

        save  FIRST, I, ICOL, INERR, IROW, IROW1, ITXTSV, KDILAB, KNT,    &
       &   LASKNT, M, MTEXT, NLINE, NSKIP, NTEXTR, NTXTSV
        save /CMESSI/, /CMESSC/
        equivalence (MTEXT(1), MTEXTR), (MTEXT(2), MTEXTC)
  !
  ! ************************** Data from common block ********************
  !
        integer, parameter :: LENBUF = 250

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
              K = abs(IVAR(M))
           end if
           OUNIT = MUNIT
           go to 120
        end if
        if (M == -MESUNI) then
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
           BUF(LBUF+1:LBUF+len(ERMSG1)) = ERMSG1
           LBUF = LBUF + len(ERMSG1)
        else
           ERMSG(30:30) = char(KS + ICHAR0)
           ERMSG(47:47) = char(KP + ICHAR0)
           BUF(LBUF+1:LBUF+len(ERMSG)) = ERMSG
           LBUF = LBUF + len(ERMSG)
        end if
        LSTRT = LBUF + 1
        call MESSFI
        LBUF = LBUF + KDI
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
        K = index(TEXT(ITEXT)(NTEXT:K2), SC)
        if (K == 0) then
  ! Want to take all that we can.
           LBUF = LSTRT + K2 - NTEXT
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
           BUF(LSTRT:LBUF1-1) = BUF(LBUF1:K-1)
           LSTRT = LBUF1
           LBUF1 = K
           go to 530
        end if
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
        go to 603
    602 continue
  !++ CODE for ~.C. is active
        FMTG='(1P,99'//TEXT(ITEXT)(NTEXT:NTEXT)
  !++ CODE for .C. is inactive
  !      FMTG(5:5) = TEXT(ITEXT)(NTEXT:NTEXT)
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
        write (BUF(LSTRT:LBUF), FMTJ) IOUT
           go to 850
        end if
  !
  !                Get format for integer output.
        call MESSFI
        LBUF = LBUF + KDI
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
  !        IWF = LENOUT
  !        lfprec = JJ
  !++ END
           if (.not. XARGOK) go to 180
           MPT = NFDAT
           NFDAT = NFDAT + KLINE
           return
  !                           Text output
    948    K1 = NTEXT + LBUF - LSTRT
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
              end if
           end if
           write(SUNIT) BUF(5:MDAT(NROW))
        end if
        if (LBUF > MDAT(NROW)) then
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
        endfile SUNIT
        KK = 1
    994 KK = KK + 1
        if (KK > NROW) go to 999
        rewind(SUNIT)
        IROW = -1
        K = KK
    995 LBUF = 5
        IROW = IROW + 1
        if (IROW /= 0) then
           write(BUF(1:4), '(I4)') mod(IROW, 10000)
        else
           BUF(1:4) = ' '
        end if
        do 996 J = 2, K
           if (J == K) LBUF = MDAT(KK)
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
        FMTR = FMTI
        KDILAB = KDI+1
        LINSTR = 2*KDILAB+2
        if (XARG) then
           if (.not. GOTFMT) go to 1150
           IWF = IWG
           FMTF = FMTG
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
        BUF(LBUF+1:K) = ' '
        LBUF = K
        go to 1110
   1085 call MESSPR
   1090 continue
  !++ CODE for ~.C. is active
        BUF = ' '
        write (BUF(1:KDILAB), FMTR) KNT+1
        BUF(KDILAB:KDILAB) = '-'
        KLINE = min(NLINE, LASKNT - KNT)
        KNT = KNT + KLINE
        write (BUF(KDILAB+1:2*KDILAB), FMTR) KNT
        BUF(2*KDILAB:LINSTR-1) = ':'
        LBUF = LINSTR
   1110 LSTRT = LBUF
        LBUF = LBUF + LENOUT * KLINE - 1
        if (XARG) return
  !                                    Integer output
  !++ CODE for ~.C. is active
        write (BUF(LSTRT:LBUF), FMTI) (IDAT(K), K = MPT,                  &
       &    MPT+INC*(KLINE-1), INC)
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
           BUF(LBUF+1:LINSTR) = ' '
           go to 1520
        end if
   1515 continue
        write (BUF(LBUF+1:LINSTR), FMTR) IROW
   1520 LSTRT = LINSTR + 1
        LBUF = LINSTR + LENOUT*KLINE
        LASTI = MPT + NDIM * KLINE - 1
        if (XARG) return
  !                                    Integer output
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

      subroutine MESSFD(IDAT)

        integer :: J !! Temporary index.
        integer :: K !! Temporary index.
        integer :: IDAT(*) !! Input array to [[MESS]]
        integer :: IMAX !! Used when computing largest integer in array.
        integer :: IMIN !! Used when computing smallest integer in array.

  ! For comments on other variables, see the listing for MESS.
        integer, parameter :: LENBUF = 250
        integer, parameter :: MEVBAS = 10
        integer, parameter :: MEVLAS = 33

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
        integer, parameter :: LENBUF = 250
        integer, parameter :: MEVBAS = 10
        integer, parameter :: MEVLAS = 33

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
        end if
        return
      end subroutine messfi
  !*************************************************************************

  !*************************************************************************
  !>
  !  Get a scratch unit assigned.
  !
  !### History
  !  * Jacob Williams : 5/29/2021 : rewrote this routine

    function messgs() result(iunit)

    integer :: iunit

    integer :: istat !! iostat code

    open (newunit=iunit, &
          status='SCRATCH', &
          access='SEQUENTIAL', &
          form='UNFORMATTED', &
          iostat=istat)

    if (istat /= 0) error stop 'Could not assign scratch unit in MESS.'

      end function messgs
  !*************************************************************************

  !*************************************************************************
  !>
  ! Processing of multiple headings.

      subroutine MESSMH(TEXT)

        character*(*) :: TEXT(*) !! Original input character vector.

        integer :: J !! Used as a temporary index.
        integer :: K !! Used as a temporary index.
        integer :: KB !! Number of initial blanks
        integer :: KK !! Used as a temporary index.
        integer :: L !! Used as a temporary index.
        integer :: LFLGDB !! 2 in usual case, 3 if see a "$ ".
        integer :: LSTRDB !! Value of LSTRT when see a "$ ".
        integer :: LTXTDB !! Value of LTEXT when see a "$ ".
        character :: C

        character, parameter :: SC = '$'

  ! For comments on other variables, see the listing for MESS.
        integer   KOLWID, LINSTR

        integer, parameter :: LENBUF = 250
        integer, parameter :: MEVBAS = 10
        integer, parameter :: MEVLAS = 33

  ! KT    Used for logic in output of headings.  Set <0 on exit if there
  !       is an input error.
  !     KT = 1 Output table headings.  (Set to -1 on fatal error.)
  !     KT = 2 Get row/column widths for matrix output.  (Set to -2 if
  !            error results in no headings.)
  !     KT = 3 Output column headings.  (Set to -1 on fatal error.)

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
           if (KB > 0) BUF(LSTRT-KB:LSTRT-1) = ' '
  !                                  Move characters
           if (NTEXT == 0) then
              BUF(LSTRT:LSTRT+3) = 'Col '
              C = '#'
              LSTRT = LSTRT+4
           else
              K = NTEXT - LTEXT - LFLGDB
              if (K <= 0) then
                 KK = max(0, 3-NTEXT)
                 BUF(LSTRT:LSTRT-K-KK)=TEXT(ITEXT-1)(LENTXT+K:LENTXT-KK)
                 LSTRT = LSTRT-K-KK+1
                 K = 1
              end if
              if (NTEXT > 3) then
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
              write (BUF(LSTRT:LBUF), FMTC) IMAG + J - 1
              if (NTEXT /= 0) NTEXT = K
              go to 300
           end if
  !                                  Set trailing blanks
           if (LSTRT <= LBUF) BUF(LSTRT:LBUF) = ' '
    300 continue
    310 return
      end subroutine messmh
  !*************************************************************************

  !*************************************************************************
  !>
  ! Prints the buffer for [[MESS]]

      subroutine MESSPR()

        integer   :: NSCRN !! Number of lines currently on CRT from messages.
        integer   :: K
        character :: SCRNAM*12
        save      NSCRN
  !
  ! For comments on other variables, see the listing for MESS.
        integer, parameter :: LENBUF = 250
        integer, parameter :: MEVBAS = 10
        integer, parameter :: MEVLAS = 33

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
                print '('' Type "Enter" to continue'')'
                read (*, *)
                NSCRN = 0
              end if
              NSCRN = NSCRN + 1
            end if
            print '(1X, A)', BUF(1:LBUF)
            if (OUNIT == 0) go to 20
          end if
  !++ CODE for ~.C. is active
          K = abs(OUNIT)
          write (K, '(A)', ERR=30) BUF(1:LBUF)
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
  !  call [[MESS]] to do the actions in MACT.  Actions in MACT can not do
  !  anything other than actions that reference MACT.
  !  This routine intended for use by library subroutines getting text in
  !  the form of a Fortran character string.

      subroutine MESSFT(MACT, FTEXT)

        integer MACT(*)
        character FTEXT*(*)
  !
        integer J, K, IDAT(1)
  !++ CODE for ~.C. is active
        character TEXT(1)*1
  !++ CODE for .C. is inactive
  !      character TEXT(1)*2
  !++ END
        integer, parameter :: MESUNI=10
        integer, parameter :: MEPRNT=21
        integer, parameter :: MECONT=50
  !
        integer, parameter :: LENBUF = 250
        integer, parameter :: MEVBAS = 10
        integer, parameter :: MEVLAS = 33

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

        do J = 1, 100, 2
           K = abs(MACT(J))
           if ((K > MEPRNT) .or. (K < MESUNI)) exit
        end do
        K = MACT(J)
        MACT(J) = MECONT
        call MESS(MACT, TEXT, IDAT)
        MACT(J) = K
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
           BUF(1:K) = FTEXT(NTEXT:NTEXT+K-1)
        end if
        ICHAR0 = ICHAR('0')
        if (MACT(J) /= MECONT) call mess(MACT(J), TEXT, IDAT)

      end subroutine messft
!*************************************************************************

!*************************************************************************
    end module messages_module
!*************************************************************************
