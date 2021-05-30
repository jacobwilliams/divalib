!*************************************************************************
!>
!  The main common block for the package.

    module divamc_module

    use diva_constants

    implicit none

    public

    double precision tg(2), tgstop(2), tmark, tmarkx, tout, tolg
    double precision alpha(kdim), beta(kdim+1)
    double precision d(maxstf+maxord,maxord), g(kdim,maxord)
    double precision v(kdim+maxord)
    double precision hc, hdec, hinc, hincc, hmax, hmaxp9, hmin
    double precision fdat(11)
    double precision ds(maxstf+maxord, maxord), gs(kdim)
    double precision sigma(kdim), rbq(kdim), dnoise
    double precision eave, eimax, eimin, emax, erep, robnd, snoise

    double precision  tmarka(2)
    integer kexit
    double precision hh
    integer ngstop(2)
    integer igflgs, izflag
    integer ivc2(65)
    double precision dvc2(7), rvc2(8)
    double precision dvc1(7)
    integer ioptc(23)

    integer icf,ics,igflg,igtype(2),igstop(2),ilgrep,ings,iop3,iop4,  &
    &   iop5,iop6,iop7,iop8,iop9,iop10,iop11,iop12,iop13,iop14,iop15,  &
    &   iop16,iop17,iop18,iop19,iop20,iop21,iop22,iop21s,itolep,iy,    &
    &   kemax,kis,kmark,kord1i,kord2i,kpred,kqdcon,kqicon,kqmaxs,      &
    &   kqmxds,kqmxil,kqmxip,kqmxis,ksc,ksout,ksstrt,kstep,lex,linc,   &
    &   lincd,lincq,lsc,maxkqd,maxkqi,method,ne,neptol,ng,ngtot,       &
    &   noiseq,noutko,ntolf,ny,idat(6)

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
    equivalence (kexit, iop17)
    equivalence (g(1, 1), hh)
    equivalence (ngstop(1), iop6)
    equivalence (izflag, iy)
    equivalence (igflgs, itolep)
    equivalence (ivc2(1), icf)
    equivalence (ioptc(3), iop3)

    end module divamc_module
!*************************************************************************
