C 
c     here a subroutine 
c     to calculate travel-time tables for hyposat.f
c
C              JOHANNES SCHWEITZER
C
C     (developed from laufps.f (version as of 28. April 1997))
c
c     11. June  1997
c
c     Nov 14, 1997   indph corrected for multiples and surface focus.
c
c     March 1999     changes due to CRUST 5.1 model input
c
c     October 2000   including converted surface reflections and
c                    converted reflections from the discontinuities 
c                    Conrad and Moho
c
c     March 2004     Usage of first order discontinuities changed
c
c     June 2004      subroutine efad included
c
c     September 2007 possible stations below surface included.
c
c     October 2010   several small changes and corrections with respect 
c                    to zero distance calculations
c
c     December 2017  CRUST 5.1 exchanged with CRUST 1.0
c
c     March 2018     changes for getting vp & vs at source depth
c
c     last changes/corrections 13. March 2018
c

      subroutine ttloc(ho,dis,czo1,nphas2,ttc,dtdd,dtdh,dpdh,dddp,
     +                 phcd,rmax,sdep,typctl,ierr,indph,emerout)

c
c     input:
c 
c             ho      = source depth in km (fixed or not, see fixho)
c
c             dis     = receiver distance in deg
c
c             czo1    = character to inform about fixed or not fixed 
c                       depth (D == no-fixed, other == fixed)
c
c             sdep    = if locgeo set and the station is below the
c                       Earth's surface, rays are only calculated
c                       until this depth.
c
c
c             typctl  = verbosity level
c
c                       only P-onsets at receiver
c
c             indph   = 00000 no travel-time calcluations
C
c             indph   = 10000 for only direct P-onsets
c                     = 11000 for P and pP
c                     = 11100 for P, pP, and PP 
c                     = 11110 for P, pP, PP, and PbP and PmP 
c                     = 11111 for P, pP, PP, PbP, PmP, and sP and SmP
c
c                       only S-onsets at receiver
c
c             indph   = 20000 for only direct S-onsets
c                     = 22000 for S and sS
c                     = 22200 for S, sS, and SS 
c                     = 22220 for S, sS, SS, and SbS and SmS 
c                     = 22222 for S, sS, SS, SbS, SmS, and pS and PmS 
c
c                       P- and S-onsets at receiver
c
c             indph   = 30000 for P and S (direct only)
c                     = 33000 for P, S, pP, and sS
c                     = 33300 for P, S, pP, sS, PP, and SS
c                     = 33330 for P, S, pP, sS, PP, SS, PbP, PmP, SbS,
c                                and Sms
c                     = 33333 for P, S, pP, sS, PP, SS, PbP, PmP, SbS.
c                                 Sms, pS, sP, PmS, and SmP
c
c             emerout = switch if emergence angles will be calculated 
c                       later
c
c	 in common MODEL :
c 
c             v0(1,i) =  P velocity in layer i
c             v0(2,i) =  S velocity in layer i
c
c             z(i)    =  depth of layer i
c
c             elev    =  topograhic elevation at modelled point 
c
c             elatc   =  latitude  to get the right CRUST 1.0 model 
c                        parameters
c             elonc   =  longitude to get the right CRUST 1.0 model 
c                        parameters
c
c             elev2   =  topograhic elevation at second modelled point 
c
c             elat2   =  latitude  to get the right CRUST 1.0 model 
c                        parameters at second modelled point
c             elon2   =  longitude to get the right CRUST 1.0 model 
c                        parameters at second modelled point
c
c             imo     <= 0
c                     <= 2 reading of model from filloc
c                 (   =  2 CRUST 1.0 only used for station
c                          corrections  )
c                     =  3 using model from CRUST 1.0
c                     =  4 dito + CRUST 1.0 used for travel-time 
c                          corrections
c
c             jmod    =  number of layers
c
c             iread   =  0 data of CRUST 1.0 model must be read in
c                     =  1 data of CRUST 1.0 model are already read
c
c             filloc  =  filename for file with local velocity model
c
c             azo(i)  =  Conrad/Moho indicator
c
c             mtyp    =  CRUST 1.0 model type
c
c             locgeo  =  false/true : if 'true', higher density of
c                        'rays' (IB) and smaller distances are used (DIS)
c            
c
c     output 
c
c             nphas2  = number of found onsets at receiver distance
c
c             phcd    = array with names of found onsets
c
c             ttc     = travel times of onsets in [sec]
c
c             dtdd    = ray parameters of onsets in [sec/deg]
c
c             dtdh    = partial derivative d(travel time)/d(ho) in 
c                       [sec/km]
c
c             dddp    = partial derivative d(ray parameter)/d(dis) in 
c                       [(sec/deg**2]
c
c             dpdh    = partial derivative d(ray parameter)/d(ho) in 
c                       [(sec/deg)/km]
c
c             ierr    = 0     everything o.k.
c                       else  some error occurred.
c
c             rmax    =  maximum distance for which this model shall be
c                        used (read in from model file).
c
c	 in common MODEL :
c 
c             v0(1,i) =  P velocity in layer i
c             v0(2,i) =  S velocity in layer i
c
c             z(i)    =  depth of layer i
c
c             rzv(1)  =  P velocity in source depth
c             rzv(2)  =  S velocity in source depth
c
c             elev    =  topograhic elevation at modelled point 
c
c             elatc   =  latitude  to get the right CRUST 1.0 model 
c                        parameters
c             elonc   =  longitude to get the right CRUST 1.0 model 
c                        parameters
c
c             elev2   =  topograhic elevation at second modelled point 
c
c             elat2   =  latitude  to get the right CRUST 1.0 model 
c                        parameters at second modelled point
c             elon2   =  longitude to get the right CRUST 1.0 model 
c                        parameters at second modelled point
c
c
c             imo     <= 0 no local model usage
c                     <= 2 reading of model from filloc
c                 (   =  2 CRUST 1.0 only used for station
c                          corrections  )
c                     =  3 using model from CRUST 1.0
c                     =  4 dito + CRUST 1.0 used for travel-time 
c                          corrections
c
c             jmod    =  number of layers
c
c             iread   =  0 data of CRUST 1.0 model must be read in
c                     =  1 data of CRUST 1.0 model are already read
c
c             filloc  =  filename for file with local velocity model
c
c             azo(i)  =  Conrad/Moho indicator
c
c             mtyp    =  CRUST 1.0 model type
c
c             zmax    =  maximum depth for which this model can be
c                        used (read in from model file).
c
c
      IMPLICIT real*8 (A-H,O-Z)

      REAL*4      ho,dis
      CHARACTER*1 czo1
      INTEGER     nphas2,ierr,indph,typctl
      logical     emerout

      include 'ttimes.h'

      dimension dpdh(mphas)

c
c     MAXLA   =    maximum number of allowed layers in model as defined in 
c
c                  model.h
c
c               (if needed then change also parameter MAXLA in subroutine reflex )
c
c
      include 'model.h'

      CHARACTER phas1*8,az(maxla)*4,line*34

c
c     NP        maximum number of calculated (defined) phases
c               is defined in phlist.h
c

      include 'phlist.h'

      integer   nphas1
      integer   indx(np)
      dimension tti(np)

      real*4    rtt(np)

      DIMENSION zo(3),pm(3,3),tm(3,3),h(maxla)

      include 'modelc.h'

      include 'ref.h'

      integer   phnum, iqq

      logical   fixho, kp, ks, surf, mul, surfc

      data phlist/'Pg','Pb','Pn','P','pPg','pPb','pPn','pP','PbP','PmP',
     +            'PgPg','PbPb','PnPn','PP','pSg','pSb','pSn','pS','PbS'
     +           ,'PmS','Sg','Sb','Sn','S','sSg','sSb','sSn','sS','SbS',
     +            'SmS','SgSg','SbSb','SnSn','SS','sPg','sPb','sPn','sP'
     +           ,'SbP','SmP',5*' '/

      SAVE      

      PHLIST (1)  = 'Pg'
      PHLIST (2)  = 'Pb'
      PHLIST (3)  = 'Pn'
      PHLIST (4)  = 'P'
      PHLIST (5)  = 'pPg'
      PHLIST (6)  = 'pPb'
      PHLIST (7)  = 'pPn'
      PHLIST (8)  = 'pP'
      PHLIST (9)  = 'PbP'
      PHLIST (10) = 'PmP'

      PHLIST (11) = 'PgPg'
      PHLIST (12) = 'PbPb'
      PHLIST (13) = 'PnPn'
      PHLIST (14) = 'PP'
      PHLIST (15) = 'pSg'
      PHLIST (16) = 'pSb'
      PHLIST (17) = 'pSn'
      PHLIST (18) = 'pS'
      PHLIST (19) = 'PbS'
      PHLIST (20) = 'PmS'

      PHLIST (21) = 'Sg'
      PHLIST (22) = 'Sb'
      PHLIST (23) = 'Sn'
      PHLIST (24) = 'S'
      PHLIST (25) = 'sSg'
      PHLIST (26) = 'sSb'
      PHLIST (27) = 'sSn'
      PHLIST (28) = 'sS'
      PHLIST (29) = 'SbS'
      PHLIST (30) = 'SmS'

      PHLIST (31) = 'SgSg'
      PHLIST (32) = 'SbSb'
      PHLIST (33) = 'SnSn'
      PHLIST (34) = 'SS'
      PHLIST (35) = 'sPg'
      PHLIST (36) = 'sPb'
      PHLIST (37) = 'sPn'
      PHLIST (38) = 'sP'
      PHLIST (39) = 'SbP'
      PHLIST (40) = 'SmP'

      PI=4.d0*DATAN(1.d0)

      PIM=PI/180.d0
      RE=6371.d0
      AA=PIM*RE

      IB = 20
      if(locgeo.or.dis.le.0.1) IB = 100
      IBN= IB*10

      ierr = 0

      if (indph.gt.0) then

         indph1 = indph  / 10000

         indphx = indph  - indph1*10000
         indph2 = indphx / 1000

         indphx = indphx - indph2*1000
         indph3 = indphx / 100

         indphx = indphx - indph3*100
         indph4 = indphx / 10

         indph5 = indphx - indph4*10

      endif

c      if (typctl.ge.8) then
c         print *,'Decoded INDPH ',indph1,indph2,indph3,indph4,indph5
c	 print *,'imo ',imo
c      endif

      if (imo.le.2) then

         if (jmod.gt.1 .and.rmax.gt.0.d0) go to 100

         OPEN (UNIT=15,FILE=trim(filloc),err=55)
C
C        Loop to get the models of P and S velocities
C

         read(15,*,err=55) rmax

         I=0
         jmod = 0

50       I=I+1
         IF(I.GT.60) THEN
             write(*,'('' Model contains too many layers (> 60)'')')
             ierr = 99
             close (15)
             GO TO 9000
         ENDIF

         READ(15,'(A)',err=55,end=56) line

         if (line(31:34).ne.'    ' .and. line(31:34).ne.'MOHO' .and.
     *       line(31:34).ne.'CONR' ) go to 55

         READ(line,'(3F10.3,A4)',err=55,end=56)
     *             Z(I),V0(1,i),V0(2,i),azo(I)

         if(Z(i) + v0(1,i) + v0(2,i).eq.0.d0) go to 56
         go to 50

55       write(*,'('' Read ERROR for file: '',a)') filloc
         ierr = 99
         close (15)
         go to 9000

56       jmod=I-1
         close (15)

         zmax  = Z(jmod)

         elev = 0.d0

         go to 100

      endif

      if (imo.ge.3)  then

         if(mtyp.ne.'C10') go to 9000

         itrue = 0
         inum = 2
         call get_mod_c10(itrue,inum,typctl,ierr)

         if(ierr.ne.0) then
            ierr = 99
            go to 9000
         endif

         rmax  = 1.5d0

         go to 100

      else


         print *,' No local/regional model defined! '
         ierr = 99
         go to 9000

      endif

100   continue

      imoh = 600
      icon = 600
      ipd  = 0
      isd  = 0

c
c     reset onset table
c
      do 110 i=1,np
      do 110 j=1,3
      do 110 k=1,3
110   ion(i,j,k) = 0

      if (indph.le.0) goto 9000

      if(ho.gt.sngl(zmax)) then
         ierr = 99
         print *,'Depth greater than maximum model depth'
         go to 9000
      endif

      if (czo1.eq.'D' .or. emerout) then
        fixho=.false.
        nzo  = 3
        jh1  = 2
        jh2  = 3
        zo(2) = dble(ho)
        if (locgeo) then
           zo(1) = zo(2) - 0.1d0
           zo(3) = zo(2) + 0.1d0
        else
           zo(1) = zo(2) - 1.d0
           zo(3) = zo(2) + 1.d0
        endif
        if(zo(1).lt.0.d0) zo(1)=0.d0
        if(zo(3).gt.zmax) zo(3)=zmax
      else
        fixho=.true.
        nzo = 1
        jh1 = 1
        jh2 = 1
        zo(1) = dble(ho)
      endif

      del(2) = dble(dis)

      if(locgeo) then
         del(1) = del(2)-0.001d0
      else
         del(1) = del(2)-0.01d0
      endif
      if(del(1).lt.0.d0) del(1)=0.d0

      if(locgeo) then
         del(3) = del(2)+0.001d0
      else
         del(3) = del(2)+0.01d0
      endif

      do 8000 iql = 1,nzo

C
C     Now model building and phase generation 
c
c     We will calculate at first P than S phases
c
C     (k-loop)
C

      ISS = -1

      DO 810 K=1,2

      ij = 0
      izo = 0

      IQQ = 0
c
c     Earth flattening approximation and interpolation of
c     source and eventually receiver layers
c
      DO 500 I=1,jmod

      i2  = i + 1
      ij = ij + 1

      az(ij) = azo(i)

      call efad (Z(I),V0(K,I),H(IJ),V(K,IJ))

      if(dabs(zo(iql)-Z(i)).lt.1.d-4 .and.izo.eq.0) then
         IQQ = IJ
         izo = 1
         if(iql.le.2) rzv(k) = v0(k,i)
         goto 499
      endif

      if(Z(i2).gt.zo(iql).and.Z(i).lt.zo(iql).and.izo.eq.0) then

         D=(zo(iql)-Z(i))*(V0(K,I2)-V0(K,I))/(Z(i2)-Z(i)) + V0(K,I)
         if(iql.le.2) rzv(k) = D

         ij = ij + 1
         IQQ = IJ
         az(ij) = ' '
         izo = 1
C
         call efad (zo(iql),D,H(IJ),V(K,IJ))

      endif

499   if(sdep.eq.Z(i).and.sdep.gt.0.d0) ISS = IJ

      if(Z(i2).gt.sdep.and.Z(i).lt.sdep.and.sdep.gt.0.d0) then

         D=(sdep-Z(i))*(V0(K,I2)-V0(K,I))/(Z(i2)-Z(i)) + V0(K,I)

         ij = ij + 1
         ISS = IJ
         az(ij) = ' '
C
         call efad (sdep,D,H(IJ),V(K,IJ))

      endif

      IF(V(1,IJ).GE.10.D0 .AND. IPD.EQ.0)  IPD  = IJ
      IF(V(2,IJ).GE.5.5D0 .AND. ISD.EQ.0)  ISD  = IJ

500   continue

      j = IJ
      m = j - 1

      DO 800 I=1,M

      I2=I+1
      V2(K,I)=V(K,I2)

      IF(AZ(I).EQ.'CONR')  ICON = I
      IF(AZ(I).EQ.'MOHO')  IMOH = I

      IF(dabs(V2(K,I)-V(K,I)).le.0.001d0) THEN
         V2(K,I)=1.0001d0*V(K,I)
         V(K,I2)=V2(K,I)
      ENDIF

      zdiff=H(I2)-H(I)
      ndisc(i) = 0
      IF(dabs(zdiff).le.0.0001d0)  then
         zdiff=0.0001d0
         H(i2)= H(i2) + zdiff
         ndisc(i) = 1
      endif

      G(K,I)=(V2(K,I)-V(K,I))/zdiff

800   continue
      
c      if(typctl.gt.5) then
c        print *,'i z h v(1) v(2) v2(1) v2(2) g(1) g(2)'
c        do 811 i=1,j
c        print*,i,az(i),z(i),h(i),v(1,i),v(2,i),v2(1,i),v2(2,i),
c     +  g(1,i),g(2,i)
c811      continue
c      endif

810   continue

      DO 7500 K=1,2

      if(k.eq.1) then
        kp=.true.
        ks=.false.
      else
        kp=.false.
        ks=.true.
      endif

      if(indph1.ne.3) then
        if(kp.and.indph1.ne.1) go to 7500
        if(ks.and.indph1.ne.2) go to 7500
      endif

      VHQ(K)=V(K,IQQ)*V(K,IQQ)
      PAD(K)=(RE-H(IQQ))*PIM/V(K,IQQ)

      conv  = .false.

      IF(IQQ.EQ.1)  GO TO 1000

C
C     direct waves ( if source deeper than 0. or
c     source at the surface and the station below)
c
c     plus defining requested direct phases
C

      if(kp) phase(1:1)='P'
      if(ks) phase(1:1)='S'

      if(iqq.le.imoh) then
         if(iqq.le.icon) phase(2:)='g     '
         if(iqq.gt.icon) phase(2:)='b     '
      else
         phase(2:)='n     '
      endif

      if((kp .and. (iqq.ge.ipd) .and. (ipd.ne.0)) .or.
     +   (ks .and. (iqq.ge.isd) .and. (isd.ne.0)) ) phase(2:)='      '

      if(iss.gt.iqq) then

         call reflex(iss,iqq,k)

      else if(iss.eq.iqq) then

         iph=phnum(phase)

         do 990 klr = 1,3

         del1 = del(klr)

         tt(2) = del1/v(k,iqq)

         ion(iph,iql,klr) = ion(iph,iql,klr)+1

         if(ion(iph,iql,klr).eq.1) then
            ttp(iph,iql,klr)=tt(2)
            ppp(iph,iql,klr)=PAD(k)
         else
            if(TT(2).lt.ttp(iph,iql,klr)) then
               ttp(iph,iql,klr)=TT(2)
               ppp(iph,iql,klr)=PAD(k)
            endif
         endif

990      continue

      else

         CALL REFLEX(IQQ,ISS,K)

      endif

C
C     body waves 
C

1000  imul=1111

      mul   = .false.
      surf  = .false.
      surfc = .false.

      IQ4=0
      MULT=0
      FFA=0.d0

      pa(1) = 0.d0
      pa(2) = 0.d0
      rr(1) = 0.d0
      rr(2) = 0.d0
      tt(1) = 0.d0
      tt(2) = 0.d0


1100  IQ5=1
      IQ6=0
      IF(IQQ.GT.1)  IQ5=IQQ
      IF(IQ4.GT.1)  IQ5=IQ4
      IF(imul.LT.M.AND.imul.GT.IQ6) IQ6=imul
      IQ5=MAX0(IQ6,IQ5)
      VMAX=V(K,IQ5)

      DO 1300 I=1,IQ5
      FA(1,I)=2.d0
      FA(2,I)=0.d0
      IF(I.GE.imul) FA(1,I)=FFA
      IF(I.LT.IQQ)  THEN
         FA(1,I)=FA(1,I)-1.d0
         if(i.lt.iss) FA(1,I)=FA(1,I)-1.d0
c        GO TO 1300
      ENDIF
      IF(I.LT.IQ4.AND.surf) FA(1,I)=FA(1,I)+1.d0
      IF(I.LT.IQ4.AND.surfc) then
         FA(2,I)=FA(2,I)+1.d0
         IF(k.eq.1 .and. VMAX.LT.V(2,I)) VMAX=V(2,I)
         IF(k.eq.2 .and. VMAX.LT.V(1,I)) VMAX=V(1,I)
      endif
      if(i.lt.iss) FA(1,I)=FA(1,I)-1.d0
      IF(fa(1,i).gt.0.9d0 .and. VMAX.LT.V(K,I)) VMAX=V(K,I)
1300  continue
C
C
      DO 3000 I=IQ5,M

      if (ndisc(i).ne.0) go to 3000

      ib2 = ib
      ib1 = 1
      ibm = 0

      FA(1,I)=2.d0
      IF(I.GE.imul) FA(1,I)=FFA
      if(i.lt.iss) FA(1,I)=FA(1,I)-1.d0

      if(kp) phase(1:1)='P'
      if(ks) phase(1:1)='S'

      if(i.lt.imoh) then
         if(i.le.icon) phase(2:)='g     '
         if(i.gt.icon) phase(2:)='b     '
      else
         phase(2:)='n     '
      endif

      if((kp .and. i.ge.ipd .and. ipd.ne.0) .or.
     +   (ks .and. i.ge.isd .and. isd.ne.0) ) phase(2:)='      '

      D=V2(K,I)
      IF(D.LE.VMAX)  GO TO    3000

      IF(imul.LT.M.AND.I.LT.imul) GO TO 2600

      C=V(K,I)
      IF(C.LT.VMAX) C=VMAX

1350  G1=DBLE(IB2-1)
      B=(D-C)/G1

      DO 2500 I2=ib1,IB2
      G2=dble(i2-1)
      VV=C+G2*B
      R=0.D0
      T=0.D0
C
      DO 2000 KK=1,I

      if(fa(1,kk).lt.0.9d0 .or. ndisc(kk).ne.0) go to 2000

      E=V(K,KK)
      G1=E/VV
      P=DSQRT(DABS(1.D0-G1*G1))
      O=1.d0/G(K,KK)

      IF(KK.GE.I)  THEN
         F=VV
         Q=0.d0
      else
         F=V2(K,KK)
         G3=F/VV
         Q=DSQRT(DABS(1.D0-G3*G3))
      ENDIF

      R=R+FA(1,KK)*(P-Q)*O
      DT=FA(1,KK)*DLOG(F*(1.D0+P)/(E*(1.D0+Q)))*O
      T=T+DT

2000  CONTINUE

c
c     extension for converted surface reflections (i.e., pS or sP)
c
      if(surfc) then

         if(k.eq.1) kc = 2
         if(k.eq.2) kc = 1

         DO 2001 KK=1,IQ4

         if(FA(2,KK).lt.0.9d0 .or. ndisc(kk).ne.0) go to 2001

         E=V(KC,KK)
         G1=E/VV
         P=DSQRT(DABS(1.D0-G1*G1))
         O=1.d0/G(KC,KK)

         F=V2(KC,KK)
         G3=F/VV
         Q=DSQRT(DABS(1.D0-G3*G3))

         R=R+FA(2,KK)*(P-Q)*O
         DT=FA(2,KK)*DLOG(F*(1.D0+P)/(E*(1.D0+Q)))*O
         T=T+DT

2001     CONTINUE

      endif

      RR(2) = R*VV/AA
      TT(2) = T

      PA(2) = AA/VV

      phas1 = phase
      iphase = len_trim(phase)

      if(mul) then
         phas1 = phase(1:iphase)//phase(1:iphase)
      endif

      if(surf) then
         if(kp) phas1 = 'p' // phase(1:iphase)
         if(ks) phas1 = 's' // phase(1:iphase)
      else if(surfc) then
         if(kp) phas1 = 's' // phase(1:iphase)
         if(ks) phas1 = 'p' // phase(1:iphase)
      endif

      iph=phnum(phas1)

      do 2400 klr = 1,3

      if(iql.ne.jh1 .and. klr.ne.2) go to 2400

      del1 = del(klr)

      IF (RR(2).EQ.del1) THEN

         if(ibm .eq. 0 ) then
            ib2 = ibn
            ib1 = (i2-1)*10
            if(ib1.le.0) ib1 = 1
            ibm = 1
            go to 1350
         endif
            
         ion(iph,iql,klr) = ion(iph,iql,klr)+1

         if(ion(iph,iql,klr).eq.1) then
            ttp(iph,iql,klr)=TT(2)
            ppp(iph,iql,klr)=PA(2)
         else
            if(TT(2).lt.ttp(iph,iql,klr)) then
               ttp(iph,iql,klr)=TT(2)
               ppp(iph,iql,klr)=PA(2)
            endif
         endif
         GO TO 2400
      ENDIF

      if(i2.le.1) go to 2400

      FCT1=DEL1-RR(1)
      FCT2=DEL1-RR(2)
      IF(FCT1*FCT2.LT.0.d0) THEN

         if(ibm .eq. 0 ) then
            ib2 = ibn
            ib1 = (i2-1)*10
            ibm = 1
            go to 1350
         endif
            
         FCT3=FCT1/(RR(2)-RR(1))
         TT1=FCT3*(TT(2)-TT(1))+TT(1)
         PA1=FCT3*(PA(2)-PA(1))+PA(1)

         ion(iph,iql,klr) = ion(iph,iql,klr)+1

         if(ion(iph,iql,klr).eq.1) then
            ttp(iph,iql,klr)=TT1
            ppp(iph,iql,klr)=PA1
         else
            if(TT1.lt.ttp(iph,iql,klr)) then
               ttp(iph,iql,klr)=TT1
               ppp(iph,iql,klr)=PA1
            endif
         endif
      ENDIF

2400  continue

      rr(1) = rr(2)
      tt(1) = tt(2)
      pa(1) = pa(2)

2500  continue
C
2600  VMAX=D
3000  CONTINUE

      if(mul) go to 3500

C
C     Now the surface reflections of the body waves will be
C     calculated (i.e. pP,sS...).
C
      IF(surf .or. iqq.eq.1) then
         if(surf) then
           IQQ=IQ4
           IQ4=1
           surf=.false.
         endif
         GO TO 3100
      endif

      if(indph2.ne.3) then
        if(kp.and.indph2.ne.1) go to 3100
        if(ks.and.indph2.ne.2) go to 3100
      endif

      surf = .true.
      IQ4=IQQ
      IQQ=1

      GO TO 1100

C
C     Now the converted surface reflections of the body waves 
C     will be calculated (i.e. pS,sP...).
C

3100  CONTINUE

      IF(surfc .or. iqq.eq.1) then
         if(surfc) then
           IQQ=IQ4
           IQ4=1
           surfc=.false.
         endif
         GO TO 3200
      endif

      if(indph5.ne.3) then
        if(kp.and.indph5.ne.1) go to 3200
        if(ks.and.indph5.ne.2) go to 3200
      endif

      surfc = .true.
      IQ4=IQQ
      IQQ=1

      GO TO 1100

C
C     Now the multiple phases will be done (e.g. PgPg, PnPn or SnSn...)
C     At the moment only single multiples can be calculated (i.e. e.g. 
C     no PPP or SnSnSn ...).
C

3200  continue

      if(indph3.ne.3) then
        if(kp.and.indph3.ne.1) go to 3500
        if(ks.and.indph3.ne.2) go to 3500
      endif

      imul = 1
      mult = 1
      mul = .true.

      FFA=dble(MULT*2+2)
      GO TO 1100
C
C     End of the body-phase and direct-wave loop
C

3500  CONTINUE

C
C     Reflections for P and S from the two possible layers:
C     the Conrad-discontinuity and the Mohorivicic-discontinuity.
C

      if(indph4.ne.3) then
        if(kp.and.indph4.ne.1) go to 6500
        if(ks.and.indph4.ne.2) go to 6500
      endif

      conv = .false.

      DO 6000 I=IQQ,J

      IF(I.NE.ICON .and. I.NE.IMOH) GO TO 6000

      if (i.eq.icon.and.iqq.lt.i) then
        if(kp) phase = 'PbP'
        if(ks) phase = 'SbS'
      else if (i.eq.imoh.and.iqq.lt.i) then
        if(kp) phase = 'PmP'
        if(ks) phase = 'SmS'
      endif

      I2=I
      CALL REFLEX(I2,iss,K)

6000  CONTINUE


6500  continue
c
C
C     Converted Reflections for P and S from the two possible layers:
C     the Conrad- the Mohorivicic-discontinuity.
C     PbS, SbP, PmS, SmP
C

      if(indph5.ne.3) then
        if(kp.and.indph5.ne.1) go to 7500
        if(ks.and.indph5.ne.2) go to 7500
      endif

      conv = .true.

      DO 7000 I=IQQ,J

      IF(I.NE.ICON .and. I.NE.IMOH) GO TO 7000

      if (i.eq.icon.and.iqq.lt.i) then
        if(kp) phase = 'SbP'
        if(ks) phase = 'PbS'
      else if (i.eq.imoh.and.iqq.lt.i) then
        if(kp) phase = 'SmP'
        if(ks) phase = 'PmS'
      endif

      I2=I
      CALL REFLEX(I2,iss,K)

7000  CONTINUE

      conv = .false.

7500  continue

c
c
c     Finally we have to do some interpolations
c
c

8000  continue

      nphas  = mphas + 1
      nphas1 = 0

      do 8800 i=1,np

      dtdh1 = 0.d0
      dtdh2 = 0.d0

      dpdh1 = 0.d0
      dpdh2 = 0.d0

      dpdd1 = 0.d0
      dpdd2 = 0.d0

      dh1   = 0.d0
      dh2   = 0.d0

      dd1   = 0.d0
      dd2   = 0.d0

      ni    = 0

      do 8700 j=1,nzo

      do 8500 k=1,3

      tm(j,k) = 0.d0
      pm(j,k) = 0.d0

      if (ion(i,j,k).eq.0) go to 8455

      tm(j,k) = ttp(i,j,k)
      pm(j,k) = ppp(i,j,k)

      if (fixho) go to 8450

      if(k.eq.2) then
        if(j.eq.1) then
          dtdh1 = tm(j,k)
          dpdh1 = pm(j,k)
          dh1   = zo(j) 
        else if(j.eq.2 .and. ion(i,j-1,k).eq.0) then
          dtdh1 = tm(j,k)
          dpdh1 = pm(j,k)
          dh1   = zo(j) 
        endif
        if(j.eq.3) then
          dtdh2 = tm(j,k)
          dpdh2 = pm(j,k)
          dh2   = zo(j)
        else if(j.eq.2 .and. ion(i,j+1,k).eq.0) then
          dtdh2 = tm(j,k)
          dpdh2 = pm(j,k)
          dh2   = zo(j)
        endif
      endif

8450  if(j.eq.jh1) then
        if(k.eq.1) then
          dpdd1 = pm(j,k)
          dd1   = del(k)
        else if(k.eq.2 .and. ion(i,j,k-1).eq.0) then
          dpdd1 = pm(j,k)
          dd1   = del(k)
        endif
        if(k.eq.3) then
          dpdd2 = pm(j,k)
          dd2   = del(k)
        else if(k.eq.2 .and. ion(i,j,k+1).eq.0) then
          dpdd2 = pm(j,k)
          dd2   = del(k)
        endif
      endif


      if(j.eq.jh1 .and. k.eq.2) then

        ni = 1
        nphas1= nphas1 + 1
        nphas = nphas  - 1

        tti(nphas1) = tm(j,k)
        dtdd(nphas) = pm(j,k)
        phcd(nphas) = phlist(i)

        if(i.ge.5) tti(nphas1) = tti(nphas1) + 1.0d-7

        if(len_trim(phcd(nphas)).ge.3) then
           tti(nphas1) = tti(nphas1) + 1.0d-6
           if(i.ge.35) tti(nphas1) = tti(nphas1) - 1.0d-6
        endif 
        if(len_trim(phcd(nphas)).ge.4) then
           tti(nphas1) = tti(nphas1) + 1.0d-6
           if(i.ge.31) tti(nphas1) = tti(nphas1) + 1.0d-6
        endif 

      endif

8455  if(j.eq.jh2 .and. k.eq.3 .and. ni.ne.0) then
        
        dtdh(nphas) = 0.d0
        dpdh(nphas) = 0.d0
        dddp(nphas) = 0.d0

        if(.not.fixho) then
           dh3 = dh2-dh1
           if(dh3.gt.0.d0) then
              dtdh(nphas) = (dtdh2-dtdh1) / dh3
              dpdh(nphas) = (dpdh2-dpdh1) / dh3
           endif
        endif

        dd3 = dd2-dd1
        if(dd3.gt.0.d0) then
           dddp(nphas) = (dpdd2-dpdd1) / dd3
        endif

        if(dabs(dddp(nphas)).ge.(dtdd(nphas)/2.d0)) then
           dddp(nphas) = 0.d0
        endif

        if(dabs(dpdh(nphas)).ge.(dtdd(nphas)/4.d0)) then
           dpdh(nphas) = 0.d0
        endif
      endif

8500  continue

8700  continue

      if(typctl.gt.8 .and. nphas1.gt.0 .and. ni.gt.0) then
         print *,'[hyposat_loc]',i,nphas,nphas1,phcd(nphas),
     *        tti(nphas1),dtdd(nphas),dtdh(nphas),dpdh(nphas),
     *        dddp(nphas)
      endif

8800  continue

      do 8850 i = 1,nphas1
8850  rtt(i) = sngl(tti(i))
      call indexx(nphas1,rtt,indx)

      do 8900 i=1,nphas1

      j       = indx(i)
      j2      = mphas + 1 - j

      ttc(i)  = tti(j)
      phcd(i) = phcd(j2)
      dtdd(i) = dtdd(j2)
      dtdh(i) = dtdh(j2)
      dddp(i) = dddp(j2)
      dpdh(i) = dpdh(j2)

      if(typctl.ge.8) then
         print *,i,j,j2,dis,phcd(i),ttc(i),dtdd(i),
     *          dtdh(i),dpdh(i),dddp(i)
      endif

8900  continue

      nphas2 = nphas1

9000  RETURN
      END
C
C
      SUBROUTINE  REFLEX(II,ISS,K)
      IMPLICIT real*8 (A-H,O-Z)
C
      PARAMETER (maxla=101,np=45)

      include 'ref.h'

      integer phnum,ii,k

      SAVE

      iph=phnum(phase)

      L=II-1

      if(conv) then

         VMAX=dmax1(V(1,IQQ),V(2,IQQ))

         if(k.eq.1) kc = 2
         if(k.eq.2) kc = 1

         DO  900  I=1,L

         FA(1,I)=1.D0
         FA(2,I)=0.D0

         IF(V(k,I) .GT.VMAX)  VMAX=V(k,I)
         IF(V2(k,I).GT.VMAX)  VMAX=V2(k,I)

         IF(I.GE.IQQ)  then
            FA(2,I)=1.D0
            IF(V(kc,I) .GT.VMAX)  VMAX=V(kc,I)
            IF(V2(kc,I).GT.VMAX)  VMAX=V2(kc,I)
         endif

         if(I.LT.ISS)  FA(1,I)=FA(1,I)-1.d0

900      CONTINUE

      else

         VMAX=V(K,IQQ)

         DO  1000  I=1,L
         FA(1,I)=2.D0
         IF(I.LT.IQQ)  FA(1,I)=FA(1,I)-1.D0
         if(I.LT.ISS)  FA(1,I)=FA(1,I)-1.d0
         IF(V(K,I) .GT.VMAX)  VMAX=V(K,I)
         IF(V2(K,I).GT.VMAX)  VMAX=V2(K,I)
1000     CONTINUE

      endif

      IBC=IB*30
      IF(IBC.GT.1200) IBC=1200

      RR(1)=0.d0
      RR(2)=0.d0
      PA(1)=0.d0
      PA(2)=0.d0
      TT(1)=0.d0
      TT(2)=0.d0

      B=PI/(2.d0*DBLE(IBC-1))

      DO  1500 I=1,IBC
      RVV=DSIN(B*DBLE(I-1))/VMAX
      T=0.D0
      R=0.D0

      DO  1100  KK=1,L
      if(FA(1,KK).lt.0.9d0 .or. ndisc(kk).ne.0) go to 1100
      E=V(K,KK)
      G1=E*RVV
      P=DSQRT(DABS(1.D0-G1*G1))
      O=1.d0/G(K,KK)
      F=V2(K,KK)
      G1=F*RVV
      Q=DSQRT(DABS(1.D0-G1*G1))
      R=R+FA(1,KK)*(P-Q)*O
      DT=FA(1,KK)*DLOG(F*(1.D0+P)/(E*(1.D0+Q)))*O
      T=T+DT
1100  CONTINUE

      if(conv) then

         DO  1110  KK=1,L
         if(FA(2,KK).lt.0.9d0 .or. ndisc(kk).ne.0) go to 1110
         E=V(KC,KK)
         G1=E*RVV
         P=DSQRT(DABS(1.D0-G1*G1))
         O=1.d0/G(KC,KK)
         F=V2(KC,KK)
         G1=F*RVV
         Q=DSQRT(DABS(1.D0-G1*G1))
         R=R+FA(2,KK)*(P-Q)*O
         DT=FA(2,KK)*DLOG(F*(1.D0+P)/(E*(1.D0+Q)))*O
         T=T+DT
1110     CONTINUE

       endif

      TT(2)=T

      IF(I.GT.1)  then
         RR(2)=R/(RVV*AA)
         P=DSQRT(DABS(1.D0/(RVV*RVV*VHQ(K))-1.D0))
         IF(P.LE.0.D0)  then
           FI=0.5d0*PI
         else
           FI=DATAN(1.D0/P)
         endif
      else
         FI    = 0.0D0
         rr(2) = 0.0d0
      endif

      IF(II.EQ.IQQ)  FI=pi-FI

      PA(2)=DSIN(FI)*PAD(K)
      if(pa(2).lt.1.d-4) pa(2)=0.0d0

      DO 1400 KLR=1,3

      if(iql.ne.jh1 .and. klr.ne.2) go to 1400

      del1 = del(klr)

      IF (dabs(RR(2)-del1).le.1.d-5) THEN

         ion(iph,iql,klr) = ion(iph,iql,klr)+1

         if(ion(iph,iql,klr).eq.1) then
            ttp(iph,iql,klr)=TT(2)
            ppp(iph,iql,klr)=PA(2)
         else
            if(TT(2).lt.ttp(iph,iql,klr)) then
               ttp(iph,iql,klr)=TT(2)
               ppp(iph,iql,klr)=PA(2)
            endif
         endif
         GO TO 1400
      ENDIF

      if(i.eq.1) go to 1400

      FCT1=del1-RR(1)
      FCT2=del1-RR(2)

      IF(FCT1*FCT2.LT.0.d0) THEN
         FCT3=FCT1/(RR(2)-RR(1))
         TT1=FCT3*(TT(2)-TT(1))+TT(1)
         PA1=FCT3*(PA(2)-PA(1))+PA(1)

         ion(iph,iql,klr) = ion(iph,iql,klr)+1

         if(ion(iph,iql,klr).eq.1) then
            ttp(iph,iql,klr)=TT1
            ppp(iph,iql,klr)=PA1
         else
            if(TT1.lt.ttp(iph,iql,klr)) then
               ttp(iph,iql,klr)=TT1
               ppp(iph,iql,klr)=PA1
            endif
         endif
      ENDIF

1400  CONTINUE

      tt(1) = tt(2)
      rr(1) = rr(2)
      pa(1) = pa(2)
 
1500  CONTINUE
C
      RETURN
      END
c
      function phnum(phase)

      include 'phlist.h'

      CHARACTER phase*8

      integer phnum

      SAVE

      phnum = 999
      do 5 i = 1,np
      if(phase.eq.phlist(i)) then
        phnum = i
        return
      endif
5     continue
      return
      end
