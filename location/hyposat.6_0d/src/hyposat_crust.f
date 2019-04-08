C
      FUNCTION CRUST(PA,PHTYP,FA,zo,TYPCTL0)
c
c     CRUST is a funtion to calculate travel-time differences 
c     for phases passing a standard crust and Model CRUST 1.0
c
C     JOHANNES SCHWEITZER
C
c     March 1999    
c
c     June 2004: subroutine EFAD included
c
c     February 2016: Common MODELG repaired
c
c     October 2017: Several adjustments needed to calculate also 
c                   corrections for Regional Models in source region and 
c                   not only Crust 5.1
c
c     December 2017: CRUST 5.1 exchanged with CRUST 1.0
c
c     input:
c 
c             pa      = ray parameter in [s/deg]
c
c             phtyp   = gives the phase type (P or S)
c
c             fa      = factor how often this phase passes the crust
c                       [i.e. : at station fa = 1, at reflection points
c                               (e.g., PP) fa = 2]
c
c             zo      = source depth
c
c             typctl0 = verbosity level
c
c             in COMMON MODEL Local Model:
c 
c             imo     <= 1 not used here
c                     =  2 CRUST 1.0 only used for travel-time 
c                           corrections at position elatc ,elonc 
c                           (station correction)
c                     =  3 not used here
c                     =  4 as (2) CRUST 1.0 used for travel-time 
c                          corrections
c
c             elatc   = latitude  to get the right CRUST 1.0 model 
c                       position for the station
c             elonc   = longitude to get the right CRUST 1.0 model 
c                       position for the station
c
c             iread   =  not used here
c
c             mtyp    =  CRUST 1.0 model type
c            
c             v0(1,i) =  P velocity in layer i
c             v0(2,i) =  S velocity in layer i
c
c             z(i)    =  depth of layer i
c
c             zmax    =  maximum model depth 
c
c             elev    =  topograhic elevation at modelled point 
c
c             jmod    =  number of layers
c
c             in COMMON MODELG parameter of standard model:
c 
c             v0g(1,i) =  P velocity in layer i
c             v0g(2,i) =  S velocity in layer i
c
c             zg(i)    =  depth of layer i
c
c             zmaxg    =  maximum model depth 
c
c             elevg    =  topograhic elevation at modelled point 
c
c             jmodg    =  number of layers
c
c
c
c     output:
c
c             crust   = travel-time correction for requested phase
c                       with respect to a standard model
c
      IMPLICIT real*8 (A-H,O-Z)
c
      save

      real*8   fa,pa,zo,crust

      INTEGER  typctl0,typctl

      include 'model.h'

      DIMENSION h(maxla),del(2),time(2),V(maxla),G(maxla),V2(maxla)

      CHARACTER phtyp*1

      include 'modelg.h'

      PI=4.d0*DATAN(1.d0)
      PIM=PI/180.d0
      re = 6371.d0
      AA=PIM*re

      typctl = typctl0

      if (mtyp.eq.'C10') then

         if (iabs(imo).ne.2 .and. imo.ne.4) go to 9998

         itrue = 1
         ierr = 0
         inum = 1
         call get_mod_c10(itrue,inum,typctl,ierr)

         if(ierr.ne.0) then
            go to 9000
         endif

      endif

      if(phtyp.eq.'P') iph = 1
      if(phtyp.eq.'S') iph = 2
C

      zma   = re - zmax - elev
      zmag  = re - zmaxg - elevg

      zmini = dmax1(zmag,zma)

c     z(jmod)   = re - zmin + elev
c     zg(jmodg) = re - zmin + elevg

      zmin  = re - zmini + elev
      zgmin = re - zmini + elevg

      do 7000 k = 1,2

      if(k.eq.1) jmodk = jmod
      if(k.eq.2) jmodk = jmodg

c
c     reset onset table
c
      del(k)  = 0.d0
      time(k) = 0.d0

      vmax = 0.d0

      i2 = 0
      iz = 0
      il = 0 

      DO 500 I=1,jmodk

      i2 = i2 + 1

      if(k.eq.1) then

         if(dabs(zo-z(i)).lt.1.d-3) iz = i

         if(i.gt.1 .and. iz.eq.0) then
            if(zo.lt.z(i).and.zo.gt.z(i-1)) then
               vnew = V0(iph,I-1) + 
     +               (V0(iph,I)-V0(iph,I-1)) * ((z(i)-zo)/(z(i)-z(i-1)))
               call efad(zo,vnew,H(I2),V(I2))
               iz = i2
               if(vmax.lt.v(i2)) vmax = v(i2)
               i2 = i2 + 1
            endif
         endif

         if(dabs(zmin-z(i)).lt.1.d-3) il = i

         if(i.gt.1 .and. il.eq.0) then
            if(zmin.lt.z(i).and.zmin.gt.z(i-1)) then
               vnew = V0(iph,I-1) + 
     +               (V0(iph,I)-V0(iph,I-1))   * 
     +               ((z(i)-zmin)/(z(i)-z(i-1)))
               call efad(zmin,vnew,H(I2),V(I2))
               il = i2
               if(vmax.lt.v(i2)) vmax = v(i2)
               i2 = i2 + 1
            endif
         endif

         call efad(z(I),V0(iph,I),H(I2),V(I2))

      else if(k.eq.2) then

         if(dabs(zo-zg(i)).lt.1.d-3) iz = i

         if(i.gt.1 .and. iz.eq.0) then
            if(zo.lt.zg(i).and.zo.gt.zg(i-1)) then
               vnew = V0g(iph,I-1) + 
     +               ((V0g(iph,I)-V0g(iph,I-1))  * 
     +              ((zg(i)-zo)/(zg(i)-zg(i-1))) )
               call efad(zo,vnew,H(I2),V(I2))
               iz = i2
               if(vmax.lt.v(i2)) vmax = v(i2)
               i2 = i2 + 1
            endif
         endif

         if(dabs(zgmin-zg(i)).lt.1.d-3) il = i

         if(i.gt.1 .and. il.eq.0) then
            if(zgmin.lt.zg(i).and.zgmin.gt.zg(i-1)) then
               vnew = V0g(iph,I-1) + 
     +               (V0g(iph,I)-V0g(iph,I-1))   * 
     +               ((zg(i)-zgmin)/(zg(i)-zg(i-1)))
               call efad(zgmin,vnew,H(I2),V(I2))
               il = i2
               if(vmax.lt.v(i2)) vmax = v(i2)
               i2 = i2 + 1
            endif
         endif

         call efad(zg(I),V0g(iph,I),H(I2),V(I2))

      endif

      if(vmax.lt.v(i2)) vmax = v(i2)

500   continue

      PA1=(re-H(i2))*PIM
      RVV = pa / PA1
      if(rvv*vmax.gt.0.985d0) go to 9998

      m = i2 - 1

      DO 800 I=1,M

      I2=I+1
      V2(I)=V(I2)

      IF(dabs(V2(I)-V(I)) .le. 0.001d0) THEN
         V2(I)=1.0001d0*V(I)
         V(I2)=V2(I)
      ENDIF

      zdiff=H(I2)-H(I)
      IF(dabs(zdiff).le.0.0001d0)  then
          zdiff=0.0001d0
          H(i2)= H(i2) + zdiff
      endif

      G(I)=(V2(I)-V(I))/zdiff

800   continue

      T=0.D0
      R=0.D0

      if(il.gt.0 .and. il.le.m) m = il-1

      if(iz.gt.0 .and. iz.le.m) m = iz-1

      DO  1100  KK=1,M
      E=V(KK)
      G1=E*RVV
      P=DSQRT(DABS(1.D0-G1*G1))
      O=1.d0/G(KK)
      F=V2(KK)
      G1=F*RVV
      Q=DSQRT(DABS(1.D0-G1*G1))
      R=R+FA*(P-Q)*O
      DT=FA*DLOG(F*(1.D0+P)/(E*(1.D0+Q)))*O
      T=T+DT
1100  CONTINUE

      del(k)  = r/(aa*rvv)
      time(k) = t

7000  continue

      ddel = del(1) - del(2)
      t1 = ddel * pa
      crust = time(1) - time(2) - t1

      if(typctl.ge.8) then
         print *,'crust-correction: ',pa,phtyp,ddel,'(',time(1),time(2),
     +     ')',t1,crust
      endif
      go to 9999

9000  continue
      print *,'crust: no CRUST 1.0 correction possible'
      ierr = 0
      
9998  crust = 0.d0
9999  RETURN
      END
c
c     funtion crustc
c
c     driver to calculate the travel-time effect 
c     for reflections at the Earth's surface for 
c     different crustal velocity structures and 
c     topography and calculates station corrections
c     for all kind of body phases.
c
c     Reflections can be corrected for the following principle 
c     phases:
c
c                 pP, sP, sS, pS 
c                 PP, SS, PS, SP
c                 P'P', S'S', P'S', S'P'
c 
c
      function crustc(phase_t0,rayp0,depth,ind,typctl)
c
c     input:
c              phase_t0 phase type (either P or S)
c
c              rayp0    ray parameter of phase in [s/deg]
c
c              depth    source depth
c
c              ind      switch for phase type
c
c                       = 1 station corrections for all phases
c
c                       = 2 surphase reflection pP, sS, PP, SS
c
c                       = 3 converted surphase reflection sP, pS, 
c                           SP & PS
c
c                       = 4 depth phase as 2, but source within the model layers
c
c                       = 5 depth pahse as 3, but source within the model layers
c
c     Several adjustments needed to calculate also corrections 
c     for Regional Models in source region and not only Crust 1.0
c
c     October 2017
c
c     some corrections spring/sommer 2018
c
c     calls FUNCTION CRUST
c

      implicit real*8 (a-h,o-z)

      real*8 rayp0, depth

      integer*4 ind, typctl

      character*1 phase_t0, phase_r

      real*8 crustc, crust

      real*8 fmult, rayp, zo

      crustc  = 0.d0
      crustc1 = 0.d0
      crustc2 = 0.d0

      fmult = 0.d0

      rayp = rayp0

      if(ind.eq.1) then

        zo      = 9999.d0
        fmult   = 1.0d0
        phase_r = phase_t0

        crustc  = crust(rayp,phase_r,fmult,zo,typctl)
        go to 9000

      else if(ind.eq.2) then

        zo      = 9999.d0
        fmult   = 2.d0
        phase_r = phase_t0

        crustc  = crust(rayp,phase_r,fmult,zo,typctl)
        go to 9000

      else if (ind.eq.3) then

        zo      = 9999.d0
        fmult   = 1.d0

        phase_r = 'P'
        crustc1 = crust(rayp,phase_r,fmult,zo,typctl)

        phase_r = 'S'
        crustc2 = crust(rayp,phase_r,fmult,zo,typctl)

        crustc = crustc1 + crustc2

        go to 9000

      else if(ind.eq.4) then

        fmult   = 2.d0
        phase_r = phase_t0

        zo      = depth
        crustc  = crust(rayp,phase_r,fmult,zo,typctl)

        go to 9000

      else if (ind.eq.5) then

        zo      = depth
        fmult   = 1.d0

        if(phase_t0 .eq. 'P') phase_r = 'P'
        if(phase_t0 .eq. 'S') phase_r = 'S'
        crustc1 = crust(rayp,phase_r,fmult,zo,typctl)

        if(phase_t0 .eq. 'P') phase_r = 'S'
        if(phase_t0 .eq. 'S') phase_r = 'P'
        crustc2 = crust(rayp,phase_r,fmult,zo,typctl)

        crustc = crustc1 + crustc2
        go to 9000

      endif

9000  return
      end
